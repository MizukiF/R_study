# =============================================================
# 01_simulate_data.R  疑似データの生成
# =============================================================
# 架空の多施設ランダム化比較試験 (RCT)
#   降圧薬 (Drug) vs プラセボ (Placebo)
#   8施設・200人・0/4/8/12週に来院・各来院で収縮期血圧 (SBP) を2回測定
#
# ポイント: 「真の値」を自分で決めてからデータを作る。
#   → 解析結果が真の値をどれだけ当てられるかを確かめられる。
#   → シミュレーションは統計手法を理解する最強の道具。
#
# データ生成モデル (患者 i, 施設 s, 週 t, 測定 r):
#   真の血圧    mu_it  = b0 + u_s + u_i + b_age*(age-65) + b_dm*dm
#                        + (b_week + b_trt*Drug + v_i) * t
#   来院日の血圧 x_it  = mu_it + e_visit        (日による生物学的変動)
#   測定値      y_itr = x_it  + e_meas          (1回ごとの測定誤差)

library(dplyr)
library(tidyr)
library(readr)

set.seed(2026)
dir.create("data", showWarnings = FALSE)

# ---- 1. 真のパラメータ ------------------------------------------
truth <- list(
  b0       = 150,   # 65歳・DMなし・0週の平均SBP (mmHg)
  b_age    = 0.3,   # 年齢 +1歳あたり
  b_dm     = 4,     # 糖尿病あり
  b_week   = -0.2,  # プラセボ群の週あたり変化 (mmHg/週)
  b_trt    = -0.8,  # 薬剤による「追加の」週あたり変化 → 12週で -9.6 mmHg
  sd_site  = 4,     # 施設間SD
  sd_pat   = 10,    # 患者間SD（切片）
  sd_slope = 0.3,   # 患者間SD（傾き）
  sd_visit = 5,     # 来院日ごとの変動SD
  sd_meas  = 6      # 測定誤差SD
)
saveRDS(truth, "data/truth.rds")

# ---- 2. 施設と患者 ----------------------------------------------
site_n <- c(32, 30, 27, 25, 24, 22, 21, 19)  # 合計200
sites <- tibble(site   = sprintf("H%02d", seq_along(site_n)),
                u_site = rnorm(length(site_n), 0, truth$sd_site))

pts <- tibble(site = rep(sites$site, site_n)) |>
  mutate(
    id      = sprintf("P%03d", row_number()),
    age     = round(rnorm(n(), 65, 10)),
    sex     = sample(c("F", "M"), n(), replace = TRUE),
    height  = round(if_else(sex == "M", rnorm(n(), 168, 6), rnorm(n(), 155, 6)), 1),
    bmi     = rnorm(n(), 24, 3.5),
    weight  = round(bmi * (height / 100)^2, 1),
    dm      = rbinom(n(), 1, plogis(-1.2 + 0.2 * (bmi - 24))),
    u_pat   = rnorm(n(), 0, truth$sd_pat),
    u_slope = rnorm(n(), 0, truth$sd_slope)
  ) |>
  group_by(site) |>
  mutate(arm = sample(rep(c("Placebo", "Drug"), length.out = n()))) |>  # 施設内で1:1割付
  ungroup() |>
  left_join(sites, by = "site")

# ---- 3. 来院ごとの血圧 ------------------------------------------
visits <- pts |>
  expand_grid(week = c(0, 4, 8, 12)) |>        # 患者 × 来院 の全組み合わせ
  mutate(
    true_sbp  = truth$b0 + u_site + u_pat +
                truth$b_age * (age - 65) + truth$b_dm * dm +
                (truth$b_week + truth$b_trt * (arm == "Drug") + u_slope) * week,
    visit_sbp = true_sbp + rnorm(n(), 0, truth$sd_visit)
  ) |>
  # 脱落: 前回の血圧が高い人ほど脱落しやすい (MAR: 観測値に依存する欠測)
  group_by(id) |>
  arrange(week, .by_group = TRUE) |>
  mutate(
    p_drop  = if_else(week == 0, 0, plogis(-3.2 + 0.06 * (lag(visit_sbp, default = 150) - 150))),
    dropped = cumsum(rbinom(n(), 1, p_drop)) > 0
  ) |>
  ungroup()

# 1来院2回測定
meas <- visits |>
  expand_grid(rep = 1:2) |>
  mutate(sbp = round(visit_sbp + rnorm(n(), 0, truth$sd_meas)),
         sbp = if_else(dropped, NA_real_, sbp))

saveRDS(meas, "data/sim_truth_long.rds")  # 真の値つき（答え合わせ用）

# ---- 4. 「現場っぽい」形式で書き出す ----------------------------
# 4a. 血圧: 横持ち (wide) — Excel で入力したような形
bp_wide <- meas |>
  mutate(key = paste0("sbp_w", week, "_", rep)) |>
  select(id, key, sbp) |>
  pivot_wider(names_from = key, values_from = sbp)
write_csv(bp_wide, "data/bp_wide.csv", na = "")

# 4b. 患者背景: 表記ゆれ・欠測コード・日付形式の混在・重複行を含む
messy_sex <- function(s) {
  n <- length(s)
  ifelse(s == "M",
         sample(c("M", "Male", "m", " M"), n, TRUE, prob = c(.7, .1, .1, .1)),
         sample(c("F", "Female", "f"),     n, TRUE, prob = c(.8, .1, .1)))
}
baseline_raw <- pts |>
  transmute(id, site, arm, age, sex = messy_sex(sex),
            height_cm = height, weight_kg = weight, dm,
            enroll_date = as.Date("2025-04-01") + sample(0:300, n(), TRUE)) |>
  mutate(enroll_date = if_else(runif(n()) < 0.3,
                               format(enroll_date, "%Y/%m/%d"),
                               format(enroll_date, "%Y-%m-%d")),
         age       = if_else(runif(n()) < 0.03, 999, age),          # 999 = 不明
         weight_kg = if_else(runif(n()) < 0.03, NA_real_, weight_kg))
baseline_raw <- bind_rows(baseline_raw, slice(baseline_raw, c(15, 120)))  # 二重入力
write_csv(baseline_raw, "data/baseline_raw.csv", na = "")

# 4c. 検査値: 縦持ち (long) — 電子カルテの抽出データのような形
labs <- pts |>
  transmute(id,
            HbA1c = round(if_else(dm == 1, rnorm(n(), 7.4, 0.9), rnorm(n(), 5.6, 0.4)), 1),
            LDL   = round(rnorm(n(), 120, 30)),
            Cr    = round(rlnorm(n(), log(0.85), 0.25), 2)) |>
  pivot_longer(-id, names_to = "test", values_to = "value") |>
  slice_sample(prop = 0.97)              # 一部の検査は未実施
write_csv(labs, "data/labs_long.csv")

# 4d. 機器比較: 0週に従来機 (A, 2回) と新しい手首式血圧計 (B, 1回)
#     B には固定バイアスと比例バイアスを仕込んである
device <- meas |>
  filter(week == 0) |>
  select(id, rep, sbp, visit_sbp) |>
  pivot_wider(names_from = rep, values_from = sbp, names_prefix = "devA_") |>
  mutate(devB = round(35 + 0.78 * visit_sbp + rnorm(n(), 0, 6))) |>
  select(id, devA_1, devA_2, devB)
write_csv(device, "data/device_compare.csv")

# 4e. アウトカム: 心血管イベントまでの時間 (週)。12週時点の血圧が高いほど起きやすい
outcome <- pts |>
  left_join(visits |> filter(week == 12) |> select(id, true_sbp12 = true_sbp), by = "id") |>
  mutate(rate    = 0.012 * exp(0.08 * (true_sbp12 - 140) + 0.6 * dm + 0.03 * (age - 65)),
         t_event = rexp(n(), rate),
         t_cens  = runif(n(), 30, 52),
         time_wk = round(pmin(t_event, t_cens), 1),
         event   = as.integer(t_event <= t_cens)) |>
  select(id, time_wk, event)
write_csv(outcome, "data/outcome.csv")

cat("data/ に疑似データを書き出しました\n")
