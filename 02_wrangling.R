# =============================================================
# 02_wrangling.R  データの整形・列/行の追加・結合
# =============================================================
# 実際の解析時間の大半はここ。「解析できる形」にするまでの作法を学ぶ。
#   tidy data の3原則: 1列=1変数, 1行=1観測, 1セル=1値

library(tidyverse)   # dplyr, tidyr, readr, stringr, lubridate, ggplot2 ...

if (!file.exists("data/baseline_raw.csv")) source("01_simulate_data.R")

# ---- 1. 読み込みと「まず眺める」 ----------------------------------
raw <- read_csv("data/baseline_raw.csv", show_col_types = FALSE)

glimpse(raw)          # 列の型を確認（日付が chr になっている！）
summary(raw)          # age の最大 999 = 欠測コード
count(raw, sex)       # 表記ゆれ: "M" "Male" "m" " M" ...
raw |> count(id) |> filter(n > 1)   # 重複ID → 二重入力

# ---- 2. クリーニング --------------------------------------------
baseline <- raw |>
  distinct() |>                                     # 完全に同じ行を削除
  mutate(
    sex = sex |> str_trim() |> str_to_upper() |> str_sub(1, 1),  # → "M" / "F"
    sex = factor(sex, levels = c("F", "M")),
    age = na_if(age, 999),                          # 欠測コードを NA に
    enroll_date = ymd(enroll_date),                 # "2025/04/01" も "2025-04-01" も読める
    arm = factor(arm, levels = c("Placebo", "Drug")),  # 先頭 = 参照群
    dm  = factor(dm, levels = 0:1, labels = c("No", "Yes"))
  )

count(baseline, sex)
stopifnot(!anyDuplicated(baseline$id))             # 想定をコードで確認する習慣

# ---- 3. 新しい列（変数）の追加 -----------------------------------
baseline <- baseline |>
  mutate(
    bmi       = weight_kg / (height_cm / 100)^2,
    bmi_cat   = cut(bmi, breaks = c(-Inf, 18.5, 25, 30, Inf), right = FALSE,
                    labels = c("低体重", "普通", "肥満1度", "肥満2度以上")),
    age_group = if_else(age >= 75, "75歳以上", "75歳未満"),
    elderly   = case_when(age >= 75 ~ "後期高齢",
                          age >= 65 ~ "前期高齢",
                          !is.na(age) ~ "65歳未満"),     # どれにも当たらない → NA
    enroll_month = floor_date(enroll_date, "month")
  )
count(baseline, bmi_cat)

# ---- 4. 新しい行の追加 ------------------------------------------
# (a) add_row(): 小さな表に1行ずつ追加（マスタ表など）
drug_master <- tibble(code = c("D01", "D02"), name = c("Amlodipine", "Losartan"),
                      class = c("CCB", "ARB"))
drug_master <- drug_master |>
  add_row(code = "D03", name = "Hydrochlorothiazide", class = "Diuretic")
drug_master

# (b) bind_rows(): 同じ形の表を縦につなぐ（追加登録・別施設データ）
#     列が足りない所は NA、.id で出所を残せる
late_entry <- tibble(id = c("P201", "P202"), site = "H09",
                     arm = factor(c("Drug", "Placebo"), levels = levels(baseline$arm)),
                     age = c(71, 58),
                     sex = factor(c("M", "F"), levels = levels(baseline$sex)))
combined <- bind_rows(original = baseline, late = late_entry, .id = "source")
tail(combined |> select(source, id, site, arm, age, sex, bmi))
# ※ 以降の解析では元の200人を使う

# (c) complete(): 「無かった組み合わせ」を明示的な行として足す
#     → 欠測来院が NA の行として見えるようになる（欠測の把握に必須）
bp_wide <- read_csv("data/bp_wide.csv", show_col_types = FALSE)

bp_long <- bp_wide |>
  pivot_longer(-id,
               names_to = c("week", "rep"),
               names_pattern = "sbp_w(\\d+)_(\\d)",   # "sbp_w12_2" → week=12, rep=2
               names_transform = list(week = as.integer, rep = as.integer),
               values_to = "sbp")

bp_observed <- bp_long |> drop_na(sbp)             # よくある「測った行だけ」のデータ
nrow(bp_observed)
bp_full <- bp_observed |> complete(id, week, rep)  # 全組み合わせに戻す
nrow(bp_full)                                      # 200人 × 4来院 × 2回 = 1600

# (d) 集計表に「合計」行を足す
bp_full |>
  group_by(week) |>
  summarise(n_obs = sum(!is.na(sbp)), mean = mean(sbp, na.rm = TRUE)) |>
  bind_rows(summarise(bp_full, n_obs = sum(!is.na(sbp)), mean = mean(sbp, na.rm = TRUE))) |>
  mutate(week = as.character(week) |> replace_na("全体"))

# ---- 5. 縦持ち ⇄ 横持ち ------------------------------------------
# 来院ごとに2回測定の平均をとる
bp_visit <- bp_full |>
  group_by(id, week) |>
  summarise(n_meas = sum(!is.na(sbp)),
            sbp    = mean(sbp, na.rm = TRUE),       # 全部 NA だと NaN になる
            .groups = "drop") |>
  mutate(sbp = if_else(n_meas == 0, NA_real_, sbp)) |>
  group_by(id) |>
  mutate(sbp0   = sbp[week == 0],                   # ベースライン値を各行に
         change = sbp - sbp0) |>
  ungroup()

# 解析用に横持ちへ（1患者1行）
bp_by_pt <- bp_visit |>
  select(id, week, sbp) |>
  pivot_wider(names_from = week, values_from = sbp, names_prefix = "sbp_w")

# 検査値: 縦持ち → 横持ち
labs_wide <- read_csv("data/labs_long.csv", show_col_types = FALSE) |>
  pivot_wider(names_from = test, values_from = value)

# ---- 6. 結合 (join) ---------------------------------------------
outcome <- read_csv("data/outcome.csv", show_col_types = FALSE)

anti_join(labs_wide, baseline, by = "id")   # 0行ならOK: 背景表に無いIDの検査値はない

patients <- baseline |>
  left_join(labs_wide, by = "id") |>
  left_join(bp_by_pt,  by = "id") |>
  left_join(outcome,   by = "id") |>
  mutate(change12 = sbp_w12 - sbp_w0)

stopifnot(nrow(patients) == nrow(baseline))   # left_join で行が増えていないか必ず確認

# 経時データ（来院レベル）に患者情報を付ける
visit_long <- bp_visit |>
  left_join(patients |> select(id, site, arm, age, sex, dm, bmi), by = "id")

rep_long <- bp_full |>
  left_join(patients |> select(id, site, arm, age, sex, dm), by = "id")

# ---- 7. 患者背景表 (Table 1) -------------------------------------
fmt_msd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_np  <- function(x) sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))

table1 <- patients |>
  mutate(arm = as.character(arm)) |>
  bind_rows(patients |> mutate(arm = "全体")) |>   # 「全体」列のために行を複製して追加
  group_by(arm) |>
  summarise(
    n                 = as.character(n()),
    `年齢 mean (SD)`  = fmt_msd(age),
    `女性 n (%)`      = fmt_np(sex == "F"),
    `BMI mean (SD)`   = fmt_msd(bmi),
    `糖尿病 n (%)`    = fmt_np(dm == "Yes"),
    `SBP 0週 mean (SD)` = fmt_msd(sbp_w0),
    `12週欠測 n (%)`  = fmt_np(is.na(sbp_w12))
  ) |>
  pivot_longer(-arm, names_to = "項目") |>
  pivot_wider(names_from = arm, values_from = value) |>
  select(項目, Placebo, Drug, 全体)                 # 列の並びを指定
print(table1)
# RCT の Table 1 に p値は付けない（ランダム化しているので差は全て偶然）

# ---- 8. 保存 ----------------------------------------------------
# rds は factor の水準や日付型を保ったまま保存できる（csv は型情報が消える）
write_rds(patients,   "data/patients_clean.rds")
write_rds(visit_long, "data/visit_long.rds")
write_rds(rep_long,   "data/rep_long.rds")
cat("整形済みデータを data/*.rds に保存しました\n")
