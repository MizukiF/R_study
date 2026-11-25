# =============================================================
# 03_measurement_error.R  測定誤差と再現性
# =============================================================
#   観測値 = 真の値 + 誤差
# 疑似データなので「真の値」を知っている → 誤差の影響を直接見られる。
# 扱うこと:
#   1. 繰り返し測定から測定誤差を推定 (within-subject SD, 再現性係数, CV)
#   2. ICC（級内相関係数）: 同日再現性と日間再現性は別物
#   3. 機器間の一致: Bland–Altman（相関係数は一致の指標ではない）
#   4. 回帰の希釈 (regression dilution): 誤差は関連を弱く見せる
#   5. 平均への回帰 (regression to the mean)
#   6. 推定値そのものの再現性: 小さい研究の ICC はどれくらいブレるか

library(tidyverse)
library(irr)
library(lme4)

if (!file.exists("data/rep_long.rds")) source("02_wrangling.R")
truth  <- read_rds("data/truth.rds")
rep_long <- read_rds("data/rep_long.rds")
device <- read_csv("data/device_compare.csv", show_col_types = FALSE)
theme_set(theme_minimal(base_size = 12))
dir.create("output", showWarnings = FALSE)

# ---- 1. 同じ来院での2回測定から誤差を推定 -------------------------
w0 <- rep_long |>
  filter(week == 0) |>
  select(id, rep, sbp) |>
  pivot_wider(names_from = rep, values_from = sbp, names_prefix = "m")

d  <- w0$m1 - w0$m2
sw <- sqrt(mean(d^2) / 2)            # within-subject SD（測定誤差のSD）
rc <- 1.96 * sqrt(2) * sw            # 再現性係数: 同じ人の2回の差の95%はこの範囲
cv <- 100 * sw / mean(c(w0$m1, w0$m2))

cat(sprintf("within-subject SD = %.2f (真値 %g)\n", sw, truth$sd_meas))
cat(sprintf("再現性係数 = %.1f mmHg / CV = %.1f%%\n", rc, cv))
# 解釈: 1回測定で「前回より10 mmHg下がった」は誤差の範囲かもしれない

# ---- 2. ICC -----------------------------------------------------
# ICC = 被験者間分散 / (被験者間分散 + 誤差分散)
#   → 「全体のばらつきのうち、人の違いで説明される割合」
icc(w0 |> select(m1, m2), model = "twoway", type = "agreement", unit = "single")
icc(w0 |> select(m1, m2), model = "twoway", type = "agreement", unit = "average")  # 2回平均の信頼性

# 混合モデルの分散成分からも同じものが出る（06章への伏線）
m_icc <- lmer(sbp ~ 1 + (1 | id), data = filter(rep_long, week == 0))
vc <- as.data.frame(VarCorr(m_icc))
cat(sprintf("ICC (lmer) = %.3f\n", vc$vcov[1] / sum(vc$vcov)))

# 日間再現性: プラセボ群の 0週 vs 4週（各2回平均）
#   来院日の変動 (sd_visit) も「誤差」側に入るので ICC は下がる
day2day <- rep_long |>
  filter(arm == "Placebo", week %in% c(0, 4)) |>
  group_by(id, week) |> summarise(sbp = mean(sbp), .groups = "drop") |>
  pivot_wider(names_from = week, values_from = sbp, names_prefix = "w") |>
  drop_na()
icc(day2day |> select(w0, w4), model = "twoway", type = "agreement", unit = "single")
# → 「同じ日に2回測ると再現性が良い」≠「日を変えても再現する」
#   何の再現性を知りたいのか（測定器? 患者の状態?）を研究デザインで決める

# ---- 3. 機器間の一致: Bland–Altman --------------------------------
ba <- device |>
  mutate(devA = (devA_1 + devA_2) / 2,
         avg  = (devA + devB) / 2,
         diff = devB - devA)

cat(sprintf("相関係数 r = %.3f  ← 高いが…\n", cor(ba$devA, ba$devB)))
bias <- mean(ba$diff)
loa  <- bias + c(-1.96, 1.96) * sd(ba$diff)
cat(sprintf("平均差 (bias) = %.1f, 一致限界 = [%.1f, %.1f]\n", bias, loa[1], loa[2]))

# 比例バイアス: 差が平均値の大きさで変わるか
summary(lm(diff ~ avg, data = ba))$coefficients

p_ba <- ggplot(ba, aes(avg, diff)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = bias, colour = "#0072B2", linewidth = 0.8) +
  geom_hline(yintercept = loa, colour = "#0072B2", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, colour = "#E69F00", linewidth = 0.8) +
  annotate("text", x = max(ba$avg), y = c(bias, loa), hjust = 1, vjust = -0.4,
           label = c("平均差", "-1.96 SD", "+1.96 SD"), size = 3.5) +
  labs(title = "Bland–Altman プロット: 手首式 (B) − 従来機 (A)",
       subtitle = "オレンジ線の傾き = 比例バイアス（高血圧ほど B が低めに出る）",
       x = "(A + B) / 2  [mmHg]", y = "B − A  [mmHg]")
print(p_ba)
ggsave("output/03_bland_altman.png", p_ba, width = 7, height = 5, dpi = 150, bg = "white")

# ---- 4. 回帰の希釈 ----------------------------------------------
# 真の血圧 x → アウトカム y の真の傾き 0.5。x を誤差つきで測ると傾きはどうなる？
set.seed(1)
n <- 5000
x_true <- rnorm(n, 150, 15)
y      <- 20 + 0.5 * x_true + rnorm(n, 0, 10)
sd_err <- 10

dilution <- map_dfr(c(1, 2, 4, 8), function(k) {
  x_obs  <- x_true + rowMeans(matrix(rnorm(n * k, 0, sd_err), n, k))  # k回測定の平均
  lambda <- 15^2 / (15^2 + sd_err^2 / k)                              # 信頼性係数
  slope  <- coef(lm(y ~ x_obs))[2]
  tibble(k = k, reliability = lambda, slope_obs = slope,
         slope_expected = 0.5 * lambda, slope_corrected = slope / lambda)
})
print(dilution)
# 測定誤差のある説明変数 → 傾きは 0 方向に引っ張られる（過小評価）
#   対策: 繰り返し測定の平均 / 信頼性で補正 / SIMEX / 誤差モデル（ベイズで自然に書ける）
# 注意: 「誤差があるのに有意 → 本当はもっと強い」は単回帰でしか成り立たない。
#       共変量に誤差があると交絡調整が不十分になり、どちら向きにも歪みうる。

# ---- 5. 平均への回帰 --------------------------------------------
# 0週の1回目が 160 以上の人だけ選ぶ（「高血圧の人を組み入れる」状況）
# → 同じ日の2回目は、何もしていないのに下がって見える
rtm <- w0 |>
  mutate(selected = m1 >= 160) |>
  group_by(selected) |>
  summarise(n = n(), mean_m1 = mean(m1), mean_m2 = mean(m2), diff = mean(m2 - m1))
print(rtm)
# 臨床研究で: 「高値の人に介入したら下がった」は介入効果とは限らない
#   → 対照群が必要 / 組み入れ測定と解析のベースラインを別にする / 複数回測定

# ---- 6. 推定値の再現性: ICC 推定のブレ ----------------------------
# 一元配置分散分析による ICC(1) を手で書く
icc1 <- function(y) {   # y: n×k 行列
  k <- ncol(y); n <- nrow(y)
  msb <- k * var(rowMeans(y))
  msw <- sum((y - rowMeans(y))^2) / (n * (k - 1))
  (msb - msw) / (msb + (k - 1) * msw)
}
sim_icc <- function(n, k = 2, icc_true = 0.8) {
  u <- rnorm(n, 0, sqrt(icc_true))
  e <- matrix(rnorm(n * k, 0, sqrt(1 - icc_true)), n, k)
  icc1(u + e)
}
set.seed(2)
icc_sims <- expand_grid(n = c(10, 30, 100), sim = 1:1000) |>
  mutate(icc = map_dbl(n, sim_icc))
icc_sims |> group_by(n) |>
  summarise(mean = mean(icc), q2.5 = quantile(icc, .025), q97.5 = quantile(icc, .975))

p_icc <- ggplot(icc_sims, aes(icc)) +
  geom_histogram(bins = 50, fill = "#0072B2", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  facet_wrap(~ n, ncol = 1, labeller = label_both) +
  labs(title = "真の ICC = 0.8 のとき、推定値はどれくらいブレるか",
       subtitle = "1000回のシミュレーション。n=10 の再現性研究は信頼できない", x = "推定 ICC", y = NULL)
print(p_icc)
ggsave("output/03_icc_sampling.png", p_icc, width = 6, height = 6, dpi = 150, bg = "white")
