# =============================================================
# 04_frequentist.R  頻度論的な統計
# =============================================================
# 臨床的な問い: 「降圧薬は12週後の血圧をプラセボより下げるか？」
#   1. 記述統計と図（検定の前に必ずデータを見る）
#   2. 2群比較: t検定 / Wilcoxon
#   3. ANCOVA（ベースライン調整）— RCTの連続アウトカムの標準
#   4. 回帰診断
#   5. 2値アウトカム: リスク差・リスク比・オッズ比, ロジスティック回帰
#   6. p値と信頼区間をシミュレーションで理解する
#   7. 多重比較・サブグループ解析の罠
#   8. 検出力とサンプルサイズ

library(tidyverse)
library(broom)
library(pwr)

if (!file.exists("data/patients_clean.rds")) source("02_wrangling.R")
pts <- read_rds("data/patients_clean.rds")
truth <- read_rds("data/truth.rds")
theme_set(theme_minimal(base_size = 12))
arm_cols <- c(Placebo = "#E69F00", Drug = "#0072B2")   # Okabe–Ito（色覚多様性に配慮）

cc <- pts |> filter(!is.na(sbp_w12))   # 12週の値がある人（欠測の扱いは 09 章）
cat("真の効果 (12週):", truth$b_trt * 12, "mmHg\n")

# ---- 1. 記述統計 ------------------------------------------------
cc |>
  group_by(arm) |>
  summarise(n = n(),
            w0 = mean(sbp_w0), w12 = mean(sbp_w12),
            change_mean = mean(change12), change_sd = sd(change12),
            change_median = median(change12))

ggplot(cc, aes(arm, change12, colour = arm)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  scale_colour_manual(values = arm_cols, guide = "none") +
  labs(x = NULL, y = "12週の変化量 (mmHg)")

# ---- 2. 2群比較 --------------------------------------------------
t.test(change12 ~ arm, data = cc)                      # Welch（等分散を仮定しない; 既定）
t.test(change12 ~ arm, data = cc, var.equal = TRUE)    # Student
wilcox.test(change12 ~ arm, data = cc, conf.int = TRUE)  # 順位ベース（位置のずれ）
# ※ factor の水準順で「Placebo − Drug」になっている点に注意（符号を必ず確認）

# ---- 3. ANCOVA ---------------------------------------------------
# 12週値 ~ 群 + ベースライン値
#   変化量の t 検定より検出力が高く、ベースラインの偶然の偏りにも頑健
fit_change <- lm(change12 ~ arm, data = cc)
fit_ancova <- lm(sbp_w12 ~ arm + sbp_w0, data = cc)
fit_adj    <- lm(sbp_w12 ~ arm + sbp_w0 + age + sex + dm + site, data = cc)

bind_rows(
  change = tidy(fit_change, conf.int = TRUE),
  ancova = tidy(fit_ancova, conf.int = TRUE),
  adjust = tidy(fit_adj,    conf.int = TRUE),
  .id = "model"
) |>
  filter(term == "armDrug") |>
  select(model, estimate, std.error, conf.low, conf.high, p.value)
# → 標準誤差 (std.error) を比べる: 変化量より ANCOVA の方が小さい（精度が高い）
#   さらに共変量を足しても、予後への影響が小さい変数（ここでは施設など）ばかりだと
#   自由度を消費するだけで SE はむしろ増えることもある。
#   調整する共変量は「強い予後因子」を解析計画で事前に決める（結果を見て選ばない）

summary(fit_ancova)

# ---- 4. 回帰診断 ------------------------------------------------
par(mfrow = c(2, 2)); plot(fit_ancova); par(mfrow = c(1, 1))
# 残差 vs 予測値: パターンがないか / QQ: 正規性 / Leverage: 影響の大きい点

# ---- 5. 2値アウトカム（追跡期間中のイベント） --------------------
tab <- table(pts$arm, pts$event, dnn = c("arm", "event"))
tab
chisq.test(tab)
fisher.test(tab)          # 期待度数が小さいとき

risk <- prop.table(tab, 1)[, "1"]
rd <- risk["Drug"] - risk["Placebo"]
rr <- risk["Drug"] / risk["Placebo"]
or <- (tab["Drug", "1"] / tab["Drug", "0"]) / (tab["Placebo", "1"] / tab["Placebo", "0"])
cat(sprintf("リスク: Placebo %.2f, Drug %.2f\nRD = %.3f, RR = %.2f, OR = %.2f, NNT = %.0f\n",
            risk["Placebo"], risk["Drug"], rd, rr, or, abs(1 / rd)))
# イベントが多いとき OR は RR より 1 から遠くなる（OR ≈ RR は稀なイベントのときだけ）

fit_logit <- glm(event ~ arm + age + dm, family = binomial, data = pts)
tidy(fit_logit, exponentiate = TRUE, conf.int = TRUE)   # exp(係数) = オッズ比
# 本当は追跡期間が人によって違うので生存時間解析が適切 → 09章
# 血圧は確実に下がったのに、イベントは「有意差なし」。
#   → 効果がない証拠？ それとも検出力不足？ p値だけでは区別できない → 05章

# ---- 6. p値と信頼区間をシミュレーションで理解する ------------------
set.seed(10)
# 6a. 帰無仮説が正しいとき p値は一様分布 → 5% は偶然 p < 0.05
p_null <- replicate(5000, t.test(rnorm(30), rnorm(30))$p.value)
p_alt  <- replicate(5000, t.test(rnorm(30), rnorm(30, 0.5))$p.value)
cat(sprintf("H0 真: p<.05 の割合 = %.3f / 効果量0.5: 検出力 = %.3f\n",
            mean(p_null < .05), mean(p_alt < .05)))
tibble(p = c(p_null, p_alt), 状況 = rep(c("効果なし", "効果あり (d=0.5)"), each = 5000)) |>
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.025), fill = "#0072B2", colour = "white") +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  facet_wrap(~ 状況, scales = "free_y") +
  labs(title = "p値の分布", x = "p値", y = NULL)

# 6b. 95%信頼区間: 「同じ研究を何度も繰り返すと 95% の区間が真値を含む」
ci_sims <- map_dfr(1:100, function(i) {
  x <- rnorm(25, mean = -9.6, sd = 15)
  ci <- t.test(x)$conf.int
  tibble(sim = i, est = mean(x), lo = ci[1], hi = ci[2])
}) |>
  mutate(covers = lo <= -9.6 & hi >= -9.6)
mean(ci_sims$covers)
ggplot(ci_sims, aes(y = sim, x = est, xmin = lo, xmax = hi, colour = covers)) +
  geom_vline(xintercept = -9.6, linetype = "dashed") +
  geom_pointrange(size = 0.15) +
  scale_colour_manual(values = c(`TRUE` = "grey60", `FALSE` = "#D55E00"),
                      labels = c(`TRUE` = "真値を含む", `FALSE` = "含まない"), name = NULL) +
  labs(title = "100回の研究の95%信頼区間", x = "推定値 (mmHg)", y = "研究")
# 1つの区間について「真値が95%の確率で入っている」とは言えない → それはベイズの信用区間（05章）

# ---- 7. 多重比較・サブグループ解析 -------------------------------
# 7a. 効果のない検査項目を20個調べると…
set.seed(14)
p20 <- replicate(20, t.test(rnorm(50), rnorm(50))$p.value)
tibble(raw = p20, bonferroni = p.adjust(p20, "bonferroni"),
       holm = p.adjust(p20, "holm"), BH = p.adjust(p20, "BH")) |>
  arrange(raw) |> head()

# 7b. サブグループ: 各群で検定して「DMでは有意, 非DMでは非有意」は根拠にならない
#     → 交互作用項で「効果の差」を直接推定する
cc |> group_by(dm) |>
  group_modify(~ tidy(lm(sbp_w12 ~ arm + sbp_w0, data = .x), conf.int = TRUE)) |>
  filter(term == "armDrug")
tidy(lm(sbp_w12 ~ arm * dm + sbp_w0, data = cc), conf.int = TRUE)
# 真のデータ生成では効果修飾はない（arm:dmYes の真値 = 0）

# ---- 8. 検出力・サンプルサイズ ------------------------------------
# 差 8 mmHg, SD 15, α=0.05, 検出力 80%
power.t.test(delta = 8, sd = 15, power = 0.8)
pwr.t.test(d = 8 / 15, power = 0.8)

# シミュレーションによる検出力（どんな解析にも使える万能法）
sim_power <- function(n_per_arm, effect = -8, sd = 15, r = 0.6, nsim = 1000) {
  mean(replicate(nsim, {
    arm <- rep(0:1, each = n_per_arm)
    base <- rnorm(2 * n_per_arm, 150, sd)
    w12 <- 150 + r * (base - 150) + effect * arm + rnorm(2 * n_per_arm, 0, sd * sqrt(1 - r^2))
    summary(lm(w12 ~ arm + base))$coefficients["arm", 4] < 0.05
  }))
}
set.seed(4)
tibble(n = c(20, 30, 40, 60)) |> mutate(power_ancova = map_dbl(n, sim_power, nsim = 500))
# ベースラインとの相関 r が高いほど ANCOVA は少ない人数で済む
