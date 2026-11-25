# =============================================================
# 05_bayes_null.R  ベイズ統計 — 特に「Null の強さ」
# =============================================================
#   事後分布 ∝ 尤度 × 事前分布
#
# 「Null の強さ」には2つの意味がある。両方扱う。
#   (A) データが Null（効果なし）をどれだけ支持するか
#       → p値ではできない。ベイズファクター BF01、ROPE で定量化できる
#   (B) Null 寄りの事前分布（懐疑的事前分布）をどれだけ強くするか
#       → 結論が事前分布にどれだけ依存するか（感度分析）
#
# 構成:
#   1. 格子近似で「事前 × 尤度 = 事後」を目で見る（Stan不要）
#   2. 有意差なし ≠ 効果なし: イベントのハザード比をベイズで読む
#   3. ベイズファクター（BayesFactor パッケージ）
#   4. サンプルサイズと BF01 / p値（シミュレーション）
#   5. brms: 事前分布を明示した回帰、ROPE、Savage–Dickey BF

library(tidyverse)
library(survival)
library(BayesFactor)
library(bayestestR)
library(brms)

if (!file.exists("data/patients_clean.rds")) source("02_wrangling.R")
pts <- read_rds("data/patients_clean.rds")
cc  <- pts |> filter(!is.na(sbp_w12)) |> as.data.frame()   # BayesFactor は data.frame を要求
theme_set(theme_minimal(base_size = 12))
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)

# Stan のバックエンド: 既定は rstan。cmdstanr が正しく設定できていれば
#   options(brms.backend = "cmdstanr") の方が速い。
# 初回はモデルのコンパイルに1分ほどかかる。file= で結果を保存するので2回目以降は一瞬。
options(mc.cores = 4, brms.backend = "rstan")

# ---- 1. 格子近似: 事前 × 尤度 = 事後 -------------------------------
# 頻度論の ANCOVA の推定値と SE を「尤度」として使う（正規近似）
fit <- lm(sbp_w12 ~ arm + sbp_w0, data = cc)
est <- unname(coef(fit)["armDrug"])      # unname: 名前が付いたままだと後で計算結果の名前が変わる
se  <- sqrt(vcov(fit)["armDrug", "armDrug"])

theta <- seq(-25, 10, by = 0.01)                # 効果の候補値 (mmHg)
priors <- list(
  `弱情報 N(0, 20)`      = dnorm(theta, 0, 20),
  `懐疑的 N(0, 3)`       = dnorm(theta, 0, 3),
  `強い懐疑的 N(0, 1)`   = dnorm(theta, 0, 1)    # 「ほぼ効かないだろう」という強い信念
)
lik <- dnorm(est, theta, se)

grid <- imap_dfr(priors, function(pr, name) {
  post <- pr * lik
  post <- post / sum(post)
  tibble(prior_name = factor(name, levels = names(priors)), theta, prior = pr / sum(pr), likelihood = lik / sum(lik), posterior = post)
})

grid |>
  group_by(prior_name) |>
  summarise(post_mean = sum(theta * posterior),
            lo95 = theta[which(cumsum(posterior) >= 0.025)[1]],
            hi95 = theta[which(cumsum(posterior) >= 0.975)[1]],
            `P(効果<0)`  = sum(posterior[theta < 0]),
            `P(効果<-5)` = sum(posterior[theta < -5]))   # 臨床的に意味のある効果の確率
# → 「P(効果 < -5 mmHg) = 0.9」のような、臨床判断に直結する言い方ができる

grid |>
  pivot_longer(c(prior, likelihood, posterior)) |>
  mutate(name = factor(name, c("prior", "likelihood", "posterior"),
                       c("事前分布", "尤度（データ）", "事後分布"))) |>
  ggplot(aes(theta, value, colour = name)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ prior_name, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#999999", "#E69F00", "#0072B2"), name = NULL) +
  labs(title = "事前分布の強さで事後分布はどう変わるか",
       subtitle = "データが十分あれば、懐疑的な事前分布も押し切られる",
       x = "治療効果 (mmHg)", y = NULL) +
  theme(axis.text.y = element_blank())
ggsave("output/05_prior_posterior.png", width = 7, height = 7, dpi = 150, bg = "white")

# 共役（解析的）に解くとこうなる: 正規事前 × 正規尤度 → 正規事後
#   事後の精度 = 事前の精度 + データの精度（精度 = 1/分散）
post_normal <- function(m0, s0, est, se) {
  prec <- 1 / s0^2 + 1 / se^2
  c(mean = (m0 / s0^2 + est / se^2) / prec, sd = sqrt(1 / prec))
}
post_normal(0, 3, est, se)

# ---- 2. 有意差なし ≠ 効果なし: イベントのハザード比 ----------------
cox <- coxph(Surv(time_wk, event) ~ arm + age + dm, data = pts)
b   <- unname(coef(cox)["armDrug"]); b_se <- sqrt(vcov(cox)["armDrug", "armDrug"])
cat(sprintf("Cox: HR = %.2f, 95%%CI [%.2f, %.2f], p = %.3f  ← 「有意差なし」\n",
            exp(b), exp(b - 1.96 * b_se), exp(b + 1.96 * b_se), 2 * pnorm(-abs(b / b_se))))

# 事前分布: log HR ~ N(0, 0.5)  （HR が 0.37〜2.7 に 95% — 穏やかに懐疑的）
pn <- post_normal(0, 0.5, b, b_se)
cat(sprintf("事後: P(HR < 1) = %.2f, P(HR < 0.8) = %.2f, P(HR > 1.1) = %.3f\n",
            pnorm(0, pn["mean"], pn["sd"]),
            pnorm(log(0.8), pn["mean"], pn["sd"]),
            1 - pnorm(log(1.1), pn["mean"], pn["sd"])))

# Savage–Dickey: 点Null (logHR = 0) の BF01 = 事後密度(0) / 事前密度(0)
bf01_hr <- dnorm(0, pn["mean"], pn["sd"]) / dnorm(0, 0, 0.5)
cat(sprintf("BF01 = %.2f  → 1 に近い = データは Null とも効果ありとも決めきれない\n", bf01_hr))
# p=0.1 は「効果なし」の証拠ではなく「データ不足」。ベイズならそれが数字で見える。

# ---- 3. ベイズファクター ------------------------------------------
# BF10 = p(データ | H1) / p(データ | H0) 。BF01 = 1 / BF10
# 目安 (Jeffreys/Lee & Wagenmakers): 1–3 弱い, 3–10 中程度, 10–30 強い, >30 非常に強い
#
# 3a. 薬剤の効果（真に効果あり）
bf_drug <- ttestBF(formula = change12 ~ arm, data = cc)
bf_drug

# 3b. 性別による変化量の違い（データ生成上、真に差なし）
bf_sex <- ttestBF(formula = change12 ~ sex, data = cc)
t.test(change12 ~ sex, data = cc)$p.value   # p値 → 「差があるとは言えない」まで
1 / bf_sex                                   # BF01 → 「差がない」方を何倍支持するか

# 3c. 回帰モデル同士の比較（ANCOVA に sex を入れるべきか）
bf_full <- lmBF(sbp_w12 ~ arm + sbp_w0 + sex, data = cc)
bf_red  <- lmBF(sbp_w12 ~ arm + sbp_w0,       data = cc)
bf_red / bf_full                             # sex なしモデルを何倍支持するか

# 3d. 事前分布の幅 (rscale) への感度 —— BF は事前分布に敏感！
#     rscale = 効果量 δ の Cauchy 事前分布の尺度。大きいほど「大きな効果」を想定
sens <- tibble(rscale = c(0.2, 0.5, 0.707, 1, 1.414)) |>
  mutate(BF01_sex = map_dbl(rscale, \(r) 1 / extractBF(ttestBF(formula = change12 ~ sex,
                                                              data = cc, rscale = r))$bf),
         BF10_drug = map_dbl(rscale, \(r) extractBF(ttestBF(formula = change12 ~ arm,
                                                           data = cc, rscale = r))$bf))
sens
# 広すぎる事前分布 → H1 が「ありえない大きな効果」に確率を割くので H0 が有利になる
#   (Jeffreys–Lindley のパラドックス)。事前分布は報告し、感度分析をする。

# ---- 4. サンプルサイズと BF01 / p値 -----------------------------------
# H0 が真のとき: p値はいくら n を増やしても一様分布のまま（Null の証拠にならない）
#                BF01 は n とともに大きくなる（Null の証拠が蓄積する）
set.seed(5)
sim_bf <- expand_grid(n = c(10, 20, 50, 100, 200, 500), sim = 1:200) |>
  mutate(res = map(n, function(n) {
    x <- rnorm(n); y <- rnorm(n)
    tibble(p = t.test(x, y)$p.value,
           bf01 = 1 / extractBF(ttestBF(x, y))$bf)
  })) |>
  unnest(res)

sim_bf |>
  group_by(n) |>
  summarise(`p<.05 の割合` = mean(p < .05),
            `BF01 中央値`  = median(bf01),
            `BF01>3 の割合` = mean(bf01 > 3),
            `BF10>3 の割合` = mean(bf01 < 1/3))

ggplot(sim_bf, aes(factor(n), bf01)) +
  geom_hline(yintercept = c(1/3, 1, 3), linetype = c("dashed", "solid", "dashed"), colour = "grey50") +
  geom_boxplot(fill = "#56B4E9", alpha = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(title = "H0 が真のときの BF01", subtitle = "n が増えるほど Null を支持する証拠が強くなる",
       x = "1群あたりの n", y = "BF01（対数軸）")
ggsave("output/05_bf01_by_n.png", width = 7, height = 4.5, dpi = 150, bg = "white")

# ---- 5. brms: 事前分布を明示した回帰 ----------------------------------
# ベースラインを中心化しておくと切片が解釈しやすく、事前分布も置きやすい
cc$sbp_w0_c <- cc$sbp_w0 - mean(cc$sbp_w0)

priors_b <- c(
  prior(normal(0, 10), class = b),          # 回帰係数: ±20 mmHg くらいまでを想定
  prior(student_t(3, 140, 20), class = Intercept),
  prior(exponential(0.1), class = sigma)
)

# まず事前予測チェック: この事前分布はありえないデータを生まないか？
fit_prior <- brm(sbp_w12 ~ arm + sbp_w0_c + sex, data = cc, prior = priors_b,
                 sample_prior = "only", seed = 1, refresh = 0,
                 file = "output/models/b05_prior")
pp_check(fit_prior, ndraws = 50) + labs(title = "事前予測チェック")

fit_b <- brm(sbp_w12 ~ arm + sbp_w0_c + sex, data = cc, prior = priors_b,
             sample_prior = "yes", seed = 1, refresh = 0,
             iter = 11000, warmup = 1000,   # BF を安定させるため多めに (4鎖 × 10000 = 40000 draws)
             file = "output/models/b05_ancova", file_refit = "on_change")
summary(fit_b)          # Rhat ≈ 1.00, ESS が十分かを必ず確認
pp_check(fit_b, ndraws = 50) + labs(title = "事後予測チェック")
plot(fit_b, variable = c("b_armDrug", "b_sexM"))   # トレースプロット

# 事後分布の要約
describe_posterior(fit_b, ci = 0.95, rope_range = c(-3, 3), test = c("p_direction", "rope"))

# 5a. ROPE（実質的同等域）: ±3 mmHg を「臨床的に無視できる差」と事前に定義
#     95%HDI が ROPE の外 → 効果あり / 完全に内 → 実質的に効果なし（Null を受容）
equivalence_test(fit_b, range = c(-3, 3), ci = 0.95)
plot(equivalence_test(fit_b, range = c(-3, 3), parameters = c("b_armDrug", "b_sexM")))
# 頻度論の対応物: TOST（2つの片側検定）による同等性検定

# 5b. Savage–Dickey 法の BF（点Null）と、区間Null（ROPE）に対する BF
bayesfactor_parameters(fit_b, parameters = c("b_armDrug", "b_sexM"), null = 0)
bayesfactor_parameters(fit_b, parameters = c("b_armDrug", "b_sexM"), null = c(-3, 3))
# 出力の BF は「Null に反する証拠」(BF10)。1 未満なら Null 寄りで、BF01 = 1 / BF。
#   点Null「ちょうど0」より区間Null「±3 以内」の方が、Null を強く支持しやすい。
#   臨床的には「ちょうど0」より「無視できるほど小さい」の方が自然な Null。

# 5c. 臨床的な問いを直接確率で答える
hypothesis(fit_b, c("armDrug < -5", "abs(sexM) < 3"))
draws <- as_draws_df(fit_b)
mean(draws$b_armDrug < -5)

# まとめ:
#   p値        : 「H0 のもとでこのデータ以上に極端な結果が出る確率」 Null の支持は言えない
#   BF01       : Null 対 対立仮説の相対的な証拠。事前分布の幅に敏感
#   ROPE/HDI   : 「臨床的に無視できる範囲」に入る確率。ROPE の幅は臨床で決める
#   事後確率   : P(効果 < -5 mmHg) のように臨床判断に直結
