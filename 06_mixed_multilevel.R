# =============================================================
# 06_mixed_multilevel.R  線形混合モデル・マルチレベルモデル
# =============================================================
# データの階層構造:
#   施設 (8) > 患者 (200) > 来院 (4回) > 測定 (2回)
# 同じ患者の測定値は似ている（独立ではない）→ 普通の回帰では SE が間違う
#
#   1. 普通の回帰 (lm) の問題
#   2. ランダム切片 → ランダム傾き → 施設も入れた3層モデル
#   3. 真のパラメータを回収できたか確認
#   4. 部分プーリング（縮約, shrinkage）を図で理解する
#   5. 測定レベルまで入れたモデル: 測定誤差と日間変動を分離（03章とつながる）
#   6. brms でベイズ版マルチレベルモデル

library(tidyverse)
library(lme4)
library(lmerTest)     # lmer に p値（Satterthwaite 自由度）を追加
library(performance)
library(emmeans)
library(brms)

if (!file.exists("data/visit_long.rds")) source("02_wrangling.R")
visit_long <- read_rds("data/visit_long.rds") |> drop_na(sbp)
rep_long   <- read_rds("data/rep_long.rds")   |> drop_na(sbp)
truth      <- read_rds("data/truth.rds")
theme_set(theme_minimal(base_size = 12))
options(mc.cores = 4, brms.backend = "rstan")
arm_cols <- c(Placebo = "#E69F00", Drug = "#0072B2")

# ---- 1. 普通の回帰 -------------------------------------------------
m_lm <- lm(sbp ~ week * arm, data = visit_long)
# 4回の来院を「別々の800人」とみなしている → 情報量を過大評価

# ---- 2. 混合モデル --------------------------------------------------
# (1 | id)     : 患者ごとに切片が違う（ランダム切片）
# (week | id)  : 患者ごとに切片も傾きも違う（ランダム傾き）
# (1 | site)   : 施設ごとに切片が違う
m1 <- lmer(sbp ~ week * arm + (1 | id),                data = visit_long)
ctrl <- lmerControl(optimizer = "bobyqa")   # 既定の最適化で収束警告が出たため変更
m2 <- lmer(sbp ~ week * arm + (week | id),             data = visit_long, control = ctrl)
m3 <- lmer(sbp ~ week * arm + (week | id) + (1 | site), data = visit_long, control = ctrl)

summary(m3)

# ランダム効果構造の比較（尤度比検定。lme4 が自動で ML に再推定する）
anova(m1, m2, m3)
# AIC/BIC も参考に。ただし分散成分の検定は境界問題で保守的になる点に注意

# 固定効果 week:armDrug（薬剤による週あたりの追加低下）の SE 比較
get_se <- function(m, term = "week:armDrug") coef(summary(m))[term, c("Estimate", "Std. Error")]
rbind(lm = get_se(m_lm), `lmer (1|id)` = get_se(m1),
      `lmer (week|id)` = get_se(m2), `+ (1|site)` = get_se(m3))
# → ここでは lm の SE の方が大きい。「時間による変化」は患者内の比較なので、
#   混合モデルは患者間のばらつき（SD 10 mmHg）を差し引いて精度よく推定できる。
#   逆に「群間差」のような患者間比較の項は、lm だと SE を過小評価しやすい（偽陽性が増える）。
#   lm の推定値自体もずれている: 血圧が高い人が脱落しやすい (MAR) ため。
#   混合モデルはその人の脱落前のデータも使うので MAR に強い（09章）。

# ICC: 全分散のうち患者間の違いの割合
icc(m1)

# ---- 3. 真値の回収 --------------------------------------------------
vc <- as.data.frame(VarCorr(m3))
tibble(
  パラメータ = c("週あたり変化 (Placebo)", "薬剤の追加効果/週", "12週での群間差",
                 "SD 患者切片", "SD 患者傾き", "SD 施設", "SD 残差"),
  真値 = c(truth$b_week, truth$b_trt, 12 * truth$b_trt,
           truth$sd_pat, truth$sd_slope, truth$sd_site,
           sqrt(truth$sd_visit^2 + truth$sd_meas^2 / 2)),   # 来院変動 + 2回平均の測定誤差
  推定値 = c(fixef(m3)["week"], fixef(m3)["week:armDrug"], 12 * fixef(m3)["week:armDrug"],
             vc$sdcor[vc$grp == "id" & vc$var1 == "(Intercept)" & is.na(vc$var2)],
             vc$sdcor[vc$grp == "id" & vc$var1 == "week" & is.na(vc$var2)],
             vc$sdcor[vc$grp == "site"],
             vc$sdcor[vc$grp == "Residual"])
)
# 患者切片の SD は年齢・DM の影響も含むので真値より大きめになる（モデルに入れていないため）
# 傾きの SD は来院4回では推定が難しく、0 に近く出やすい（境界付近 → 収束警告の原因にもなる）

# 12週時点の群間差（推定周辺平均）
emm <- emmeans(m3, ~ arm | week, at = list(week = 12))
emm
contrast(emm, "revpairwise")
confint(m3, parm = "week:armDrug", method = "Wald")

# 残差の診断
plot(m3)                                # 残差 vs 予測値
qqnorm(resid(m3)); qqline(resid(m3))
qqnorm(ranef(m3)$id[, 1]); qqline(ranef(m3)$id[, 1])   # ランダム効果の正規性

# ---- 4. 部分プーリング（縮約） ----------------------------------------
# 3つの考え方:
#   完全プーリング: 全員同じ傾き（個人差を無視）
#   プーリングなし: 患者ごとに別々に回帰（データの少ない人は不安定）
#   部分プーリング: 混合モデル。データの少ない人ほど全体平均に引き寄せられる
no_pool <- visit_long |>
  group_by(id) |>
  filter(n() >= 2) |>
  summarise(n_visits = n(), arm = first(arm),
            slope_nopool = coef(lm(sbp ~ week))[2])

partial <- coef(m2)$id |>
  as_tibble(rownames = "id") |>
  mutate(slope_partial = week + if_else(id %in% visit_long$id[visit_long$arm == "Drug"],
                                        `week:armDrug`, 0)) |>
  select(id, slope_partial)

shrink <- inner_join(no_pool, partial, by = "id")

p_shrink <- shrink |>
  pivot_longer(c(slope_nopool, slope_partial), names_to = "method", values_to = "slope") |>
  mutate(method = recode(method, slope_nopool = "患者ごとの回帰", slope_partial = "混合モデル")) |>
  ggplot(aes(method, slope, group = id, colour = factor(n_visits))) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.5) +
  facet_wrap(~ arm) +
  scale_colour_manual(values = c(`2` = "#D55E00", `3` = "#E69F00", `4` = "grey60"),
                      name = "来院回数") +
  labs(title = "部分プーリング: 推定値は全体平均に引き寄せられる",
       subtitle = "来院回数が少ない（情報が少ない）患者ほど大きく縮む",
       x = NULL, y = "傾き (mmHg/週)")
print(p_shrink)
ggsave("output/06_shrinkage.png", p_shrink, width = 8, height = 5, dpi = 150, bg = "white")

# ---- 5. 測定レベルまで入れた4層モデル --------------------------------
# (1 | id:week): 来院ごとの変動（その日の血圧） / 残差: 同じ来院内の測定誤差
m_rep <- lmer(sbp ~ week * arm + (week | id) + (1 | id:week) + (1 | site), data = rep_long,
              control = ctrl)
VarCorr(m_rep)
cat(sprintf("来院変動 SD: 推定 %.2f (真値 %g) / 測定誤差 SD: 推定 %.2f (真値 %g)\n",
            attr(VarCorr(m_rep)$`id:week`, "stddev"), truth$sd_visit,
            sigma(m_rep), truth$sd_meas))
# 反復測定があれば「測定誤差」と「生物学的変動」を分けて推定できる
# → 何回測れば十分か、どちらを減らすべきか（測定法の改善 vs 測定日を増やす）の判断材料

# ---- 6. brms: ベイズ版マルチレベルモデル ------------------------------
# 利点: 施設が8つしかない → 施設SDの推定が不安定。事前分布で穏やかに正則化できる
#       不確実性がすべて事後分布として得られる（ランダム効果の SD も）
priors_mm <- c(
  prior(normal(0, 10), class = b),                 # 固定効果
  prior(student_t(3, 0, 10), class = sd),          # ランダム効果の SD（半t分布）
  prior(exponential(0.1), class = sigma),
  prior(lkj(2), class = cor)                       # 切片と傾きの相関
)
fit_mm <- brm(sbp ~ week * arm + (week | id) + (1 | site),
              data = visit_long, prior = priors_mm,
              seed = 1, refresh = 0, control = list(adapt_delta = 0.95),
              iter = 4000, warmup = 1000,   # 2000 だと sd(week) の ESS 不足の警告が出た
              file = "output/models/b06_mixed", file_refit = "on_change")
summary(fit_mm)
# 必ず確認: Rhat < 1.01, Bulk/Tail_ESS > 400, divergent transitions なし
#   警告が出たら → iter を増やす / adapt_delta を上げる / 事前分布やモデルを見直す

# lmer と比べる
fixef(fit_mm)
fixef(m3)

# 臨床的な問い: 「12週で群間差が 5 mmHg を超える確率は？」
hypothesis(fit_mm, "12 * week:armDrug < -5")

# 施設SDの事後分布（lmer では点推定だけ）
as_draws_df(fit_mm) |>
  ggplot(aes(sd_site__Intercept)) +
  geom_density(fill = "#56B4E9", alpha = 0.6) +
  geom_vline(xintercept = truth$sd_site, linetype = "dashed") +
  labs(title = "施設間SD の事後分布", subtitle = "破線 = 真値。施設が8つだけなので不確実性が大きい",
       x = "SD (mmHg)", y = NULL)

# 予測の図
plot(conditional_effects(fit_mm, effects = "week:arm"), plot = FALSE)[[1]] +
  scale_colour_manual(values = arm_cols) + scale_fill_manual(values = arm_cols) +
  labs(x = "週", y = "SBP (mmHg)")

pp_check(fit_mm, ndraws = 50)

# 発展: 2値アウトカムの混合モデル
#   glmer(y ~ x + (1 | id), family = binomial)   /  brm(y ~ x + (1 | id), family = bernoulli())
