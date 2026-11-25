# =============================================================
# 07_bootstrap.R  ブートストラップ
# =============================================================
# 考え方: 手元のデータを「母集団の代わり」とみなし、復元抽出を繰り返して
#         推定値のばらつき（標準誤差・信頼区間）を作る。
#   → 公式が無い統計量（中央値の差, ICC, 比, AUC…）でも CI が出せる
#
#   1. 手作りブートストラップ（中身を理解する）
#   2. boot パッケージ: パーセンタイル法と BCa 法
#   3. 再現性指標（ICC・一致限界）の CI
#   4. クラスターブートストラップ: 反復測定データでは「患者」を単位に再抽出
#   5. パラメトリックブートストラップ (lme4::bootMer)
#   6. ベイジアンブートストラップ（05章とのつながり）
#   7. 注意点

library(tidyverse)
library(boot)
library(lme4)

if (!file.exists("data/patients_clean.rds")) source("02_wrangling.R")
pts <- read_rds("data/patients_clean.rds")
cc  <- pts |> filter(!is.na(sbp_w12))
visit_long <- read_rds("data/visit_long.rds") |> drop_na(sbp)
rep_long   <- read_rds("data/rep_long.rds")
theme_set(theme_minimal(base_size = 12))

# ---- 1. 手作りブートストラップ: 変化量の「中央値の差」 ------------------
med_diff <- function(d) {
  median(d$change12[d$arm == "Drug"]) - median(d$change12[d$arm == "Placebo"])
}
med_diff(cc)

set.seed(1)
B <- 2000
boot_dist <- replicate(B, {
  i <- sample(nrow(cc), replace = TRUE)   # 同じ人が2回選ばれてもよい（復元抽出）
  med_diff(cc[i, ])
})
sd(boot_dist)                              # ブートストラップ標準誤差
quantile(boot_dist, c(0.025, 0.975))       # パーセンタイル 95% CI

ggplot(tibble(x = boot_dist), aes(x)) +
  geom_histogram(bins = 40, fill = "#0072B2", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = med_diff(cc), linetype = "dashed") +
  labs(title = "中央値の差のブートストラップ分布", x = "Drug − Placebo (mmHg)", y = NULL)
# 中央値はデータの値をとるのでガタガタした分布になる（統計量による）

# ---- 2. boot パッケージ ---------------------------------------------
# 関数は (データ, 添字) を受け取る形で書く
f_med <- function(d, i) med_diff(d[i, ])
set.seed(2)
b_med <- boot(cc, f_med, R = 2000, strata = cc$arm)   # strata: 群ごとの人数を固定して再抽出
b_med
boot.ci(b_med, type = c("perc", "bca"))
# BCa: 偏りと歪みを補正。一般にパーセンタイル法より正確（R を多めに: 2000以上）

# 回帰係数のブートストラップ（残差の正規性などに頼らない CI）
f_ancova <- function(d, i) coef(lm(sbp_w12 ~ arm + sbp_w0, data = d[i, ]))["armDrug"]
set.seed(3)
b_ancova <- boot(cc, f_ancova, R = 2000, strata = cc$arm)
boot.ci(b_ancova, type = c("perc", "bca"))
confint(lm(sbp_w12 ~ arm + sbp_w0, data = cc))["armDrug", ]   # 解析的 CI とほぼ一致するはず

# ---- 3. 再現性指標の CI ---------------------------------------------
w0 <- rep_long |> filter(week == 0) |> select(id, rep, sbp) |>
  pivot_wider(names_from = rep, values_from = sbp, names_prefix = "m")

f_repro <- function(d, i) {
  d <- d[i, ]
  msb <- 2 * var((d$m1 + d$m2) / 2)
  msw <- sum((d$m1 - d$m2)^2 / 2) / nrow(d)
  c(icc = (msb - msw) / (msb + msw),
    sw  = sqrt(msw),                     # within-subject SD
    loa_upper = mean(d$m1 - d$m2) + 1.96 * sd(d$m1 - d$m2))
}
set.seed(4)
b_repro <- boot(w0, f_repro, R = 2000)
boot.ci(b_repro, index = 1, type = "bca")   # ICC
boot.ci(b_repro, index = 2, type = "bca")   # within-subject SD
boot.ci(b_repro, index = 3, type = "bca")   # 一致限界の上限

# ---- 4. クラスターブートストラップ ------------------------------------
# 経時データの各行は独立でない。「行」を再抽出すると同じ患者の情報を独立とみなしてしまう。
# → 「患者」を丸ごと再抽出する
by_id <- split(visit_long, visit_long$id)
ids <- names(by_id)
# 2つの統計量を見る:
#   slope   : 週あたりの低下の群間差（患者「内」の変化を比べる）
#   level   : 全来院を通した平均血圧の群間差（患者「間」の比較）
two_stats <- function(d) {
  c(slope = coef(lm(sbp ~ week * arm, data = d))[["week:armDrug"]],
    level = coef(lm(sbp ~ arm, data = d))[["armDrug"]])
}

set.seed(5)
B <- 1000
boot_row     <- replicate(B, two_stats(visit_long[sample(nrow(visit_long), replace = TRUE), ]))
boot_cluster <- replicate(B, two_stats(bind_rows(by_id[sample(ids, replace = TRUE)])))

tibble(統計量 = rownames(boot_row),
       `SE 行を再抽出（誤り）`   = apply(boot_row, 1, sd),
       `SE 患者を再抽出`         = apply(boot_cluster, 1, sd))
# level: 行の再抽出は SE を過小評価（同じ人の4回を「4人分の情報」と数えてしまう）→ 偽陽性
# slope: 逆に行の再抽出は SE を過大評価（患者間の差がノイズ扱いになる）
#        患者単位の再抽出は 06章の混合モデルの SE (≈ 0.11) とほぼ一致する
# どちらにしても「データの独立な単位」で再抽出しないと答えが狂う
# 施設も考えるなら「施設 → 患者」の2段階で再抽出する方法もある（施設数が少ないと不安定）

# ---- 5. パラメトリックブートストラップ（混合モデル） ---------------------
# 推定したモデルから新しいデータを生成 → 再推定 を繰り返す
m <- lmer(sbp ~ week * arm + (1 | id), data = visit_long)
sd_fun <- function(fit) {
  c(sd_id = attr(VarCorr(fit)$id, "stddev")[[1]], sd_resid = sigma(fit),
    trt = fixef(fit)[["week:armDrug"]])
}
set.seed(6)
b_mer <- bootMer(m, sd_fun, nsim = 300, type = "parametric")   # 本番は 1000 以上
apply(b_mer$t, 2, quantile, probs = c(0.025, 0.975))
# confint(m, method = "boot", nsim = 500) でも同じことができる

# ---- 6. ベイジアンブートストラップ (Rubin, 1981) ------------------------
# 再抽出の代わりに、各人に Dirichlet 分布の重みをつける
#   → 「事前情報をほぼ入れないベイズ推定」の事後分布とみなせる
set.seed(7)
bayes_boot <- replicate(4000, {
  w <- rexp(nrow(cc)); w <- w / sum(w)      # Dirichlet(1,...,1) の重み
  drug <- cc$arm == "Drug"
  weighted.mean(cc$change12[drug], w[drug]) - weighted.mean(cc$change12[!drug], w[!drug])
})
quantile(bayes_boot, c(0.025, 0.5, 0.975))
mean(bayes_boot < -5)                        # P(差 < -5 mmHg) のような事後確率も出せる

# ---- 7. 注意点 ------------------------------------------------------
# - n が小さいと手元のデータが母集団を代表しない → ブートストラップも当てにならない
# - 最大値・最小値など「端」の統計量は苦手
# - データの構造（クラスター, 層, 時系列）を再抽出で再現すること
# - 乱数の種 set.seed() を固定して再現可能にする。R は結果が安定するまで増やす
#   (SE なら 200〜1000, CI なら 2000 以上が目安)
