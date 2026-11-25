# =============================================================
# 09_extras.R  その他の重要事項
# =============================================================
#   1. 生存時間解析: Kaplan–Meier, log-rank, Cox 回帰, 比例ハザードの確認
#   2. 診断・予測: ROC 曲線と AUC（+ ブートストラップ CI）
#   3. 欠測データ: 完全ケース解析 / LOCF / 混合モデル / 多重代入 を真値と比較
#   4. 解析の再現性: 乱数の種, パッケージのバージョン, 解析計画

library(tidyverse)
library(survival)
library(pROC)
library(lme4)
library(mice)

if (!file.exists("data/patients_clean.rds")) source("02_wrangling.R")
pts <- read_rds("data/patients_clean.rds")
visit_long <- read_rds("data/visit_long.rds")
truth <- read_rds("data/truth.rds")
theme_set(theme_minimal(base_size = 12))

# ---- 1. 生存時間解析 -------------------------------------------------
# 追跡期間が人によって違う（打ち切り）→ 「イベントあり/なし」の割合比較では情報を捨てている
km <- survfit(Surv(time_wk, event) ~ arm, data = pts)
km                                   # 中央生存時間
summary(km, times = c(12, 26, 52))   # 指定時点のイベント非発生割合
survdiff(Surv(time_wk, event) ~ arm, data = pts)   # log-rank 検定

cox <- coxph(Surv(time_wk, event) ~ arm + age + dm, data = pts)
summary(cox)                         # exp(coef) = ハザード比
cox.zph(cox)                         # 比例ハザード性の検定（p が小さいと仮定に疑い）
# ※ 05章: この HR の p値は 0.1 前後。ベイズで読むと P(HR<1) ≈ 0.93

# ---- 2. ROC 曲線 ----------------------------------------------------
# 「ベースラインの HbA1c で糖尿病を判別できるか」（簡単な例）
roc_a1c <- roc(pts$dm, pts$HbA1c, levels = c("No", "Yes"), direction = "<")
auc(roc_a1c)
ci.auc(roc_a1c, method = "bootstrap", boot.n = 1000)     # 07章のブートストラップ
coords(roc_a1c, "best", best.method = "youden",
       ret = c("threshold", "sensitivity", "specificity"))
ggroc(roc_a1c, colour = "#0072B2", linewidth = 1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey60") +
  labs(title = sprintf("HbA1c による糖尿病の判別 (AUC = %.2f)", auc(roc_a1c)),
       x = "特異度", y = "感度")

# 予測モデルを作る場合の注意:
#   - 同じデータで作って評価すると性能を過大評価（楽観バイアス）→ ブートストラップ/交差検証で補正
#   - AUC（識別能）だけでなく較正 (calibration) も見る
#   - 報告は TRIPOD に従う

# ---- 3. 欠測データ ---------------------------------------------------
# このデータの脱落は「前回の血圧が高い人ほど脱落」= MAR（観測データで説明できる欠測）
#   MCAR: 完全にランダム / MAR: 観測値に依存 / MNAR: 観測されない値そのものに依存
cat("真の効果 (12週):", 12 * truth$b_trt, "\n")
md.pattern(pts |> select(sbp_w0, sbp_w4, sbp_w8, sbp_w12), plot = FALSE)

# (a) 完全ケース解析
est_cc <- coef(lm(sbp_w12 ~ arm + sbp_w0, data = pts))["armDrug"]

# (b) LOCF（最終観測値の繰り越し）— 古くからあるが一般に推奨されない
locf <- pts |>
  mutate(sbp_w12_locf = coalesce(sbp_w12, sbp_w8, sbp_w4, sbp_w0))
est_locf <- coef(lm(sbp_w12_locf ~ arm + sbp_w0, data = locf))["armDrug"]

# (c) 混合モデル（脱落前のデータも全部使う。MAR なら妥当）
m <- lmer(sbp ~ week * arm + (week | id), data = drop_na(visit_long, sbp),
          control = lmerControl(optimizer = "bobyqa"))
est_mm <- 12 * fixef(m)["week:armDrug"]

# (d) 多重代入 (MICE): 欠測値を「もっともらしい値」で複数回埋めて、結果を統合する
imp_data <- pts |> select(arm, age, sex, dm, sbp_w0, sbp_w4, sbp_w8, sbp_w12)
imp <- mice(imp_data, m = 20, method = "pmm", seed = 1, printFlag = FALSE)
pooled <- pool(with(imp, lm(sbp_w12 ~ arm + sbp_w0)))
summary(pooled, conf.int = TRUE)
est_mi <- summary(pooled) |> filter(term == "armDrug") |> pull(estimate)

tibble(方法 = c("真値", "完全ケース", "LOCF", "混合モデル", "多重代入"),
       推定値 = c(12 * truth$b_trt, est_cc, est_locf, est_mm, est_mi))
# LOCF は「脱落時点で効果が止まる」と仮定するので、効果を過小評価しやすい
# 1回のデータでは偶然のズレも大きい → 01章の seed を変えて何度か比べてみよう（練習問題）
# MNAR が疑われるなら感度分析（delta 調整法など）を行う

# ---- 4. 解析の再現性 -------------------------------------------------
# - 乱数を使う処理の前に set.seed()
# - 生データは書き換えない。整形はすべてコード（02章）で行い、生データ → 結果 を再実行可能に
# - パッケージのバージョンを記録: sessionInfo()、プロジェクトごとに固定するなら renv
#     install.packages("renv"); renv::init(); renv::snapshot()
# - 解析計画（主要評価項目・解析方法・共変量・欠測の扱い）はデータを見る前に決める
#   報告ガイドライン: CONSORT（RCT）, STROBE（観察研究）, STARD（診断精度）, TRIPOD（予測モデル）
# - Quarto / R Markdown でコードと結果と文章を1つの文書に（レポートの再現性）
sessionInfo()
