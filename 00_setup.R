# =============================================================
# 00_setup.R  環境準備（最初に1回だけ実行）
# =============================================================
# 作業ディレクトリはこのフォルダ (R/) にしておくこと。
#   RStudio なら R.Rproj を開けば自動でそうなる。
#   確認: getwd()

pkgs <- c(
  "tidyverse",   # データ整形 + ggplot2
  "broom",       # モデル結果を data.frame に
  "patchwork",   # 複数の図を並べる
  "irr",         # ICC など一致度
  "lme4", "lmerTest", "performance", "emmeans", "pbkrtest",  # 混合モデル
  "BayesFactor", "bayestestR", "brms", "see", "logspline",  # ベイズ
  "boot",        # ブートストラップ
  "pwr",         # 検出力
  "survival", "pROC", "mice"                     # 生存時間, ROC, 多重代入
)
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)

# brms は裏で Stan を使う。Windows では Rtools が必要。
#   pkgbuild::has_build_tools(debug = TRUE)   # TRUE ならOK
# 速いバックエンド cmdstanr を使う場合（任意）:
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
#   cmdstanr::install_cmdstan()

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("output/models", showWarnings = FALSE)

cat("準備OK\n")
