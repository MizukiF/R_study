# =============================================================
# 08_visualization.R  ビジュアライズ (ggplot2)
# =============================================================
# ggplot2 の文法: データ + aes(対応づけ) + geom(図形) + scale/facet/theme
#
#   1. 分布を見る: ヒストグラム, 密度
#   2. 群比較: 棒グラフ+エラーバー（ダイナマイトプロット）はやめて、点を見せる
#   3. 経時データ: スパゲッティプロット + 平均と95%CI
#   4. フォレストプロット（サブグループ解析）
#   5. Kaplan–Meier 曲線
#   6. 事後分布の図
#   7. 複数の図を並べて論文用に保存
#
# 色: Okabe–Ito パレット（色覚の多様性に配慮した定番）
#   日本語が文字化けする場合（主に Mac）:
#   theme_set(theme_minimal(base_family = "HiraKakuProN-W3"))

library(tidyverse)
library(survival)
library(broom)
library(patchwork)

if (!file.exists("data/patients_clean.rds")) source("02_wrangling.R")
pts <- read_rds("data/patients_clean.rds")
visit_long <- read_rds("data/visit_long.rds")
cc <- pts |> filter(!is.na(sbp_w12))

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))
arm_cols <- c(Placebo = "#E69F00", Drug = "#0072B2")
dir.create("output", showWarnings = FALSE)

# ---- 1. 分布 --------------------------------------------------------
p_hist <- ggplot(pts, aes(sbp_w0, fill = arm)) +
  geom_histogram(binwidth = 4, colour = "white", linewidth = 0.2, alpha = 0.85) +
  facet_wrap(~ arm, ncol = 1) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  labs(title = "ベースライン SBP の分布", x = "SBP (mmHg)", y = "人数")
p_hist

ggplot(pts, aes(HbA1c, colour = dm)) +
  geom_density(linewidth = 0.9) +
  scale_colour_manual(values = c(No = "#999999", Yes = "#D55E00"), name = "糖尿病") +
  labs(title = "HbA1c の分布", x = "HbA1c (%)", y = "密度")

# ---- 2. 群比較: ダイナマイトプロット vs 点を見せる図 ---------------------
summ <- cc |> group_by(arm) |>
  summarise(m = mean(change12), se = sd(change12) / sqrt(n()))

p_bad <- ggplot(summ, aes(arm, m, fill = arm)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  labs(title = "× 棒 + SE", subtitle = "分布も n も外れ値も見えない", x = NULL, y = "変化量 (mmHg)")

p_good <- ggplot(cc, aes(arm, change12, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_violin(fill = NA, colour = "grey80") +
  geom_jitter(width = 0.1, alpha = 0.45, size = 1.6) +
  stat_summary(fun.data = mean_cl_normal, colour = "black", linewidth = 0.8, size = 0.5) +
  scale_colour_manual(values = arm_cols, guide = "none") +
  labs(title = "○ 個々の点 + 平均と95%CI", subtitle = "分布とばらつきがそのまま見える",
       x = NULL, y = "変化量 (mmHg)")
p_bad + p_good

# ---- 3. 経時データ ---------------------------------------------------
traj <- visit_long |>
  drop_na(sbp) |>
  group_by(arm, week) |>
  summarise(n = n(), m = mean(sbp), se = sd(sbp) / sqrt(n),
            lo = m - qt(.975, n - 1) * se, hi = m + qt(.975, n - 1) * se, .groups = "drop")

p_traj <- ggplot() +
  geom_line(data = drop_na(visit_long, sbp), aes(week, sbp, group = id, colour = arm),
            alpha = 0.12) +
  geom_ribbon(data = traj, aes(week, ymin = lo, ymax = hi, fill = arm), alpha = 0.25) +
  geom_line(data = traj, aes(week, m, colour = arm), linewidth = 1.2) +
  geom_point(data = traj, aes(week, m, colour = arm), size = 2.5) +
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 4, 8, 12)) +
  coord_cartesian(ylim = c(110, 190)) +
  labs(title = "SBP の推移", subtitle = "細線 = 個々の患者、太線と帯 = 群平均と95%CI",
       x = "週", y = "SBP (mmHg)") +
  theme(legend.position = "top")
p_traj

# 患者を小分けにして見る（facet）: 個人差・脱落のパターンが見える
visit_long |>
  filter(id %in% sprintf("P%03d", 1:16)) |>
  ggplot(aes(week, sbp, colour = arm)) +
  geom_line() + geom_point() +
  facet_wrap(~ id, ncol = 4) +
  scale_colour_manual(values = arm_cols) +
  scale_x_continuous(breaks = c(0, 4, 8, 12)) +
  labs(x = "週", y = "SBP")

# ---- 4. フォレストプロット -------------------------------------------
sub_vars <- c(性別 = "sex", 糖尿病 = "dm", 年齢 = "age_group", BMI = "bmi25")
cc_sub <- cc |> mutate(bmi25 = if_else(bmi >= 25, "25以上", "25未満"))

forest <- imap_dfr(sub_vars, function(v, label) {
  d <- cc_sub |> mutate(sub = .data[[v]]) |> filter(!is.na(sub))   # 変数を "sub" 列にコピー
  p_int <- anova(lm(sbp_w12 ~ arm + sbp_w0 + sub, data = d),
                 lm(sbp_w12 ~ arm * sub + sbp_w0, data = d))$`Pr(>F)`[2]
  d |>
    group_by(level = as.character(sub)) |>
    group_modify(~ tidy(lm(sbp_w12 ~ arm + sbp_w0, data = .x), conf.int = TRUE) |>
                   filter(term == "armDrug") |> mutate(n = nrow(.x))) |>
    ungroup() |>
    mutate(group = label, p_int = p_int)
})
overall <- tidy(lm(sbp_w12 ~ arm + sbp_w0, data = cc), conf.int = TRUE) |>
  filter(term == "armDrug") |> mutate(group = "全体", level = "全体", n = nrow(cc), p_int = NA)

forest_df <- bind_rows(overall, forest) |>
  mutate(label = if_else(group == "全体", "全体", paste0(group, ": ", level)),
         label = sprintf("%s (n=%d)", label, n),
         label = fct_rev(fct_inorder(label)))

p_forest <- ggplot(forest_df, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_vline(xintercept = overall$estimate, linetype = "dashed", colour = "#0072B2") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, orientation = "y") +
  geom_point(aes(size = n), shape = 15, colour = "#0072B2") +
  geom_text(aes(x = 6, label = if_else(is.na(p_int), "", sprintf("p交互作用=%.2f", p_int))),
            hjust = 0, size = 3, colour = "grey40",
            data = forest_df |> group_by(group) |> slice(1) |> ungroup()) +
  scale_size_continuous(range = c(2, 4), guide = "none") +
  coord_cartesian(xlim = c(-22, 12)) +
  labs(title = "サブグループ別の治療効果 (ANCOVA)",
       subtitle = "左ほど薬剤が有利。サブグループ間の差は交互作用の p値で判断する",
       x = "群間差 Drug − Placebo (mmHg) と 95%CI", y = NULL)
p_forest

# ---- 5. Kaplan–Meier 曲線 ---------------------------------------------
km <- survfit(Surv(time_wk, event) ~ arm, data = pts)
km_df <- tidy(km) |>
  mutate(arm = str_remove(strata, "arm=")) |>
  group_by(arm) |>
  group_modify(~ bind_rows(tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1), .x)) |>
  mutate(next_time = lead(time, default = last(time))) |>   # 群の中で次の時点（群をまたがない）
  ungroup() |>
  mutate(arm = factor(arm, levels = names(arm_cols)))

p_km <- ggplot(km_df, aes(time, estimate, colour = arm, fill = arm)) +
  geom_step(linewidth = 1) +
  geom_rect(aes(xmin = time, xmax = next_time, ymin = conf.low, ymax = conf.high),
            alpha = 0.12, colour = NA) +
  geom_point(data = filter(km_df, n.censor > 0), shape = 3, size = 1.5) +  # 打ち切り
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "イベント非発生割合 (Kaplan–Meier)", subtitle = "+ = 打ち切り",
       x = "週", y = "イベント非発生") +
  theme(legend.position = "top")
p_km
# Number at risk 表も付けるなら survminer::ggsurvplot() が便利

# ---- 6. 事後分布の図（05章のモデルがあれば） -------------------------------
if (file.exists("output/models/b05_ancova.rds")) {
  fit_b <- read_rds("output/models/b05_ancova.rds")
  draws <- brms::as_draws_df(fit_b)
  p_post <- ggplot(tibble(x = draws$b_armDrug), aes(x)) +
    annotate("rect", xmin = -3, xmax = 3, ymin = -Inf, ymax = Inf, fill = "grey85") +
    geom_density(fill = "#56B4E9", alpha = 0.6, colour = "#0072B2") +
    geom_vline(xintercept = -5, linetype = "dashed") +
    annotate("text", x = 0, y = 0, label = "ROPE", vjust = -0.5, colour = "grey40") +
    labs(title = "治療効果の事後分布",
         subtitle = sprintf("P(効果 < -5 mmHg) = %.3f", mean(draws$b_armDrug < -5)),
         x = "Drug − Placebo (mmHg)", y = NULL)
  print(p_post)
}

# ---- 7. 並べて保存 ----------------------------------------------------
fig <- (p_traj | p_km) / p_forest +
  plot_annotation(tag_levels = "A", title = "Figure 1")
fig
ggsave("output/08_figure1.png", fig, width = 12, height = 10, dpi = 300, bg = "white")
# 論文: 300 dpi 以上 / ベクタ形式なら ggsave("x.pdf", device = cairo_pdf)
# 図のルール: 軸ラベルに単位 / 0基準の要否を考える / 二重軸は使わない / 色だけに頼らない
