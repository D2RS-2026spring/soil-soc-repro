# 加载环境配置
source(here("code/setup.R"))

# ----------------------
# 1. 图1
# ----------------------
cat("正在生成图1...\n")
source(here("data/raw/NEON_NMR_boxplots_fig1_and_figED4.R"))
# ----------------------
# 2. 图2
# ----------------------
cat("正在生成图2...\n")
source(here("data/raw/NEON_NMR_ordination_fig2.R"))
# ----------------------
# 3. 图3
# ----------------------
cat("正在生成图3...\n")
source(here("data/raw/NEON_NMR_boxplot_figED2.R"))
# ----------------------
# 4. 图4
# ----------------------
cat("正在生成图4...\n")
source(here("data/raw/NEON_NMR_heatmap_figED5.R"))
