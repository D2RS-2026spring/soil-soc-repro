# 安装并加载所有依赖包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  here,
  ggplot2,
  multcomp,
  gridExtra,
  cowplot,
  reshape,
  viridis,
  psych,
  Hmisc,
  GGally,
  vegan,
  pheatmap
)

# 锁定项目根目录
here::i_am("code/setup.R")

cat("✅ 所有依赖包已加载完成\n")
