# soil-soc-repro
课程作业：Molecular trade-offs in soil organic carbon composition at continental scale论文数据可复现性研究项目
---
title: "Molecular trade-offs in soil organic carbon composition at continental scale论文数据复现报告"
author: [姜宇洋, 岳诗文, 吕智彤]
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```
1. 项目概述
本报告用于复现一篇土壤生物与生物地球化学相关论文中的部分图表。复现内容包括：
图 1：土壤有机碳组分箱线图
图 2：土壤有机碳分子的旋转主成分分析（PCA）双图
扩展数据图 ED2：土壤与植物生物地球化学特征箱线图
扩展数据图 ED4：不同植被类型与火干扰下的碳组分差异箱线图
扩展数据图 ED5：土壤有机碳分子组分相关性热图
复现目标：
完全复现论文原图样式
保证代码可运行、可重复
解决路径、绘图、窗口显示等常见问题
形成标准化、可交付的科研项目
2. 项目环境
本项目使用 R + RStudio + Quarto 完成。
```{r}
sessionInfo()
```
3. 文件结构
项目采用标准化结构，确保路径稳定、可移植：
```plaintext
project/
├── code/                # 绘图脚本
├── data/raw/            # 原始数据
├── output/figures/      # 复现结果图片
├── .Rproj               # RStudio 项目文件
├── README.txt           # 环境说明
└── report.qmd           # 本报告
```
4. 复现内容与结果
本项目成功复现以下图表：
4.1 图 1：土壤有机碳组分箱线图
使用 ggplot2 绘制土壤有机碳官能团与分子类别的箱线图，展示不同碳组分的相对丰度分布，图表样式、标签与论文原图完全一致。
4.2 图 2：PCA 排序图（双图）
使用 psych::biplot() 绘制三组合成图，展示土壤碳组分与主成分关系。
4.3 扩展数据图 ED2：环境因子分布箱线图
使用 ggplot2 + facet_wrap() 绘制 20 个环境变量的分布。
4.4 扩展数据图 ED4：碳组分与环境因子关系图
成功复现分组箱线图，标签、布局、样式与论文一致。
4.5 扩展数据图 ED5：有机碳分子组分相关性热图
成功复现碳组分相关性热图，包含相关系数与显著性标记，图表样式与论文完全一致。
5. 遇到的问题及解决方案
5.1 路径错误
原因：原始代码使用绝对路径解决：使用 here::here() 实现跨平台路径
5.2 Mac / Windows 不兼容
原因：quartz() 是 Mac 专用绘图函数解决：删除或替换为 windows()
5.3 灰色空窗口
原因：ggplot 未被 print()解决：print(ggplot(...))
5.4 基础绘图无法保存
原因：biplot 与 ggsave 不兼容解决：使用 png() 设备保存
6. 项目分工
姜宇洋：负责整体流程设计、代码调试、图表复现、问题解决与报告撰写，承担主要工作内容。岳诗文：参与数据整理与部分图表代码初步运行，协助完成基础核对工作及报告撰写。吕智彤：参与图表复现、协助项目文件整理、结果核对与环境备份工作。
7. 复现结论
✅ 所有图表 100% 复现论文原图✅ 代码无报错、可重复运行✅ 项目结构规范、可移植、可交付✅ 所有环境、路径、绘图问题已解决
本项目已达到科研论文可复现性标准。
