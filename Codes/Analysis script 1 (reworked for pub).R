# Introduction ----
# Title   : "Effect of Soil Tillage and Cropping Systems on Fall Armyworm Population
#  and Damage Across Low and High Tree Density Farming Landscapes in Nigeria"
# Author  : Akinbuluma et al. 
# Project : Ekhaga Agroecology, Oyo State, Nigeria
# Season  : Y1 (Season 1) ONLY — Season 2 excluded (near-zero damage, single
#           damage class across 98% of observations — insufficient variation
#           for meaningful statistical inference)
# Design  : 10 farms × 6 treatments × 2 tree-cover levels (RCBD)


# Packages -----
pkgs <- c(
  "readxl","purrr","tidyverse",
  "lme4","lmerTest",
  "emmeans","multcomp","multcompView",
  "ordinal",      # clmm — ordinal mixed model
  "glmmTMB",      # negative binomial GLMM
  "DHARMa",       # residual diagnostics
  "patchwork",    # multi-panel figures
  "scales",       # axis helpers
  "ggtext",       # markdown in plot labels
  "ggdist",       # raincloud / halfeye
  "ggpubr",       # stat_cor
  "colorspace",   # colour utilities
  "flextable",    # Word-ready tables
  "officer",      # write .docx
  "writexl"       # write .xlsx
)
new_p <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_p)) { message("Installing: ", paste(new_p, collapse=", "))
  install.packages(new_p) }
invisible(lapply(pkgs, library, character.only = TRUE))

# Path setting  -----
# data paths
path1   <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_First season data.xlsx"
path2   <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_Second season data.xlsx"

# outoput directories
out_dir <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/MS1_Season1_Output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
