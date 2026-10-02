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
library(readxl)
library(tidyverse)
library(lmerTest)
library(lme4)
library(multcompView)
library(multcomp)
library(emmeans)
library(ordinal)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(scales)
library(ggtext)
library(ggdist)
library(ggpubr)
library(colorspace)
library(flaxtable)
library(officer)
library(writexl)

# Path setting  -----
# data paths
path1   <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_First season data.xlsx"
path2   <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_Second season data.xlsx"

# outoput directories
out_dir <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/MS1_Season1_Output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# data wrangling starts here ----
