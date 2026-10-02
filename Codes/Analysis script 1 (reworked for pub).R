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

# output directories
out_dir <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/MS1_Season1_Output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# data wrangling starts here ----

TRT_LEVELS <- c("NICT","NIMTM","CPCT","CPMTM","AICT","AIMTM")
TRT_SHORT  <- c(NICT="NI-CT", NIMTM="NI-MTM", CPCT="CP-CT",
                CPMTM="CP-MTM", AICT="AI-CT",  AIMTM="AI-MTM")
TRT_FULL   <- c(
  NICT  = "Natural Intercrop + Conventional Tillage",
  NIMTM = "Natural Intercrop + Minimum Tillage",
  CPCT  = "Cowpea Intercrop + Conventional Tillage",
  CPMTM = "Cowpea Intercrop + Minimum Tillage",
  AICT  = "Agroforestry Intercrop + Conventional Tillage",
  AIMTM = "Agroforestry Intercrop + Minimum Tillage"
)
TREE_LEVELS <- c("High","Low")
WEEK_LEVELS <- c(3, 6, 9)   # WAS — ordered earliest first

# selected color pallete for visualization 
trt_col <- c(
  NICT  = "#009E73",  # green
  NIMTM = "#E69F00",  # amber
  CPCT  = "#56B4E9",  # sky blue
  CPMTM = "#CC79A7",  # pink
  AICT  = "#D55E00",  # vermillion
  AIMTM = "#0072B2"   # blue
)
tree_col <- c(High = "#1B7837", Low = "#A6DBA0")

# tree cover label functions 
tc_lab <- labeller(
  Tree_cover = c(High = "High Tree Cover", Low = "Low Tree Cover")
)


# load data ----
sn <- excel_sheets(path1)

rd <- function(path, sheets, yr) {
  map(sheets, ~ read_excel(path, sheet = .x, na = "NA")) %>%
    set_names(sheets) %>%
    map(~ mutate(.x, Year = yr))
}

dfs1 <- rd(path1, sn, "Y1")
dfs2 <- rd(path2, sn, "Y2")   # loaded but only Y1 used below

FAW_raw  <- bind_rows(dfs1$FAW_DAMAGE, dfs2$FAW_DAMAGE)
HARV_raw <- bind_rows(dfs1$HARVEST_DATA, dfs2$HARVEST_DATA)
EL_raw   <- bind_rows(dfs1$EGG_LARVAE_BORDER_CORE_PLOTS,
                      dfs2$EGG_LARVAE_BORDER_CORE_PLOTS)