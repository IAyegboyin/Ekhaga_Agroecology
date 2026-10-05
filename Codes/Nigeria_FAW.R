# Author of the full script for preliminary analysis was Rhett Harrison 
# Few tweaks were made to make it work for me including some libraries and data paths,
# Others were just check and discussion with other contributors of the project 
# Including Dr Akinbuluma, M.D (University of Ibadan); Dr Akinyemi (University of Osun)

# library(gitcreds)
library(readxl)
library(purrr)
library(tidyverse)
library(performance)
library(see)
library(multcomp)
library(multcompView)

# gitcreds_set()

sheet_names <- excel_sheets(
  "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_First season data.xlsx"
)
dfs1 <- map(
  sheet_names,
  ~ read_excel(
    "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_First season data.xlsx",
    sheet = .x,
    na = "NA"
  )
)

dfs2 <- map(
  sheet_names,
  ~ read_excel(
    "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_Second season data.xlsx",
    sheet = .x,
    na = "NA"
  )
)

dfs1 <- map(sheet_names, ~read_excel("/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_First season data.xlsx", sheet = .x, na = "NA"))
names(dfs1) <- sheet_names

#file_path <- "/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data"
#sheet_names <- excel_sheets(file_path)

dfs2 <- map(sheet_names, ~read_excel("/Users/Esmael/Desktop/Data Science Library/Data for play/Ekhaga_Agroecology/Data/Ekhaga_Second season data.xlsx", sheet = .x, na = "NA"))
names(dfs2) <- sheet_names

dfs1 <- map(dfs1, ~ .x %>% mutate(Year = "Y1"))
dfs2 <- map(dfs2, ~ .x %>% mutate(Year = "Y2"))

FAW_damage <- list(dfs1$FAW_DAMAGE, dfs2$FAW_DAMAGE)
FAW_damage <- bind_rows(FAW_damage[[1]],FAW_damage[[2]])

Harvest_data <- list(dfs1$HARVEST_DATA, dfs2$HARVEST_DATA)
Harvest_data <- bind_rows(Harvest_data[[1]],Harvest_data[[2]])

Pheromone_data <- list(dfs1$PHEROMONES_DATA, dfs2$PHEROMONES_DATA)
Pheromone_data <- bind_rows(Pheromone_data[[1]],Pheromone_data[[2]])

Arthropods_data <- list(dfs1$SOIL_ARTHRPODS, dfs2$SOIL_ARTHROPODS)
Arthropods_data <- bind_rows(Arthropods_data[[1]],Arthropods_data[[2]])

Egg_Larvae_data <- list(dfs1$EGG_LARVAE_BORDER_CORE_PLOTS, dfs2$EGG_LARVAE_BORDER_CORE_PLOTS)
Egg_Larvae_data <- bind_rows(Egg_Larvae_data[[1]],Egg_Larvae_data[[2]])

Germination_data <- list(dfs1$CROP_GERMINATION_SUCCESS, dfs2$CROP_GERMINATION_SUCCESS)
Germination_data <- bind_rows(Germination_data[[1]],Germination_data[[2]])


with(FAW_damage, table(Farm, Name_of_Farmers)) # Farm number consistent with name
with(FAW_damage, table(Tree_cover))
with(FAW_damage, table(Treatment))
with(FAW_damage, table(Weeks))
with(FAW_damage, table(Year))

with(FAW_damage, table(is.na(Tree_cover)))
with(FAW_damage, table(is.na(Treatment)))
with(FAW_damage, table(is.na(Weeks)))
with(FAW_damage, table(is.na(Year)))
# Dataset is balanced with no missing data

library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ggplot2)
library(ggforce)
library(ggeffects)
library(glmmTMB)

with(FAW_damage, hist(FAW_Damage))
with(FAW_damage, table(FAW_Damage,Year)) # only 40 observations of damage in Yr2
str(FAW_damage)

FAW_damage <- FAW_damage %>%
        mutate(
                Farm = as.factor(Farm),
                Tree_cover = as.factor(Tree_cover),
                Treatment = factor(Treatment, levels = c('NICT','NIMTM','CPCT','CPMTM','AICT','AIMTM')),
                W_points = as.factor(W_points),
                Name_of_Farmers = as.factor(Name_of_Farmers),
                FAW_Damage = FAW_Damage - 1
        )

Year1 <- FAW_damage %>% filter(Year == "Y1")
with(Year1, hist(sqrt(FAW_Damage)))
with(Year1, table(Treatment))

# Damage model for Year1 only

M1 <-  lmer(sqrt(Value_pct) ~
                    Tree_cover +
                    Weeks +
                    Treatment +
                    #Tree_cover:Treatment +    # comment out because of high collinearity
                    (1 | Name_of_Farmers),
            data = dmg_yr1)
summary(M1)
M1.emm <- emmeans(M1, ~ Treatment | Tree_cover, pbkrtest.limit = 8000)
contrast(M1.emm, 'eff')
contrast(M1.emm, "pairwise")
check_model(M1)

M1 <- lmer(
  sqrt(Value_pct) ~
    Treatment +
    Tree_cover +
    Weeks +
    (1 + Weeks | Name_of_Farmers),
  data = dmg_yr1
)

M1 <- lmer(
  sqrt(Value_pct) ~
    Treatment +
    Tree_cover +
    Weeks +
    (1 | Name_of_Farmers) +
    (0 + Weeks | Name_of_Farmers),
  data = dmg_yr1
)

M1.emm <- emmeans(M1, ~ Treatment)
cld(M1.emm, adjust = "tukey", Letters = letters)

AIC(M1)

resids <- resid(M1, type='pearson')
plot(resids~fitted(M1))
lines(lowess(resids~fitted(M1)), col='red') # OK problem of integer values
plot(sqrt(abs(resids))~ fitted(M1))
lines(lowess(sqrt(abs(resids))~
                     fitted(M1)), col='red') # obvious curve at lower values but not bad



Tree_cover_labels = c("High" = "HIGH TREE COVER", "Low" = "LOW TREE COVER")
sig_difference = c(0,-1,-1,-1,0,0,
                   0,0,0,0,0,0)

Cover_Trt <- paste(Year1$Tree_cover,Year1$Treatment, sep = "_")


Year1 <- Year1 %>%
        mutate(
                SIGDIFF = factor(Cover_Trt, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Year1$SIGDIFF <- as.numeric(as.character(Year1$SIGDIFF))

Year1 %>%
        ggplot(aes(x = Treatment, y = FAW_Damage, fill = sign(SIGDIFF))) +
        geom_violin(kernel = "r") +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        stat_summary(fun=mean, geom="point", size=3) +
        geom_hline(
                data = . %>%
                        group_by(Tree_cover) %>%
                        summarise(line = mean(FAW_Damage, na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylim(0,6) +
        xlab("") +
        ylab("Fall armyworm damage") +
        facet_col("Tree_cover", labeller = labeller(Tree_cover = Tree_cover_labels))


### Model for FAW damage Year2
Year2 <- FAW_damage %>% filter(Year == "Y2")
with(Year2, hist(sqrt(FAW_Damage)))
with(Year2, table(Treatment))

# Damage model for Year1 only

M2 <-  glmer(sign(FAW_Damage) ~
                    #Tree_cover +
                    #as.factor(Weeks) +
                    Treatment +
                    #Tree_cover:Treatment +
                    (1 | Name_of_Farmers),
             data = Year2,
            family = binomial(link = "logit")
)

summary(M2) # None of the variables are significant

#### Look at Pheromone data

Pheromone_data
Pheromone_data <- Pheromone_data %>%
        mutate(
                Farm = factor(Farm),
                Tree_cover = factor(Tree_cover)
        )

with(Pheromone_data, hist(Male_FAW))
with(Pheromone_data, hist(log(Male_FAW+1)))

# Problem of zero inflated distribution
# Try Binomial model first

M3 <- glmer(sign(Male_FAW) ~ 
                    Weeks +
                    Tree_cover +
                    (1 | Farmer_Name),
            data = Pheromone_data,
            family = binomial(link = "logit")
                    )
summary(M3) # Best performing model - slightly overdispersed but acceptible


view(Pheromone_data)
M4 <- glmmTMB(Male_FAW ~
                       #Tree_cover +
                       Weeks +
                       (1 | Farmer_Name),
              data = Pheromone_data,
              family = nbinom2(link = "log"),
              ziformula = ~ .
                       )

summary(M4) # Zero inflation model is not significant and neg.binomial model is overdispersed
# AIC and BIC also much higher

M5 <- glmer(Male_FAW ~
                      Tree_cover +
                      Weeks +
                      (1 | Farmer_Name),
              data = Pheromone_data,
              family = negative.binomial(theta = 1)
)

summary(M5) # Performs even worse than M4

# Look at harvest data
print(Harvest_data, width = Inf) # Need to deal with NAs
with(Harvest_data, table(Treatments))

Harvest_data <- Harvest_data %>%
        mutate(
                Farm = factor(Farm),
                Tree_Cover = factor(Tree_Cover),
                Year = factor(Year),
                Name = factor(Name_of_Farmers),
                Maize_Cob_damage = Maize_Cob_damage - 1,
                Treatments = factor(Treatments, levels = c('NICT','NIMTM','CPCT','CPMTM','AICT','AIMTM')),
        )

# Model of cob damage

M6 <- glmer(sign(Maize_Cob_damage) ~
                    Year +
                    Tree_Cover +
                    Treatments +
                    #Treatments:Year +
                    #Treatments:Tree_Cover +
                    (1 | Name),
            data = Harvest_data,
            family = binomial(link = "logit")
)

summary(M6) # Not over dispersed
M6.emm <- emmeans(M6, ~ Treatments, pbkrtest.limit = 8000)
contrast(M6.emm, 'eff')
contrast(M6.emm, "pairwise")

sig_difference = c(0,-1,-1,-1,0,0)

Harvest_data <- Harvest_data %>%
        mutate(
                SIGDIFF = factor(Treatments, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Harvest_data$SIGDIFF <- as.numeric(as.character(Harvest_data$SIGDIFF))

view(Harvest_data)
Harvest_data %>%
        ggplot(aes(x = Treatments, y = Maize_Cob_damage, fill = sign(SIGDIFF))) +
        geom_violin(kernel = "r") +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        stat_summary(fun=mean, geom="point", size=3) +
        geom_hline(
                data = . %>%
                        summarise(line = mean(Maize_Cob_damage, na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylim(0,4) +
        xlab("") +
        ylab("Cob damage")


print(Harvest_data, width = Inf)

Harvest_data <- Harvest_data %>%
        group_by(Year,Tree_Cover,Treatments,Name) %>%
        summarise(
              Maize_harvest_weight = mean(Maize_harvest_weight, na.rm = TRUE),
              Cowpea_weight = mean(Cowpea_weight, na.rm = TRUE),
              Vegetable_plant_weight = mean(Vegetable_plant_weight, na.rm = TRUE)
        ) %>% ungroup()

print(Harvest_data, n = Inf)

with(Harvest_data, hist(Maize_harvest_weight))

M7 <- lmer(Maize_harvest_weight ~
                   Year +
                   Tree_Cover +
                   Treatments +
                 #Treatments:Year,
                 #Treatments:Tree_Cover +
           (1 | Name),
         data = Harvest_data
)

summary(M7)

resids <- resid(M7, type='pearson')
plot(resids~fitted(M7))
lines(lowess(resids~fitted(M7)), col='red') # Good
plot(sqrt(abs(resids))~ fitted(M7))
lines(lowess(sqrt(abs(resids))~
                     fitted(M7)), col='red') # obvious curve but not bad


sig_difference = c(0,0,1,0,0,0)

Harvest_data <- Harvest_data %>%
        mutate(
                SIGDIFF = factor(Treatments, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Harvest_data$SIGDIFF <- as.numeric(as.character(Harvest_data$SIGDIFF))


Harvest_data %>%
        ggplot(aes(x = Treatments, y = Maize_harvest_weight, fill = sign(SIGDIFF))) +
        geom_violin() +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        stat_summary(fun=mean, geom="point", size=3) +
        geom_hline(
                data = . %>%
                        summarise(line = mean(Maize_harvest_weight, na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylim(0,4) +
        xlab("") +
        ylab("Maize yield")

## Do models for Cowpea and veg
with(Harvest_data, hist(Cowpea_weight))

M8 <- lmer(Cowpea_weight ~
                 Year +
                 #Tree_Cover +
                 #Treatments +
         #Treatments:Year +
         #Treatments:Tree_Cover +
                 (1 | Name),
         data = Harvest_data
)

summary(M8)

resids <- resid(M8, type='pearson')
plot(resids~fitted(M8))
lines(lowess(resids~fitted(M8)), col='red') # Some trend but not worth worrying about

Harvest_data %>%
        ggplot(aes(x = Year, y = Cowpea_weight)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=3) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Cowpea yield")


M9 <- lmer(Vegetable_plant_weight ~
                 Year +
         Tree_Cover +
         Treatments +
         Treatments:Year +
         Treatments:Tree_Cover +
                 (1 | Name),
         data = Harvest_data
)

summary(M9)

resids <- resid(M9, type='pearson')
plot(resids~fitted(M9))
lines(lowess(resids~fitted(M9)), col='red') # Some kink but not worth worrying about

Harvest_data <- Harvest_data %>%
        filter(Treatments == "AICT" | Treatments == "AIMTM") %>%
        droplevels()


Harvest_data %>%
        ggplot(aes(x = Treatments, y = Vegetable_plant_weight)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=3) +
        
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Vegetable weight") +
        facet_grid(rows = vars(Year), cols = vars(Tree_Cover), scales = 'free') 


## Egg / larvae data
print(Egg_Larvae_data, width = Inf)
Egg_Larvae_data <- Egg_Larvae_data %>%
        mutate(
                Farm = factor(Farm),
                Tree_cover = factor(`Tree cover`),
                Year = factor(Year),
                Name = factor(`Name of farmer`),
                Treatment = factor(Treatment,
                                   levels = c('NICT', 'NIMTM', 'CPCT', 'CPMTM', 'AICT', 'AIMTM')),
                Larvae = `Larvae No in border row`,
                Eggs = `Egg masses in border row`,
                Egg_parasitoid = ifelse(`Egg parasitoid` == "YES", 1, 0),
                Larvae_parasitoid = ifelse(`Larvae Parasitoid` == "YES", 1, 0)
        )

Egg_Larvae_data

with(Egg_Larvae_data, table(Treatment, Larvae))
with(Egg_Larvae_data, table(Treatment, Eggs))
with(Egg_Larvae_data, hist(Larvae))
with(Egg_Larvae_data, hist(Eggs))
with(Egg_Larvae_data, hist(Larvae_parasitoid))
with(Egg_Larvae_data, hist(Egg_parasitoid))

with(Egg_Larvae_data, table(Tree_cover, Name)) # Spelling mistakes with names

Egg_Larvae_data <- Egg_Larvae_data %>%
        mutate(
                Name = fct_collapse(Name, "Elijah Ominyi" = c("Elijah Ominyi", "Elijaj Ominyi")),
                Name = fct_collapse(Name, "Mr Rahman" = c("Mr Rahman", "Mr_Rahman")),
                Name = fct_collapse(Name, "Mr and Mrs Onaolapo" = c("Mr and Mrs Onaolapo", "Mr_and_Mrs_Onaolapo")),
                Name = fct_collapse(Name, "Abdullahi Abubakar" = c("Abdullahi Abubakar", "Mr_Abdullahi_Abubakar")),
                Name = fct_collapse(Name, "Godwin Tachanta" = c("Godwin Tachanta", "Mr_Godwin_Tachanta")),
                Name = fct_collapse(Name, "Kademla Tachanta" = c("Kademla Tachanta", "Mr_Kademla_Tachanta")),
                Name = fct_collapse(Name, "Bamidele Seyi" = c("Bamidele Seyi", "Mr_Bamidele_Seyi"))
        )

M10 <- glmer(sign(Eggs) ~ 
                    Year +
                    Week +
                    Tree_cover +
                    Treatment +
                    #Treatment:Week +
                    #Treatment:Tree_cover +
                    #Treatment:Year +
                    (1 | Name),
            data = Egg_Larvae_data,
            family = binomial(link = "logit")
)

summary(M10) # 

M11 <- glmmTMB(Eggs ~
                       Year +
                       Week +
                       Tree_cover +
                       Treatment +
                       #Treatment:Week +
                       #Treatment:Tree_cover +
                       (1 | Name),
               data = Egg_Larvae_data,
               family = nbinom2(link = "log"),
               ziformula = ~ .
)

summary(M11) #not giving any better result than the binomial model and convergence problems

sig_difference = c(0,-1,-1,-1,-1,-1)

Egg_Larvae_data <- Egg_Larvae_data %>%
        mutate(
                SIGDIFF = factor(Treatment, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Egg_Larvae_data$SIGDIFF <- as.numeric(as.character(Egg_Larvae_data$SIGDIFF))

Egg_Larvae_data %>%
        ggplot(aes(x = Treatment, y = sign(Eggs), fill = SIGDIFF)) +
        geom_violin(kernel = "r") +
        stat_summary(fun=mean, geom="point", size=3) +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        geom_hline(
                data = . %>%
                        summarise(line = mean(sign(Eggs), na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Proportion of plots with egg masses") 

M12 <- glmer(sign(Larvae) ~ 
                     #Year +
                     Week +
                     Tree_cover +
                     Treatment +
                     #Treatment:Week +
                     #Treatment:Tree_cover +
                     #Treatment:Year +
                     (1 | Name),
             data = Egg_Larvae_data,
             family = binomial(link = "logit")
)

summary(M12) # 

M13 <- glmmTMB(Larvae ~
                       Year +
                       Week +
                       Tree_cover +
                       Treatment +
                       #Treatment:Week +
                       #Treatment:Tree_cover +
                       #Treatment:Year +
                       (1 | Name),
               data = Egg_Larvae_data,
               family = nbinom2(link = "log"),
               ziformula = ~ .
)

summary(M13) #

sig_difference = c(0,-1,-1,-1,-1,-1)

Egg_Larvae_data <- Egg_Larvae_data %>%
        mutate(
                SIGDIFF = factor(Treatment, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Egg_Larvae_data$SIGDIFF <- as.numeric(as.character(Egg_Larvae_data$SIGDIFF))

Egg_Larvae_data %>%
        ggplot(aes(x = Treatment, y = Larvae, fill = SIGDIFF)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=3) +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        geom_hline(
                data = . %>%
                        summarise(line = mean(Larvae, na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Abundance of FAW Larvae") 


## Next arthropod data
print(Arthropods_data, width = Inf)

Arthropods_data <- Arthropods_data %>%
        mutate(
                Farm = factor(Farm),
                Tree_Cover = factor(Tree_Cover),
                Treatments = factor(Treatments, 
                                    levels = c('NICT','NIMTM','CPCT','CPMTM','AICT','AIMTM')),
                Famer_Name = factor(Farmer_Name),
                Year = factor(Year)
        )

apply(Arthropods_data[,5:18], 2, sum, na.rm=TRUE)
# Only Hymenoptera (mostly ants), Orthoptera, and Areanae are abundant
# Focus on Hymenoptera and Areanae as these are predators

with(Arthropods_data, hist(log(Hymenoptera+1)))
with(Arthropods_data, hist(log(Areanae_Other_Arthropods)))
with(Arthropods_data, table(sign(Hymenoptera),Treatments, Year))
# Almost no variance in Year2 and variance in Year1 is small


M14 <- glmer(sign(Hymenoptera) ~ 
                     Year +
                     #weeks +
                     #Tree_Cover +
                     Treatments +
                     #Treatments:weeks +
                     #Treatments:Tree_Cover +
                     #Treatments:Year +
                     (1 | Farmer_Name),
             data = Arthropods_data,
             family = binomial(link = "logit")
) # 

summary(M14)

Arthropods_data %>%
        group_by(Year,Treatments) %>%
        summarise(
                Hymenoptera = sum(Hymenoptera, na.rm = TRUE),
                Areanae = sum(Areanae_Other_Arthropods, na.rm = TRUE)
        ) %>% ungroup()

print(Arthropods_data, n = Inf) %>%
        group_by(Tree_Cover,Farmer_Name) %>%
        summarise(
                Hymenoptera = sum(Hymenoptera, na.rm = TRUE),
                Areanae = sum(Areanae_Other_Arthropods, na.rm = TRUE)
        ) %>% ungroup()

#High variance across Farms


M15 <- glmmTMB(log(Hymenoptera+1) ~
                       #Year +
                       #weeks +
                       #Tree_Cover +
                       Treatments +
                       #Treatments:weeks +
                       #Treatments:Tree_Cover +
                       #Treatments:Year +
                       (1 | Farmer_Name),
               data = Arthropods_data,
               family = 'gaussian',
               ziformula = ~ .
)

summary(M15) # Binomial model better

sig_difference = c(0,0,-1,0,-1,-1)

Arthropods_data <- Arthropods_data %>%
        mutate(
                SIGDIFF = factor(Treatments, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Arthropods_data$SIGDIFF <- as.numeric(as.character(Arthropods_data$SIGDIFF))

Arthropods_data %>%
        ggplot(aes(x = Treatments, y = sign(Hymenoptera), fill = SIGDIFF)) +
        geom_violin(kernel = "r") +
        stat_summary(fun=mean, geom="point", size=3) +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        geom_hline(
                data = . %>%
                        summarise(line = mean(sign(Hymenoptera), na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Proportion of plots with ants") 


M16 <- glmer(sign(Areanae_Other_Arthropods) ~
                     #Year +
                     #weeks +
                     #Tree_Cover +
                     #Treatments +
                     #Treatments:weeks +
                     #Treatments:Tree_Cover +
                     #Treatments:Year +
                     (1 | Farmer_Name),
             data = Arthropods_data,
             family = binomial(link = "logit")
) 

summary(M16) # Singular


M17 <- glmmTMB(log(Areanae_Other_Arthropods+1) ~
                       Year +
                       weeks +
                       Tree_Cover +
                       Treatments +
                       #Treatments:weeks +
                       #Treatments:Tree_Cover +
                       #Treatments:Year +
                       (1 | Farmer_Name),
               data = Arthropods_data,
               family = 'gaussian',
               ziformula = ~ .
)

summary(M17) # 

sig_difference = c(0,1,0,1,0,0)

Arthropods_data <- Arthropods_data %>%
        mutate(
                SIGDIFF = factor(Treatments, labels = sig_difference),
                SIGDIFF = factor(SIGDIFF, levels =  c(-1,0,1))
        )
Arthropods_data$SIGDIFF <- as.numeric(as.character(Arthropods_data$SIGDIFF))

Arthropods_data %>%
        ggplot(aes(x = Treatments, y = log(Areanae_Other_Arthropods+1), fill = SIGDIFF)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=3) +
        scale_fill_gradient2(low='red', mid='snow3', high='blue', space='Lab') +
        geom_hline(
                data = . %>%
                        summarise(line = mean(log(Areanae_Other_Arthropods+1), na.rm=TRUE)),
                mapping = aes(yintercept = line), linewidth = 0.5, linetype =5
        ) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        ylab("Log(Abundance of Spiders") 



