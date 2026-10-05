# Statistical modelling ----
# 
# helper functions
clean_cld <- function(df)
  mutate(as.data.frame(df), .group = trimws(.group))

smry <- function(data, resp, grp_vars) {
  data %>%
    group_by(across(all_of(grp_vars))) %>%
    summarise(
      Mean = mean(.data[[resp]], na.rm = TRUE),
      SE   = sd(.data[[resp]],   na.rm = TRUE) / sqrt(n()),
      N    = n(),
      .groups = "drop"
    )
}
auto_nudge <- function(df, frac = 0.10, floor_v = 0.5)
  max(max(df$Mean + df$SE, na.rm = TRUE) * frac, floor_v)


# save plot functions 
save_fig <- function(p, nm, w = 170, h = 120, units = "mm", dpi_tiff = 600) {
  ggsave(file.path(out_dir, paste0(nm, ".png")), p,
         width = w, height = h, units = units, dpi = 600, bg = "white")
  cat("  ✓", nm, "\n")
}

## ── 8a. FAW Foliar Damage — Ordinal CLMM ─────────────────────────────────────
# Season 1 has full score spread (0–6) — ordinal model is appropriate
# Fixed: Treatment × Tree_cover + Weeks (temporal covariate)
# Random: (1 | Farm) — blocking factor

cat("\n── 8a. Ordinal CLMM: FAW Foliar Damage ──\n")
mod_dmg <- clmm(
  Damage_ord ~ Treatment * Tree_cover + Weeks + (1 | Farm),
  data = FAW,
  link = "logit"
)
print(summary(mod_dmg))

# Treatment × Tree_cover contrasts
emm_dmg_trt <- emmeans(mod_dmg,
                       ~ Treatment | Tree_cover,
                       mode = "linear.predictor")
cld_dmg_trt <- clean_cld(as.data.frame(
  cld(emm_dmg_trt, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — FAW Damage by Treatment | Tree Cover:\n")
print(cld_dmg_trt)

# Tree cover main effect
emm_dmg_tree <- emmeans(mod_dmg, ~ Tree_cover,
                        mode = "linear.predictor")
cld_dmg_tree <- clean_cld(as.data.frame(
  cld(emm_dmg_tree, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — FAW Damage by Tree Cover:\n")
print(cld_dmg_tree)

# Weeks main effect (ordered earliest first)
emm_dmg_wk <- emmeans(mod_dmg, ~ Weeks,
                      mode = "linear.predictor")
cld_dmg_wk <- clean_cld(as.data.frame(
  cld(emm_dmg_wk, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — FAW Damage by Sampling Week:\n")
print(cld_dmg_wk)

# Odds ratios for fixed effects
cat("\nOdds Ratios (exp[coef]) — FAW Damage:\n")
cf_d <- coef(summary(mod_dmg))
or_d <- data.frame(
  Term    = rownames(cf_d),
  OR      = round(exp(cf_d[,"Estimate"]), 3),
  CI_low  = round(exp(cf_d[,"Estimate"] - 1.96*cf_d[,"Std. Error"]), 3),
  CI_high = round(exp(cf_d[,"Estimate"] + 1.96*cf_d[,"Std. Error"]), 3),
  p       = round(cf_d[,"Pr(>|z|)"], 4)
)
print(or_d)


## ── 8b. FAW Larval Count (Core Plot) — NB GLMM ───────────────────────────────
cat("\n── 8b. NB GLMM: FAW Larval Count (Core Plot) ──\n")

EL_clean <- droplevels(EL)

mod_larv <- glmmTMB(
  Larvae_core ~ Treatment * Tree_cover + Week + (1 | Farm),
  data   = EL_clean,
  family = nbinom2()
)
print(summary(mod_larv))

# DHARMa diagnostics
sim_larv <- simulateResiduals(mod_larv)
png(file.path(out_dir, "DHARMa_Larvae.png"), width = 900, height = 500)
plot(sim_larv, main = "DHARMa: Larval Count (NB GLMM)")
dev.off()
cat("  DHARMa plot saved.\n")

emm_larv_trt <- emmeans(mod_larv,
                        ~ Treatment | Tree_cover,
                        type = "response")
cld_larv_trt <- clean_cld(as.data.frame(
  cld(emm_larv_trt, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Larval Count by Treatment | Tree Cover:\n")
print(cld_larv_trt)

emm_larv_tree <- emmeans(mod_larv, ~ Tree_cover, type = "response")
cld_larv_tree <- clean_cld(as.data.frame(
  cld(emm_larv_tree, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Larval Count by Tree Cover:\n")
print(cld_larv_tree)

emm_larv_wk <- emmeans(mod_larv, ~ Week, type = "response")
cld_larv_wk <- clean_cld(as.data.frame(
  cld(emm_larv_wk, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Larval Count by Week:\n")
print(cld_larv_wk)

# Incidence Rate Ratios
cat("\nIRR (exp[coef]) — Larval Count:\n")
cf_l <- as.data.frame(coef(summary(mod_larv))$cond)
irr_l <- data.frame(
  Term    = rownames(cf_l),
  IRR     = round(exp(cf_l$Estimate), 3),
  CI_low  = round(exp(cf_l$Estimate - 1.96*cf_l$`Std. Error`), 3),
  CI_high = round(exp(cf_l$Estimate + 1.96*cf_l$`Std. Error`), 3),
  p       = round(cf_l$`Pr(>|z|)`, 4)
)
print(irr_l)


## ── 8c. Egg Masses (Border Row) — NB GLMM ────────────────────────────────────
cat("\n── 8c. NB GLMM: Egg Masses (Border Row) ──\n")

mod_eggs <- glmmTMB(
  Eggs_border ~ Treatment * Tree_cover + Week + (1 | Farm),
  data   = EL_clean,
  family = nbinom2()
)
print(summary(mod_eggs))

sim_eggs <- simulateResiduals(mod_eggs)
png(file.path(out_dir, "DHARMa_Eggs.png"), width = 900, height = 500)
plot(sim_eggs, main = "DHARMa: Egg Masses (NB GLMM)")
dev.off()

emm_eggs_trt <- emmeans(mod_eggs,
                        ~ Treatment | Tree_cover,
                        type = "response")
cld_eggs_trt <- clean_cld(as.data.frame(
  cld(emm_eggs_trt, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Egg Masses by Treatment | Tree Cover:\n")
print(cld_eggs_trt)

emm_eggs_tree <- emmeans(mod_eggs, ~ Tree_cover, type = "response")
cld_eggs_tree <- clean_cld(as.data.frame(
  cld(emm_eggs_tree, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Egg Masses by Tree Cover:\n")
print(cld_eggs_tree)

emm_eggs_wk <- emmeans(mod_eggs, ~ Week, type = "response")
cld_eggs_wk <- clean_cld(as.data.frame(
  cld(emm_eggs_wk, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Egg Masses by Week:\n")
print(cld_eggs_wk)

cat("\nIRR (exp[coef]) — Egg Masses:\n")
cf_e <- as.data.frame(coef(summary(mod_eggs))$cond)
irr_e <- data.frame(
  Term    = rownames(cf_e),
  IRR     = round(exp(cf_e$Estimate), 3),
  CI_low  = round(exp(cf_e$Estimate - 1.96*cf_e$`Std. Error`), 3),
  CI_high = round(exp(cf_e$Estimate + 1.96*cf_e$`Std. Error`), 3),
  p       = round(cf_e$`Pr(>|z|)`, 4)
)
print(irr_e)


## ── 8d. Maize Yield — Linear Mixed Model ─────────────────────────────────────
cat("\n── 8d. LMM: Maize Harvest Weight ──\n")

mod_yld <- lmer(
  Maize_harvest_weight ~ Treatment * Tree_cover + (1 | Farm),
  data = Yield,
  REML = TRUE
)
print(summary(mod_yld))
cat("\nType III ANOVA (Satterthwaite df):\n")
print(anova(mod_yld, type = "III"))

# Variance components + ICC
vc     <- as.data.frame(VarCorr(mod_yld))
v_farm <- vc$vcov[vc$grp == "Farm"]
v_res  <- vc$vcov[vc$grp == "Residual"]
icc    <- v_farm / (v_farm + v_res)
cat(sprintf("\nRandom effects: Farm variance = %.4f, Residual = %.4f\n",
            v_farm, v_res))
cat(sprintf("ICC = %.4f  (%.1f%% of variance attributable to Farm blocking)\n",
            icc, icc * 100))

emm_yld_trt <- emmeans(mod_yld, ~ Treatment | Tree_cover)
cld_yld_trt <- clean_cld(as.data.frame(
  cld(emm_yld_trt, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Yield by Treatment | Tree Cover:\n")
print(cld_yld_trt)

emm_yld_tree <- emmeans(mod_yld, ~ Tree_cover)
cld_yld_tree <- clean_cld(as.data.frame(
  cld(emm_yld_tree, Letters = letters, adjust = "tukey", reversed = TRUE)
))
cat("\nCLD — Yield by Tree Cover:\n")
print(cld_yld_tree)

cat("\nFixed Effect Coefficients ±95% CI:\n")
fe <- as.data.frame(coef(summary(mod_yld)))
fe$LCL <- fe$Estimate - 1.96 * fe$`Std. Error`
fe$UCL <- fe$Estimate + 1.96 * fe$`Std. Error`
print(round(fe, 4))


# ── 9. RAW MEAN SUMMARIES joined to CLD letters ───────────────────────────────
raw_dmg <- smry(FAW,   "Damage_pct",          c("Treatment","Tree_cover")) %>%
  left_join(cld_dmg_trt %>%
              dplyr::select(Treatment, Tree_cover, .group),
            by = c("Treatment","Tree_cover"))

raw_dmg_wk <- smry(FAW, "Damage_pct", c("Treatment","Tree_cover","Weeks"))

raw_lar <- smry(EL, "Larvae_core", c("Treatment","Tree_cover")) %>%
  left_join(cld_larv_trt %>%
              dplyr::select(Treatment, Tree_cover, .group),
            by = c("Treatment","Tree_cover"))

raw_lar_wk <- smry(EL, "Larvae_core", c("Treatment","Tree_cover","Week"))

raw_egg <- smry(EL, "Eggs_border", c("Treatment","Tree_cover")) %>%
  left_join(cld_eggs_trt %>%
              dplyr::select(Treatment, Tree_cover, .group),
            by = c("Treatment","Tree_cover"))

raw_yld <- smry(Yield, "Maize_harvest_weight", c("Treatment","Tree_cover")) %>%
  left_join(cld_yld_trt %>%
              dplyr::select(Treatment, Tree_cover, .group),
            by = c("Treatment","Tree_cover"))

# NA diagnostic
cat("\n── CLD join NA check ──\n")
for (nm in c("raw_dmg","raw_lar","raw_egg","raw_yld")) {
  df  <- get(nm)
  nas <- sum(is.na(df$.group))
  cat(sprintf("  %-12s  NA .group: %d / %d  %s\n",
              nm, nas, nrow(df),
              if (nas > 0) "⚠ check join" else "✓"))
}


# ── 10. EXPORT CLD TABLES ─────────────────────────────────────────────────────
cat("\n── CLD table exports ──\n")
write_csv(cld_dmg_trt,  file.path(out_dir, "CLD_Damage_Treatment.csv"))
write_csv(cld_dmg_tree, file.path(out_dir, "CLD_Damage_TreeCover.csv"))
write_csv(cld_dmg_wk,   file.path(out_dir, "CLD_Damage_Week.csv"))
write_csv(cld_larv_trt, file.path(out_dir, "CLD_Larvae_Treatment.csv"))
write_csv(cld_larv_tree,file.path(out_dir, "CLD_Larvae_TreeCover.csv"))
write_csv(cld_larv_wk,  file.path(out_dir, "CLD_Larvae_Week.csv"))
write_csv(cld_eggs_trt, file.path(out_dir, "CLD_Eggs_Treatment.csv"))
write_csv(cld_eggs_tree,file.path(out_dir, "CLD_Eggs_TreeCover.csv"))
write_csv(cld_eggs_wk,  file.path(out_dir, "CLD_Eggs_Week.csv"))
write_csv(cld_yld_trt,  file.path(out_dir, "CLD_Yield_Treatment.csv"))
write_csv(cld_yld_tree, file.path(out_dir, "CLD_Yield_TreeCover.csv"))
cat("  ✓ 11 CLD CSV files exported\n")


# ── 11. PUBLICATION FIGURES ───────────────────────────────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  SECTION 11: FIGURES\n")
cat(strrep("═", 60), "\n")


## ── Figure 1. FAW Damage Over Time (line + SE ribbon) ────────────────────────
f1 <- ggplot(raw_dmg_wk,
             aes(x = Weeks, y = Mean,
                 colour = Treatment, group = Treatment)) +
  geom_ribbon(aes(ymin = pmax(Mean - SE, 0),
                  ymax = Mean + SE,
                  fill = Treatment),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.3) +
  scale_colour_manual(values = trt_col, labels = TRT_SHORT,
                      name = "Treatment") +
  scale_fill_manual(  values = trt_col, guide  = "none") +
  scale_y_continuous(labels = label_percent(scale = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 1.  FAW Foliar Damage Over Time",
    subtitle = "Mean (±SE) damage severity (Davis scale → % leaf area);<br>
weeks ordered earliest (3 WAS) to latest (9 WAS)",
    x = "Weeks After Sowing (WAS)",
    y = "Mean FAW Damage (%)"
  ) +
  theme_jrnl() +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.2,"cm")) +
  guides(colour = guide_legend(nrow = 2,
                               override.aes = list(linewidth = 1.4)))

save_fig(f1, "Fig1_Damage_Over_Time", w = 170, h = 115)


## ── Figure 2. FAW Damage CLD Bar Chart ───────────────────────────────────────
nd2 <- auto_nudge(raw_dmg, frac = 0.10)

f2 <- ggplot(raw_dmg,
             aes(x = Treatment, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE),
                width = 0.22, linewidth = 0.65, colour = "grey25") +
  geom_text(aes(y = Mean + SE + nd2, label = .group),
            size = 3.5, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_discrete(labels = TRT_SHORT) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.16))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 2.  FAW Foliar Damage by Treatment",
    subtitle = "Mean ±SE; bars with the same letter do not differ significantly<br>
(Tukey-adjusted pairwise comparisons, α = 0.05; ordinal CLMM)",
    x = "Treatment",
    y = "Mean FAW Damage (%)"
  ) +
  theme_jrnl() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1))

save_fig(f2, "Fig2_Damage_CLD", w = 170, h = 115)


## ── Figure 3. FAW Larval Count CLD Bar Chart ─────────────────────────────────
nd3 <- auto_nudge(raw_lar)

f3 <- ggplot(raw_lar,
             aes(x = Treatment, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE),
                width = 0.22, linewidth = 0.65, colour = "grey25") +
  geom_text(aes(y = Mean + SE + nd3, label = .group),
            size = 3.5, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_discrete(labels = TRT_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 3.  FAW Larval Density in Core Plots",
    subtitle = "Mean (±SE) larval count; NB GLMM; Tukey-adjusted CLD (α = 0.05)",
    x = "Treatment",
    y = "Mean Larval Count per Core Plot"
  ) +
  theme_jrnl() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1))

save_fig(f3, "Fig3_Larvae_CLD", w = 170, h = 115)


## ── Figure 4. Larval Heatmap (Treatment × Week) ───────────────────────────────
lheat <- raw_lar_wk %>%
  mutate(Treatment = fct_rev(Treatment))

f4 <- ggplot(lheat, aes(x = Week, y = Treatment, fill = Mean)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(Mean >= 0.5, round(Mean, 1), "")),
            size = 2.8, colour = "white", fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#FFFFD4","#FED98E","#FE9929","#D94701","#7F2704"),
    name    = "Mean\nlarvae",
    limits  = c(0, NA)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, " WAS")) +
  scale_y_discrete(labels = rev(TRT_SHORT)) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 4.  Temporal Dynamics of FAW Larval Density",
    subtitle = "Mean larval count per core plot; sampling weeks ordered earliest→latest;<br>
darker cells = greater infestation pressure",
    x = "Sampling Week (WAS)",
    y = "Treatment"
  ) +
  theme_jrnl(base_size = 10) +
  theme(panel.grid  = element_blank(),
        axis.line   = element_blank(),
        axis.ticks  = element_blank())

save_fig(f4, "Fig4_Larval_Heatmap", w = 170, h = 105)


## ── Figure 5. Egg Masses CLD Bar Chart ───────────────────────────────────────
nd5 <- auto_nudge(raw_egg, frac = 0.10, floor_v = 0.05)

f5 <- ggplot(raw_egg,
             aes(x = Treatment, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE),
                width = 0.22, linewidth = 0.65, colour = "grey25") +
  geom_text(aes(y = Mean + SE + nd5, label = .group),
            size = 3.5, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_discrete(labels = TRT_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 5.  FAW Egg Masses in Border Rows",
    subtitle = "Mean (±SE) egg mass count; NB GLMM; Tukey-adjusted CLD (α = 0.05)",
    x = "Treatment",
    y = "Mean Egg Masses per Border-Row Transect"
  ) +
  theme_jrnl() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1))

save_fig(f5, "Fig5_EggMasses_CLD", w = 170, h = 115)


## ── Figure 6. Border vs Core Larval Movement ─────────────────────────────────
el_mv <- EL %>%
  group_by(Treatment, Tree_cover, Week) %>%
  summarise(
    Border = mean(Larvae_border, na.rm = TRUE),
    Core   = mean(Larvae_core,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(Border:Core,
               names_to  = "Position",
               values_to = "Mean") %>%
  mutate(Week_num = as.numeric(as.character(Week)))

f6 <- ggplot(el_mv,
             aes(x = Week_num, y = Mean,
                 colour = Position, linetype = Position,
                 group  = Position)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_colour_manual(values  = c(Border = "#D55E00", Core = "#0072B2"),
                      name    = "Plot position") +
  scale_linetype_manual(values = c(Border = "dashed", Core = "solid"),
                        name   = "Plot position") +
  scale_x_continuous(
    breaks = sort(unique(as.numeric(as.character(EL$Week)))),
    labels = function(x) paste0(x, " WAS")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_grid(Treatment ~ Tree_cover,
             labeller = labeller(
               Tree_cover = c(High="High Cover", Low="Low Cover"),
               Treatment  = TRT_SHORT
             )) +
  labs(
    title    = "Figure 6.  FAW Larval Dispersal: Border Row vs Core Plot",
    subtitle = "Mean larval count across sampling weeks (earliest→latest);<br>
dashed = border row (entry point), solid = core plot (infestation pressure)",
    x = "Sampling Week (WAS)",
    y = "Mean Larval Count"
  ) +
  theme_jrnl(base_size = 9) +
  theme(legend.position = "top")

save_fig(f6, "Fig6_Border_vs_Core", w = 170, h = 200)


## ── Figure 7. Maize Yield CLD Bar Chart ──────────────────────────────────────
nd7 <- auto_nudge(raw_yld, frac = 0.07, floor_v = 0.05)

f7 <- ggplot(raw_yld,
             aes(x = Treatment, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE),
                width = 0.22, linewidth = 0.65, colour = "grey25") +
  geom_text(aes(y = Mean + SE + nd7, label = .group),
            size = 3.5, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_discrete(labels = TRT_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 7.  Maize Harvest Weight by Treatment",
    subtitle = "Mean (±SE) kg; linear mixed model; Tukey-adjusted CLD (α = 0.05)",
    x = "Treatment",
    y = "Maize Harvest Weight (kg)"
  ) +
  theme_jrnl() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1))

save_fig(f7, "Fig7_Yield_CLD", w = 170, h = 115)


## ── Figure 8. Damage–Yield Scatter + Marginal Densities ──────────────────────
dm_f <- FAW %>%
  group_by(Farm, Tree_cover, Treatment) %>%
  summarise(Damage = mean(Damage_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(Farm = as.character(Farm))

yl_f <- Yield %>%
  group_by(Farm, Tree_cover, Treatment) %>%
  summarise(Yield_kg = mean(Maize_harvest_weight, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Farm = as.character(Farm))

dy <- left_join(dm_f, yl_f,
                by = c("Farm","Tree_cover","Treatment")) %>%
  drop_na(Damage, Yield_kg)

p_sc <- ggplot(dy,
               aes(x = Damage, y = Yield_kg,
                   colour = Treatment, fill = Treatment,
                   shape  = Tree_cover)) +
  stat_ellipse(aes(group = Treatment), type = "norm", level = 0.68,
               geom = "polygon", alpha = 0.08, colour = NA) +
  geom_point(size = 3.0, alpha = 0.88, stroke = 0.5) +
  geom_smooth(
    data    = dy,
    mapping = aes(x = Damage, y = Yield_kg,
                  colour = Tree_cover, fill = Tree_cover),
    method = "lm", se = TRUE, linewidth = 0.85,
    alpha  = 0.15, inherit.aes = FALSE
  ) +
  stat_cor(
    data    = dy,
    mapping = aes(x = Damage, y = Yield_kg, colour = Tree_cover),
    inherit.aes = FALSE,
    size = 2.8, label.x.npc = "left", label.y.npc = "top"
  ) +
  scale_colour_manual(
    values = c(trt_col, tree_col),
    breaks = TRT_LEVELS,
    labels = TRT_SHORT,
    name   = "Treatment"
  ) +
  scale_fill_manual(values = c(trt_col, tree_col), guide = "none") +
  scale_shape_manual(values = c(High = 21, Low = 24),
                     name   = "Tree cover") +
  labs(x = "Mean FAW Damage (% leaf area)",
       y = "Maize Harvest Weight (kg)") +
  theme_jrnl(base_size = 10) +
  theme(legend.position = "right")

p_top <- ggplot(dy, aes(x = Damage, fill = Treatment)) +
  geom_density(alpha = 0.45, colour = NA) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_jrnl(base_size = 8) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        plot.margin = margin(4,4,0,10))

p_rt <- ggplot(dy, aes(x = Yield_kg, fill = Treatment)) +
  geom_density(alpha = 0.45, colour = NA) +
  scale_fill_manual(values = trt_col, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_jrnl(base_size = 8) +
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(4,4,4,0))

f8 <- (p_top + plot_spacer() + p_sc + p_rt) +
  plot_layout(ncol = 2, nrow = 2, widths = c(4,1), heights = c(1,4)) +
  plot_annotation(
    title    = "Figure 8.  FAW Damage–Yield Trade-off",
    subtitle = "Each point = farm-level mean; shaded ellipses = ±1 SD treatment cluster\nTrend lines fitted per tree-cover level (±95% CI); Pearson r shown per panel",
    theme = theme(
      plot.title    = element_text(face="bold", size=11, colour="grey10"),
      plot.subtitle = element_text(colour="grey35", size=9),
      plot.background = element_rect(fill="white", colour=NA)
    )
  )

save_fig(f8, "Fig8_Damage_Yield_Scatter", w = 170, h = 155)


## ── Figure 9. Tillage × Cropping System Interaction ──────────────────────────
ti_smry <- FAW %>%
  group_by(Tillage, Cropping, Tree_cover) %>%
  summarise(
    Mean = mean(Damage_pct, na.rm = TRUE),
    SE   = sd(Damage_pct,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

f9 <- ggplot(ti_smry,
             aes(x = Cropping, y = Mean,
                 colour = Tillage, group = Tillage)) +
  geom_line(linewidth = 1.2, position = position_dodge(0.12)) +
  geom_point(size = 3.5, position = position_dodge(0.12)) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0.14, linewidth = 0.75,
                position = position_dodge(0.12)) +
  scale_colour_manual(
    values = c("Conv-Till" = "#555555", "Min-Till" = "#D55E00"),
    name   = "Tillage system"
  ) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0.02, 0.12))) +
  facet_wrap(~ Tree_cover, labeller = tc_lab) +
  labs(
    title    = "Figure 9.  Tillage × Cropping System Interaction on FAW Damage",
    subtitle = "Mean ±SE; crossing lines indicate a treatment interaction between factors",
    x = "Cropping System",
    y = "Mean FAW Damage (%)"
  ) +
  theme_jrnl()

save_fig(f9, "Fig9_Tillage_Cropping_Interaction", w = 170, h = 110)


## ── Figure 10. Raincloud: Full Damage Distribution ───────────────────────────
f10 <- ggplot(FAW,
              aes(x      = Treatment,
                  y      = Damage_pct,
                  fill   = Tree_cover,
                  colour = Tree_cover)) +
  ggdist::stat_halfeye(
    adjust = 0.7, width = 0.55, .width = 0,
    point_colour = NA, alpha = 0.50,
    position = position_dodge(width = 0.75)
  ) +
  geom_boxplot(
    width = 0.18, outlier.shape = NA,
    colour = "grey20", linewidth = 0.5,
    position = position_dodge(width = 0.75)
  ) +
  stat_summary(
    fun = mean, geom = "point",
    shape = 23, size = 2.8,
    fill = "white", colour = "grey20",
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(  values = tree_col, name = "Tree Cover") +
  scale_colour_manual(values = tree_col, name = "Tree Cover") +
  scale_x_discrete(labels = TRT_SHORT) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     limits = c(-2, NA)) +
  facet_wrap(~ Weeks,
             labeller = labeller(
               Weeks = function(x) paste0(x, " WAS")
             )) +
  labs(
    title    = "Figure 10.  FAW Damage Distribution by Treatment, Tree Cover and Week",
    subtitle = "Raincloud: density cloud + IQR box + mean (◇);<br>
panels ordered by sampling week (earliest left, latest right)",
    x = "Treatment",
    y = "FAW Damage (% leaf area)"
  ) +
  theme_jrnl() +
  theme(legend.position = "top",
        axis.text.x     = element_text(angle = 35, hjust = 1))

save_fig(f10, "Fig10_Raincloud_Damage", w = 200, h = 120)


# ── 12. PUBLICATION TABLES ────────────────────────────────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  SECTION 12: PUBLICATION TABLES\n")
cat(strrep("═", 60), "\n")

## Table 1. Trial Design -------------------------------------------------------
tbl1 <- data.frame(
  Feature = c(
    "Study location",
    "Agro-ecology",
    "Cropping season reported",
    "Number of farms (replicates)",
    "Tree cover levels",
    "Treatments (n = 6)",
    "   NI-CT",
    "   NI-MTM",
    "   CP-CT",
    "   CP-MTM",
    "   AI-CT",
    "   AI-MTM",
    "Sampling points per plot",
    "FAW assessment method",
    "Damage scale",
    "Observations (FAW damage)",
    "Observations (larvae/eggs)",
    "Observations (yield)"
  ),
  Details = c(
    "Oyo State, southwest Nigeria",
    "Guinea savannah; smallholder subsistence maize farming",
    "Season 1 (Y1)",
    "10 (5 high tree cover · 5 low tree cover)",
    "High (≥30% canopy cover) · Low (<10% canopy cover)",
    "",
    "Natural intercrop + Conventional tillage",
    "Natural intercrop + Minimum tillage",
    "Cowpea intercrop + Conventional tillage",
    "Cowpea intercrop + Minimum tillage",
    "Agroforestry intercrop + Conventional tillage",
    "Agroforestry intercrop + Minimum tillage",
    "6 per plot (W-scouting pattern at 3, 6, 9 WAS)",
    "Visual scouting; Davis scale (1–7) per 10-plant cluster",
    "Converted to % leaf area: 1=0%, 2=10%, 3=20%, 4=40%, 5=60%, 6=90%, 7=100%",
    paste0(nrow(FAW)),
    paste0(nrow(EL)),
    paste0(nrow(Yield))
  ),
  stringsAsFactors = FALSE
)

ft1 <- flextable(tbl1) %>%
  set_header_labels(Feature = "Feature / Variable", Details = "Description") %>%
  bold(j = 1) %>%
  bg(i = 6, bg = "#E8E8E8") %>%
  fontsize(size = 10) %>%
  font(fontname = "Times New Roman") %>%
  width(j = 1, width = 2.4) %>%
  width(j = 2, width = 4.2) %>%
  set_caption("Table 1. Summary of trial design, site characteristics, sampling scheme and treatment descriptions (Season 1).") %>%
  theme_vanilla()

## Table 2. FAW Foliar Damage --------------------------------------------------
tbl2 <- raw_dmg %>%
  mutate(
    Treatment  = TRT_SHORT[as.character(Treatment)],
    Tree_cover = as.character(Tree_cover),
    Mean       = round(Mean, 1),
    SE         = round(SE,   1),
    .group     = ifelse(is.na(.group), "—", .group)
  ) %>%
  rename(`Treatment` = Treatment, `Tree Cover` = Tree_cover,
         `Mean Damage (%)` = Mean, `±SE` = SE, `n` = N, `CLD` = .group)

ft2 <- flextable(tbl2) %>%
  bold(part = "header") %>%
  fontsize(size = 10) %>%
  font(fontname = "Times New Roman") %>%
  bg(i = which(tbl2$`Tree Cover` == "High"), bg = "#F5F5F5") %>%
  autofit() %>%
  set_caption("Table 2. Mean (±SE) FAW foliar damage (% leaf area) per treatment and tree-cover level. Letters denote Tukey-adjusted pairwise group membership (CLD; ordinal CLMM, α = 0.05).")

## Table 3. Larval Density + Egg Masses ----------------------------------------
tbl3 <- bind_rows(
  raw_lar %>% mutate(Variable = "Larval count (core plot)"),
  raw_egg %>% mutate(Variable = "Egg masses (border row)")
) %>%
  mutate(
    Treatment  = TRT_SHORT[as.character(Treatment)],
    Tree_cover = as.character(Tree_cover),
    Mean       = round(Mean, 2),
    SE         = round(SE,   2),
    .group     = ifelse(is.na(.group), "—", .group)
  ) %>%
  dplyr::select(Variable, Treatment, `Tree Cover` = Tree_cover,
                Mean, `±SE` = SE, n = N, CLD = .group)

ft3 <- flextable(tbl3) %>%
  bold(part = "header") %>%
  fontsize(size = 10) %>%
  font(fontname = "Times New Roman") %>%
  merge_v(j = "Variable") %>%
  autofit() %>%
  set_caption("Table 3. Mean (±SE) FAW larval density (core plots) and egg mass deposition (border rows) per treatment and tree-cover level. Tukey-adjusted CLD from negative binomial GLMM (α = 0.05).")

## Table 4. Maize Yield --------------------------------------------------------
tbl4 <- raw_yld %>%
  mutate(
    Treatment  = TRT_SHORT[as.character(Treatment)],
    Tree_cover = as.character(Tree_cover),
    Mean       = round(Mean, 2),
    SE         = round(SE,   2),
    .group     = ifelse(is.na(.group), "—", .group)
  ) %>%
  rename(`Treatment` = Treatment, `Tree Cover` = Tree_cover,
         `Mean Yield (kg)` = Mean, `±SE` = SE, `n` = N, `CLD` = .group)

ft4 <- flextable(tbl4) %>%
  bold(part = "header") %>%
  fontsize(size = 10) %>%
  font(fontname = "Times New Roman") %>%
  bg(i = which(tbl4$`Tree Cover` == "High"), bg = "#F5F5F5") %>%
  autofit() %>%
  set_caption("Table 4. Mean (±SE) maize harvest weight (kg) per treatment and tree-cover level. Tukey-adjusted CLD from linear mixed model (α = 0.05).")

## Table 5. Model Reference ----------------------------------------------------
tbl5 <- data.frame(
  `Response variable`  = c("FAW foliar damage","FAW larval count (core)",
                           "Egg masses (border)","Maize yield"),
  `Distribution`       = c("Ordinal (logistic)","Negative binomial",
                           "Negative binomial","Gaussian"),
  `Model`              = c("clmm","glmmTMB","glmmTMB","lmer"),
  `Fixed effects`      = c(
    "Treatment × Tree_cover + Weeks",
    "Treatment × Tree_cover + Week",
    "Treatment × Tree_cover + Week",
    "Treatment × Tree_cover"),
  `Random effects`     = rep("(1 | Farm)", 4),
  `R package`          = c("ordinal","glmmTMB","glmmTMB","lme4 + lmerTest"),
  `Post-hoc`           = rep("emmeans; Tukey CLD", 4),
  check.names = FALSE
)

ft5 <- flextable(tbl5) %>%
  bold(part = "header") %>%
  fontsize(size = 9.5) %>%
  font(fontname = "Times New Roman") %>%
  autofit() %>%
  set_caption("Table 5. Summary of statistical models fitted for each response variable (Season 1). All models include Farm as a random intercept (blocking factor).")

## Save Word document ----------------------------------------------------------
doc <- read_docx() %>%
  body_add_par("MANUSCRIPT 1 — PUBLICATION TABLES (Season 1)", style = "heading 1") %>%
  body_add_par("Ekhaga Agroecology — Oyo State Nigeria", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft1) %>% body_add_break() %>%
  body_add_flextable(ft2) %>% body_add_break() %>%
  body_add_flextable(ft3) %>% body_add_break() %>%
  body_add_flextable(ft4) %>% body_add_break() %>%
  body_add_flextable(ft5)
print(doc, target = file.path(out_dir, "MS1_Tables.docx"))
cat("  ✓ MS1_Tables.docx\n")

## Excel -----------------------------------------------------------------------
writexl::write_xlsx(
  list(
    Table1_Design   = tbl1,
    Table2_Damage   = tbl2,
    Table3_LarvEggs = tbl3,
    Table4_Yield    = tbl4,
    Table5_Models   = tbl5
  ),
  path = file.path(out_dir, "MS1_Tables.xlsx")
)
cat("  ✓ MS1_Tables.xlsx\n")


# ── 13. RESULTS TEXT: KEY NUMBERS ─────────────────────────────────────────────
sink(file.path(out_dir, "Results_Numbers.txt"), split = TRUE)

cat("═══════════════════════════════════════════════════════\n")
cat("  RESULTS SECTION — KEY NUMBERS  (Season 1 only)\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("── FAW FOLIAR DAMAGE ──────────────────────────────────\n")
cat("Grand mean damage (%):\n")
print(FAW %>% group_by(Tree_cover) %>%
        summarise(Mean=round(mean(Damage_pct),1),
                  SE=round(sd(Damage_pct)/sqrt(n()),1)))
cat("\nDamage by treatment:\n")
print(raw_dmg %>% arrange(Tree_cover, Mean) %>%
        mutate(Treatment=TRT_SHORT[as.character(Treatment)]))
cat("\nCLMM model summary:\n")
print(summary(mod_dmg))
cat("\nCLD — Treatment | Tree Cover:\n"); print(cld_dmg_trt)
cat("\nCLD — Tree Cover:\n");             print(cld_dmg_tree)
cat("\nCLD — Week:\n");                   print(cld_dmg_wk)

cat("\n── FAW LARVAL DENSITY ──────────────────────────────────\n")
print(raw_lar %>% arrange(Tree_cover, Mean) %>%
        mutate(Treatment=TRT_SHORT[as.character(Treatment)]))
cat("\nNB GLMM summary:\n"); print(summary(mod_larv))
cat("\nCLD — Treatment | Tree Cover:\n"); print(cld_larv_trt)
cat("\nCLD — Tree Cover:\n");             print(cld_larv_tree)
cat("\nCLD — Week:\n");                   print(cld_larv_wk)

cat("\n── EGG MASSES ──────────────────────────────────────────\n")
print(raw_egg %>% arrange(Tree_cover, Mean) %>%
        mutate(Treatment=TRT_SHORT[as.character(Treatment)]))
cat("\nNB GLMM summary:\n"); print(summary(mod_eggs))
cat("\nCLD — Treatment | Tree Cover:\n"); print(cld_eggs_trt)
cat("\nCLD — Tree Cover:\n");             print(cld_eggs_tree)
cat("\nCLD — Week:\n");                   print(cld_eggs_wk)

cat("\n── MAIZE YIELD ─────────────────────────────────────────\n")
print(raw_yld %>% arrange(Tree_cover, Mean) %>%
        mutate(Treatment=TRT_SHORT[as.character(Treatment)]))
cat("\nLMM Type III ANOVA:\n"); print(anova(mod_yld, type="III"))
cat("\nCLD — Treatment | Tree Cover:\n"); print(cld_yld_trt)
cat("\nCLD — Tree Cover:\n");             print(cld_yld_tree)
cat(sprintf("\nICC = %.4f (%.1f%% Farm blocking)\n", icc, icc*100))

sink()
cat("  ✓ Results_Numbers.txt\n")


# ── 14. COMPLETION SUMMARY ────────────────────────────────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  ANALYSIS COMPLETE — SEASON 1 ONLY\n")
cat("  Output folder:", out_dir, "\n\n")

cat("  FIGURES  (TIFF 600 dpi + PNG 300 dpi)\n")
figs <- c(
  "Fig1_Damage_Over_Time",
  "Fig2_Damage_CLD",
  "Fig3_Larvae_CLD",
  "Fig4_Larval_Heatmap",
  "Fig5_EggMasses_CLD",
  "Fig6_Border_vs_Core",
  "Fig7_Yield_CLD",
  "Fig8_Damage_Yield_Scatter",
  "Fig9_Tillage_Cropping_Interaction",
  "Fig10_Raincloud_Damage"
)
for (f in figs) cat("  ·", f, ".tiff / .png\n", sep="")

cat("\n  TABLES\n")
cat("  · MS1_Tables.docx  (Tables 1–5, Word format)\n")
cat("  · MS1_Tables.xlsx  (5 sheets)\n")

cat("\n  CLD CSV FILES  (11 files)\n")
cat("  · CLD_Damage_Treatment/TreeCover/Week\n")
cat("  · CLD_Larvae_Treatment/TreeCover/Week\n")
cat("  · CLD_Eggs_Treatment/TreeCover/Week\n")
cat("  · CLD_Yield_Treatment/TreeCover\n")

cat("\n  DIAGNOSTICS\n")
cat("  · DHARMa_Larvae.png\n")
cat("  · DHARMa_Eggs.png\n")
cat("  · Results_Numbers.txt\n")

cat("\n  MODELS IN ENVIRONMENT\n")
cat("  · mod_dmg   — Ordinal CLMM  (FAW Damage)\n")
cat("  · mod_larv  — NB GLMM       (Larval Count)\n")
cat("  · mod_eggs  — NB GLMM       (Egg Masses)\n")
cat("  · mod_yld   — LMM           (Maize Yield)\n")
cat(strrep("═", 60), "\n")