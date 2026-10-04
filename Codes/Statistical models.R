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