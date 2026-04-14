library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(data.table)
library(ggplot2)

base_dir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_4"
file_ids <- 1:13
source("/projectnb/dietzelab/guYANG/pecan/runners/test10/diag_functions.R")

###### SDA Validation: Comparison between forecast and analysis
res <- build_sda_compare_data(
  base_dir = base_dir,
  file_ids = file_ids
)
### forecast vs analysis
compare_dt <- res$compare_dt
compare_dt[variable == "NEE", mean_value_analysis := nee_model_to_obs(mean_value_analysis)]
compare_dt[variable == "NEE", sd_value_analysis   := nee_model_to_obs(sd_value_analysis)]

### Stats validation for forecast and analysis
stats_dt <- calc_error_stats(compare_dt)

### Total comparison plot
p1 <- plot_obs_vs_fc_an(compare_dt)
ggsave(
  filename = file.path(base_dir, "obs_vs_forecast_vs_analysis_all.png"),
  plot = p1, width = 16, height = 10, dpi = 300
)

### Selected sites comparison plot
selected_sites <- c("4421", "4480", "6805")
p2 <- plot_obs_vs_fc_an(compare_dt, sites = selected_sites)
ggsave(
  filename = file.path(base_dir, "obs_vs_forecast_vs_analysis_selected_sites.png"),
  plot = p2, width = 14, height = 9, dpi = 300
)

### Selected variable comparison plot
selected_vars <- c("NEE", "LAI")
p3 <- plot_obs_vs_fc_an(compare_dt, variables = selected_vars)
ggsave(
  filename = file.path(base_dir, "obs_vs_forecast_vs_analysis_selected_vars.png"),
  plot = p3, width = 14, height = 8, dpi = 300
)

### Scatterplot for all
p4 <- plot_scatter_obs_vs_model(compare_dt)
ggsave(
  filename = file.path(base_dir, "obs_vs_model_scatter.png"),
  plot = p4, width = 12, height = 8, dpi = 300
)

###### MCMC diagnose: Extract all diagnostics
all_diag <- map_dfr(file_ids, extract_from_one_file) %>%
  left_join(site_map_all, by = c("file_id", "year", "block_id")) %>%
  dplyr::select(
    par_name,
    rhat,
    ess,
    param_type,
    file_id,
    year,
    block_id,
    site_id,
    state_index,
    state_var,
    used_in_H
  )

###### Fix this bug
res <- diagnose_mcmc_convergence(all_diag,
                                 rhat_bad = 1.05,
                                 ess_bad = 400)

df_diag <- as.data.frame(res$by_param_type)

mapping_table <- all_diag %>%
  dplyr::select(block_id, site_id) %>%
  distinct(block_id, .keep_all = TRUE)
df_diag_extended <- df_diag %>%
  left_join(mapping_table, by = "block_id") %>%
  relocate(site_id, .after = block_id)
mapping_lc <- tiff %>%
  dplyr::select(index, LC) %>%
  distinct(index, .keep_all = TRUE) %>%
  mutate(site_id = as.character(index)) %>% 
  dplyr::select(-index)
df_diag_final <- df_diag_extended %>%
  mutate(site_id = as.character(site_id)) %>% 
  left_join(mapping_lc, by = "site_id") %>%
  relocate(LC, .after = site_id)
df_lc_summary <- df_diag_final %>%
  group_by(LC) %>%
  summarise(
    avg_frac_bad = mean(frac_bad, na.rm = TRUE),
    n_blocks = n()
  )
ggplot(df_lc_summary, aes(x = factor(LC), y = avg_frac_bad, fill = factor(LC))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0("n=", n_blocks)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "MCMC Failure Rate by Land Cover Type",
       x = "Land Cover Type (LC)",
       y = "Mean Fraction of Bad Parameters",
       fill = "LC Type")




