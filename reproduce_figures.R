source("figures_heatmaps.R")
source("data_gen_vizualization.R")

set.seed(0)

figure1 <- show_treatment_effect_and_variance_difference(n_control_neurons = 5, n_intervention_neurons = 5)

dat <- gen_gridded_dataset_draws(n_control_neurons = 6, n_intervention_neurons = 6, treatment_effect = 0)
figure3 <- plot_gridded_dataset_draws(dat, FALSE)
figure4 <- plot_gridded_dataset_draws(dat, TRUE, TRUE)


effect_size_interv_ctrl <- plot_heatmaps(x_axis = "effect_size", y_axis = "intervention_variance_ratio",
                                         y_axis_min=0.5, y_axis_max=1.5, proportion_explainable_variance = 0)

Nn_var_ratio <- plot_heatmaps(x_axis = "number_of_neurons", y_axis = "proportion_explainable_variance", total_sample_size = 100,
                              x_axis_min = 2, x_axis_max = 20)

figure2 <- effect_size_interv_ctrl$ks_only

figure5 <- Nn_var_ratio$ks_only

figure6 <- effect_size_interv_ctrl$coloredraster

figure7 <- Nn_var_ratio$coloredraster


figures = list(figure1, figure2, figure3, figure4, figure5, figure6, figure7)
fig_width <- 7

fig_heights = c(2, 4.5, 6, 6, 4.5, 7, 7)

fignames <- sapply(1:length(figures), function(x) paste("figure", x))

dir.create("vector_figures", showWarnings = FALSE)
for (fig in 1:length(figures)){
  for (format in c("pdf", "svg")){
    fig_object <- figures[[fig]]
    filename   <- paste("vector_figures/",fignames[[fig]], ".", format, sep="")
    fig_height <- fig_heights[[fig]]
    ggsave(filename, fig_object, device = format, width = fig_width, height = fig_height, units = "in")
  }
}

save(figures, dat, fig_width, fig_heights, effect_size_interv_ctrl, Nn_var_ratio, "figures.RData")