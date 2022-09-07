source("figures_heatmaps.R")
source("data_gen_vizualization.R")

set.seed(0)

figure1 <- show_treatment_effect_and_variance_difference(n_control_neurons = 5, n_intervention_neurons = 5)

dat <- gen_gridded_dataset_draws(n_control_neurons = 6, n_intervention_neurons = 6, treatment_effect = 0)
figure3 <- plot_gridded_dataset_draws(dat, FALSE)
figure4 <- plot_gridded_dataset_draws(dat, TRUE, TRUE)

#Since these plots are only used for illustrative purposes, set the seed
set.seed(0)
figure6 <- show_ideal_vs_realistic_paired_datasets()
figure7 <- show_sources_of_variance_in_paired_datasets()


stop()

effect_size_interv_ctrl <- plot_heatmaps(x_axis = "effect_size", y_axis = "intervention_variance_ratio",
                                         y_axis_min=0.5, y_axis_max=1.5, proportion_explainable_variance = 0)

Nn_var_ratio_100samples <- plot_heatmaps(x_axis = "number_of_neurons", y_axis = "proportion_explainable_variance", total_sample_size = 100,
                              x_axis_min = 2, x_axis_max = 20)

Nn_var_ratio_1000samples <- plot_heatmaps(x_axis = "number_of_neurons", y_axis = "proportion_explainable_variance", total_sample_size = 1000,
                                         x_axis_min = 2, x_axis_max = 100)

paired_Nn<- plot_heatmaps(x_axis = "number_of_neurons",
                          x_axis_min = 2,
                          x_axis_max = 100,
                          y_axis = "prop_explainable_variance_between_pairs",
                          y_axis_min = 0.05,
                          y_axis_max = 0.95,
                          total_sample_size = 1000,
                          paired=TRUE,
                          prop_explainable_variance_between_pairs=0.30)


figure2 <- effect_size_interv_ctrl$ks_only

figure5 <- Nn_var_ratio_1000samples$ks_only

figure6 <- effect_size_interv_ctrl$coloredraster

figure7 <- Nn_var_ratio$coloredraster


figures = list(figure1, figure2, figure3, figure4, figure5, figure6, figure7)
fig_width <- 7

fig_heights = c(2, 4.5, 6, 6, 4.5, 7, 7)

fignames <- sapply(1:length(figures), function(x) paste("figure", x))

dir.create("vector_figures", showWarnings = FALSE)
format = "pdf"
for (fig in 1:length(figures)){
    fig_object <- figures[[fig]]
    filename   <- paste("vector_figures/",fignames[[fig]], ".", format, sep="")
    fig_height <- fig_heights[[fig]]
    ggsave(filename, fig_object, device = format, width = fig_width, height = fig_height, units = "in")
}

supplimentry_table <- print(variance_decomp, digits = 2)

save(figures, dat, fig_width, fig_heights, effect_size_interv_ctrl, Nn_var_ratio_100samples, Nn_var_ratio_1000samples, paired_Nn, file="figures.RData")
