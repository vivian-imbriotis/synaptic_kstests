source("data_generation.r")


gen_gridded_dataset_draws <- function(nrows=3, ncols=3){
  rows = 1:nrows
  cols = 1:ncols
  
  between_neuron_vars <- seq(0.1,1,length=nrows)
  within_neuron_vars <- c(0.1,1,length=ncols)
  
  datasets <- list()
  
  for (row in rows){
    for (col in cols){
      between <- between_neuron_vars[[row]]
      within <- within_neuron_vars[[col]]
      
      dat<- gen_unpaired_data(control_group_interneuron_sd = sqrt(between), intervention_group_interneuron_sd = sqrt(between),
                              within_neuron_sd = sqrt(within), treatment_effect = 1)
      
      dat$between_neuron_variance <- between
      dat$within_neuron_variance  <- within
      
      datasets[[(row-1)*length(cols) + col]] <- dat
    }
  }
  
  
  longform_datasets <- Reduce(rbind, datasets)
}

plot_gridded_dataset_draws <- function(longform_datasets, marginalize_neurons = FALSE){
  
  
  
  library(ggplot2)
  plt <- ggplot(data = longform_datasets, aes(x = dependant, color = group))
  
  if (marginalize_neurons){
    plt <- plt + geom_density()
  }else{
    for (neuron in unique(longform_datasets$neuron_id)){
      plt <- plt + geom_density(data = longform_datasets[longform_datasets$neuron_id == neuron,])
      
    }
  }
  
  plt <- plt + facet_grid(rows = vars(between_neuron_variance), cols = vars(within_neuron_variance), labeller = label_both)
  return(plt)
}

dat <- gen_gridded_dataset_draws()
unmarginalized <- plot_gridded_dataset_draws(dat, FALSE)
marginalized   <- plot_gridded_dataset_draws(dat, TRUE)
print(marginalized)