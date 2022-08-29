source("data_generation.r")



gen_gridded_dataset_draws <- function(nrows=3, ncols=3, paired = T, ...){
  rows = 1:nrows
  cols = 1:ncols
  
  between_neuron_vars <- seq(0.1,1,length=nrows)
  within_neuron_vars <- seq(0.1,1,length=ncols)
  
  datasets <- list()
  
  for (row in rows){
    for (col in cols){
      between <- between_neuron_vars[[row]]
      within <- within_neuron_vars[[col]]
      if(paired){
        dat<- gen_paired_data(interneuron_sd = sqrt(between), within_neuron_sd = sqrt(within), ...)
      }else{      
        dat<- gen_unpaired_data(control_group_interneuron_sd = sqrt(between), intervention_group_interneuron_sd = sqrt(between),
                                          within_neuron_sd = sqrt(within), ...)
      }

      
      dat$between_neuron_variance <- between
      dat$within_neuron_variance  <- within
      
      datasets[[(row-1)*length(cols) + col]] <- dat
    }
  }
  
  
  longform_datasets <- Reduce(rbind, datasets)
}

plot_gridded_dataset_draws <- function(longform_datasets, marginalize_neurons = FALSE, cumulative = FALSE, paired = F){
  
  
  
  library(ggplot2)
  plt <- ggplot(data = longform_datasets, aes(x = dependant, color = group))
  
  if (marginalize_neurons && !cumulative){
    plt <- plt + geom_density()
  }else if (marginalize_neurons && cumulative){
    plt <- plt + stat_ecdf(geom="step")
  }else{
    if(paired){
      plt <- plt + geom_density(aes(color=neuron_id, linetype = group))
    }else{
      #Keep all the neurons in a group the same color
      for (neuron in unique(longform_datasets$neuron_id)){
        plt <- plt + geom_density(data = longform_datasets[longform_datasets$neuron_id == neuron,])
        
      }
    }
  }
  
  plt <- plt + facet_grid(rows = vars(between_neuron_variance), cols = vars(within_neuron_variance), labeller = label_both)
  return(plt)
}

dat <- gen_gridded_dataset_draws(treatment_effect = 0)
unmarginalized <- plot_gridded_dataset_draws(dat, FALSE)
cumulative     <- plot_gridded_dataset_draws(dat, TRUE, TRUE)
marginalized   <- plot_gridded_dataset_draws(dat, TRUE, FALSE)

paired_dat <- gen_gridded_dataset_draws(n_neurons = 3, treatment_effect = 0, paired = T, residual_interneuron_sd = 0.2)
paired_unmarginalized <- plot_gridded_dataset_draws(paired_dat, FALSE , paired = T)
paired_cumulative     <- plot_gridded_dataset_draws(paired_dat, TRUE, TRUE, paired = T)
paired_marginalized   <- plot_gridded_dataset_draws(paired_dat, TRUE, FALSE, paired = T)
