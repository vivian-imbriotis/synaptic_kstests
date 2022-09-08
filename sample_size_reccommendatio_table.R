source("data_generation.R")

produce_sample_size_rec_table <- function(){
  set.seed(0)
  N = 500
  
  TOTAL_VAR = 1
  BETWEEN_NEURON_VAR = 0.1
  WITHIN_NEURON_VAR = TOTAL_VAR - BETWEEN_NEURON_VAR
  
  BETWEEN_NEURON_SD <- sqrt(BETWEEN_NEURON_VAR)
  WITHIN_NEURON_SD  <- sqrt(WITHIN_NEURON_VAR)
  
  TREATMENT_EFFECT <- 0.5 * sqrt(TOTAL_VAR)  #Inverting the formula for cohen's D == 0.5
  
  record <- data.frame()
  args   <- list()
  
  args$intervention_between_neuron_sd <- BETWEEN_NEURON_SD
  args$control_between_neuron_sd      <- BETWEEN_NEURON_SD
  args$control_within_neuron_sd          <- WITHIN_NEURON_SD
  args$intervention_within_neuron_sd     <- WITHIN_NEURON_SD
  args$treatment_effect                  <- TREATMENT_EFFECT
  args$suppress <- T
  for (n_neurons in c(10,25)){
    for (samples_per_neuron in c(10,40,100)){
      for(effect_size in c(0.2,0.5)){
  
        args$n_samples = N
        args$paired    = F
        
        args$n_control_neurons                 <- n_neurons %/% 2 + n_neurons %% 2
        args$n_intervention_neurons            <- n_neurons %/%2
        args$samples_per_neuron                <- samples_per_neuron
      
        pow <- do.call(get_power, args)$lmer
        fpr <- do.call(get_fpr,   args)$lmer
        
        
        
        row <- list()
        row$n_neurons <- n_neurons
        row$samples_per_neuron <- samples_per_neuron
        row$effect_size <- effect_size
        row$power <- pow[[1]]
        row$power_lwr <- pow[[2]]
        row$power_upr <- pow[[3]]
        row$fpr <- fpr[[1]]
        row$fpr_lwr <- fpr[[2]]
        row$fpr_upr <- fpr[[3]]
        
        record <- rbind(record, row)
        
      }
    }
  }
  
  sample_size_recs_table <- subset(record, !(n_neurons==10 & samples_per_neuron==40) & !(n_neurons==25 & samples_per_neuron==100))
  return(sample_size_recs_table)
}