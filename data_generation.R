
#' gen_unpaired_data
#' Generates a set of mock miniature synaptic event data.
#' The data consists of 2 groups, a control and intervention group.
#' Each group contains a number of neurons.
#' Each neuron is associated with several samples.
#'
#' @param control_group_mean The baseline mean of the control group
#' @param treatment_effect  The mean difference between control and intervention groups
#' @param control_group_interneuron_sd  the SD between the means of individual neurons in the control group
#' @param intervention_group_interneuron_sd the SD between the means of individual neurons in the intervention group
#' @param within_neuron_sd the SD between samples taken from a single neuron
#' @param n_control_neurons number of neurons in the control group
#' @param n_intervention_neurons number of neurons in the intervention group
#' @param samples_per_neuron number of samples per neuron
#'
#' @return a dataframe with the columns (group, neuron_id, dependant), where each row is an observation
#' 

HALF_UNIT_VARIANCE_SD <- sqrt(2)/2


gen_unpaired_data <- function(control_group_mean=0, treatment_effect=1, control_group_interneuron_sd = 2, 
                                   intervention_group_interneuron_sd=2, within_neuron_sd=2, 
                                   n_control_neurons=10, n_intervention_neurons=10, samples_per_neuron=100) {
  

  #Making some basic calculations explicit
  intervention_group_mean <- control_group_mean + treatment_effect
  n_total_neurons <- n_control_neurons + n_intervention_neurons
  
  #Generate control group neuron means about 0, repeat each mean for all the samples
  control_grp_neuron_means <- rnorm(n_control_neurons, mean=0, sd = control_group_interneuron_sd)
  neuron_effect <- rep(control_grp_neuron_means, each = samples_per_neuron)
  
  #All the neurons have the same within-neuron variance, so it suffices to add a residual noise term
  sample_effect <- rnorm(n_control_neurons*samples_per_neuron, mean=0, sd = within_neuron_sd)
  
  #Then the samples are the sum of the group mean, the neuron's mean, and the sample effect
  control_group_samples <- neuron_effect + sample_effect + control_group_mean
  
  
  #Repeat for the intervention group
  intervention_grp_neuron_means <- rnorm(n_intervention_neurons, mean=0, sd = intervention_group_interneuron_sd)
  neuron_effect <- rep(intervention_grp_neuron_means, each = samples_per_neuron)
  sample_effect <- rnorm(n_intervention_neurons*samples_per_neuron, mean=0, sd = within_neuron_sd)
  
  intervention_group_samples <- neuron_effect + sample_effect + intervention_group_mean
  
  #Glue control and intervention groups together
  all_samples <- c(control_group_samples, intervention_group_samples)
  
  #Add the independant vars
  neuron_id <- rep(1:n_total_neurons, each = samples_per_neuron)
  neuron_id <- as.factor(neuron_id)
  
  group <- c(rep("Control",n_control_neurons*samples_per_neuron),rep("Intervention",n_intervention_neurons*samples_per_neuron))
  
  #Return as dataframe
  dat <- data.frame(neuron_id = neuron_id, group = group, dependant = all_samples)
  
  return(dat)
}


plot_data <- function(dat){
  #Plot each neuron's empirical PDF (unfilled histogram) colored by that neuron's group
  plt <- ggplot(data=dat, aes(x = dependant, color = group))
  for (id in unique(dat$neuron_id)){
    plt <- plt + geom_freqpoly(data = dat[dat$neuron_id==id,], binwidth = 0.2)
  }
  
  return(plt)
}


is_ks_test_positive <- function(dat, alpha = 0.05){
  #Perform the KS test
  res <- ks.test(dat$dependant[dat$group=="Control"],dat$dependant[dat$group=="Intervention"])
  ks_positive <- (res$p.value < alpha)
  return(ks_positive)
}


is_lmer_test_positive <- function(dat, alpha = 0.05, suppress = FALSE){
  #Create a mixed model, then create a mixed model with a dropped term, then compare
  if (suppress){
    return(suppressMessages(is_lmer_test_positive(dat,alpha)))
  }
  model_h1 <- lme4::lmer(dependant ~ group + (1|neuron_id), data = dat, REML = FALSE)
  model_h0 <- lme4::lmer(dependant ~ (1|neuron_id),         data = dat, REML = FALSE)
  
  res <- anova(model_h1, model_h0)
  #get the p-value and compare to alpha
  lmer_positive <- (res$`Pr(>Chisq)`[2] < alpha)
  return(lmer_positive)
}

