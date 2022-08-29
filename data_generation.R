
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

DEFAULT_INTERNEURON_SD <- 0.34
DEFAULT_WITHIN_NEURON_SD <- 0.61


gen_unpaired_data <- function(control_group_mean=0, treatment_effect=1, 
                              control_group_interneuron_sd      = DEFAULT_INTERNEURON_SD, 
                              intervention_group_interneuron_sd = DEFAULT_INTERNEURON_SD, 
                              control_within_neuron_sd          = DEFAULT_WITHIN_NEURON_SD,
                              intervention_within_neuron_sd     = DEFAULT_WITHIN_NEURON_SD, 
                              n_control_neurons = 10, 
                              n_intervention_neurons = 10, 
                              samples_per_neuron = 100) {
  

  #Making some basic calculations explicit
  intervention_group_mean <- control_group_mean + treatment_effect
  n_total_neurons <- n_control_neurons + n_intervention_neurons
  
  #Generate control group neuron means about 0, repeat each mean for all the samples
  control_grp_neuron_means <- rnorm(n_control_neurons, mean=0, sd = control_group_interneuron_sd)
  neuron_effect <- rep(control_grp_neuron_means, each = samples_per_neuron)
  
  #All the neurons have the same within-neuron variance, so it suffices to add a residual noise term
  sample_effect <- rnorm(n_control_neurons*samples_per_neuron, mean=0, sd = control_within_neuron_sd)
  
  #Then the samples are the sum of the group mean, the neuron's mean, and the sample effect
  control_group_samples <- neuron_effect + sample_effect + control_group_mean
  
  
  #Repeat for the intervention group
  intervention_grp_neuron_means <- rnorm(n_intervention_neurons, mean=0, sd = intervention_group_interneuron_sd)
  neuron_effect <- rep(intervention_grp_neuron_means, each = samples_per_neuron)
  sample_effect <- rnorm(n_intervention_neurons*samples_per_neuron, mean=0, sd = intervention_within_neuron_sd)
  
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

gen_paired_data <- function(control_group_mean=0, treatment_effect=1, interneuron_sd = 2, within_neuron_sd=2, residual_interneuron_sd = 0.5,
                              n_neurons=10, samples_per_neuron=100) {
  
  intervention_group_mean <- control_group_mean + treatment_effect
  
  #Generate neuron means about 0, then add some gaussian noise to the observed means in each group
  fixed_neuron_means            <- rnorm(n_neurons, mean=0, sd = interneuron_sd)
  control_grp_neuron_means      <- fixed_neuron_means + rnorm(n_neurons, mean=0, sd = residual_interneuron_sd)
  intervention_grp_neuron_means <- fixed_neuron_means + rnorm(n_neurons, mean=0, sd = residual_interneuron_sd)
  
  #Repeat the neuron means
  neuron_effect <- rep(control_grp_neuron_means, each = samples_per_neuron)
  
  #All the neurons have the same within-neuron variance, so it suffices to add a residual noise term
  sample_effect <- rnorm(n_neurons*samples_per_neuron, mean=0, sd = within_neuron_sd)
  
  #Then the samples are the sum of the group mean, the neuron's mean, and the sample effect
  control_group_samples <- neuron_effect + sample_effect + control_group_mean
  
  
  #Repeat for the intervention group
  neuron_effect <- rep(intervention_grp_neuron_means, each = samples_per_neuron)
  sample_effect <- rnorm(n_neurons*samples_per_neuron, mean=0, sd = within_neuron_sd)
  intervention_group_samples <- neuron_effect + sample_effect + intervention_group_mean
  
  #Glue control and intervention groups together
  all_samples <- c(control_group_samples, intervention_group_samples)
  
  #Build the exogenous variables. Each neuron was sampled N times, first in the control, then in the intervention group
  neuron_id <- rep(1:n_neurons, each = samples_per_neuron)
  neuron_id <- rep(neuron_id, 2)
  neuron_id <- as.factor(neuron_id)
  
  group <- c(rep("Control",n_neurons*samples_per_neuron),rep("Intervention",n_neurons*samples_per_neuron))
  
  #Return as dataframe
  dat <- data.frame(neuron_id = neuron_id, group = group, dependant = all_samples)
  
  return(dat)
}

plot_data <- function(dat, colorby = "group"){
  #Plot each neuron's empirical PDF (unfilled histogram) colored by that neuron's group
  if (colorby=="group"){
  plt <- ggplot(data=dat, aes(x = dependant, color = group))
  }else if(colorby=="neuron id"){
    plt <- ggplot(data=dat, aes(x = dependant, color = neuron_id, linetype = group))
  }else{
    stop("in function plot_data, argument colorby must be one of {'group','neuron_id'}")
  }
  for (grp in unique(dat$group)){
    grp_data <- subset(dat, group==grp)
    for (id in unique(grp_data$neuron_id)){
      plt <- plt + geom_freqpoly(data = subset(grp_data,neuron_id==id), binwidth = 0.2)
    }
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
  #Suppress model fitting warnings (when fitting many models)
  if (suppress){
    return(suppressMessages(is_lmer_test_positive(dat,alpha)))
  }
  #Create a mixed model, then create a mixed model with a dropped term, then compare
  model_h1 <- lme4::lmer(dependant ~ group + (1|neuron_id), data = dat, REML = FALSE)
  model_h0 <- lme4::lmer(dependant ~ (1|neuron_id),         data = dat, REML = FALSE)
  
  res <- anova(model_h1, model_h0)
  #get the p-value and compare to alpha
  lmer_positive <- (res$`Pr(>Chisq)`[2] < alpha)
  return(lmer_positive)
}



get_fpr <- function(n_samples, saveplot = FALSE, suppress=FALSE, 
                    paired = FALSE, filename="dataset", ...){
  ks_positives   <- 0
  lmer_positives <- 0
  dots <- list(...)
  dots["treatment_effect"] <- 0
  for(i in 1:n_samples){
    #Generate a mock dataset with the requested properties
    if(paired){
      dat <- do.call(gen_paired_data, dots)
    }else{
      dat <- do.call(gen_unpaired_data, dots)
    }
    if(saveplot&&i==1){
      plt<-plot_data(dat, colorby="neuron id"); ggsave(filename=paste("datasets/",filename,".png"),
                                  plot=plt,device="png", w = 7, h = 7)
    }
    
    if (is_ks_test_positive(dat)){
      ks_positives <- ks_positives + 1
    }
    if (is_lmer_test_positive(dat, suppress=suppress)){
      lmer_positives <- lmer_positives + 1
    }
  }
  #Get point estimates and CIs
  ks_fpr   <- c(ks_positives / n_samples, c(prop.test(ks_positives,n_samples)$conf.int))
  lmer_fpr <- c(lmer_positives / n_samples, c(prop.test(lmer_positives,n_samples)$conf.int))
  return(data.frame(ks=ks_fpr,lmer=lmer_fpr, row.names = c("value", "lower_CI", "upper_CI")))
}

get_power <- function(n_samples, saveplot = FALSE, paired = FALSE, filename = "dataset", suppress = FALSE, ...){
  ks_positives   <- 0
  lmer_positives <- 0
  for(i in 1:n_samples){
    #Generate a mock dataset with the requested properties
    if (paired){
      dat <- gen_paired_data(...)
    }else{
      dat <- gen_unpaired_data(...)
    }
    if(saveplot&&i==1){
      plt<-plot_data(dat,colorby="neuron id"); ggsave(filename=paste("datasets/",filename,".png"),
                                  plot=plt,device="png", w = 7, h = 7)
    }
    if (is_ks_test_positive(dat)){
      ks_positives <- ks_positives + 1
    }
    if (is_lmer_test_positive(dat, suppress=suppress)){
      lmer_positives <- lmer_positives + 1
    }
  }
  ks_pow   <- c(ks_positives / n_samples, c(prop.test(ks_positives,n_samples)$conf.int))
  lmer_pow <- c(lmer_positives / n_samples, c(prop.test(lmer_positives,n_samples)$conf.int))
  return(data.frame(ks=ks_pow,lmer=lmer_pow, row.names = c("value", "lower_CI", "upper_CI")))
}
