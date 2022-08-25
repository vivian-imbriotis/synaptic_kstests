source("data_generation.R")
library(ggplot2)
library(patchwork)
require(beepr)

ALPHA = 0.05

HALF_UNIT_VARIANCE_SD <- sqrt(2)/2



get_fpr <- function(n_samples, saveplot = FALSE, ...){
  ks_positives   <- 0
  lmer_positives <- 0
  for(i in 1:n_samples){
    #Generate a mock dataset with the requested properties
    dat <- gen_unpaired_data(treatment_effect = 0, ...)
    
    if (is_ks_test_positive(dat)){
      ks_positives <- ks_positives + 1
    }
    if (is_lmer_test_positive(dat)){
      lmer_positives <- lmer_positives + 1
    }
  }
  ks_fpr   <- ks_positives / n_samples
  lmer_fpr <- lmer_positives / n_samples
  return(c(ks_fpr,lmer_fpr))
}

get_power <- function(n_samples, treatment_effect, saveplot = FALSE, ...){
  if (treatment_effect==0){
    stop("Attempted to find power with no treatment effect")
  }
  ks_postives   <- 0
  lmer_postives <- 0
  for(i in 1:n_samples){
    #Generate a mock dataset with the requested properties
    dat <- gen_unpaired_data(treatment_effect = treatment_effect, ...)
    
    if (is_ks_test_positive(dat)){
      ks_positives <- ks_positives + 1
    }
    if (is_lmer_test_positive(dat)){
      lmer_positives <- lmer_positives + 1
    }
  }
  ks_power   <- ks_positives / n_samples
  lmer_power <- lmer_positives / n_samples
  return(c(ks_fpr,lmer_fpr))
}


check_test_performances_while_varying_data_parameter <- function(varying_param = "interneuron variance", 
                                                                 param_values = seq(from=0,to=1,by=0.1),
                                                                 alpha = ALPHA,
                                                                 n_samples = 50,
                                                                 reserve_plots=FALSE,
                                                                 power_only = F,
                                                                 fpr_only = F,
                                                                 n_neurons = 200,
                                                                 ...){

  if(fpr_only && power_only){stop()}
  
  if(reserve_plots){
    dir.create("datasets", showWarnings = FALSE)
  }
  
  #This record is going to store the FPR/power/both for our tests
  record = data.frame(parameter_value = numeric(0),               #The data-generating parameter value
                      test = character(0),                        #{"KS Test", "Mixed Effects"}
                      test_stat = character(0),                   #{"False Positive Rate", "Power"}
                      test_stat_value = numeric(0))               #The FPR/power value
  

  
  #Half our neurons are controls, the other half are intervention
  n_control_neurons = n_neurons%/%2 + n_neurons%%2
  n_intervention_neurons = n_neurons%/%2

  #For each value of our data parameter, we'll generate a bunch of datasets and test each one
  for (parameter_idx in 1:length(param_values)){
    
    cat("Simulating with parameter set to",
        (as.character(parameter_idx)),"of",
        as.character(length(param_values)),
        "possible values\n")
    
    param = param_values[[parameter_idx]]

    if (varying_param == "interneuron variance"){
      pow <- get_power(control_group_interneuron_sd = sqrt(param), intervention_group_interneuron_sd = sqrt(param),
                             treatment_effect = effect, ...)
      fpr <- get_fpr(control_group_interneuron_sd = sqrt(param), intervention_group_interneuron_sd = sqrt(param),
                     treatment_effect = effect, ...)
    }else if(varying_param=="intraneuron variance"){
      pow <- get_power(within_neuron_sd = sqrt(param),treatment_effect = effect, ...)
      fpr <- get_fpr(within_neuron_sd = sqrt(param),treatment_effect = effect, ...)
    }else if(varying_param=="effect size"){
      pow <- get_power(treatment_effect = param, ...)
      fpr <- get_fpr(...)
    }else{
      stop("Argument 'mode' must be one of interneuron variance, intraneuron variance, or effect size.")
    }
    
    # #This stores the test results for the datasets generated with a particular value of the data parameter
    # tmp <- data.frame(ground_truth = TRUE, ks_positive = TRUE, lmer_positive = TRUE)
    # 
    # for (i in 1:samples_per_param_value){
    #   if (power_only){
    #     ground_truth = T
    #   }else if (fpr_only){
    #     ground_truth = F
    #   }else{
    #   #A treatment effect is present for the first N/2 sim runs, then absent for the remaining runs
    #   ground_truth = i<(samples_per_param_value / 2)
    #   }
    #   
    #   #The effect size is constant (if present) unless its the varying parameter...
    #   if (varying_param == "effect size"){
    #     effect<-as.integer(ground_truth) * param
    #   }else{
    #     effect <- as.integer(ground_truth) * treatment_effect
    #   }
    #   
    #   #Generate a mock dataset
    #   if (varying_param == "interneuron variance"){
    #     dat <- gen_unpaired_data(control_group_interneuron_sd = sqrt(param), intervention_group_interneuron_sd = sqrt(param),
    #                            treatment_effect = effect, ...)
    #   }else if(varying_param=="intraneuron variance"){
    #     dat <- gen_unpaired_data(within_neuron_sd = sqrt(param),treatment_effect = effect, ...)
    #   }else if(varying_param=="effect size"){
    #     dat <- gen_unpaired_data(treatment_effect = effect, ...)
    #   }else{
    #     stop("Argument 'mode' must be one of interneuron variance, intraneuron variance, or effect size.")
    #   }
    #   
    #   #if the reserve_plots flag is set, save a vizualization of the first example of each kind of dataset  
    #   if(reserve_plots){
    #     if(i==1 || (i==samples_per_param_value/2 && !fpr_only && !power_only)){
    #       plt<-plot_data(dat); ggsave(filename=paste("datasets/",varying_param,"=",as.character(param),as.character(ground_truth),".png"),
    #                                   plot=plt,device="png", w = 7, h = 7)}
    #   }
    #   
    #   #Test the dataset
    #   ks_positive <- is_ks_test_positive(dat)
    #   lmer_positive <- is_lmer_test_positive(dat, suppress = TRUE)
    #   
    #   #Store the results for this particular dataset
    #   tmp[i,"ground_truth"]  <- ground_truth
    #   tmp[i,"ks_positive"]   <- ks_positive
    #   tmp[i,"lmer_positive"] <- lmer_positive
    # }
    
    #Record the FPR/power/both for the group of datasets with this particular data-generating parameter
    if(!power_only){
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("KS Test", "False Positive Rate", sum((!tmp$ground_truth) & tmp$ks_positive)   / sum(!tmp$ground_truth))
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("Mixed Effects", "False Positive Rate", sum((!tmp$ground_truth) & tmp$lmer_positive) / sum(!tmp$ground_truth))
    }
    
    if(!fpr_only){
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("KS Test", "Power", 1 - (sum((tmp$ground_truth) & !tmp$ks_positive)   / sum(tmp$ground_truth)))
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("Mixed Effects", "Power", 1 - (sum((tmp$ground_truth) & !tmp$lmer_positive) / sum(tmp$ground_truth)))
    }
  }
  #Add a column to our record indicating the parameter values used. This dataframe is in long format, hence the repeats.
  record$parameter_value <- rep(param_values, each=if(power_only||fpr_only) 2 else 4)

  return(record)
}



plot_test_performances_against_varying_data_parameter <- function(varying_param = "interneuron variance",treatment_effect=1,
                                                                  control_group_interneuron_sd = 1,
                                                                  within_neuron_sd = 1,...){
  
  record <- check_test_performances_while_varying_data_parameter(varying_param = varying_param, treatment_effect = treatment_effect,
                                                                 control_group_interneuron_sd = control_group_interneuron_sd,
                                                                 within_neuron_sd = within_neuron_sd,
                                                                 ...)
  record$test_stat_value <- as.numeric(record$test_stat_value)
  base <- ggplot(data = record, aes(x = parameter_value, y = as.numeric(test_stat_value), color = test))
  base <- base + scale_x_continuous(labels = scales::percent)
  
  if(varying_param=="interneuron variance"){
    xlabel<- "Within-group, between-neuron variance (%ATE)"
    base <- base + scale_y_continuous(labels = scales::percent)
    record$test_stat_value <- record$test_stat_value / treatment_effect
  }else if(varying_param=="intraneuron variance (%ATE)"){
    xlabel<- "Within-neuron variance (%total variance)"
    base <- base + scale_y_continuous(labels = scales::percent)
    record$test_stat_value <- record$test_stat_value / treatment_effect
  }else if(varying_param=="effect size"){
    xlabel<- "Effect size (Cohen's D)"
    record$test_stat_value <- record$test_stat_value / sqrt(control_group_interneuron_sd**2 + within_neuron_sd**2)
    
  }else{
    stop("Argument 'varying_param' must be one of interneuron, intraneuron, or effect size")
  }
  
  labels <- labs(x = xlabel, y = "Rate", color = "Test")

  
  plt <- base + geom_line() + geom_point() + facet_grid(cols=vars(test_stat)) + labels
  

  return(plt) 
}

# interneuron <- plot_test_performances_against_varying_data_parameter("interneuron variance", param_values = seq(from=0.1,to=2.1,by=0.2),
#                                                                      samples_per_param_value = 1000, n_neurons = 20, reserve_plots=TRUE)
# intraneuron <- plot_test_performances_against_varying_data_parameter("intraneuron variance", samples_per_param_value = 100)


# plt <- plot_test_performances_against_varying_data_parameter(varying_param = "interneuron variance", 
#                                                              param_values = seq(0,0.6,0.1), 
#                                                              samples_per_param_value = 300)

effect_size <- plot_test_performances_against_varying_data_parameter("effect size", samples_per_param_value = 1000, reserve_plots = F,
                                                                     power_only = F, param_values = seq(from=0,to=1,by=0.2),
                                                                     n_neurons = 500, samples_per_neuron = 10,
                                                                     control_group_interneuron_sd = 0.711,
                                                                     intervention_group_interneuron_sd = 0.711,
                                                                     within_neuron_sd = 0.711)
print(effect_size)

# beepr::beep()
