source("data_generation.R")
library(ggplot2)
library(patchwork)
require(beepr)

ALPHA = 0.05

HALF_UNIT_VARIANCE_SD <- sqrt(2)/2


#' Check test performance while varying a single parameter
#' 
#'
#' @param varying_param One of 'interneuron variance', 'intraneuron variance', 'effect size', 'intervention group variance'. See details.
#' @param param_values The values varying_param should take on
#' @param alpha The critical value, defaults to 0.05
#' @param n_samples The number of datasets to generate for each value varying_param takes on
#' @param reserve_plots If True, save a figure for a representative example of each dataset
#' @param power_only If True, only determine the power (not the false positive rate)
#' @param fpr_only If True, only determine the false positive rate (not the power)
#' @param n_neurons The number of neurons in each generated dataset (if paired=FALSE, they will be evenly split between control and intervention groups)
#' @param paired Whether the datasets should be paired or unpaired. See ?data_generation::gen_paired_data and ?data_generation::gen_unpaired_data
#' @param ... Further parameters are passed on to gen_paired_data if paired=TRUE and gen_unpaired_data otherwise
#'
#' @details This figure holds all but one data-generating parameter constant, varying that one parameter.
#' If varying_param='interneuron variance', then the variance between neuronal means is varied.
#' If varying_param='intraneuron variance', then the variance between different observations from the same neuron is varied
#' If varying_param='effect size', then the effect size is varied, holding all variances in the data constant.
#' If varying_param='intervention group variance' is varied, then the intervention group's variance is varied, holding the control group's variance constant.
#' 
#' @return data.frame object
check_test_performances_while_varying_data_parameter <- function(varying_param = "interneuron variance", 
                                                                 param_values = seq(from=0,to=1,by=0.1),
                                                                 alpha = ALPHA,
                                                                 n_samples = 50,
                                                                 reserve_plots=FALSE,
                                                                 power_only = F,
                                                                 fpr_only = F,
                                                                 n_neurons = 200,
                                                                 paired = F,
                                                                 ...){

  if(fpr_only && power_only){stop()}
  
  if(reserve_plots){
    dir.create("datasets", showWarnings = FALSE)
  }
  
  n_control_neurons <- n_neurons %/% 2 + n_neurons %% 2
  n_intervention_neurons <- n_neurons %/% 2
  
  #This record is going to store the FPR/power/both for our tests
  record = data.frame(parameter_value = numeric(0),               #The data-generating parameter value
                      test = character(0),                        #{"KS Test", "Mixed Effects"}
                      test_stat = character(0),                   #{"False Positive Rate", "Power"}
                      test_stat_value = numeric(0))               #The FPR/power value

  #Set constant parameters to our data generating function
  kwargs <- list(...)
  kwargs["n_samples"] <- n_samples
  kwargs["paired"] <- paired
  if(paired){
    kwargs["n_neurons"] <- n_neurons
  }else{  
    kwargs["n_control_neurons"] <- n_control_neurons
    kwargs["n_intervention_neurons"] <- n_intervention_neurons
  }
  
  #For each value our parameter takes on...
  for (parameter_idx in 1:length(param_values)){
    

    #Do a basic printout
    cat("Simulating with parameter set to",
        (as.character(parameter_idx)),"of",
        as.character(length(param_values)),
        "possible values\n")
    
    #Fetch the relevant variable parameter value
    param = param_values[[parameter_idx]]
    
    #Set the variable parameter(s)
    if (varying_param == "interneuron variance"){
      if(paired){
        kwargs["interneuron_sd"] <- sqrt(param)
      }else{        
        kwargs["control_group_interneuron_sd"] <- sqrt(param)
        kwargs["intervention_group_interneuron_sd"] <- sqrt(param)}
    }else if(varying_param=="intraneuron variance"){
        kwargs["within_neuron_sd"] <- sqrt(param)
    }else if(varying_param=="effect size"){
      kwargs["treatment_effect"] <- param
    }else if(varying_param=="intervention group variance"){
      kwargs["intervention_within_neuron_sd"] <- sqrt(param * kwargs$control_within_neuron_sd**2)
    }else{
      stop("Argument 'mode' must be one of interneuron variance, intraneuron variance, or effect size.")
    }

    
    #Calculate the FDR, or Power, or both, calling those functions with our named list of arguments
    if(!fpr_only) {pow <- do.call(get_power, kwargs)}
    if(!power_only){fpr <- do.call(get_fpr,  kwargs)}

    if(!power_only){
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("KS Test", "False Positive Rate", fpr$ks)
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("Mixed Effects", "False Positive Rate", fpr$lmer)
    }
    
    if(!fpr_only){
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("KS Test", "Power", pow$ks)
      record[nrow(record) + 1,c("test", "test_stat", "test_stat_value")] <- c("Mixed Effects", "Power", pow$lmer)
    }
  }
  #Add a column to our record indicating the parameter values used. This dataframe is in long format, hence the repeats.
  record$parameter_value <- rep(param_values, each=if(power_only||fpr_only) 2 else 4)

  return(record)
}



#' Plot test performances against a varying data parameter
#'
#' Create a figure showing how the performance (FPR and statistical power) of the tests (KS tests and a mixed model comparison test) vary as 
#' a single parameter used to generate the dataset takes on different values. 
#' 
#' @inheritParams varying_param 
#' @inheritDotParams varying_param
#' 
#' @inherit varying_params details
#' 
#' @return ggplot object

plot_test_performances_against_varying_data_parameter <- function(varying_param = "interneuron variance",...){
  
  record <- check_test_performances_while_varying_data_parameter(varying_param = varying_param, ...)
  

  kwargs <- data.frame(...)
  
  record$test_stat_value <- as.numeric(record$test_stat_value)
  base <- ggplot(data = record, aes(x = parameter_value, y = as.numeric(test_stat_value), color = test))
  base <- base + scale_x_continuous(labels = scales::percent)
  
  if(varying_param=="interneuron variance"){
    xlabel<- "Within-group, between-neuron variance (%ATE)"
    base <- base + scale_y_continuous(labels = scales::percent)
    record$test_stat_value <- record$test_stat_value / kwargs$treatment_effect
  }else if(varying_param=="intraneuron variance"){
    xlabel<- "Within-neuron variance (%total variance)"
    base <- base + scale_y_continuous(labels = scales::percent)
    record$test_stat_value <- record$test_stat_value / kwargs$treatment_effect
  }else if(varying_param=="effect size"){
    xlabel<- "Effect size (Cohen's D)"
    record$test_stat_value <- record$test_stat_value / sqrt(kwargs$control_group_interneuron_sd**2 + kwargs$within_neuron_sd**2)
    
  }else if(varying_param=="intervention group variance"){
    xlabel <- "Variance between intervention group neurons (%Variance between controls)"
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

# effect_size <- plot_test_performances_against_varying_data_parameter("effect size", n_samples = 50, reserve_plots = F,
#                                                                      power_only = F, param_values = seq(from=0,to=1,by=0.2),
#                                                                      n_neurons = 500, samples_per_neuron = 10,
#                                                                      control_group_interneuron_sd = 0.711,
#                                                                      intervention_group_interneuron_sd = 0.711,
#                                                                      within_neuron_sd = 0.711)
# print(effect_size)

# r <- plot_test_performances_against_varying_data_parameter("intraneuron variance", param_values = seq(from=0.05, to=0.95, by = 0.1), n_samples = 100, intervention_group_interneuron_sd = sqrt(0.05), control_group_interneuron_sd= sqrt(0.05),
#                                                       n_neurons = 2, samples_per_neuron = 60, treatment_effect = 1, suppress = TRUE)

# p <- plot_test_performances_against_varying_data_parameter("intervention group variance", param_values = seq(0.8,1.2,by=0.1), treatment_effect = 1, control_within_neuron_sd = 1,
#                                                            n_neurons = 200, samples_per_neuron = 100, suppress = TRUE)

# beepr::beep()
