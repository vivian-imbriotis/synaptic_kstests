source("data_generation.R")
library(ggplot2)
library(patchwork)


check_test_performances_while_varying_data_parameter <- function(varying_param = "interneuron variance", 
                                                                 param_values = seq(from=0,to=1,by=0.1),
                                                                 alpha = 0.05,
                                                                 samples_per_param_value = 50,
                                                                 n_neurons = 20,
                                                                 treatment_effect = 1){
  #Make a dataframe to record the false positive and negative rates for each test, for each value of the data param
  record = data.frame(parameter_values = param_values, 
                      ks_false_positive_rate = 0, 
                      lmer_false_positive_rate = 0,
                      ks_false_negative_rate = 0,
                      lmer_false_negative_rate = 0)
  
  #Half our neurons are controls, the other half are intervention
  n_control_neurons = n_neurons%/%2 + n_neurons%%2
  n_intervention_neurons = n_neurons%/%2

  #For each value of our data parameter, we'll generate a bunch of datasets and test each one
  for (parameter_idx in 1:length(param_values)){
    
    param = param_values[[parameter_idx]]
    
    #This stores the test results for the datasets generated with a particular value of the data parameter
    tmp <- data.frame(ground_truth = TRUE, ks_positive = TRUE, lmer_positive = TRUE)
    
    for (i in 1:samples_per_param_value){
      #A treatment effect is present for the first N/2 sim runs, then absent for the remaining runs
      ground_truth = i<(samples_per_param_value / 2)
      
      #The effect size is constant (if present) unless its the varying parameter...
      if (varying_param == "effect_size"){
        effect<-as.integer(ground_truth) * param
      }else{
        effect <- as.integer(ground_truth) * treatment_effect
      }
      
      #Generate a mock dataset
      if (varying_param == "interneuron variance"){
        dat <- gen_unpaired_data(control_group_interneuron_sd = sqrt(param), intervention_group_interneuron_sd = sqrt(param),
                               treatment_effect = effect, n_control_neurons = n_control_neurons, n_intervention_neurons = n_intervention_neurons)
      }else if(varying_param=="intraneuron variance"){
        dat <- gen_unpaired_data(within_neuron_sd = sqrt(param),treatment_effect = effect, n_control_neurons = n_control_neurons, 
                                 n_intervention_neurons = n_intervention_neurons)
      }else if(varying_param=="effect size"){
        dat <- gen_unpaired_data(treatment_effect = effect, n_control_neurons = n_control_neurons, 
                                 n_intervention_neurons = n_intervention_neurons)
      }else{
        stop("Argument 'mode' must be one of interneuron, intraneuron, or effect size.")
      }
      
      #Calculate the KS test and drop1 mixed model test
      res <- ks.test(dat$dependant[dat$group=="Control"],dat$dependant[dat$group=="Intervention"])
      ks_positive <- res$p.value < alpha
      
      model <- lmerTest::lmer(dependant ~ group + (1|neuron_id), data = dat)
      res <- drop1(model)
      lmer_positive <- res$`Pr(>F)` < alpha
      
      tmp[i,"ground_truth"] <- ground_truth
      tmp[i,"ks_positive"] <- ks_positive
      tmp[i,"lmer_positive"] <- lmer_positive
    }
    
    #FPR = false positives / total negatives
    record[parameter_idx,"ks_false_positive_rate"]   = sum((!tmp$ground_truth) & tmp$ks_positive)   / sum(!tmp$ground_truth)
    record[parameter_idx,"lmer_false_positive_rate"] = sum((!tmp$ground_truth) & tmp$lmer_positive) / sum(!tmp$ground_truth)
  
    #FNR = false negatives / total positives
    record[parameter_idx,"ks_false_negative_rate"]   = sum((tmp$ground_truth) & !tmp$ks_positive)   / sum(tmp$ground_truth)
    record[parameter_idx,"lmer_false_negative_rate"] = sum((tmp$ground_truth) & !tmp$lmer_positive) / sum(tmp$ground_truth)
  }
  
  return(record)
}



plot_test_performances_against_varying_data_parameter <- function(varying_param, samples_per_param_value, 
                                                                  param_values = seq(from=0,to=1,by=0.1)){
  
  record <- check_test_performances_while_varying_data_parameter(varying_param = varying_param, param_values = param_values,
                                                              samples_per_param_value = samples_per_param_value)
  
  base <- ggplot(data = record, aes(x = parameter_values)) + scale_x_continuous(labels = scales::percent)
  
  if(varying_param=="interneuron variance"){
    xlabel<- "Within-group, between-neuron variance (%ATE)"
  }else if(varying_param=="intraneuron variance"){
    xlabel<- "Within-neuron variance (%ATE)"
  }else if(varying_param=="effect size"){
    xlabel<- "Average treatment effect (%between-neuron variance)"
  }else{
    stop("Argument 'varying_param' bust be one of interneuron, intraneuron, or effect_size")
  }
  
  labels <- labs(x = xlabel, y = "Rate", color = "Test")
  legend_scale <- scale_color_discrete(limits = c("blue","orange"), labels = c("KS Test", "Mixed Model"))
  false_negatives <- base + geom_line(aes(y = ks_false_negative_rate, color = "blue")) + geom_line(aes(y = lmer_false_negative_rate, color="orange"))
  false_negatives <- false_negatives + scale_y_continuous(limits = c(0,1)) + labels + labs(y = "False Negative Rate")
  false_negatives <- false_negatives + legend_scale
  
  false_positives <- base + geom_line(aes(y = ks_false_positive_rate, color = "blue")) + geom_line(aes(y = lmer_false_positive_rate, color="orange"))
  false_positives <- false_positives + scale_y_continuous(limits = c(0,1)) + labels + labs(y = "False Positive Rate")
  false_positives <- false_positives + theme(legend.position = "None")
  plt = false_positives + false_negatives
  return(plt) 
}

interneuron <- plot_test_performances_against_varying_data_parameter("interneuron variance", samples_per_param_value = 100)
intraneuron <- plot_test_performances_against_varying_data_parameter("intraneuron variance", samples_per_param_value = 100)
effect_size <- plot_test_performances_against_varying_data_parameter("effect size",          samples_per_param_value = 100)