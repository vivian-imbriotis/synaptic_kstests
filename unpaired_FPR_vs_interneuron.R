source("data_generation.R")

ALPHA <- 0.05
interneuron_variances <- seq(from=0,to=1,by=0.1)
n_neurons = 50

alpha = ALPHA
samples_per_variance_setting <- 100
record = data.frame(interneuron_variance = interneuron_variances, ks_false_positive_rate = 0, lmer_false_positive_rate = 0)
n_control_neurons = n_neurons%/%2 + n_neurons%%2
c_intervention_neurons = n_neurons%/%2

for (variance_idx in 1:length(interneuron_variances)){
  
  interneuron_var = interneuron_variances[[variance_idx]]
  
  tmp <- data.frame(ground_truth = TRUE, ks_positive = TRUE, lmer_positive = TRUE)
  
  for (i in 1:samples_per_variance_setting){
    ground_truth = i<(samples_per_variance_setting / 2)
    
    effect <- as.integer(ground_truth)
    
    dat <- gen_unpaired_data(control_group_interneuron_sd = sqrt(interneuron_var), intervention_group_interneuron_sd = sqrt(interneuron_var),
                             treatment_effect = effect, n_control_neurons = n_control_neurons, n_intervention_neurons = n_intervention_neurons)
    
    res <- ks.test(dat$dependant[dat$neuron_group=="Control"],dat$dependant[dat$neuron_group=="Intervention"])
    ks_positive <- res$p.value < alpha
    
    model <- lmerTest::lmer(dependant ~ neuron_group + (1|neuron_id), data = dat)
    res <- drop1(model)
    lmer_positive <- res$`Pr(>F)` < alpha
    
    tmp[i,"ground_truth"] <- ground_truth
    tmp[i,"ks_positive"] <- ks_positive
    tmp[i,"lmer_positive"] <- lmer_positive
  }
  
  #FPR = false positives / total negatives
  record[variance_idx,"ks_false_positive_rate"]   = sum((!tmp$ground_truth) & tmp$ks_positive)   / sum(!tmp$ground_truth)
  record[variance_idx,"lmer_false_positive_rate"] = sum((!tmp$ground_truth) & tmp$lmer_positive) / sum(!tmp$ground_truth)

  #FNR = false negatives / total positives
  record[variance_idx,"ks_false_negative_rate"]   = sum((tmp$ground_truth) & !tmp$ks_positive)   / sum(tmp$ground_truth)
  record[variance_idx,"lmer_false_negative_rate"] = sum((tmp$ground_truth) & !tmp$lmer_positive) / sum(tmp$ground_truth)
}

#Make the plot...

library(ggplot2)
library(patchwork)

base <- ggplot(data = record, aes(x = interneuron_variance)) + scale_x_continuous(labels = scales::percent)
labels <- labs(x = "Within-group, between-neuron variance (%ATE)", y = "Rate", color = "Test")
legend_scale <- scale_color_discrete(limits = c("blue","orange"), labels = c("KS Test", "Mixed Model"))
false_negatives <- base + geom_line(aes(y = ks_false_negative_rate, color = "blue")) + geom_line(aes(y = lmer_false_negative_rate, color="orange"))
false_negatives <- false_negatives + scale_y_continuous(limits = c(0,1)) + labels + labs(y = "False Negative Rate")
false_negatives <- false_negatives + legend_scale

false_positives <- base + geom_line(aes(y = ks_false_positive_rate, color = "blue")) + geom_line(aes(y = lmer_false_positive_rate, color="orange"))
false_positives <- false_positives + scale_y_continuous(limits = c(0,1)) + labels + labs(y = "False Positive Rate")
false_positives <- false_positives + theme(legend.position = "None")
plt = false_positives + false_negatives
print(plt)
