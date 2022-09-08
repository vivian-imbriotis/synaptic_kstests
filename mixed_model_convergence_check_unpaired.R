source("data_generation.R")
require(lme4)
require(beepr)
library(ggplot2)

unpaired_convirgence_check_figure <- function(){
  
  N_BOOTSTRAP = 250
  ATE = 0
  VAR_BETWEEN_NEURONS = 1.5
  VAR_WITHIN_NEURONS = 0.8
  
  record = data.frame(overall_sd_control = numeric(N_BOOTSTRAP), overall_sd_interv = numeric(N_BOOTSTRAP), control_group_mean_estimate = numeric(N_BOOTSTRAP),
                      ATE_estimate = numeric(N_BOOTSTRAP), between_neuron_sd = numeric(N_BOOTSTRAP), residual_sd = numeric(N_BOOTSTRAP), pval = numeric(N_BOOTSTRAP))
  
  record$run <- 1:N_BOOTSTRAP
  
  for (i in 1:N_BOOTSTRAP){
    
    dat <- gen_unpaired_data(treatment_effect = ATE, 
                             control_between_neuron_sd = VAR_BETWEEN_NEURONS, 
                             intervention_between_neuron_sd = VAR_BETWEEN_NEURONS, 
                             control_within_neuron_sd = VAR_WITHIN_NEURONS, 
                             intervention_within_neuron_sd = VAR_WITHIN_NEURONS,
                             n_control_neurons = 20, 
                             n_intervention_neurons = 20)
    
    overall_sd_control <- sd(dat[dat$group=="Control"     ,"dependant"])
    overall_sd_interv  <- sd(dat[dat$group=="Intervention","dependant"])
    
    
    model <- lme4::lmer(dependant ~ group + (1|neuron_id), data = dat, REML = F)
    model2<- lme4::lmer(dependant ~ (1|neuron_id), data = dat, REML = F)
    a <- anova(model,model2)
    
    
    s <- summary(model)
    fixed  <- s$coefficients
    random <- as.data.frame(s$varcor)$sdcor
    
    record$overall_sd_control[[i]] <- overall_sd_control
    record$overall_sd_interv[[i]]  <- overall_sd_interv
    record$control_group_mean_estimate[[i]] <- fixed[[1,1]]
    record$ATE_estimate[[i]] <- fixed[[2,1]]
    record$between_neuron_sd[[i]] <- random[[1]]
    record$residual_sd[[i]] <- random[[2]]
    record$pval[[i]] <- a$`Pr(>Chisq)`[[2]]
  
  }
  
  long <- tidyr::gather(record,statistic,value,control_group_mean_estimate,ATE_estimate,between_neuron_sd,residual_sd)
  long$is_model <- long$statistic %in% colnames(record)[3:7]
  
  long$statistic <- factor(long$statistic, levels = c("control_group_mean_estimate","ATE_estimate","between_neuron_sd", "residual_sd"))
  
  ground_truths <- data.frame(statistic = character(6), value = 0)
  ground_truths[,"statistic"]<-colnames(record)[1:6]
  values <- c(expected_overall_sd,expected_overall_sd,0,ATE,VAR_BETWEEN_NEURONS,VAR_WITHIN_NEURONS)
  ground_truths[,"value"] <- values
  ground_truths[,"is_model"] <- ground_truths$statistic %in% colnames(record)[3:6]
  
  
  plt <- ggplot(data = subset(long, is_model==T), aes(x = statistic, y=value)) + geom_violin(fill = "#56B4E9", color="black", scale = "width")
  
  
  plt <- plt + geom_point(data = subset(ground_truths, is_model==T), aes(x = statistic, y = value))
  
  plt <- plt + scale_x_discrete(labels = c(expression(Intercept),expression(ATE), expression(Var[Between~Neurons]), expression(Var[VAR_WITHIN_NEURONS~Neurons])))
  
  plt <- plt + labs(x="Parameter", y="Value")
  
  plt <- plt + theme(axis.text = element_text(size = 11))
  
  # ggsave("supplimentary1.pdf", width = 7, height=2, units="in", plt, device = "pdf")
  
  # fpr = sum(record$pval < 0.01) / length(record$pval)
  # print(fpr)
  return(plt)
}