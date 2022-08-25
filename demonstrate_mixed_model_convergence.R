source("data_generation.R")
source("unpaired_FPR_vs_interneuron.R")
require(lme4)
require(beepr)

N = 250
ATE = 0
INTER = 1.5
WITHIN = 0.8

expected_overall_sd = sqrt(INTER^2 + WITHIN^2)
record = data.frame(overall_sd_control = numeric(N), overall_sd_interv = numeric(N), control_group_mean_estimate = numeric(N),
                    ATE_estimate = numeric(N), between_neuron_sd = numeric(N), residual_sd = numeric(N), pval = numeric(N))

record$run <- 1:N

for (i in 1:N){
  
  dat <- gen_unpaired_data_fast(treatment_effect = ATE, control_group_interneuron_sd = INTER, intervention_group_interneuron_sd = INTER, 
                           within_neuron_sd = WITHIN, n_control_neurons = 20, n_intervention_neurons = 20)
  
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

long <- tidyr::gather(record,statistic,value,overall_sd_control,overall_sd_interv,control_group_mean_estimate,ATE_estimate,between_neuron_sd,residual_sd)
long$is_model <- long$statistic %in% colnames(record)[3:7]

ground_truths <- data.frame(statistic = character(6), value = 0)
ground_truths[,"statistic"]<-colnames(record)[1:6]
values <- c(expected_overall_sd,expected_overall_sd,0,ATE,INTER,WITHIN)
ground_truths[,"value"] <- values
ground_truths[,"is_model"] <- ground_truths$statistic %in% colnames(record)[3:6]



plt <- ggplot(data = long, aes(x = statistic, y=value)) + geom_violin(aes(color = is_model))


plt <- plt + geom_point(data = ground_truths, aes(x = statistic, y = value))


print(plt)

fpr = sum(record$pval < 0.01) / length(record$pval)
print(fpr)

beepr::beep()