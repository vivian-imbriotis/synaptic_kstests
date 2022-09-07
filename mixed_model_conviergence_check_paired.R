source("data_generation.R")
source("figures_linegraphs.R")
require(lme4)
require(beepr)


N = 25
ATE = 0

#Define the variances we use to generate the data
BETWEEN_NEURONS = 1.5
BETWEEN_PAIRS= 2
WITHIN = 0.8

TOTAL_VAR = BETWEEN_NEURONS + BETWEEN_PAIRS + WITHIN

#The models don't fit these variances exactly, instead they fit these ones...
INTERCEPT_VAR = BETWEEN_NEURONS + BETWEEN_PAIRS
TE_VAR = 2 * BETWEEN_PAIRS
E_COVAR = -1 * BETWEEN_PAIRS


expected_residual_var = (BETWEEN_PAIRS**2 + WITHIN**2)

expected_overall_var = (BETWEEN_NEURONS**2 + BETWEEN_PAIRS**2 + WITHIN**2)

#For each sample, record the statistics the random slope and random intercept models extract
#Also record the sample standard deviation in the control and intervention groups, so 
#we can check that our data generating process is right (i.e. using the law of total variance)
model_record = data.frame(control_group_mean_estimate = numeric(N),
                          ATE_estimate = numeric(N), 
                          between_neuron_var = numeric(N), 
                          residual_var = numeric(N), 
                          pval = numeric(N))

slope_model_record = data.frame(control_group_mean_estimate = numeric(N),
                                ATE_estimate = numeric(N),
                                between_control_neuron_means = numeric(N),
                                between_treatment_effects = numeric(N),
                                residual_var = numeric(N),
                                covar = numeric(N),
                                pval = numeric(N))



var_record <- data.frame(control_grp = numeric(N), 
                         interve_grp = numeric(N),
                         intercepts  = numeric(N),
                         slopes      = numeric(N),
                         cov         = numeric(N))

slope_model_record$run <- model_record$run <- rep(1:N)

for (i in 1:N){
  
  dat <- gen_paired_data(treatment_effect = ATE, interneuron_sd = sqrt(BETWEEN_NEURONS), residual_interneuron_sd = sqrt(BETWEEN_PAIRS),
                           within_neuron_sd = sqrt(WITHIN), n_neurons = 40, samples_per_neuron = 20000%/%40)
  
  #Calculate the sample variances
  cntrl <- subset(dat, group=="Control")
  inter <- subset(dat, group=="Intervention")
  
  overall_var_control <- var(cntrl$dependant)
  overall_var_interv  <- var(inter$dependant)
  
  neuron_intercepts <- tapply(cntrl$dependant, cntrl$neuron_id, mean)
  overall_var_intercepts <- var(neuron_intercepts)
  
  neuron_treatment_effects <- tapply(inter$dependant, inter$neuron_id, mean) - neuron_intercepts
  
  covar <- cov(neuron_intercepts, neuron_treatment_effects)
  
  overall_var_slopes <- var(neuron_treatment_effects)
  var_record[i,] <- c(overall_var_control, overall_var_interv, overall_var_intercepts, overall_var_slopes, covar)
  
  #Fit a model with random intercepts for each neuron
  intercept_model <- lme4::lmer(dependant ~ group + (1|neuron_id), data = dat, REML = F)
  intercept_model_h0<- lme4::lmer(dependant ~ (1|neuron_id), data = dat, REML = F)
  a_intercept <- anova(intercept_model,intercept_model_h0)
  s_intercept <- summary(intercept_model)
  
  #Fit a model with random (correlated) intercepts and slopes for each neuron
  slope_model   <- lme4::lmer(dependant ~ group + (1+group|neuron_id), data = dat, REML = F)
  slope_model_h0<- lme4::lmer(dependant ~ (1+group|neuron_id), data = dat, REML = F)
  a_slope <- anova(slope_model, slope_model_h0)
  s_slope <- summary(slope_model)

  #Get fixed effect estimates
  intercept_fixed  <- s_intercept$coefficients
  #Get the random effect variance and the residual variance
  intercept_random <- as.data.frame(s_intercept$varcor)$vcov
  
  #Get fixed effects estimates
  slope_fixed <- s_slope$coefficients
  #Get the variance-covariance matrix and the residual variance
  slope_random <- as.data.frame(s_slope$varcor)
  
  
  slope_residual_var <- subset(slope_random, grp=="Residual")$vcov
  slope_between_neuron_var <- subset(slope_random,var1=="(Intercept)" & is.na(var2))$vcov
  slope_between_pairs_var <- subset(slope_random, var1=="groupIntervention" & is.na(var2))$vcov
  slope_covariance <- subset(slope_random, var1=="(Intercept)" & var2=="groupIntervention")$vcov
  
  
  model_record$control_group_mean_estimate[[i]] <- intercept_fixed[[1,1]]
  model_record$ATE_estimate[[i]] <- intercept_fixed[[2,1]]
  model_record$between_neuron_var[[i]] <- intercept_random[[1]]
  model_record$residual_var[[i]] <- intercept_random[[2]]
  model_record$pval[[i]] <- a_intercept$`Pr(>Chisq)`[[2]]
  
  
  slope_model_record$control_group_mean_estimate[[i]] <- slope_fixed[['(Intercept)','Estimate']]
  slope_model_record$ATE_estimate[[i]]                <- slope_fixed[['groupIntervention','Estimate']]
  slope_model_record$between_control_neuron_means[[i]] <- slope_between_neuron_var
  slope_model_record$between_treatment_effects[[i]] <- slope_between_pairs_var
  slope_model_record$residual_var[[i]]  <- slope_residual_var
  slope_model_record$covar[[i]] <- slope_covariance
  
  slope_model_record$pval[[i]] <- a_slope$`Pr(>Chisq)`[[2]]
  
}

long_inter <- tidyr::gather(model_record,statistic,value,control_group_mean_estimate,ATE_estimate,between_neuron_var,residual_var, pval)
long_slope <- tidyr::gather(slope_model_record,statistic,value,control_group_mean_estimate,ATE_estimate,between_control_neuron_means,
                            between_treatment_effects,
                            residual_var,
                            covar,
                            pval)
long_inter$model_type <- "Random Intercept"
long_slope$model_type <- "Random Slope"

long <- rbind(long_inter,long_slope)

long <- subset(long, statistic!="pval")

# long$statistic <- factor(long$statistic, c("ATE_estimate", "control_group_mean_estimate", "between_neuron_var", "pair_var", "covar", "residual_var"))

ground_truths <- data.frame(statistic = character(7), value = numeric(7))
ground_truths[,"statistic"]<-unique(long$statistic)
values <- c(0, ATE, BETWEEN_NEURONS, WITHIN, INTERCEPT_VAR, TE_VAR, E_COVAR)
ground_truths[,"value"] <- values

plt <- ggplot(data = long, aes(x = statistic, y=value)) + facet_grid(rows = vars(model_type))

plt <- plt + geom_violin(scale = "width")


plt <- plt + geom_point(data = ground_truths, aes(x = statistic, y = value), color = "red", shape = "square")

plt <- plt + labs(x="Parameter", y="Value") 

# plt <- plt + scale_x_discrete(labels = c("ATE",bquote(mu[Control]),
#                                                                           bquote(Var[Between]), bquote(Var[TE]),
#                                                                           bquote(cov[list(Between,TE)]), bquote(Var[Residual])))

plt <- plt + theme(axis.text = element_text(size = 10))

fpr <- data.frame(Model = character(2), FPR_upper_CI = numeric(2), FPR_lower_CI = numeric(2))

fpr$Model <- c("Intercept", "Slope")
fpr[1,2:3] <- binom.test(sum(model_record$pval < 0.05), N)$conf.int
fpr[2,2:3] <- binom.test(sum(slope_model_record$pval < 0.05),     N)$conf.int

require(patchwork)
require(gridExtra)

tablegrob <- gridExtra::tableGrob(head(fpr), cols = c("Model", "FPR (95%CI lower)", "FPR (95% CI upper)"), rows = c("",""))

figure <- plt / tablegrob + plot_layout(heights = c(3,1))

print(figure)

vplt <- ggplot(data = tidyr::gather(var_record, statistic, value, intercepts, slopes, cov), aes(x=statistic, y=value)) + geom_violin()


slope_randcoefs <- coef(slope_model)$neuron_id
slope_randcoefs$neuron_id <- as.factor(1:nrow(slope_randcoefs))
slopes_dat <- tidyr::gather(slope_randcoefs, coef, value, `(Intercept)`, groupIntervention)

slopeplt <- ggplot(slopes_dat, aes(x=coef, y=value)) + geom_line()

beepr::beep()