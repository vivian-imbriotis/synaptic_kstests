source("data_generation.R")

TOTAL_VAR <- 1
TOTAL_SS <- 120
BOOTSTRAP_SAMPLES <- 40
GRID_FREQ <- 5
PAIRED <- FALSE
PROP_RESIDUAL_VAR_BETWEEN_PAIRS <- 0.5

record <- data.frame(prop_res_var=numeric(), n_neurons_on_ss=numeric(), test=character(), test_stat=character(), test_stat_value=numeric())

proportion_residual_variance <- seq(from=0.05,to=0.95,length=GRID_FREQ)

min_neurons <- 2
max_neurons <- (TOTAL_SS %/% 2)
by <- (max_neurons - min_neurons) %/% GRID_FREQ
n_neurons <- seq(from=min_neurons, to = max_neurons, by = by)
if (any(n_neurons != as.integer(n_neurons))) {stop("Noninteger neurons requested")}

for (prop_res_var in proportion_residual_variance){
  cat("setting %within-neuron variance to ", prop_res_var, "\n")
  for (Nn in n_neurons){
    residual_var <- TOTAL_VAR * prop_res_var
    between_neuron_var <- TOTAL_VAR - residual_var
    
    events_per_neuron <- TOTAL_SS %/% Nn
    if(Nn < 2){
      stop("Must be >2 neurons")
    }
    if (events_per_neuron < 1){
      stop("Must be at least 1 event per neuron")
    }
    args <- list()
    args$n_samples <- BOOTSTRAP_SAMPLES
    args$samples_per_neuron <- events_per_neuron
    args$suppress <- T
    args$paired <- PAIRED
    if(PAIRED){
      args$n_neurons <- Nn
      
      args$residual_interneuron_sd <- sqrt(residual_var * PROP_RESIDUAL_VAR_BETWEEN_PAIRS)
      args$within_neuron_sd <- sqrt(residual_var * (1 - PROP_RESIDUAL_VAR_BETWEEN_PAIRS))
      args$interneuron_sd <- sqrt(between_neuron_var)
    }else{    
      args$n_control_neurons <- Nn
      args$n_intervention_neurons <- Nn
      
      args$control_group_interneuron_sd <- sqrt(between_neuron_var)
      args$intervention_group_interneuron_sd <- sqrt(between_neuron_var)
      args$within_neuron_sd <- sqrt(residual_var)
    }
    
    args$saveplot <- T
    args$filename <- paste("neurons=",Nn,"prop_res_var=",prop_res_var,"FALSE")
    fpr <- do.call(get_fpr, args)
    args$filename <- paste("neurons=",Nn,"prop_res_var=",prop_res_var, "TRUE")
    pow <- do.call(get_power,args)
    
    record[nrow(record) + 1,] <- list(prop_res_var, Nn / TOTAL_SS, "KS Test", "False Positive Rate", fpr["lower_CI","ks"])
    record[nrow(record) + 1,] <- list(prop_res_var, Nn / TOTAL_SS, "KS Test", "Power", pow["upper_CI","ks"])
    record[nrow(record) + 1,] <- list(prop_res_var, Nn / TOTAL_SS, "Mixed Model", "False Positive Rate",  fpr["lower_CI","lmer"])
    record[nrow(record) + 1,] <- list(prop_res_var, Nn / TOTAL_SS, "Mixed Model", "Power", pow["upper_CI","lmer"])
    
    
  }
}

clip_vector <- function(x, a, b) {
  ifelse(x <= a,  a, ifelse(x >= b, b, x))
}


library(ggplot2)
require(interp)

# res <- 10
# dfs <- list()
# for(test in unique(record$test)){
#   for(test_stat in unique(record$test_stat)){
#     grid <- with(record[record$test==test & record$test_stat==test_stat,], interp::interp(n_neurons_on_ss, prop_explainable_var, test_stat_value, nx=res,ny=res,
#                                                                                           method = "akima", extrap = T))
# 
#     griddf <- subset(data.frame(n_neurons_on_ss = rep(grid$x, nrow(grid$z)),
#                                 prop_explainable_var = rep(grid$y, each = ncol(grid$z)),
#                                 test_stat_value = as.numeric(grid$z)))
#     griddf$test <- test
#     griddf$test_stat <- test_stat
#     dfs[[length(dfs) + 1]] <- griddf
#   }
# }
# 
# df <- do.call(rbind,dfs)
record$prop_explainable_var <- 1 - record$prop_res_var
# df$test_stat_value <- clip_vector(df$test_stat_value, 0, 1)

df <- record



#OLD PLOT
df$acceptable <- (df$test_stat=="False Positive Rate" & df$test_stat_value <= 0.1) | (df$test_stat=="Power" & df$test_stat_value >= 0.8)
oldplt <- ggplot(data = df, aes(x=n_neurons_on_ss, y=prop_explainable_var)) + geom_raster(aes(fill=acceptable), interpolate = F) + facet_grid(rows = vars(test), cols = vars(test_stat))
labels <- labs(x = "Number of neurons / sample size", y = "Between-neuron variance (% total variance)")
oldplt <- oldplt + labels + scale_fill_discrete(name = "Acceptable")

#NEW PLOT
plt <- ggplot(data = df, aes(x=n_neurons_on_ss, y=prop_explainable_var))
plt <- plt + geom_contour_filled(aes(z=test_stat_value), breaks = c(0,0.05,0.2,0.4,0.6,0.8,Inf)) + facet_grid(rows = vars(test), cols = vars(test_stat))
plt <- plt + scale_fill_viridis_d(drop=F, name = "Rate", labels = c(expression(""<=0.05), expression(""<=0.2), expression(""<=0.4), expression(""<=0.6), expression(""<=0.8), expression(""<=1)))
plt <- plt + labels

#RASTER VERSION
mybreaks <- c(0,0.05,0.2,0.4,0.6,0.8,Inf)
df$cuts <- cut(df$test_stat_value, mybreaks)
raster <- ggplot(data = df, aes(x=n_neurons_on_ss, y=prop_explainable_var))
raster <- raster + geom_raster(aes(fill=cuts), breaks = mybreaks) + facet_grid(rows = vars(test), cols = vars(test_stat))
raster <- raster + scale_fill_viridis_d(drop=F, name = "Rate", labels = c(expression(""<=0.05), expression(""<=0.2), expression(""<=0.4), expression(""<=0.6), expression(""<=0.8), expression(""<=1)))
raster <- raster + labels

print(raster)
