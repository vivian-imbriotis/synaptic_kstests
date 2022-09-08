source("data_generation.R")
require(gridExtra)
require(lme4)


get_real_data <- function(){
  dgg_amp <- read.csv("real_data/dggcipscamp.csv")
  dgg_fre <-read.csv("real_data/dggcipscfreq.csv")
  pvi_amp <- read.csv("real_data/pvinipscamp.csv")
  pvi_fre <-read.csv("real_data/pvinipscfreq.csv")
  
  dgg_amp$area <- "Dentate gyrus granule cells"
  dgg_fre$area <- "Dentate gyrus granule cells"
  
  
  pvi_amp$area <- "Parvalbumin+ interneurons"
  pvi_fre$area <- "Parvalbumin+ interneurons"
  
  dgg_amp$dat_type <- "Amplitude"
  pvi_amp$dat_type <- "Amplitude"
  
  dgg_fre$dat_type <- "Frequency"
  pvi_fre$dat_type <- "Frequency"
  
  dat <- rbind(dgg_amp,dgg_fre,pvi_amp,pvi_fre)
  
  return(list(dgg_amp, dgg_fre, pvi_amp, pvi_fre, dat))
  
}

plot_real_data <- function(data, marginalize = T){
  plt <- ggplot(data=dat, aes(x=ln, y=after_stat(density), color=gt)) + facet_grid(cols = vars(dat_type), rows = vars(area))
  
  if(!marginalize){
    for(id in unique(dat$cellid)){
      cell <- subset(dat, cellid==id)
      plt <- plt + geom_freqpoly(data = cell, bins = 25)
    }
  }else{
    plt <- plt + geom_freqpoly(data = dat, bins = 50)
  }
  return(plt)
}

get_var_decomp_from_model <- function(model){
  
  variances <- VarCorr(model)
  variances <- as.data.frame(variances)
  
  between_cell_var <- variances$vcov[[1]]
  residual_var     <- variances$vcov[[2]]
  total_var        <- between_cell_var + residual_var
  
  return(list(Var_between_neurons = between_cell_var,
              Var_residual = residual_var,
              Var_proportion = between_cell_var / total_var))
}


get_models_for_real_data <- function(){
  data <- get_real_data()
  amp_model1 <- lme4::lmer(ln ~ gt + (1|cellid), data = pvi_amp)
  amp_model2 <- lme4::lmer(ln ~ gt + (1|cellid), data = dgg_amp)
  
  freq_model1 <- lme4::lmer(ln~gt + (1|cellid), data = pvi_fre)
  freq_model2 <- lme4::lmer(ln~gt + (1|cellid), data = dgg_fre)
}

generate_variance_decomp_table <- function(){
  amp_model1_vars <- get_var_decomp_from_model(amp_model1)
  amp_model1_vars$area <- "Parvalbumin+ interneurons"
  
  amp_model2_vars <- get_var_decomp_from_model(amp_model2)
  amp_model2_vars$area <- "Dentate Gyrus Granule Cells"
  
  amp_model1_vars$dat_type <- "ln(Amplitude)"
  amp_model2_vars$dat_type <- "ln(Amplitude)"
  
  freq_model1_vars <- get_var_decomp_from_model(freq_model1)
  freq_model1_vars$area <- "Parvalbumin+ interneurons"
  
  freq_model2_vars <- get_var_decomp_from_model(freq_model2)
  freq_model2_vars$area <- "Dentate Gyrus Granule Cells"
  
  freq_model1_vars$dat_type <- "ln(Frequency)"
  freq_model2_vars$dat_type <- "ln(Frequency)"
  
  variance_decomp <- rbind(amp_model1_vars, amp_model2_vars, freq_model1_vars, freq_model2_vars)
  rownames(variance_decomp) <- 1:4
  variance_decomp <- variance_decomp[,c(4,5,1,2,3)]
  
  return (variance_decomp)
}




