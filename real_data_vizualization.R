source("data_generation.R")
require(lme4)



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


amp_model1 <- lme4::lmer(ln ~ gt + (1|cellid), data = pvi_amp)
amp_model2 <- lme4::lmer(ln ~ gt + (1|cellid), data = dgg_amp)



print(VarCorr(amp_model1))
