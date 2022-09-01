source("data_generation.R")
library(ggplot2)
library(patchwork)

set.seed(0)

DEFAULT_NEURONS = 20
DEFAULT_RESIDUAL_VARIANCE = 0.34
DEFAULT_EXPLAINABLE_VARIANCE = 0.61
DEFAULT_EFFECT_SIZE = 0.5

#' Check test performances while varying two data parameters
#'
#' @param x_axis One of "number_of_neurons" or "effect_size"
#' @param y_axis One of 'proportion_explainable_variance' or 'intervention_variance_ratio'
#' @param x_axis_min The minimum value of the x-axis variate
#' @param x_axis_max The maximum value of the y-axis variate
#' @param y_axis_min  The minimum value of the x-axis variate
#' @param y_axis_max  The maximum value of the y=axis variate
#' @param total_variance The total amount of variance in the control group (In the unpaired case,
#' this is composed of between_neuron/explainable and within-neuron/unexplainable variance, in
#' the paired case it is composed of between-neuron/explainable, within-neuron, and between-pair variances)
#' @param total_sample_size The total sample size (number of neuron * samples per neuron) of each generated dataset
#' @param bootstrap_samples The number of datasets to generate for each (x,y) variate pair
#' @param grid_frequency How many different values for x and y to take on
#' @param paired Whether the dataset should be paired or unpaired
#' @param prop_residual_variance_between_pairs Only has an effect when paired=TRUE. How 
#' much of the residual (unexplained) variance should be between the paired means of a
#' neuron in the intervention vs control states; the remaining residual variance is
#' simply the variance between observations from a particular neuron.
#' @param default_neurons If number of neurons does not vary, how many neurons should be in each dataset?
#' @param default_residual_variance If the residual/explainable variance ratio does not vary, set the residual variance
#' @param default_explainable_variance If the residual/explainable variance ratio does not vary, set the explainable variance
#' @param default_treatment_effect If the effect size does not vary, set the treatment effect
#'
#' @return
#' @export
#'
#' @examples
check_test_performances_while_varying_two_data_parameters <- function(
                                 x_axis  = "number_of_neurons", 
                                 y_axis  = "proportion_explainable_variance",
                                 x_axis_min = 0,
                                 x_axis_max = 0.8,
                                 y_axis_min = 0.05,
                                 y_axis_max = 0.95,
                                 total_variance = 1,
                                 total_sample_size = 1000,
                                 bootstrap_samples = 200,
                                 grid_frequency = 10,
                                 paired = FALSE,
                                 prop_residual_variance_between_pairs = 0.5,
                                 number_of_neurons = DEFAULT_NEURONS,
                                 proportion_explainable_variance = DEFAULT_EXPLAINABLE_VARIANCE / (DEFAULT_EXPLAINABLE_VARIANCE+DEFAULT_RESIDUAL_VARIANCE),
                                 effect_size = DEFAULT_EFFECT_SIZE,
                                 intervention_variance_ratio = 1){
  record <- data.frame(x_axis=numeric(), y_axis=numeric(), test=character(), test_stat=character(), test_stat_value=numeric())
  
  #Start with some error checking. Both the x and y axis have to be from the list of data-generating parameters
  #that are allowed to vary. intervention_variance_ratio cannot currently vary in the paired case (within a neuron
  #the mean and variance are equal, regardless of if its in the intervention or control state, only the mean can change).
  #Similarly, the prop_residual_variance_between_pairs can't change in the unpaired case!
  
  if(!(x_axis %in% names(axislabels))){
    stop(paste("x_axis must be set to one of the following:", toString(names(axislabels))))
  }
  if(!(y_axis %in% names(axislabels))){
    stop(paste("y_axis must be set to one of the following:", toString(names(axislabels))))
  }
  if("prop_residual_variance_between_pairs"%in%c(x_axis,y_axis) && (!paired)){
    stop(paste("Cannot vary variance between pairs in unpaired data. Did you mean to set paired=TRUE?"))
  }
  if("intervention_variance_ratio"%in%c(x_axis,y_axis) && paired){
    stop(paste("Cannot vary intervention group variance in paired data. Did you mean to set paired=FALSE?"))
  }
  
  names(record)[names(record)=="x_axis"] <- x_axis
  names(record)[names(record)=="y_axis"] <- y_axis
  
  y_values <- seq(from=y_axis_min,to=y_axis_max,length=grid_frequency)
  x_values <- seq(from=x_axis_min,to=x_axis_max,length=grid_frequency)
  #Even through we just determined the xvalues above based on args, since number_of_neurons needs to be both an INTEGER and 
  #needs to be evenly spaced, we're going to do the best we can to match the user's request if that's on the xaxis
  if(x_axis=="number_of_neurons" || y_axis=="number_of_neurons"){
    if(x_axis=="number_of_neurons"){
      min_neurons <- max(2, x_axis_min)
      max_neurons <- min(x_axis_max, (total_sample_size %/% 2))
    }else{
      min_neurons <- max(2, y_axis_min)
      max_neurons <- min(y_axis_max, (total_sample_size %/% 2))
    }
    by <- (max_neurons - min_neurons) %/% grid_frequency
    vals <- seq(from=min_neurons, to = max_neurons, by = by)
    if (any(vals != as.integer(vals))) {stop("Noninteger neurons requested")}
    if(x_axis=="number_of_neurons") x_values<-vals else y_values<-vals
  }
  

  for (y in y_values){
    cat("setting", y_axis, "to", y, "while varying", x_axis, "\n")
    
    for (x in x_values){
      
      #Assign our varying parameters
      assign(y_axis, y)
      assign(x_axis, x)
      
      
      #Derive the variances in the dataset
      between_neuron_var <- total_variance * proportion_explainable_variance
      residual_var       <- total_variance  - between_neuron_var
      #If data are paired, residual variance is split up into unexplainable between_pair variance and unexplainable within_neuron variance
      if(paired){
        paired_var <- residual_var * prop_residual_variance_between_pairs
        within_neuron_var <- residual_var - paired_var
      }
      
      #Derive the number of events per neuron
      events_per_neuron <- total_sample_size %/% number_of_neurons
      if(number_of_neurons < 2){
        stop("Must be >2 neurons")
      }
      if (events_per_neuron < 1){
        stop("Must be at least 1 event per neuron")
      }
      
      #Derive the treatment effect from the effect size and the variances
      control_group_total_variance <- residual_var + between_neuron_var
      treatment_effect <- effect_size * sqrt(control_group_total_variance) #Inverting the formula for glass's delta

      #Build our list of arguments to our data-generating function
      args <- list()
      args$n_samples <- bootstrap_samples
      args$samples_per_neuron <- events_per_neuron
      args$suppress <- T
      args$paired <- paired
      args$treatment_effect <- treatment_effect
      

      if(paired){
        args$n_neurons <- number_of_neurons
        args$residual_interneuron_sd <- sqrt(paired_var)
        args$within_neuron_sd <- sqrt(within_neuron_var)
        args$interneuron_sd <- sqrt(between_neuron_var)
      }else{    #The unpaired case
        args$n_control_neurons      <- number_of_neurons
        args$n_intervention_neurons <- number_of_neurons
        
        args$control_group_interneuron_sd      <- sqrt(between_neuron_var)
        args$intervention_group_interneuron_sd <- sqrt(between_neuron_var)
        args$control_within_neuron_sd <- sqrt(residual_var)
        args$intervention_within_neuron_sd <- intervention_variance_ratio * args$control_within_neuron_sd
      }
      
      #Generate a lot of datasets with the given parameters!!
      # args$saveplot <- T
      # args$filename <- paste(x_axis,"=",x,y_axis,"=",y,"FALSE")
      fpr <- do.call(get_fpr, args)
      # args$filename <- paste(x_axis,"=",x,y_axis,"=",y, "TRUE")
      pow <- do.call(get_power,args)
      
      record[nrow(record) + 1,] <- list(x, y, "KS Test", "False Positive Rate", fpr["lower_CI","ks"])
      record[nrow(record) + 1,] <- list(x, y, "KS Test", "Power", pow["value","ks"])
      record[nrow(record) + 1,] <- list(x, y, "Mixed Model", "False Positive Rate",  fpr["lower_CI","lmer"])
      record[nrow(record) + 1,] <- list(x, y, "Mixed Model", "Power", pow["value","lmer"])
    }
  }
  return(record)
}



#' Plot Heatmaps
#'
#' @inheritDotParams check_test_performance_while_varying_two_data_parameters
#'
#' @return list of ggplot objects, with names 'acceptibility', 'contour', 'raster'

axislabels <- list()
axislabels$effect_size <- expression(Effect~size~("Glass's"~Delta))
axislabels$number_of_neurons <- expression(N[Neurons]~(Constant~sample~size))
axislabels$proportion_explainable_variance <- expression(frac(Var[Between~neurons],Var[Total]))
axislabels$intervention_variance_ratio <- expression(frac(Var[Intervention],Var[Control]))
axislabels$prop_residual_variance_between_pairs <- expression(frac(Var[pair],Var[residual]))

plot_heatmaps <- function(...){
  
  df <- check_test_performances_while_varying_two_data_parameters(...)
  
  x_axis <- colnames(df)[1]
  y_axis <- colnames(df)[2]
  
  xlabel <- axislabels[[x_axis]]
  ylabel <- axislabels[[y_axis]]
  

  labels <- labs(x = xlabel, y = ylabel)
  breaks <- c(-Inf,0.05,0.2,0.4,0.6,0.8,Inf)
  facets <- facet_grid(rows = vars(test), cols = vars(test_stat))
  colors <- scale_fill_viridis_d(drop=F, name = "Rate", 
                                 labels = c(expression(""<=0.05), expression(""<=0.2), expression(""<=0.4), 
                                            expression(""<=0.6), expression(""<=0.8), expression(""<=1))
                                 )
  
  #Acceptability plot
  df$acceptable <- (df$test_stat=="False Positive Rate" & df$test_stat_value <= 0.05) | (df$test_stat=="Power" & df$test_stat_value >= 0.8)
  accept <- ggplot(data = df, aes_string(x=x_axis, y=y_axis)) + geom_raster(aes(fill=acceptable), interpolate = F) + facets
  accept <- accept + labels + scale_fill_discrete(name = "Acceptable")
  
  #Contour plot
  contour <- ggplot(data = df, aes_string(x=x_axis, y=y_axis))
  contour <- contour + geom_contour_filled(aes(z=test_stat_value), breaks = breaks) + facets
  contour <- contour + colors
  contour <- contour + labels
  
  #Raster plot
  df$cuts <- cut(df$test_stat_value, breaks)
  raster <- ggplot(data = df, aes_string(x=x_axis, y=y_axis))
  raster <- raster + geom_raster(aes(fill=cuts)) + facet_grid(rows = vars(test), cols = vars(test_stat))
  raster <- raster + colors
  raster <- raster + labels
  
  #Raster plot, different color scales
  rasters <- list()
  for(ts in unique(df$test_stat)){
    options(repr.plot.width = 10, repr.plot.height =2)
    df_subset <- subset(df,test_stat==ts)
    clabels <- sapply(breaks, (function(x) (if(x<1) bquote(""<=.(x)) else bquote(""<=1))))
    #drop the first element, which is the <=-inf level
    clabels <- clabels[2:length(clabels)]

    colors <- scale_fill_brewer(palette = "RdYlBu", drop=F, name = ts, labels = clabels, direction=if(ts=="Power") 1 else -1)
    plt <- ggplot(data = df_subset, aes_string(x=x_axis, y=y_axis))  + scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0))
    plt <- plt + geom_raster(aes(fill=cuts)) + facet_grid(rows=vars(test), cols=vars(test_stat))
    plt <- plt + colors
    plt <- plt + theme(legend.position = "bottom", legend.spacing.x = unit(0,'cm'), legend.title.align=0.5) + guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1))
    # plt <- plt + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank)
    plt <- plt + labels
    if(ts=="Power"){
      plt <- plt + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
    }else{
      plt <- plt + theme(strip.text.y = element_blank())
    }
    rasters[[length(rasters) + 1]] <- plt
    
  }
  
  coloredrasters <- rasters[[1]] + rasters[[2]]
  
  #Raster plot, different color scales, KS test only
  rasters <- list()
  for(ts in unique(df$test_stat)){
    options(repr.plot.width = 10, repr.plot.height =2)
    df_subset <- subset(df,test_stat==ts)
    df_subset <- subset(df_subset, test=="KS Test")
    clabels <- sapply(breaks, (function(x) (if(x<1) bquote(""<=.(x)) else bquote(""<=1))))
    #drop the first element, which is the <=-inf level
    clabels <- clabels[2:length(clabels)]
    
    colors <- scale_fill_brewer(palette = "RdYlBu", drop=F, name = ts, labels = clabels, direction=if(ts=="Power") 1 else -1)
    plt <- ggplot(data = df_subset, aes_string(x=x_axis, y=y_axis))  + scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0))
    plt <- plt + geom_raster(aes(fill=cuts)) + facet_grid(cols=vars(test_stat))
    plt <- plt + colors
    plt <- plt + theme(legend.position = "bottom", legend.spacing.x = unit(0,'cm'), legend.title.align=0.5) + guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1))
    # plt <- plt + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank)
    plt <- plt + labels
    if(ts=="Power"){
      plt <- plt + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
    }else{
      plt <- plt + theme(strip.text.y = element_blank())
    }
    rasters[[length(rasters) + 1]] <- plt
    
  }
  ks_only <- rasters[[1]] + rasters[[2]]
  plots <- list()
  plots$data <- df
  plots$acceptable <- accept
  plots$contour <- contour
  plots$raster<-raster
  plots$ks_only <- ks_only
  plots$coloredraster <- coloredrasters
  return(plots)
}

# paired_effectsize <- plot_heatmaps(x_axis = "effect_size",
#                                    x_axis_min = 0.2,
#                                    x_axis_max = 0.8,
#                                    y_axis = "prop_residual_variance_between_pairs", 
#                                    y_axis_min = 0,
#                                    y_axis_max = 0.5,
#                                    number_of_neurons = 20, 
#                                    total_sample_size = 1000, 
#                                    bootstrap_samples = 10, 
#                                    grid_frequency = 3, 
#                                    paired=TRUE)
# 
# paired_variances <- plot_heatmaps(x_axis = "proportion_explainable_variance",
#                                   x_axis_min = 0.05,
#                                   x_axis_max = 0.95,
#                                   y_axis = "prop_residual_variance_between_pairs", 
#                                   y_axis_min = 0,
#                                   y_axis_max = 0.5,
#                                   number_of_neurons = 20, 
#                                   total_sample_size = 1000, 
#                                   bootstrap_samples = 10, 
#                                   grid_frequency = 3, 
#                                   paired=TRUE)
# 
# paired_var <- plot_heatmaps(x_axis = "effect_size", 
#                                x_axis_min = 0.05,
#                                x_axis_max = 0.95,
#                                y_axis = "proportion_explainable_variance", 
#                                y_axis_min = 0.05,
#                                y_axis_max = 0.95,
#                                number_of_neurons = 20, 
#                                total_sample_size = 1000, 
#                                bootstrap_samples = 10, 
#                                grid_frequency = 3, paired=TRUE, 
#                                prop_residual_variance_between_pairs=0.05)


if(TRUE){
  #MAKE SOME UNPAIRED FIGURES
  # effect_size_var_ratio <- plot_heatmaps(x_axis = "effect_size", y_axis = "proportion_explainable_variance", number_of_neurons = 100)
  effect_size_interv_ctrl <- plot_heatmaps(x_axis = "effect_size", y_axis = "intervention_variance_ratio", 
                                           y_axis_min=0.5, y_axis_max=1.5, proportion_explainable_variance = 0)

  Nn_var_ratio <- plot_heatmaps(x_axis = "number_of_neurons", y_axis = "proportion_explainable_variance", total_sample_size = 100,
                                x_axis_min = 2, x_axis_max = 20)
  # Nn_interv_ctrl <- plot_heatmaps(x_axis = "number_of_neurons", y_axis = "intervention_variance_ratio", y_axis_min=0.5, y_axis_max=1.5)
  
  #MAKE 
}