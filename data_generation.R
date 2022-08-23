
#' gen_unpaired_data
#' Generates a set of mock miniature synaptic event data.
#' The data consists of 2 groups, a control and intervention group.
#' Each group contains a number of neurons.
#' Each neuron is associated with several samples.
#'
#' @param control_group_mean The baseline mean of the control group
#' @param treatment_effect  The mean difference between control and intervention groups
#' @param control_group_interneuron_sd  the SD between the means of individual neurons in the control group
#' @param intervention_group_interneuron_sd the SD between the means of individual neurons in the intervention group
#' @param within_neuron_sd the SD between samples taken from a single neuron
#' @param n_control_neurons number of neurons in the control group
#' @param n_intervention_neurons number of neurons in the intervention group
#' @param samples_per_neuron number of samples per neuron
#'
#' @return a dataframe with the columns (group, neuron_id, dependant), where each row is an observation
gen_unpaired_data <- function(control_group_mean=0, treatment_effect=1, control_group_interneuron_sd = 0.5, 
                             intervention_group_interneuron_sd=0.5, within_neuron_sd=0.5, n_control_neurons=5,
                             n_intervention_neurons=5, samples_per_neuron=100) {

  intervention_group_mean <- control_group_mean + treatment_effect
  
  #Set up an empty (correctly sized) dataframe to hold our control group data
  rows = n_control_neurons * samples_per_neuron
  dat_control <- data.frame(neuron_id = numeric(rows), group=character(rows),dependant=numeric(rows))
  
  #For each neuron...
  for (i in 1:n_control_neurons){
    #How far away from the group (control or intervention) mean is this neuron's mean?
    neuron_deviation <- rnorm(1, mean=0, sd = control_group_interneuron_sd)
    #Randomly generate all that neuron's observations...
    row_range = seq((i-1)*samples_per_neuron + 1, i*samples_per_neuron)
    dat_control[row_range,"dependant"] <- rnorm(samples_per_neuron, mean = control_group_mean + neuron_deviation, sd = within_neuron_sd)
    
    #Fill out the neuron's exogenous variables
    dat_control[row_range,"neuron_id"] <- i
    dat_control[row_range,"group"] <- "Control"
  }
  
  #Repeat for the intervention group
  rows = n_intervention_neurons * samples_per_neuron
  dat_intervention <- data.frame(neuron_id = numeric(rows), group=character(rows),dependant=numeric(rows))
  for (i in 1:n_intervention_neurons){
    neuron_deviation <- rnorm(1, mean=0, sd = intervention_group_interneuron_sd)
    row_range = seq((i-1)*samples_per_neuron + 1, i*samples_per_neuron)
    dat_intervention[row_range,"dependant"] <- rnorm(samples_per_neuron, mean = control_group_mean + neuron_deviation, sd = within_neuron_sd)
    dat_intervention[row_range,"neuron_id"] <- i
    dat_intervention[row_range,"group"] <- "Intervention"
  }
  
  #Glue the control and intervention groups together and return them
  dat <- rbind(dat_control,dat_intervention)
  
  return(dat)
}


plot_data <- function(dat){
  #Plot each neuron's empirical PDF (unfilled histogram) colored by that neuron's group
  plt <- ggplot(data=dat, aes(x = dependant, color = group))
  for (id in unique(dat$neuron_id)){
    plt <- plt + geom_freqpoly(data = dat[dat$neuron_id==id,])
  }
  
  return(plt)
}
