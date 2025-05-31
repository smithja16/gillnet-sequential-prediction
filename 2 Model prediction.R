###############################################################
###             ESTUARY MESH NET ANALYSIS                   ###
###  Predicting unrecorded fishing method and gear details  ### 
###   to aid prediction of fishing effort and discards      ###
###         J.A.Smith - NSW DPIRD - Jan 2025+               ###
###############################################################

## This script uses the fitted models to predict new data
## In our case, the new data were commercial logbooks

library(brms)

source("0 Prediction functions.R")  # load main prediction functions
source("0 Additional functions.R")  # load additional helper functions


######################
##     STEP 1       ##
## Predict New Data ##
######################

## To predict gear characteristics and discards for small amount of new data (i.e. logbooks)
# If more than ~ 500 new data and 1000 posterior samples use the parallel version
predictions <- make_sequential_predictions(
  method_model = method_model,
  mesh_model = mesh_model,
  net_model = net_model,
  discards_model = discards_model,
  new_data = logbook_data,
  scaled_data = trip_catch_data,
  n_samples_use = 1000 )  #how many posterior samples to predict from


## To predict for many new samples
# The new_data is split into chunks and predicted in parallel

t1 <- Sys.time()
predictions <- run_parallel_predictions(
  method_model = method_model,
  mesh_model = mesh_model,
  net_model = net_model,
  discards_model = discards_model,
  new_data = logbook_data,
  scaled_data = trip_catch_data,
  n_samples_use = 500,
  chunk_size = 250,  # rows of new_data per chunk (keep < 500)
  n_cores = 8 )
t2 <- Sys.time()
t2 - t1  #time elapsed


##########################
##       STEP 2         ##
## Evaluate Predictions ##
##########################

## Basic summary
predictions$summary

# Trip discard rates
mean(predictions$summary$discards$mean)
quantile(predictions$summary$discards$mean, probs=c(0.025,0.5,0.975))

## Access categorical intermediate predictions
# Get method, mesh, net categories
method_levels <- levels(method_model$data$Fishing_method)
mesh_levels <- levels(mesh_model$data$Mesh_size)
net_levels <- levels(net_model$data$Net_length)

## Summary of categorical predictions
# e.g. can be used to find the % of overnight sets that use large mesh
cat_sum <- get_categorical_predictions(predictions = predictions)

## Trip level probability of each fishing method
trip_probs <- get_method_probabilities(predictions, trip_index = 1)

## Aggregation of fishing method with uncertainty
# Monthly totals
monthly_summary <- get_method_totals_by_group(
  predictions = predictions,
  new_data = logbook_data,
  group_var = "Month",
  method_levels = method_levels )
# Estuary totals
estuary_summary <- get_method_totals_by_group(
  predictions =predictions,
  new_data = logbook_data,
  group_var = "Estuary",
  method_levels = method_levels )

## Process discard predictions by estuary and mesh size
summary_stats <- process_predictions(predictions, new_data)

## Aggregate discards into month totals
agg_pred <- process_predictions_temporal(
  predictions = predictions,
  new_data = logbook_data )
discard_series <- agg_pred$summary

