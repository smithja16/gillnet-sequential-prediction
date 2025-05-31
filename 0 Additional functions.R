###############################################################
###             ESTUARY MESH NET ANALYSIS                   ###
###  Predicting unrecorded fishing method and gear details  ### 
###   to aid prediction of fishing effort and discards      ###
###         J.A.Smith - NSW DPIRD - Jan 2025+               ###
###############################################################

## This script loads some helper functions for processing data


## Function to convert predictions to categories

get_categorical_predictions <- function(predictions) {
  # For fishing method (3 categories)
  method_probs <- predictions$intermediate$method
  method_preds <- array(NA, dim = c(dim(method_probs)[1], dim(method_probs)[2]))
  
  for(i in 1:dim(method_probs)[1]) {  # for each sample
    for(j in 1:dim(method_probs)[2]) {  # for each trip
      # Get category with highest probability
      method_preds[i,j] <- method_levels[which.max(method_probs[i,j,])]
    }
  }
  
  # For mesh size (binary)
  mesh_preds <- matrix(mesh_levels[1 + (predictions$intermediate$mesh > 0.5)], 
                       nrow = dim(predictions$intermediate$mesh)[1])
  
  # For net length (3 categories)
  net_probs <- predictions$intermediate$net
  net_preds <- array(NA, dim = c(dim(net_probs)[1], dim(net_probs)[2]))
  
  for(i in 1:dim(net_probs)[1]) {  # for each sample
    for(j in 1:dim(net_probs)[2]) {  # for each trip
      # Get category with highest probability
      net_preds[i,j] <- net_levels[which.max(net_probs[i,j,])]
    }
  }
  
  return(list(
    method = method_preds,
    mesh = mesh_preds,
    net = net_preds ))
}


## Trip-level probability of each fishing method

get_method_probabilities <- function(
    predictions,
    trip_index) {
  
  # Get probabilities across all samples for this trip
  method_probs <- predictions$intermediate$method[, trip_index,]
  
  # Calculate summary statistics
  prob_summary <- data.frame(
    Method = method_levels,
    Mean_prob = colMeans(method_probs),
    Lower_CI = apply(method_probs, 2, quantile, 0.025),
    Upper_CI = apply(method_probs, 2, quantile, 0.975) )
  
  return(prob_summary)
}


## Aggregation of fishing method with uncertainty

get_method_totals_by_group <- function(
    predictions,
    new_data,
    group_var,
    method_levels) {
  
  # Input validation
  if(!group_var %in% names(new_data)) {
    stop(paste("Grouping variable", group_var, "not found in new_data"))
  }
  
  n_samples <- dim(predictions$intermediate$method)[1]
  
  # Print method levels for verification
  print("Method levels:")
  print(method_levels)
  
  # Get unique groups and sort them
  groups <- sort(unique(new_data[[group_var]]))
  print(paste("Processing", length(groups), "groups:"))
  print("All groups:")
  print(groups)
  
  # For each posterior sample, count methods by group
  group_counts <- list()  # Initialize empty list
  
  for(g in groups) {
    # Get trips for this group
    trip_idx <- which(new_data[[group_var]] == g)
    print(paste("Processing group:", g, "- has", length(trip_idx), "trips"))
    
    # Initialize counts matrix for this group
    sample_counts <- matrix(0, nrow = n_samples, ncol = length(method_levels),
                            dimnames = list(NULL, method_levels))
    
    if(length(trip_idx) > 0) {
      for(i in 1:n_samples) {
        tryCatch({
          # Get predicted methods for this sample
          pred_matrix <- predictions$intermediate$method[i, trip_idx, ]
          if(length(trip_idx) == 1) {
            # Handle single trip case
            methods <- method_levels[which.max(pred_matrix)]
          } else {
            methods <- apply(pred_matrix, 1, function(x) method_levels[which.max(x)])
          }
          sample_counts[i,] <- table(factor(methods, levels = method_levels))
        }, error = function(e) {
          print(paste("Error processing group:", g, "sample:", i))
          print(e)
        })
      }
    }
    
    group_counts[[as.character(g)]] <- sample_counts
  }
  
  # Create summary by group
  results <- lapply(names(group_counts), function(g) {
    counts <- group_counts[[g]]
    data.frame(
      Group = g,
      Method = method_levels,
      Mean_count = colMeans(counts),
      Lower_CI = apply(counts, 2, quantile, 0.025),
      Upper_CI = apply(counts, 2, quantile, 0.975) )
  })
  
  # Combine results and rename Group column
  results_df <- do.call(rbind, results)
  names(results_df)[names(results_df) == "Group"] <- group_var
  
  return(results_df)
}


## Process discard predictions by estuary and mesh size

process_predictions <- function(
    predictions,
    new_data) {
  
  # Get number of samples and observations
  n_samples <- nrow(predictions$predictions)
  n_obs <- ncol(predictions$predictions)
  
  # Get predicted fishing methods for each observation
  method_predictions <- apply(predictions$intermediate$method, 2, function(x) {
    # For each observation, get the most common predicted method across samples
    method_counts <- apply(matrix(x, nrow = n_samples, byrow = FALSE), 1, which.max)
    # Return the most common method index
    as.numeric(names(sort(table(method_counts), decreasing = TRUE)[1]))
  })
  
  # Get predicted mesh sizes for each observation
  mesh_predictions <- apply(predictions$intermediate$mesh, 2, function(x) {
    # For each observation, calculate the proportion of "Small" predictions
    mean_prob <- mean(x > 0.5)
    ifelse(mean_prob > 0.5, "Small", "Large")
  })
  
  # Create a data frame of predictions with grouping variables
  results_df <- data.frame(
    Type_of_fishing_method = rep(method_predictions, each = n_samples),
    Mesh_size_f = rep(mesh_predictions, each = n_samples),
    Discards = as.vector(predictions$predictions) )
  
  # Ensure method levels are correct
  results_df$Fishing_method <- factor(
    results_df$Fishing_method,
    levels = 1:3,
    labels = method_levels )
  
  # Ensure mesh size levels are correct
  results_df$Mesh_size <- factor(
    results_df$Mesh_size,
    levels = mesh_levels )
  
  # Calculate summary statistics by group
  summary_by_group <- results_df %>%
    group_by(Fishing_method, Mesh_size) %>%
    summarise(
      mean_discards = mean(Discards),
      median_discards = median(Discards),
      lower_ci = quantile(Discards, 0.025),
      upper_ci = quantile(Discards, 0.975),
      n_trips = n()/n_samples  # Divide by n_samples to get actual number of trips
    ) %>%
    ungroup()
  
  return(summary_by_group)
}


## Aggregate discards into month totals

process_predictions_temporal <- function(
    predictions,
    new_data) {
  
  # Get dimensions
  n_samples <- nrow(predictions$predictions)
  n_trips <- ncol(predictions$predictions)
  
  # Create array of predictions: samples x trips
  pred_array <- predictions$predictions
  
  # Create unique month-year combinations
  month_year_combos <- unique(data.frame(
    Month = new_data$Monthf,
    Year = new_data$cal_year
  ))
  
  # For each posterior sample, sum discards by month and year
  monthly_yearly_sums <- lapply(1:nrow(month_year_combos), function(i) {
    # Get current month and year
    m <- month_year_combos$Month[i]
    y <- month_year_combos$Year[i]
    
    # Get trips for this month-year combination
    trip_idx <- which(new_data$Monthf == m & new_data$cal_year == y)
    
    # For each posterior sample, sum those trips
    month_year_totals <- rowSums(pred_array[, trip_idx, drop = FALSE])
    
    return(month_year_totals)  # Returns n_samples sums for this month-year
  })
  
  # Convert to data frame with uncertainty
  results <- data.frame(
    Month = month_year_combos$Month,
    Year = month_year_combos$Year,
    Mean_total = sapply(monthly_yearly_sums, mean),
    Lower_CI = sapply(monthly_yearly_sums, quantile, 0.025),
    Upper_CI = sapply(monthly_yearly_sums, quantile, 0.975) )
  
  return(list(
    summary = results,
    sample_sums = monthly_yearly_sums  # Keep all posterior samples
  ))
}