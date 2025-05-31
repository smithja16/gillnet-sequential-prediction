###############################################################
###             ESTUARY MESH NET ANALYSIS                   ###
###  Predicting unrecorded fishing method and gear details  ### 
###   to aid prediction of fishing effort and discards      ###
###         J.A.Smith - NSW DPIRD - Jan 2025+               ###
###############################################################

## This script provides the functions for sequential prediction with 
## full uncertainty propagation.
## The first function will run in a sequential loop for a set of new data,
## but this will become very slow once n > 1000. 
## The second function uses the first and runs in parallel for larger data sets.
## The parallel function saves text files to track progress in a dedicated folder.


## Sequential loop

make_sequential_predictions <- function(
    method_model,          # the 'brms' model objects
    mesh_model,
    net_model,
    discards_model,
    new_data,              # the new data for prediction, in our case the logbook data
    scaled_data,           # the data used for model fitting with scaled species predictors
    n_samples_use = 1000,  # the samples from the posterior of the fitted model
    progress_file = NULL,  # only used when running in parallel
    chunk_num = NULL) {    # only used when running in parallel
  
  # Get available number of posterior samples from all models
  n_samples_method <- nrow(as_draws_matrix(method_model))
  n_samples_mesh <- nrow(as_draws_matrix(mesh_model))
  n_samples_net <- nrow(as_draws_matrix(net_model))
  n_samples_discards <- nrow(as_draws_matrix(discards_model))
  
  # minimum samples available
  n_samples_available <- min(n_samples_method, n_samples_mesh, n_samples_net, n_samples_discards)
  n_new_data <- nrow(new_data)
  
  # Check if requested samples exceeds available
  if(n_samples_use > n_samples_available) {
    warning(paste("Requested samples exceed available. Using", n_samples_available, "samples"))
    n_samples_use <- n_samples_available
  }
  
  # Randomly select sample indices to use
  sample_indices <- sample(1:n_samples_available, n_samples_use)
  
  # Initialize storage
  predictions_discarded_full <- matrix(NA, nrow = n_samples_use, ncol = n_new_data)
  predictions_method_all <- array(NA, dim = c(n_samples_use, n_new_data, 3))  # 3 fishing methods
  predictions_mesh_all <- matrix(NA, nrow = n_samples_use, ncol = n_new_data)  # matrix for binary
  predictions_net_all <- array(NA, dim = c(n_samples_use, n_new_data, 3))  # 3 net length categories
  
  scale_params <- list()  #for storing species scaling parameters
  
  # Species that require scaling
  species_vars <- c("Bream_Yellowfin", "Mullet_Sea", "Mulloway", "Luderick", 
                    "Flathead_Dusky", "Whiting_Sand", "Shark_Bull")
  
  # Add species scaling parameters for reference
  for(species in species_vars) {
    if(species %in% names(scaled_data)) {
      attr_center <- attr(scaled_data[[species]], "scaled:center")
      attr_scale <- attr(scaled_data[[species]], "scaled:scale")
      
      if(!is.null(attr_center) && !is.null(attr_scale)) {
        scale_params[[paste0(species, "_mean")]] <- attr_center
        scale_params[[paste0(species, "_sd")]] <- attr_scale
      }
    }
  }
  
  # Predict in loop
  for(i in 1:n_samples_use){
    if(i %% 50 == 0) {
      # Write progress to file if provided
      if(!is.null(progress_file) && !is.null(chunk_num)) {
        cat(sprintf("Chunk %d: Processing sample %d of %d\n", 
                    chunk_num, i, n_samples_use),
            file = progress_file,
            append = TRUE)
      } else {
        message(sprintf("Processing sample %d of %d", i, n_samples_use))
      }
    }
    
    s_idx <- sample_indices[i]  # Get the actual sample index to use
    temp_new_data <- new_data %>% as.data.frame()
    
    # Scale species variables in place
    for(species in species_vars) {
      if(species %in% names(temp_new_data) && species %in% names(scaled_data)) {
        attr_center <- attr(scaled_data[[species]], "scaled:center")
        attr_scale <- attr(scaled_data[[species]], "scaled:scale")
        
        if(!is.null(attr_center) && !is.null(attr_scale)) {
          temp_new_data[[species]] <- (temp_new_data[[species]] - attr_center) / attr_scale
        }
      }
    }
    
    # Predict method
    predictions_method_temp <- posterior_epred(
      method_model,
      newdata = temp_new_data,
      allow_new_levels = T,  # * when Estuary is a random effect
      summary = FALSE)[s_idx,,]
    predicted_method_class <- apply(predictions_method_temp, 1, function(x){
      colnames(predictions_method_temp)[which.max(x)]
    })
    temp_new_data$Fishing_method <- predicted_method_class
    predictions_method_all[i,,] <- predictions_method_temp
    
    # Predict Mesh (binary outcome)
    predictions_mesh_temp <- posterior_epred(
      mesh_model, newdata = temp_new_data,
      allow_new_levels = T,
      summary = FALSE)[s_idx,]
    predicted_mesh_class <- ifelse(predictions_mesh_temp > 0.5, "Small", "Large")
    temp_new_data$Mesh_size <- predicted_mesh_class
    predictions_mesh_all[i,] <- predictions_mesh_temp
    
    # Predict categorical net length
    predictions_net_temp <- posterior_epred(
      net_model, newdata = temp_new_data,
      allow_new_levels = T,
      summary = FALSE)[s_idx,,]
    predicted_net_class <- apply(predictions_net_temp, 1, function(x){
      colnames(predictions_net_temp)[which.max(x)]
    })
    temp_new_data$Net_length <- predicted_net_class
    predictions_net_all[i,,] <- predictions_net_temp
    
    # Predict discards
    predictions_discarded_full[i,] <- posterior_predict(
      discards_model, 
      newdata = temp_new_data,
      allow_new_levels = T,
      summary = FALSE)[s_idx,]
  }
  
  # Calculate summary statistics
  summary_stats <- list(
    discards = list(
      mean = colMeans(predictions_discarded_full),
      median = apply(predictions_discarded_full, 2, median),
      q025 = apply(predictions_discarded_full, 2, quantile, probs = 0.025),
      q975 = apply(predictions_discarded_full, 2, quantile, probs = 0.975) ) )
  
  return(list(
    predictions = predictions_discarded_full,
    intermediate = list(
      method = predictions_method_all,
      mesh = predictions_mesh_all,
      net = predictions_net_all ),
    summary = summary_stats,
    n_samples_use = n_samples_use,
    scale_params = scale_params ))
}



## When running in parallel

library(future)
library(furrr)
library(dplyr)

run_parallel_predictions <- function(
    method_model,
    mesh_model,
    net_model,
    discards_model,
    new_data,
    scaled_data,
    n_samples_use,
    chunk_size,        # the number of new observations per chunk
    n_cores = parallelly::availableCores() - 2) {
  
  # Clean up memory
  gc()
  
  # Create progress directory if it doesn't exist
  progress_dir <- "progress"
  if (!dir.exists(progress_dir)) {
    dir.create(progress_dir)
  } else {
    # Clean up any existing progress files
    unlink(file.path(progress_dir, "*.txt"))
  }
  
  # Set up parallel processing
  plan(multisession, workers = n_cores)
  
  # Split data into chunks and create chunk info
  n_rows <- nrow(new_data)
  chunk_indices <- split(1:n_rows, ceiling(seq_along(1:n_rows)/chunk_size))
  n_chunks <- length(chunk_indices)
  
  # Create a list with chunk numbers and indices
  chunk_list <- Map(function(idx, num) list(indices = idx, chunk_num = num),
                    chunk_indices,
                    seq_along(chunk_indices))
  
  # Show initial message
  message(sprintf("Processing %d chunks of data in parallel using %d cores", n_chunks, n_cores))
  
  # Process chunks in parallel
  chunk_results <- future_map(
    chunk_list,
    function(chunk_info) {
      
      # Garbage collection
      gc()
      
      # Suppress package startup messages
      suppressPackageStartupMessages({
        library(posterior)
        library(brms)
        library(dplyr)
      })
      
      # Get chunk information
      chunk_num <- chunk_info$chunk_num
      indices <- chunk_info$indices
      
      progress_file <- file.path(progress_dir, sprintf("chunk_%d_progress.txt", chunk_num))
      
      # Write initial progress
      cat(sprintf("Starting chunk %d with %d rows\n", chunk_num, length(indices)), 
          file = progress_file)
      
      # Extract chunk data
      chunk_data <- new_data[indices, ]
      
      tryCatch({
        # Run make_sequential_predictions_full on this chunk
        results <- make_sequential_predictions(
          method_model = method_model,
          mesh_model = mesh_model,
          net_model = net_model,
          discards_model = discards_model,
          new_data = chunk_data,
          scaled_data = scaled_data,
          n_samples_use = n_samples_use,
          progress_file = progress_file,
          chunk_num = chunk_num )
        
        # Add chunk information
        results$chunk_info <- list(
          start_idx = min(indices),
          end_idx = max(indices),
          n_rows = length(indices) )
        
        # Write completion to progress file
        cat(sprintf("Completed chunk %d\n", chunk_num),
            file = progress_file,
            append = TRUE)
        
        # Garbage collection
        gc()
        
        results
        
      }, error = function(e) {
        # Write error to progress file
        cat(sprintf("Error in chunk %d: %s\n", chunk_num, as.character(e)),
            file = progress_file,
            append = TRUE)
        stop(e)
      })
    },
    .options = furrr_options(seed = TRUE) )
  
  # Read and display final progress
  progress_files <- list.files(progress_dir, full.names = TRUE)
  message("\nProgress summary:")
  for(file in progress_files) {
    cat("\n", readLines(file), sep = "\n")
  }
  
  # Combine results
  message("\nCombining results from all chunks...")
  
  # Initialize combined storage
  n_methods <- dim(chunk_results[[1]]$intermediate$method)[3]
  n_nets <- dim(chunk_results[[1]]$intermediate$net)[3]
  combined_results <- list(
    predictions = matrix(NA_real_, 
                         nrow = n_samples_use, 
                         ncol = n_rows),
    intermediate = list(
      method = array(NA_real_, 
                     dim = c(n_samples_use, n_rows, n_methods)),
      mesh = matrix(NA_real_, 
                    nrow = n_samples_use, 
                    ncol = n_rows),
      net = array(NA_real_, 
                  dim = c(n_samples_use, n_rows, n_nets)) ) )
  
  # Fill in results from each chunk
  for(i in seq_along(chunk_results)) {
    chunk <- chunk_results[[i]]
    idx_start <- chunk$chunk_info$start_idx
    idx_end <- chunk$chunk_info$end_idx
    
    combined_results$predictions[, idx_start:idx_end] <- chunk$predictions
    combined_results$intermediate$method[, idx_start:idx_end, ] <- chunk$intermediate$method
    combined_results$intermediate$mesh[, idx_start:idx_end] <- chunk$intermediate$mesh
    combined_results$intermediate$net[, idx_start:idx_end, ] <- chunk$intermediate$net
  }
  
  # Calculate summary statistics
  message("Calculating final summary statistics...")
  combined_results$summary <- list(
    discards = list(
      mean = colMeans(combined_results$predictions),
      median = apply(combined_results$predictions, 2, median),
      q025 = apply(combined_results$predictions, 2, quantile, probs = 0.025),
      q975 = apply(combined_results$predictions, 2, quantile, probs = 0.975) ) )

  message("Processing complete!")
  return(combined_results)
}

