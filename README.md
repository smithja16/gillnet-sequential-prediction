# Code to run sequential prediction for gillnet analysis

This R code demonstrates the model fitting and sequential prediction process used for the NSW Estuary General gillnet fishery.
It accompanies the paper by Smith, Tyler, and Johnson (in review), NSW DPIRD.
This code cannot be run because the data is confidential. Instead, it is used to show in detail the structure of the analysis.
This code is also not developed as software that could easily be used with different data, because the code contains things
like variables names that are specific to our data. However, the general 'Bayesian sequential prediction' approach is generally
useful for similar studies that require the prediction of multiple linked variables.

First the models are fit (1 Model fitting.R) and these are then used for prediction using new data (2 Model prediction.R).
For insight into the sequential prediction part, examine the 'make_sequential_predictions' and 'process_predictions_temporal'
functions.
