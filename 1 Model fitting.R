###############################################################
###             ESTUARY MESH NET ANALYSIS                   ###
###  Predicting unrecorded fishing method and gear details  ### 
###   to aid prediction of fishing effort and discards      ###
###         J.A.Smith - NSW DPIRD - Jan 2025+               ###
###############################################################

## This script fits the Bayesian GLMs
## In our case, the data for fitting were those from an observer survey

library(brms)
library(bayesplot)


##################
##   STEP 1     ##
## Fit Models   ##
##################

## Ensure the continuous predictors are centered and scaled
## Month and Estuary are factors

## Model 1: fishing method
method_model <- brm(
  bf(Fishing_method ~ Mullet_Sea + Mulloway + Shark_Bull + 
       Bream_Yellowfin + Flathead_Dusky + Luderick + Whiting_Sand +
       Estuary + Month),
  family = categorical(link = "logit"),
  data = trip_catch_data,
  prior = c(
    # Priors for intercepts
    prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "muSetmorethan3hr"),  #names of non-reference method levels
    prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "muSplash"),
    # Priors for coefficients
    prior(normal(0, 3), class = "b", dpar = "muSetmorethan3hr"),
    prior(normal(0, 3), class = "b", dpar = "muSplash") ),
  control = list(adapt_delta = 0.98),
  iter = 8000, warmup = 4000 )


## Model 2: gillnet mesh size
mesh_model <- brm(
  bf(Mesh_size ~ Mullet_Sea + Mulloway + Shark_Bull + 
       Bream_Yellowfin + Flathead_Dusky + Luderick + Whiting_Sand +
       Estuary + Month),
  family = bernoulli(link = "logit"),  #only two levels so model is binary
  data = trip_catch_data,
  prior = c(
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(normal(0, 3), class = "b") ),
  control = list(adapt_delta = 0.98),
  iter = 8000, warmup = 4000 )


## Model 3: net length
net_model <- brm(
  bf(Net_length ~ Estuary + Month + Mesh_size),
  family = categorical(link = "logit"),
  data = trip_catch_data,
  prior = c(
    # Priors for intercepts
    prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "muless200"),  #names of non-reference length levels
    prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "muover600"),
    # Priors for coefficients
    prior(normal(0, 3), class = "b", dpar = "muless200"),
    prior(normal(0, 3), class = "b", dpar = "muover600") ),
  control = list(adapt_delta = 0.98),
  iter = 8000, warmup = 4000 )


## Model 4: discards for chosen species; hurdle gamma
discards_model <- brm(
  bf(discards ~ Mulloway + Mullet_Sea + Fishing_method + 
       Net_length + Mesh_size + Estuary + Month,  # positive part
     hu ~ Mulloway + Mullet_Sea + Fishing_method + 
       Net_length + Mesh_size + Estuary + Month),  # presence part
  family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
  data = trip_catch_data,
  chains = 4,
  iter = 8000,
  warmup = 4000,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  prior = c(
    # Priors for the gamma part (fairly weak priors)
    prior(normal(0, 5), class = "b"),
    prior(student_t(3, 0, 5), class = "Intercept"),
    prior(exponential(1), class = "shape"),  # favours right-skewed data
    
    # Priors for hurdle part (fairly weak priors)
    prior(normal(0, 5), class = "b", dpar = "hu"),
    prior(normal(0, 5), class = "Intercept", dpar = "hu") ) )


#####################
##     STEP 2      ##
## Evaluate Models ##
#####################

## Basic evaluation

summary(discards_model)

# Chain convergence and diagnostics
rhat <- brms::rhat(discards_model)
neff <- neff_ratio(discards_model)
print(data.frame(
  Rhat_max = max(rhat),
  ESS_min = min(neff) ))
#Rhat should be ~ 1; ESS should be > 0.1

# Trace plots; check mixing and normal estimates
plot(discards_model, ask=F)

# Posterior predictive checks
pp_check(discards_model, ndraws = 100)
pp_check(discards_model, type = "stat", stat = "mean", bins=50)  # Check mean
pp_check(discards_model, type = "stat", stat = "sd", bins=50)    # Check standard deviation

# Look at conditional effects of predictors
plot(conditional_effects(discards_model), ask=F)

# For the hurdle part specifically
plot(conditional_effects(discards_model, dpar = "hu"), ask=F)

# Model performance metrics
loo_model <- loo(discards_model, moment_match = T)  #takes a minute
print(loo_model)

# Explore accuracy
predicted_values <- posterior_predict(discards_model)
actual <- trip_catch_data$discards
predicted <- colMeans(predicted_values)
rmse <- sqrt(mean((actual - predicted)^2))
mae <- mean(abs(actual - predicted))
cor_val <- cor(actual, predicted)

## Explore estimates

# Plot coefficient estimates with uncertainty (both hurdle components)
posterior <- as.matrix(discards_model)

main_pars <- grep("^b_[^h]", colnames(posterior), value = TRUE)
mcmc_intervals(posterior, pars = main_pars)

hurdle_pars <- grep("^b_hu", colnames(posterior), value = TRUE)
mcmc_intervals(posterior, pars = hurdle_pars)


