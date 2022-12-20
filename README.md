![Check package](https://github.com/Craig44/SpatialSablefishAssessment/actions/workflows/r.yml/badge.svg)
# SpatialSablefishAssessment
An R package that contains the TMB model and auxillary functions for running a spatially explicit model for Alaskan sablefish.


# Install 
To install this R package run the following command. Before installing, it is advised to check that the R package is passing all unit-tests.

```r
devtools::install_github("Craig44/SpatialSablefishAssessment")
```

# Using `SpatialSablefishAssessment`

## Configuring and checking model inputs (data and parameters)
This function contains TMB models and summarizing and plotting functions for the spatial model, users are responsible for configuring the `data` and `parameter` structures that are passed to the TMB `MakeADFun` function which compiles the model. The best place to look for what dimensions and names of expected elements for `data` and `parameter` it is best to look at an example or the source code for the model ([see here](https://github.com/Craig44/SpatialSablefishAssessment/blob/master/src/TMB/TagIntegrated.hpp)). To view an example data run 


```r
load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
names(data)
names(parameters)
region_key
```


Once you have built the `data` and `parameter` objects then you can check that the dimensions are consistent with the model by using,

```r
validate_input_data_and_parameters(data, parameters)
```

This function will report a message telling you either success or if one of the objects dimensions are off, then you will need to fix this because it is likely if you pass this to the TMBs `MakeADFun` it will crash your R session.


Once you have checked the `data` and `parameter` objects you can build the TMB model following,

```r
my_model = TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
```

### Check for zero gradients
It is possible to configure a model that has zero gradients for parameters, and you want to check this is not the case using the 
```r
check_gradients(my_model)
```


## Optimisation

```r
mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))

# Try and improve the optimsation running the model for two additional Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(my_model$gr(mle_optim$par))
                           h = optimHess(mle_spatial$par, fn = my_model$fn, gr = my_model$gr)
                           mle_optim$par = mle_spatial$par - solve(h,g)
                           mle_optim$objective = my_model$fn(mle_optim$par)
                         }
                       , error = function(e){e})

try_improve
```


### Check convergence


```r
## largest gradient 
names(mle_optim$par)[which.max(abs(my_model$gr(mle_optim$par)))]
```


### Fixing estimable parameters at input values


### Sharing parameter estimates among estimable variables i.e., male and female having a common selectivity



## Model Summary
Once you are satisfied that the model has converged at global minimum, then you can get the model to report model quantities

```r
## get MLE outputs
mle_report = my_model$report(mle_optim$par)



## movement assumption
plot_movement(mle_report, region_key = region_key)
# starting values
plot_movement(OM_report, region_key)
# plot selectivities
plot_selectivities(mle_report)
# starting values
plot_selectivities(OM_report)

## Regional SSBs
ssb_plt = plot_SSB(mle_report, region_key = region_key)
ssb_plt
## sum SSB over all regions
ggplot(ssb_plt$data %>% group_by(Year) %>% summarise(total_ssb = sum(SSB)), aes(x = Year, y = total_ssb)) +
  geom_line(linewidth = 1.1) +
  ylim(0,NA)

## plot F's
plot_fishing_mortalities(MLE_report = mle_report, region_key = region_key)

## plot fits to observations

# catch and index
plot_catch_fit(MLE_report = mle_report, region_key = region_key) + facet_wrap(label~Region, ncol = n_regions) + ylab("Catch (mt)")
plot_index_fit(MLE_report = mle_report, region_key = region_key)

# survey AFs
first_year_set = c(1981, seq(from = 1985, to = 1993, by = 2), 1996:1999)
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "male") +
  ggtitle("Male survey AF") +
  guides(shape = "none", linetype = "none") +
  ylim(0,20)
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "female") +
  ggtitle("Female survey AF") +
  guides(shape = "none", linetype = "none")
  
# fishery LFs
# Fixed gear
plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1999, sex = "male") +
  ggtitle("Male fixed LF") +
  guides(shape = "none", linetype = "none")
plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1999, sex = "female") +
  ggtitle("Female fixed LF") +
  guides(shape = "none", linetype = "none")

# plot by tag releases event
tag_fits_by_age_sex = get_tag_recovery_obs_fitted_values(mle_report, region_key = region_key)
tag_fits = tag_fits_by_age_sex %>% group_by(release_event, unique_recovery_id) %>% summarise(observed = sum(observed), predicted = sum(predicted), recovery_year = unique(recovery_year), recovery_region = unique(recovery_region), release_region = unique(release_region), release_year = unique(release_year))
tag_fits$resid = tag_fits$observed - tag_fits$predicted
tag_fits$stand_resid = tag_fits$resid / tag_fits$predicted
tag_fits$resid_sign = ifelse(tag_fits$resid < 0, "negative", "positive")
unique_release_events = unique(tag_fits$release_event)
gplt = ggplot(tag_fits %>% filter(release_event %in% unique_release_events[1:16]), aes(x = recovery_year, y = recovery_region)) +
  geom_point(aes(size = abs(stand_resid), col = resid_sign)) +
  labs(x = "Recovery year", y = "Recovery region", size = "Pearson residuals") +
  theme_bw() +
  scale_size_area() +
  facet_wrap(~release_event) +
  scale_x_discrete(breaks = every_nth(n = 5))

```



## Standard errors for derived quantities
