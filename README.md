![Check package](https://github.com/Craig44/SpatialSablefishAssessment/actions/workflows/r.yml/badge.svg)
# About
`SpatialSablefishAssessment` is a R package that contains TMB models that can be used as operating models (OMs) or Estimation Models (EMs), it also contains a bunch of utility R functions for checking inputs and visualizing/summarizing model outputs/fits to observations. The models are developed for the Alaskan sablefish stock region, but the TMB models are fairly general and can be applied to any assessment that has two fisheries, two sexes and assumes an annual cycle (time-step). There are two online Gitbook resources that are related to this package. The first is the [online documentation](https://craig44.github.io/SpatialSablefishAssessment/) which describes the TMB models, expected inputs and functions to investigate and summaries model fits. There is an additional online Gitbook that relates to my post-doctoral research which uses this R package [can be found here](https://craig44.github.io/SableFishResearch/).

There are currently three TMB models contained in this package

- `TagIntegrated` A generalized spatially disaggregated model that assumes two fisheries, two sexes and an annual cycle
- `TagIntegratedValidate` This model is used to unit-test `TagIntegrated`. Any change to `TagIntegrated` should be incorporated into `TagIntegratedValidate`
- `Assessment` The closest version to the current ADMB model that is used for the current assessment (needs further testing and comparisons to current ADMB assessment model before being used for management advice).


# Install 
Before installing, it is advised to check that the R package is passing all unit-tests and you have Rtools correctly installed [(see here)](https://cran.r-project.org/bin/windows/Rtools/). It is advised that you have TMB installed before attempting to install this package. In addition to installing TMB, make sure you can compile a TMB example model i.e., `TMB::compile(file = system.file("examples", "simple.cpp",package = "TMB"))`. Once TMB is installed and you can successfully compile a TMB model, you should be able to install this package following

```r
devtools::install_github("Craig44/SpatialSablefishAssessment")
```

# Using `SpatialSablefishAssessment`

## Configuring and checking model inputs (data and parameters)
This package contains TMB models and summarizing and plotting functions for a range of models used during my sablefish research. Users are responsible for configuring the `data` and `parameter` structures that are passed to the `TMB::MakeADFun` function which compiles the model. Each model type (defined by `data$model`) will expect slightly different input data and parameters structures. For

- `TagIntegrated` has described inputs [found here](https://craig44.github.io/SpatialSablefishAssessment/TagIntegrated.html)
- `Assessment` I haven't documented this model yet, the best place to look is the source TMB code [found here](https://github.com/Craig44/SpatialSablefishAssessment/blob/master/src/TMB/CurrentAssessment.hpp)


An alternative place to look for expected names and dimensions for elements of `data` and `parameter` is an example or the source code for the model ([see here](https://github.com/Craig44/SpatialSablefishAssessment/blob/master/src/TMB/TagIntegrated.hpp)). To view an example `data` and `parameter` for the `TagIntegrated` model run the following code 


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


The following list are some useful utility functions contained in this package to visualize data inputs. All functions should be documented and queried using the usual R `?` method

- `plot_input_observations`
- `plot_input_catches`
- `plot_mean_weight`
- `plot_age_length_matrix`


Once you have checked the `data` and `parameter` objects you can build the TMB model following,

```r
my_model = TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
```

### Simple sanity checks
It is possible to configure a model that has zero gradients for parameters, as well as NaN or Inf in the likelihood evaluations. These usually occur when users ask the model to calculate a predicted value in a region or year that has no fish which can cause undefined behavior in some likelihood functions. To check this is not the case use the following function
```r
pre_optim_sanity_checks(my_model)
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

After a model has successfully converged during optimization, you should check there are no parameters at bounds, gradients are small and you have a positive definite hessian matrix. This can be achieved by using the following function.

```r
post_optim_sanity_checks(my_model)
```


### Fixing estimable parameters at input values
```r
?fix_pars()
```

### Sharing parameter estimates among estimable variables i.e., male and female having a common selectivity or survey catchability the same among spatial regions
```r
?set_pars_to_be_the_same()
```

If you are using the `TagIntegrated` model then there is a special function called `set_up_parameters` that will set a range of parameters as fixed or shared depending on user inputs.


## Model Summary
Once you are satisfied that the model has converged at global minimum, then you can get the model to report model quantities

```r
## get MLE quantities
mle_report = my_model$report(mle_optim$par)

#########
## Plot interesting
## Model quantities i.e., recruitment SSBs model fits etc
#########

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
