#
# Bookdown functions for summarising a model run with automated R-markdown code
#

#' summarise_individual_models
#' A bookdown function that will summarise a suite of models in individual Sections of a single Rmarkdown file
#' @param model_dir a vector of absolute paths were we can find model output to summarise
#' @param bookdown_dir an absolute path where we will build the Bookdown
#' @param bookdown_labels a vector of model labels
#' @param model_description A string that describes all the models that are pointed at by model_dir
#' @return Creates a suite of Rmarkdown files in bookdown_dir and compiles an R markdown file
#' @details this function expects the following RDS objects to be in each `model_dir` for each model, data.RDS region_key.RDS, parameters.RDS, mle_report.RDS, sd_report.RDS, mle_optim.RDS, map_fixed_pars.RDS
#' @export

summarise_individual_models <- function(model_dir, bookdown_dir, bookdown_labels, model_description = "") {
  if(length(model_dir) != length(bookdown_labels))
    stop("number of model_dir provided was not the same as bookdown_labels. These need to be the same length")

  ## set up Bookdown directory
  verbose = T
  output_dir = normalizePath(file.path(bookdown_dir), winslash = "/")
  if(dir.exists(output_dir)) {
    unlink(output_dir, recursive = T, force = T)
    if(verbose)
      print("deleting 'output_dir'")
  }

  dir.create(output_dir)
  if(verbose)
    print(paste0("creating 'output_dir' ", output_dir))
  ## build a skeleton Bookdown
  #bookdown:::bookdown_skeleton(output_dir)
  bookdown::create_bs4_book(output_dir)

  ## get rid of these skeleton changes. This may need tweaking over time
  files_to_remove = c("01-intro.Rmd", "02-cross-refs.Rmd", "03-parts.Rmd", "04-citations.Rmd","05-blocks.Rmd","06-share.Rmd","07-references.Rmd")
  for(i in 1:length(files_to_remove)) {
    file.remove(file.path(output_dir, files_to_remove[i]))
  }

  ## change the CSS
  css_file = file.path(output_dir, "style.css")
  #write("\n.book .book-body .page-wrapper .page-inner {", file = css_file, append = T)
  #write("  max-width: 1200px !important;", file = css_file, append = T)
  #write("}", file = css_file, append = T)

  ## change the _output.yml
  output_yml_file = file.path(output_dir, "_output.yml")
  first_chunk =
    'bookdown::bs4_book:
  css: style.css
  theme:
    primary: "#096B72"
  repo: https://github.com/rstudio/bookdown-demo
bookdown::pdf_book:
  includes:
    in_header: preamble.tex
  latex_engine: xelatex
  citation_package: natbib
  keep_tex: yes
bookdown::epub_book: default
'
  write(first_chunk, file = output_yml_file, append = F)

  ### change index.md

  index_md =
    "---
title: 'Sablefish models'
author: 'C Marsh'
date: '`r Sys.Date()`'
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
#cover-image: path to the social sharing image like images/cover.jpg
description: A detailed summary of sablefish model fits and estimated quantities
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---
"

  tail_index_chunk = "
```{r include=FALSE}
# automatically create a bib database for R packages
knitr::write_bib(c(
  .packages(), 'bookdown', 'knitr', 'rmarkdown'
), 'packages.bib')
```


```{r setup, include=FALSE}
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```
"
  write(index_md, file = file.path(output_dir, "index.Rmd"), append = F)
  write('# Model descriptions', file = file.path(output_dir, "index.Rmd"), append = T)
  write(model_description, file = file.path(output_dir, "index.Rmd"), append = T)
  write(tail_index_chunk, file = file.path(output_dir, "index.Rmd"), append = T)

  for(j in 1:length(bookdown_labels)) {

    model_run_label = bookdown_labels[j]
    ## header
    model_input_file = file.path(output_dir, paste0("0", j + 1, "-ModelInterrogate.Rmd"))
    file.create(model_input_file)
    header = paste0("# ", model_run_label, " {#model_",j,"}")
    model_ref = paste0("#model_",j,"")
    write(header, file = model_input_file, append = T)


    ## hide the path setting so people can't see your directory structure
    write(paste0("```{r ", j,"_set_path, eval = T, echo = F, message=FALSE, warning=FALSE, comment=FALSE}"), file = model_input_file, append = T)
    header = '
    library(SpatialSablefishAssessment)
    library(tidyverse)
    library(data.table)
    library(ggplot2)
    library(RColorBrewer)
    library(TMB)
    library(reshape2)
    obs_pallete <- c("#E69F00", "#56B4E9","#009E73")
    names(obs_pallete) = c("Observed", "Predicted","Residuals")
    '
    write(header, file = model_input_file, append = T)
    write(paste0('this_model_dir = "', normalizePath(model_dir[j], winslash = "/"),'"'), file = model_input_file, append = T)
    write("```", file = model_input_file, append = T)

    ## read in model
    write(paste0("```{r ", j,"_read_in_info, eval = T, echo = T, warning = F, results = F}"), file = model_input_file, append = T)
    read_in = '
    data = readRDS(file.path(this_model_dir, "data.RDS"))
    region_key = readRDS(file.path(this_model_dir, "region_key.RDS"))
    parameters = readRDS(file.path(this_model_dir, "parameters.RDS"))
    mle_report = readRDS(file.path(this_model_dir, "mle_report.RDS"))
    sd_report = readRDS(file.path(this_model_dir, "sd_report.RDS"))
    mle_spatial = readRDS(file.path(this_model_dir, "mle_optim.RDS"))
    map_fixed_pars = readRDS(file.path(this_model_dir, "map_fixed_pars.RDS"))
    mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", silent = T)
```\n\n
    '
    write(read_in, file = model_input_file, append = T)

    ## Plot model inputs
    write("## Inputs", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_inputs, eval = T, echo = T, warning = F, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    inputs =
      '
    plot_input_observations(data, region_key = region_key)
    plot_input_timeblocks(data)
    plot_input_catches(data, region_key = region_key)
    '
    write(inputs, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Log likelihoods", file = model_input_file, append = T)
    write(paste0("```{r ", j,"log-likelihoods, eval = T, warning = F, echo = T, results = T}"), file = model_input_file, append = T)
    log_like =
      '
    log_like = get_negloglike(mle_report)
    knitr::kable(log_like, digits = 2, caption = "Negative log likelihood values and corresponding distrbutions")
    '
    write(log_like, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)



    write("## Derived quantities", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_init_age, eval = T, echo = T, warning = F, results = T}"), file = model_input_file, append = T)
    init_age =
      '
    plot_init_nage(mle_report, region_key = region_key)
    '
    write(init_age, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_movement, eval = T, echo = T, results = T, warning = F, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    movement =
      '
    plot_movement(mle_report, region_key = region_key)
    '
    write(movement, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_select, eval = T, echo = T, results = T, warning = F, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    select =
      '
    plot_selectivities(mle_report)
    '
    write(select, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_ssb, eval = T, echo = T, results = T, warning = F, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    ssb =
      '
    plot_SSB(mle_report, region_key = region_key, depletion = F)
    '
    write(ssb, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write(paste0("```{r ", j,"_depletion, eval = T, echo = T, warning = F, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    depletion =
      '
    plot_SSB(mle_report, region_key = region_key, depletion = T) +
    ylim(0,NA) +
    geom_hline(yintercept = 40 , col = "darkgreen", linetype = "dashed") +
    geom_hline(yintercept = 20 , col = "orange", linetype = "dashed") +
    geom_hline(yintercept = 10 , col = "red", linetype = "dashed")

    '
    write(depletion, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_recruitment, eval = T, echo = T, warning = F, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    recruitment =
      '
    plot_recruitment(MLE_report = mle_report, region_key = region_key)
    '
    write(recruitment, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_Fs, eval = T, echo = T, results = T, warning = F, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    fish_mort =
      '
    plot_fishing_mortalities(MLE_report = mle_report, region_key = region_key)
    '
    write(fish_mort, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    ## Check MPD convergence
    write("## MPD convergence", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_convergence, eval = T, echo = T, warning = F, results = T}"), file = model_input_file, append = T)
    chol_covar =
      '
    try_decompose = tryCatch(chol(sd_report$cov.fixed), error = function(e){e})
    if(inherits(try_decompose, "error")) {
      message("Failed to decompose the covariance, suggesting non-converged model")
    } else if(any(class(try_decompose) %in% "matrix")) {
       message("Successfully decomposed the covariance, suggesting convergence")
    }
    '
    write(chol_covar, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Other Parameters", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_other_params, eval = T, echo = T, warning = F, results = T}"), file = model_input_file, append = T)
    param_code =
      '
      other_params = get_other_derived_quantities(MLE_report = mle_report, data = data, region_key = region_key)
      other_params$catchabilities$q = other_params$catchabilities$q
      knitr::kable(other_params$catchabilities, row.names = NA, digits = 2, caption = "Estimated survey catchability coeffecients")
      knitr::kable(cbind(names(other_params$scalar_quants), Reduce(c, as.vector(other_params$scalar_quants))), colnames = c("Parameter", "Value"),row.names = NA, digits = 2, caption = "Estimable scalar parameters")
      knitr::kable(other_params$spatial_params, colnames = NA,row.names = NA, digits = 2, caption = "Estimable spatial parameters/derived quantities")
      knitr::kable(other_params$tag_reporting_rate, colnames = NA,row.names = NA, digits = 2, caption = "Tag reporting rate parameters")
    '
    write(param_code, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Reference Points", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_reference_code, eval = T, echo = T, warning = F, results = T}"), file = model_input_file, append = T)
    ref_code =
      '
    n_proj_years = 90
    regional_spr = find_regional_Fspr(data = data, MLE_report = mle_report, n_years_for_fleet_ratio = 5,
                   percent_Bzero = 40, n_future_years = n_proj_years, verbose = F)
    F40 = round(regional_spr$Fspr,3)
    terminal_ssb = mle_report$SSB_yr[length(mle_report$years),]
    terminal_depletion = terminal_ssb / mle_report$Bzero * 100
    B40 = regional_spr$terminal_ssb$terminal_ssb
    Bzero = mle_report$Bzero
    ref_df = data.frame(Region = region_key$area, F40 = F40, B40 = B40, terminal_ssb = terminal_ssb, Bzero = Bzero, terminal_depletion = terminal_depletion)
    knitr::kable(ref_df, col.names = c("Region","F40%","B40%", "Bcurrent", "Bzero","depletion terminal year"), digits = 2)
    '
    write(ref_code, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    ################
    write("## Abundance Fits", file = model_input_file, append = T)
    ## plot Abundance index
    write(paste0("```{r ", j,"_plot_abund, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    plot_abun = '
    plot_index_fit(mle_report, region_key = region_key) + ylim(0,NA) +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            title=element_blank(),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            legend.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 12)) +
      scale_color_manual(values = obs_pallete[1:2])
      '
    write(plot_abun, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_plot-abund-resid, eval = T, echo = F, warning = F, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    plot_abun_resid = '
    index_df = get_index(mle_report, region_key = region_key) %>% rename(Residuals = Pearsons_residuals)
    ggplot(index_df, aes(x = Year)) +
      geom_point(aes(y = Residuals, col = "Residuals"), size = 3) +
      geom_smooth(aes(y = Residuals, col = "Residuals", fill = "Residuals"), size = 1, alpha = 0.4, linetype = "dashed") +
      labs(x = "Year", y = "Pearsons Residuals", fill = "", col = "") +
      geom_hline(yintercept = 0, col = "#999999", linetype = "dashed", size = 1) +
      geom_hline(yintercept = c(-2,2), col = "red", linetype = "dashed", size = 1, alpha = 0.3) +
      ylim(-6,6) +
      theme_classic() +
      facet_wrap(~Region, ncol = 2) +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            title=element_blank(),
            strip.text = element_text(size=14),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[3]) +
      scale_fill_manual(values = obs_pallete[3])
      '
    write(plot_abun_resid, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Aggregated tag fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_tag_agg, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 8}"), file = model_input_file, append = T)

    tag_fit = '
    if(sum(mle_report$tag_recovery_indicator_by_year) != 0 & sum(mle_report$tag_recovery_indicator) != 0) {
      tag_pred = get_tag_recovery_obs_fitted_values(MLE_report = mle_report, region_key = region_key)
      tag_aggregated = tag_pred %>% group_by(release_region, recovery_region) %>% summarise(observed = sum(observed), predicted = sum(predicted))
      tag_aggregated$resid = tag_aggregated$observed - tag_aggregated$predicted
      tag_aggregated$resid_sign = ifelse(tag_aggregated$resid < 0, "negative", "positive")
      tag_aggregated$release_region = factor(tag_aggregated$release_region, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
      tag_aggregated$recovery_region = factor(tag_aggregated$recovery_region, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)

      gplt = ggplot(tag_aggregated, aes(x = release_region, y = recovery_region, fill = resid)) +
        geom_tile() +
        scale_fill_gradient(low = "red", high = "steelblue") +
        geom_text(aes(x = release_region, y = recovery_region, label = round(resid,1)), color = "black", size = 4) +
        labs(x = "Recovery", y = "Release", fill = "Aggregated\nresiduals\n (O - E)")
      gplt
    }
    '

    write(tag_fit, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Catch fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_catch, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 6}"), file = model_input_file, append = T)
    catch_fits = '
    plot_catch_fit(MLE_report = mle_report, region_key = region_key)
    '
    write(catch_fits, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Mean age plots", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_mean_age, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 4.5}"), file = model_input_file, append = T)
    mean_age = '
    plot_mean_age(MLE_report = mle_report, region_key = region_key,observation = "srv_dom_ll",sex = "male")+
      ggtitle("Male mean age survey")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,observation = "srv_dom_ll",sex = "female")+
      ggtitle("Female mean age survey")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,observation = "fixed",sex = "male")+
      ggtitle("Male mean age fixed")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,observation = "fixed",sex = "female") +
    ggtitle("Female mean age fixed")
    '
    write(mean_age, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Mean length plots", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_mean_length, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 4.5}"), file = model_input_file, append = T)
    mean_len = '
    plot_mean_length(MLE_report = mle_report, region_key = region_key,observation = "fixed",sex = "male") +
      ggtitle("Male mean length fixed")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,observation = "fixed",sex = "female") +
      ggtitle("Female mean length fixed")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,observation = "trwl",sex = "male") +
      ggtitle("Male mean length trawl")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,observation = "trwl",sex = "female") +
      ggtitle("Female mean length trawl")
    '
    write(mean_len, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    if(file.exists(file.path(model_dir[j], "sim_obs.RDS"))) {
      write("## Quantile residual plots", file = model_input_file, append = T)
      write(paste0("```{r ", j,"_quantile_resids, eval = T, echo = F, warning = F, results = T, out.width =  '100%', fig.height = 9}"), file = model_input_file, append = T)
      create_resids = '
      sim_obs = readRDS(file.path(this_model_dir, "sim_obs.RDS"))
      full_sim_resid_srv_abundance = calculate_simulated_residuals(sim_ob = sim_obs$sim_srv_bio, type = "abundance")
      full_sim_resid_srv_AF = calculate_simulated_residuals(sim_ob = sim_obs$sim_srv_AF, type = "AF")
      full_sim_resid_fixed_AF = calculate_simulated_residuals(sim_ob = sim_obs$sim_fixed_AF, type = "AF")
      full_sim_resid_fixed_LF = calculate_simulated_residuals(sim_ob = sim_obs$sim_fixed_LF, type = "LF")
      full_sim_resid_trwl_LF = calculate_simulated_residuals(sim_ob = sim_obs$sim_trwl_LF, type = "LF")
      if(sum(mle_report$tag_recovery_indicator_by_year) != 0 & sum(mle_report$tag_recovery_indicator) != 0) {
        full_sim_resid_tag_recovery = calculate_simulated_residuals(sim_ob = sim_obs$sim_tag_recovery, type = "tag-recovery")
      }
      '
      write(create_resids, file = model_input_file, append = T)
      ## plot them
      plot_quant_resids = '
      summarise_AF_quant_resids(AF_sim_resids = full_sim_resid_srv_AF, sex = "M", obs_label = "Survey AF")
      summarise_AF_quant_resids(AF_sim_resids = full_sim_resid_srv_AF, sex = "F", obs_label = "Survey AF")
      summarise_AF_quant_resids(AF_sim_resids = full_sim_resid_fixed_AF, sex = "M", obs_label = "Fixed AF")
      summarise_AF_quant_resids(AF_sim_resids = full_sim_resid_fixed_AF, sex = "F", obs_label = "Fixed AF")
      summarise_LF_quant_resids(LF_sim_resids = full_sim_resid_trwl_LF, sex = "M", obs_label = "Trawl LF")
      summarise_LF_quant_resids(LF_sim_resids = full_sim_resid_trwl_LF, sex = "F", obs_label = "Trawl LF")
      summarise_LF_quant_resids(LF_sim_resids = full_sim_resid_fixed_LF, sex = "M", obs_label = "Fixed LF")
      summarise_LF_quant_resids(LF_sim_resids = full_sim_resid_fixed_LF, sex = "F", obs_label = "Fixed LF")
      if(sum(mle_report$tag_recovery_indicator_by_year) != 0 & sum(mle_report$tag_recovery_indicator) != 0) {
        summarise_tag_quant_resids(tag_sim_resids = full_sim_resid_tag_recovery, T)
        summarise_tag_quant_resids(tag_sim_resids = full_sim_resid_tag_recovery, F)
      }
      '
      write(plot_quant_resids, file = model_input_file, append = T)

      write("```\n\n", file = model_input_file, append = T)
    }

    write("## Age compositional Fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_age_comp, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 10}"), file = model_input_file, append = T)

    age_comp = '
    first_year_set = c(1981, seq(from = 1985, to = 1993, by = 2), 1996:1999)
    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = first_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = first_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")

    second_year_set = 2000:2010
    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = second_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = second_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")

    third_year_set = 2011:2021
    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = third_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "srv_dom_ll", subset_years = third_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")


    second_year_set = 1999:2010
    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = second_year_set, sex = "male") +
      ggtitle("Male fixed AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = second_year_set, sex = "female") +
      ggtitle("Female fixed AF") +
      guides(shape = "none", linetype = "none")

    third_year_set = 2011:2021
    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = third_year_set, sex = "male") +
      ggtitle("Male fixed AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = third_year_set, sex = "female") +
      ggtitle("Female fixed AF") +
      guides(shape = "none", linetype = "none")
      '
    write(age_comp, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Length compositional Fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_len_comp, eval = T, echo = F, results = T, warning = F, out.width =  '100%', fig.height = 10}"), file = model_input_file, append = T)

    len_comp = '

    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = 1991:1998, sex = "male") +
      ggtitle("Male fixed LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "fixed", subset_years = 1991:1998, sex = "female") +
      ggtitle("Female fixed LF") +
      guides(shape = "none", linetype = "none")

    # Trawl gear
    trwl_LF_years = data$years[which(colSums(data$trwl_catchatlgth_indicator) > 0)]
    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "trwl", subset_years = trwl_LF_years[1:12], sex = "male") +
      ggtitle("Male trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "trwl", subset_years = trwl_LF_years[1:12], sex = "female") +
      ggtitle("Female trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "trwl", subset_years = trwl_LF_years[13:length(trwl_LF_years)], sex = "male") +
      ggtitle("Male trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, observation = "trwl", subset_years = trwl_LF_years[13:length(trwl_LF_years)], sex = "female") +
      ggtitle("Female trawl LF") +
      guides(shape = "none", linetype = "none")
    '

    write(len_comp, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)
  }
  ## render or compile the bookdown
  bookdown::render_book(input = output_dir)
}


#' summarise_multiple_models
#' A bookdown function that will summaries a suite of models in individual Sections of a single Rmarkdown file
#' @param model_dir a vector of absolute paths were we can find model output to summarise
#' @param bookdown_dir an absolute path where we will build the Bookdown
#' @param bookdown_labels a vector of model labels
#' @param model_description A string must be single quotes '' that describes all the models preferably in bullet points that are pointed at by model_dir
#' @return Creates a suite of Rmarkdown files in bookdown_dir and compiles an R markdown file
#' @export
#' @details this function expects the following RDS objects to be in each `model_dir` for each model, data.RDS region_key.RDS, parameters.RDS, mle_report.RDS, sd_report.RDS, mle_optim.RDS, map_fixed_pars.RDS

summarise_multiple_models <- function(model_dir, bookdown_dir, bookdown_labels, model_description = "") {
  if(length(model_dir) != length(bookdown_labels))
    stop("number of model_dir provided was not the same as bookdown_labels. These need to be the same length")

  ## set up Bookdown directory
  verbose = T
  output_dir = normalizePath(file.path(bookdown_dir), winslash = "/")
  if(dir.exists(output_dir)) {
    unlink(output_dir, recursive = T, force = T)
    if(verbose)
      print("deleting 'output_dir'")
  }

  dir.create(output_dir)
  if(verbose)
    print(paste0("creating 'output_dir' ", output_dir))
  ## build a skeleton Bookdown
  #bookdown:::bookdown_skeleton(output_dir)
  bookdown::create_bs4_book(output_dir)

  ## get rid of these skeleton changes. This may need tweaking over time
  files_to_remove = c("01-intro.Rmd", "02-cross-refs.Rmd", "03-parts.Rmd", "04-citations.Rmd","05-blocks.Rmd","06-share.Rmd","07-references.Rmd")
  for(i in 1:length(files_to_remove)) {
    file.remove(file.path(output_dir, files_to_remove[i]))
  }

  ## change the CSS
  css_file = file.path(output_dir, "style.css")
  #write("\n.book .book-body .page-wrapper .page-inner {", file = css_file, append = T)
  #write("  max-width: 1200px !important;", file = css_file, append = T)
  #write("}", file = css_file, append = T)

  ## change the _output.yml
  output_yml_file = file.path(output_dir, "_output.yml")
  first_chunk =
    'bookdown::bs4_book:
  css: style.css
  theme:
    primary: "#096B72"
  repo: https://github.com/rstudio/bookdown-demo
bookdown::pdf_book:
  includes:
    in_header: preamble.tex
  latex_engine: xelatex
  citation_package: natbib
  keep_tex: yes
bookdown::epub_book: default
'
  write(first_chunk, file = output_yml_file, append = F)

  ### change index.md

  index_md =
    "---
title: 'Sablefish models'
author: 'C Marsh'
date: '`r Sys.Date()`'
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
#cover-image: path to the social sharing image like images/cover.jpg
description: A detailed summary of sablefish model fits and estimated quantities
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---
"

  tail_index_chunk = "
```{r include=FALSE}
# automatically create a bib database for R packages
knitr::write_bib(c(
  .packages(), 'bookdown', 'knitr', 'rmarkdown'
), 'packages.bib')
```
"
  write(index_md, file = file.path(output_dir, "index.Rmd"), append = F)
  write('# Model descriptions', file = file.path(output_dir, "index.Rmd"), append = T)
  write(model_description, file = file.path(output_dir, "index.Rmd"), append = T)
  write(tail_index_chunk, file = file.path(output_dir, "index.Rmd"), append = T)


  #'
  #' Start comparing model fits
  #'

  model_input_file = file.path(output_dir, "01-Fits.Rmd")
  write("# Model Fits", file = model_input_file, append = F)
  write("```{r read_in_models, eval = T, echo = F, message=FALSE, warning=FALSE, comment=FALSE}", file = model_input_file, append = T)

  header = '
    library(SpatialSablefishAssessment)
    library(tidyverse)
    library(data.table)
    library(ggplot2)
    library(RColorBrewer)
    library(TMB)
    library(reshape2)
    library(viridis)
    mle_ls = list()
    n_pars = nll = vector();

  '
  write(header, file = model_input_file, append = T)
  write("\n\n", file = model_input_file, append = T)
  for(i in 1:length(bookdown_labels)) {
    write(paste0('this_model_dir = "', normalizePath(model_dir[i], winslash = "/"),'"'), file = model_input_file, append = T)
    write('mle_report = readRDS(file.path(this_model_dir, "mle_report.RDS"))', file = model_input_file, append = T)
    write(paste0('mle_ls[["',bookdown_labels[i],'"]] = mle_report'), file = model_input_file, append = T)

    write('mle_optim = readRDS(file.path(this_model_dir, "mle_optim.RDS"))', file = model_input_file, append = T)
    write(paste0('n_pars["',bookdown_labels[i],'"] = length(mle_optim$par)'), file = model_input_file, append = T)
    write(paste0('nll["',bookdown_labels[i],'"] = mle_optim$objective'), file = model_input_file, append = T)

  }
  write('region_key = readRDS(file.path(this_model_dir, "region_key.RDS"))', file = model_input_file, append = T)

  pallette_rmd = '
  n_pars_df = data.frame(label = names(n_pars), n_pars = n_pars)
  obs_pallete <- c("black", brewer.pal(length(mle_ls), "Set1"))
  names(obs_pallete) = c("Observed",names(mle_ls))
  '
  write(pallette_rmd, file = model_input_file, append = T)
  write("```", file = model_input_file, append = T)

  write("## Catch fits", file = model_input_file, append = T)
  write("```{r catch_fits, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  catch_fits =
    '
    catch_fit_df = get_multiple_catch_fits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)

    ggplot() +
      geom_line(data = catch_fit_df %>% dplyr::filter(type == "Predicted"), aes(x = Year, y = Catch, col = label, linetype = label), linewidth = 1) +
      geom_point(data = catch_fit_df %>% dplyr::filter(type == "Observed"), aes(x = Year, y = Catch, col = "Observed")) +
      labs(x = "Year", y = "Catch", col = "", linetype= "") +
      facet_grid(Region~Fishery) +
      ggtitle("Catch fits") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete) +
      guides(linetype = "none")

    '
  write(catch_fits, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)


  write("## Abundance fits", file = model_input_file, append = T)
  write("```{r abund_fits, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  abundance_fits =
    '
    index_fit_df = get_multiple_index_fits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)

    ggplot() +
      geom_line(data = index_fit_df, aes(x = Year, y = Predicted, col = label, linetype = label), linewidth = 1) +
      geom_point(data = index_fit_df, aes(x = Year, y = Observed, col = "Observed")) +
      geom_errorbar(data = index_fit_df, aes(x = Year, ymin=L_CI, ymax=U_CI, col = "Observed"), width=.2, position=position_dodge(.9)) +
      labs(x = "Year", y = "Survey", col = "", linetype= "") +
      facet_wrap(~Region) +
      ggtitle("Survey fits") +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 12),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete) +
      guides(linetype = "none")

    '
  write(abundance_fits, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)



  write("## Mean Age", file = model_input_file, append = T)
  write("```{r get_comp_obs, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  get_comps =
    '
        mean_age_df = get_multiple_mean_age_fits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
    '
  write(get_comps, file = model_input_file, append = T)

  plt_fixed_age =
    '
ggplot(data = mean_age_df %>% dplyr::filter(observation == "fixed")) +
  geom_point(aes(x = Year, y = Oy, col = "Observed"), size = 1.6) +
  geom_errorbar(aes(x = Year, ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
  geom_line(aes(x = Year, y = Ey, col = label, linetype = label)) +
  labs(x = "Year", y = "Mean age", col = "", linetype= "", title = "Fixed gear") +
  facet_grid(Region~Sex) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=14)) +
  scale_color_manual(values = obs_pallete)
  '
  write(plt_fixed_age, file = model_input_file, append = T)

  plt_srv_age =
    '
ggplot(data = mean_age_df %>% dplyr::filter(observation == "srv_dom_ll")) +
  geom_point(aes(x = Year, y = Oy, col = "Observed"), size = 1.6) +
  geom_errorbar(aes(x = Year, ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
  geom_line(aes(x = Year, y = Ey, col = label, linetype = label)) +
  labs(x = "Year", y = "Mean age", col = "", linetype= "", title = "Survey") +
  facet_grid(Region~Sex) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=14)) +
  scale_color_manual(values = obs_pallete)
  '
  write(plt_srv_age, file = model_input_file, append = T)

  write("```\n\n", file = model_input_file, append = T)


  write("## Mean Length", file = model_input_file, append = T)
  write("```{r get_len_comp_obs, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  get_comps =
    '
        mean_length_df = get_multiple_mean_length_fits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
    '
  write(get_comps, file = model_input_file, append = T)

  plt_fixed_len =
    '
ggplot(data = mean_length_df %>% dplyr::filter(observation == "fixed")) +
  geom_point(aes(x = Year, y = Oy, col = "Observed"), size = 1.6) +
  geom_errorbar(aes(x = Year, ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
  geom_point(aes(x = Year, y = Ey, col = label, shape = "Predicted"), size = 1.6) +
  geom_line(aes(x = Year, y = Ey, col = label, linetype = label)) +
  labs(x = "Year", y = "Mean length", col = "", linetype= "", title = "Fixed gear") +
  facet_grid(Region~Sex) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=14)) +
  scale_color_manual(values = obs_pallete)
  '
  write(plt_fixed_len, file = model_input_file, append = T)

  plt_trwl_len =
    '
ggplot(data = mean_length_df %>% dplyr::filter(observation == "trwl")) +
  geom_point(aes(x = Year, y = Oy, col = "Observed", shape = "Observed"), size = 1.6) +
  geom_point(aes(x = Year, y = Ey, col = label, shape = "Predicted"), size = 1.6) +
  geom_errorbar(aes(x = Year, ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
  geom_line(aes(x = Year, y = Ey, col = label, linetype = label)) +
  labs(x = "Year", y = "Mean length", col = "", linetype= "", title = "Trawl") +
  facet_grid(Region~Sex) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=14)) +
  scale_color_manual(values = obs_pallete)
  '
  write(plt_trwl_len, file = model_input_file, append = T)

  write("```\n\n", file = model_input_file, append = T)

  #'
  #' Start comparing model quantities
  #'

  model_input_file = file.path(output_dir, "02-DerivedQuantities.Rmd")
  write("# Derived Quantities", file = model_input_file, append = F)
  write("```{r read_in_models_dqs, eval = T, echo = F, message=FALSE, warning=FALSE, comment=FALSE}", file = model_input_file, append = T)

  header = '
    library(SpatialSablefishAssessment)
    library(tidyverse)
    library(data.table)
    library(ggplot2)
    library(RColorBrewer)
    library(TMB)
    library(reshape2)
    library(viridis)
    mle_ls = list()
    n_pars = nll = vector();
    data_ls = list();

  '
  write(header, file = model_input_file, append = T)
  write("\n\n", file = model_input_file, append = T)

  for(i in 1:length(bookdown_labels)) {
    write(paste0('this_model_dir = "', normalizePath(model_dir[i], winslash = "/"),'"'), file = model_input_file, append = T)
    write('mle_report = readRDS(file.path(this_model_dir, "mle_report.RDS"))', file = model_input_file, append = T)
    write(paste0('mle_ls[["',bookdown_labels[i],'"]] = mle_report'), file = model_input_file, append = T)

    write('data_df = readRDS(file.path(this_model_dir, "data.RDS"))', file = model_input_file, append = T)
    write(paste0('data_ls[["',bookdown_labels[i],'"]] = data_df'), file = model_input_file, append = T)


    write('mle_optim = readRDS(file.path(this_model_dir, "mle_optim.RDS"))', file = model_input_file, append = T)
    write(paste0('n_pars["',bookdown_labels[i],'"] = length(mle_optim$par)'), file = model_input_file, append = T)
    write(paste0('nll["',bookdown_labels[i],'"] = mle_optim$objective'), file = model_input_file, append = T)

  }
  write('region_key = readRDS(file.path(this_model_dir, "region_key.RDS"))', file = model_input_file, append = T)

  pallette_rmd = '
    n_pars_df = data.frame(label = names(n_pars), n_pars = n_pars)
    nll_df = data.frame(label = names(n_pars), n_pars = n_pars)

  obs_pallete <- c("black", brewer.pal(length(mle_ls), "Set1"))
  names(obs_pallete) = c("Observed",names(mle_ls))
  '
  write(pallette_rmd, file = model_input_file, append = T)
  write("```", file = model_input_file, append = T)

  write("## SSBs", file = model_input_file, append = T)
  write("```{r ssbs, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  ssbs =
    '
    ssb_df = get_multiple_ssbs(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key, depletion = F)
    depletion_df = get_multiple_ssbs(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key, depletion = T)

    ggplot() +
      geom_line(data = ssb_df, aes(x = Year, y = SSB, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "SSB (kt)", col = "", linetype= "") +
      facet_wrap(~Region) +
      ylim(0,NA)+
      ggtitle("Female spawning stock biomass") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(ssbs, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)


  write("```{r depletion, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  deplet =
    '

    ggplot() +
      geom_line(data = depletion_df, aes(x = Year, y = Depletion, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "Depletion (SSB/B0)", col = "", linetype= "") +
      facet_wrap(~Region) +
      ylim(0,NA)+
      ggtitle("SSB depletion") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(deplet, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)


  write("## Recruitment", file = model_input_file, append = T)

  write("```{r recuitment, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  recruit =
    '
    recruit_df = get_multiple_recruits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
    ggplot() +
      geom_line(data = recruit_df, aes(x = Year, y = Recruitment, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "Recruitment (millions)", col = "", linetype= "") +
      facet_wrap(~Region) +
      ylim(0,NA)+
      ggtitle("Recruitment") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(recruit, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)

  write("## Recruitment YCS", file = model_input_file, append = T)

  write("```{r recuitmentYCS, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  recruit_ycs =
    '
    recruit_df = get_multiple_recruits(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
    ggplot() +
      geom_line(data = recruit_df, aes(x = Year, y = YCS, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "Year class strength multipliers (YCS)", col = "", linetype= "") +
      facet_wrap(~Region) +
      ylim(0,NA)+
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.1) +
      ggtitle("Year class strength multipliers (YCS)") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(recruit_ycs, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)

  write("## Initial numbers at age", file = model_input_file, append = T)
  write("```{r Initnage, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  init_nage =
    '
    init_nage = get_multiple_init_nage(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key);
    ggplot() +
      geom_line(data = init_nage, aes(x = Age, y = Numbers, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "Initial numbers at age", col = "", linetype= "") +
      ylim(0,NA)+
      facet_wrap(~Region) +
      ggtitle("Initial numbers at age") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")
  '
  write(init_nage, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)

  write("## Fishing Mortality", file = model_input_file, append = T)

  write("```{r fishMort, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  fish_mort =
    '
    fish_mort_df = get_multiple_Fs(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
    ggplot() +
      geom_line(data = fish_mort_df, aes(x = Year, y = F, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "Annual fishing mortality (F)", col = "", linetype= "") +
      ylim(0,NA)+
      facet_wrap(Fishery~Region) +
      ggtitle("Fishing Mortality") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(fish_mort, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)



  write("## Movement", file = model_input_file, append = T)
  write("```{r move, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 10}", file = model_input_file, append = T)
  move_rmd =
    '
    move_df = get_multiple_movements(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)

    ggplot(move_df, aes(x = To, y = From, fill = Proportion)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      geom_text(aes(x = To, y = From, label = round(Proportion,2)), color = "black", size = 4) +
      labs(x = "To", y = "From") +
      facet_wrap(~label, ncol = 2)
    '
  write(move_rmd, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)

  write("## Selectivities", file = model_input_file, append = T)
  write("```{r select, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
  select =
    '
    select_df = get_multiple_selectivities(mle_ls = mle_ls, run_labels = names(mle_ls))
    ggplot() +
      geom_line(data = select_df, aes(x = age, y = value, col = label, linetype = label), linewidth = 1) +
      labs(x = "Year", y = "", col = "", linetype= "") +
      facet_wrap(sex~gear) +
      ggtitle("Selectivities") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size=14),
            plot.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size=14)) +
      scale_color_manual(values = obs_pallete[-1]) +
      guides(linetype = "none")

    '
  write(select, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)


    write("## Tag reporting rates", file = model_input_file, append = T)
    write("```{r tagreportrates, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}", file = model_input_file, append = T)
    report_rate =
      '
      report_df = get_multiple_tag_reporting_rates(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = region_key)
      if(!is.null(report_df)) {
        ggplot() +
          geom_line(data = report_df, aes(x = Year, y = ReportingRate, col = label, linetype = label), linewidth = 1) +
          labs(x = "Tag recovery year", y = "", col = "", linetype= "") +
          facet_wrap(~Region) +
          ggtitle("Tag reporting rates") +
          ylim(0,NA)+
          theme_bw() +
          theme(legend.position = "bottom",
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 14),
                strip.text = element_text(size=14),
                plot.title = element_text(size = 20, face = "bold"),
                legend.text = element_text(size=14)) +
          scale_color_manual(values = obs_pallete[-1]) +
          guides(linetype = "none")
      }

      '
    write(report_rate, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


  write("# Other Parameters", file = model_input_file, append = T)
  write("```{r otherParams, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8, message = F}", file = model_input_file, append = T)
  param_code =
    '
      q_df = scalar_df = spatial_quants_df = tag_reporting_df = NULL
      for(i in 1:length(mle_ls)) {
        other_params = get_other_derived_quantities(MLE_report = mle_ls[[i]], data = data_ls[[i]], region_key = region_key)
        other_params$catchabilities$q = other_params$catchabilities$q
        other_params$catchabilities$model = names(mle_ls)[i]
        temp_scalar_df = data.frame(values = Reduce(c, as.vector(other_params$scalar_quants[1,])), parameter = colnames(other_params$scalar_quants), model = names(mle_ls)[i])
        temp_tag_report_df = data.frame(time_block = paste0("block-", 1:length(other_params$tag_reporting_rate)), value = other_params$tag_reporting_rate, model = names(mle_ls)[i])

        tag_reporting_df = rbind(tag_reporting_df, temp_tag_report_df)
        other_params$spatial_params$model = names(mle_ls)[i]
        spatial_quants_df = rbind(spatial_quants_df, other_params$spatial_params)
        scalar_df = rbind(scalar_df, temp_scalar_df)
        q_df = rbind(q_df, other_params$catchabilities)
      }
      knitr::kable(q_df %>% pivot_wider(id_cols = c("time-block","Region"), values_from = q, names_from = model), row.names = NA, digits = 2, caption = "Estimated survey catchability coefficients")
      knitr::kable(scalar_df %>% pivot_wider(id_cols = c("parameter"), values_from = values, names_from= model), row.names = NA, digits = 2, caption = "Estimable scalar parameters")
      knitr::kable(tag_reporting_df %>% pivot_wider(id_cols = time_block, values_from = value, names_from = model), colnames = NA,row.names = NA, digits = 2, caption = "Estimated tag-reporting rates")
      knitr::kable(spatial_quants_df %>% select(Region, Binit, model) %>% pivot_wider(id_cols = Region, values_from = Binit, names_from = model), colnames = NA,row.names = NA, digits = 2, caption = "Binitial (after accounting for Finit and initial age deviations)")

    '
  write(param_code, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)

  write("# Reference Points", file = model_input_file, append = T)
  write("```{r RefPoints, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8, message = F}", file = model_input_file, append = T)
  ref_rmd =
    '
    n_proj_years = 90
    full_ref_df = NULL
    for(i in 1:length(mle_ls)) {
      regional_spr = find_regional_Fspr(data = data_ls[[i]], MLE_report = mle_ls[[i]], n_years_for_fleet_ratio = 5,
                     percent_Bzero = 40, n_future_years = n_proj_years, verbose = F)
      F40 = round(regional_spr$Fspr,3)
      terminal_ssb = mle_ls[[i]]$SSB_yr[length(mle_ls[[i]]$years),]
      terminal_depletion = terminal_ssb / mle_ls[[i]]$Bzero * 100
      B40 = regional_spr$terminal_ssb$terminal_ssb
      Bzero = mle_ls[[i]]$Bzero
      ref_df = data.frame(Region = region_key$area, F40 = F40, B40 = B40, terminal_ssb = terminal_ssb, Bzero = Bzero, terminal_depletion = terminal_depletion)
      ref_df$model = names(mle_ls)[i]
      full_ref_df = rbind(full_ref_df, ref_df)
    }
    knitr::kable(full_ref_df %>% select(Region, Bzero, model) %>% pivot_wider(id_cols = Region, names_from = model, values_from = Bzero), digits = 2, caption  = "Bzero estimates by region and model")
    knitr::kable(full_ref_df %>% select(Region, B40, model) %>% pivot_wider(id_cols = Region, names_from = model, values_from = B40), digits = 2, caption  = "B40% estimates by region and model")
    knitr::kable(full_ref_df %>% select(Region, F40, model) %>% pivot_wider(id_cols = Region, names_from = model, values_from = F40), digits = 2, caption  = "F40% estimates by region and model")
    knitr::kable(full_ref_df %>% select(Region, terminal_depletion, model) %>% pivot_wider(id_cols = Region, names_from = model, values_from = terminal_depletion), digits = 2, caption  = "Depletion estimates by region and model for the terminal year")
    knitr::kable(full_ref_df %>% select(Region, Bzero, model) %>% group_by(model) %>% summarise(global_Bzero = sum(Bzero)), digits = 2, caption  = "Global Bzero (summed over all regions) estimates by model")
  '
  write(ref_rmd, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)



  write("# Likelihoods", file = model_input_file, append = T)
  write("```{r likelihoods, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8, message = F}", file = model_input_file, append = T)
  nll_rmd =
    '
    n_datasets = get_multiple_input_datasets(data_ls, run_labels = names(data_ls), region_key) %>% group_by(label) %>% summarise(n_obs = sum(indicator, na.rm =T))
    nll_df = get_multiple_nlls(mle_ls = mle_ls, run_labels = names(mle_ls), region_key = n_datasets)
    nll_df = nll_df %>% dplyr::inner_join(n_pars_df )  %>% dplyr::inner_join(n_datasets )
    summarised_nll_df = nll_df %>% filter(observations == "Total") %>% group_by(label) %>% mutate(negloglike  = negloglike, AIC = 2*n_pars + 2*negloglike + 2*n_pars*(n_pars+1)/(n_obs-n_pars-1), BIC = 2 * negloglike + n_pars * log(n_obs))
    knitr::kable((summarised_nll_df %>% select(!c("observations", "distribution")))[order(summarised_nll_df$AIC),], digits = 2, caption = "Model summary comparison among models")
    knitr::kable(nll_df %>% pivot_wider(, id_cols = label, values_from = negloglike, names_from = observations), digits = 2, caption = "Negative loglikelihoods by observation and model.")
  '
  write(nll_rmd, file = model_input_file, append = T)
  write("```\n\n", file = model_input_file, append = T)


  bookdown::render_book(input = output_dir)

}



