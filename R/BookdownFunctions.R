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
#' @export

summarise_individual_models <- function(model_dir, bookdown_dir, bookdown_labels, model_description = "") {
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
    write(paste0("```{r ", j,"_read_in_info, eval = T, echo = T, results = F}"), file = model_input_file, append = T)
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
    write(paste0("```{r ", j,"_inputs, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    inputs =
      '
    plot_input_observations(data, region_key = region_key)
    plot_input_timeblocks(data)
    plot_input_catches(data, region_key = region_key)
    '
    write(inputs, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Derived quantities", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_init_age, eval = T, echo = T, results = T}"), file = model_input_file, append = T)
    init_age =
      '
    plot_init_nage(mle_report, region_key = region_key)
    '
    write(init_age, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_movement, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    movement =
      '
    plot_movement(mle_report, region_key = region_key)
    '
    write(movement, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_select, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    select =
      '
    plot_selectivities(mle_report)
    '
    write(select, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_ssb, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    ssb =
      '
    plot_SSB(mle_report, region_key = region_key, depletion = F)
    '
    write(ssb, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write(paste0("```{r ", j,"_depletion, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
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

    write(paste0("```{r ", j,"_recruitment, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    recruitment =
      '
    plot_recruitment(MLE_report = mle_report, region_key = region_key)
    '
    write(recruitment, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write(paste0("```{r ", j,"_Fs, eval = T, echo = T, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
    fish_mort =
      '
    plot_fishing_mortalities(MLE_report = mle_report, region_key = region_key)
    '
    write(fish_mort, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    ## Check MPD convergence
    write("## MPD convergence", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_convergence, eval = T, echo = T, results = T}"), file = model_input_file, append = T)
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



    write("## Reference Points", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_reference_code, eval = T, echo = T, results = T}"), file = model_input_file, append = T)
    ref_code =
      '
    n_proj_years = 90
    regional_spr = find_regional_Fspr(data = data, MLE_report = mle_report, n_years_for_fleet_ratio = 5,
                   percent_Bzero = 40, n_future_years = n_proj_years, verbose = F)
    F40 = round(regional_spr$Fspr,3)
    names(F40) = region_key$area
    knitr::kable(F40, col.names = "F40%")
    '
    write(ref_code, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    ################
    write("## Abundance Fits", file = model_input_file, append = T)
    ## plot Abundance index
    write(paste0("```{r ", j,"_plot_abund, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 7}"), file = model_input_file, append = T)
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

    write(paste0("```{r ", j,"_plot-abund-resid, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}"), file = model_input_file, append = T)
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


    write("## Tag Fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_tag_agg, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 10}"), file = model_input_file, append = T)

    tag_fit = '
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
    '

    write(tag_fit, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Catch fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_catch, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 8}"), file = model_input_file, append = T)
    catch_fits = '
    plot_catch_fit(MLE_report = mle_report, region_key = region_key)
    '
    write(catch_fits, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Mean age plots", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_mean_age, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 6}"), file = model_input_file, append = T)
    mean_age = '
    plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "srv_dom_ll",sex = "male")+
      ggtitle("Male mean age survey")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "srv_dom_ll",sex = "female")+
      ggtitle("Female mean age survey")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male")+
      ggtitle("Male mean age fixed")
    plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female") +
    ggtitle("Female mean age fixed")
    '
    write(mean_age, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Mean length plots", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_mean_length, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 6}"), file = model_input_file, append = T)
    mean_len = '
    plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male") +
      ggtitle("Male mean length fixed")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female") +
      ggtitle("Female mean length fixed")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "male") +
      ggtitle("Male mean length trawl")
    plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "female") +
      ggtitle("Female mean length trawl")
    '
    write(mean_len, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)

    write("## Age compositional Fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_age_comp, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 10}"), file = model_input_file, append = T)

    age_comp = '
    first_year_set = c(1981, seq(from = 1985, to = 1993, by = 2), 1996:1999)
    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")

    second_year_set = 2000:2010
    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = second_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = second_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")

    third_year_set = 2011:2021
    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = third_year_set, sex = "male") +
      ggtitle("Male survey AF") +
      guides(shape = "none", linetype = "none")

    plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = third_year_set, sex = "female") +
      ggtitle("Female survey AF") +
      guides(shape = "none", linetype = "none")
      '
    write(age_comp, file = model_input_file, append = T)
    write("```\n\n", file = model_input_file, append = T)


    write("## Length compositional Fits", file = model_input_file, append = T)
    write(paste0("```{r ", j,"_len_comp, eval = T, echo = F, results = T, out.width =  '100%', fig.height = 10}"), file = model_input_file, append = T)

    len_comp = '

    plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1998, sex = "male") +
      ggtitle("Male fixed LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1998, sex = "female") +
      ggtitle("Female fixed LF") +
      guides(shape = "none", linetype = "none")

    # Trawl gear
    trwl_LF_years = data$years[which(colSums(data$trwl_catchatlgth_indicator) > 0)]
    plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = trwl_LF_years[1:12], sex = "male") +
      ggtitle("Male trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = trwl_LF_years[1:12], sex = "female") +
      ggtitle("Female trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = trwl_LF_years[13:length(trwl_LF_years)], sex = "male") +
      ggtitle("Male trawl LF") +
      guides(shape = "none", linetype = "none")

    plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = trwl_LF_years[13:length(trwl_LF_years)], sex = "female") +
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
#' @param model_description A string that describes all the models that are pointed at by model_dir
#' @return Creates a suite of Rmarkdown files in bookdown_dir and compiles an R markdown file
#' @export

summarise_multiple_models <- function(model_dir, bookdown_dir, bookdown_labels, model_description = "") {
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

  read_all_mods = '
  mle_ls = list()
  for(j in 1:length(bookdown_labels)) {
    mle_ls[[bookdown_labels[j]]] = readRDS(file.path(normalizePath(model_dir[j], winslash = "/"), "mle_report.RDS"))
  }
'

}

