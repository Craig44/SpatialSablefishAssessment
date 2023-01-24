#' record_grooming_rule
#' @details a utility function to record the effect of grooming rules on data sets. It will note remove the index from the dataset. TODO: As well as recording the effect of a groomong rule on catch and number of records summarise spatail distribution.
#'
#' @param rule the label for the grooming rule (string)
#' @param df the df we are manipulating
#' @param index a conditional index length(index) == nrow(df). TRUE means it stays, FALSE means it is being excluded from
#' @param catch.col column label for the catch variable
#' @param year.col column label for the year variable. Can be ommited, if supplied the Data.frame will break down catch and number of events removed from each grooming rule by year.
#' @param record a data.frame to append the summaries to.
#' @param keep_NA if index has NA's and this is true then they stay in the data
#' @param attribute if by year specify the attribute to record
#' \itemize{
#'   \item catch absoulte catch
#'   \item events number of events
#'   \item relative_catch is the percent change in catch from applying the rule
#'   \item relative_events is the percent change in number of events from applying the rule
#'   \item cumulative_catch track the catch left in the data set after each grooming rule
#'   \item cumulative_events track the events left in the data set after each grooming rule
#' }
#' @export
#' @return the record data frame with new entries for the grooming rule
#' @examples
#'\dontrun{
#'## example of plotting this output
#'## plot records
#'melt_total = melt(total_grooming_record, id.vars = "rule")
#'melt_total$rule = factor(melt_total$rule, ordered = T, levels = total_grooming_record$rule)
#' ggplot(melt_total %>% filter(variable %in%  c("events", "catch")),
#'       aes(y = value /1000, x = rule, group = 1)) +
#'  geom_line(size = 2, linetype = "dotted") +
#'  geom_point(size = 4, aes(col = rule)) +
#'  theme(axis.text.x = element_blank()) +
#'  ylab("") +
#'  ylim(0,1600) +
#'  facet_wrap(~variable) +
#'  xlab("")
#'}
record_grooming_rule = function(df, index, catch.col, record = NULL, rule, year.col = NULL, attribute = "catch", keep_NA = T) {
  use_year = F # record values by year
  first_record = F # first record
  years = NULL;
  attributes_allowed = c("catch", "events", "relative_catch", "relative_events", "cumulative_catch","cumulative_events")
  attr_ndx = which(attribute == attributes_allowed)
  ## business rules
  if(!attribute %in% attributes_allowed)
    stop(paste0("attribute must be one of: ", paste(attributes_allowed, collapse = ", ")))
  if(length(index) != nrow(df))
    stop(paste0("length of index must equal rows in dataframe"))
  if(!catch.col %in% colnames(df)) {
    stop(paste0("Could not find the catch.col colname ", catch.col, " in colnames of dataframe"))
  }
  if(sum(is.na(index)) >= 1) {
    cat(paste0("Found ", sum(is.na(index)), " NA's in your index "))
    if(keep_NA) {
      index[is.na(index)] = TRUE
      cat("converting them to True\n")
    } else {
      index[is.na(index)] = FALSE
      cat("converting them to False\n")
    }
  }


  if(!is.null(year.col)) {
    if(!year.col %in% colnames(df)) {
      stop(paste0("Could not find the year.col colname ", year.col, " in colnames of dataframe"))
    }
    use_year = T
    if(is.null(record)) {
      years = sort(unique(df[, year.col]))
    } else {
      years = colnames(record)
    }
  }
  if(is.null(record)) {
    first_record = T
  } else {
    if(!is.data.frame(record))
      stop("record needs to be a data.frame if supplied")
    ## check it has expected format
    if(use_year) {
      if(ncol(record) != (length(record)))
        stop(paste0("If you supply the parameter 'year.col' and record, there needs to be a column in record for each year. There were ", ncol(record), " columns in record and ", length(years), " years in the dataframe"))
    } else {
      if(all(colnames(record) == c("rule", attributes_allowed)))
        stop(paste0("If you supply record, it needs 6 columns one for: ", paste(c("rule", attributes_allowed), collapse = ", ")))
    }
  }
  if(use_year) {
    #init calcs
    catch_df = df %>% group_by(!!sym(year.col)) %>% summarise(catch = sum(!!sym(catch.col), na.rm = T))
    event_df = df %>% group_by(!!sym(year.col)) %>% summarise(events = n())
    catch_df = catch_df %>% mutate_if(is.character, factor, levels = years)
    event_df = event_df %>% mutate_if(is.character, factor, levels = years)

    catch_df = catch_df %>% complete(!!sym(year.col), fill = list(catch = 0.0)) %>%
      pivot_wider(names_from = year.col, values_from = !year.col)
    event_df = event_df %>% complete(!!sym(year.col), fill = list(events = 0.0)) %>%
      pivot_wider(names_from = year.col, values_from = !year.col)

    # apply rule
    new_df = subset(df, index)
    sub_catch_df = new_df %>% group_by(!!sym(year.col)) %>% summarise(catch = sum(!!sym(catch.col), na.rm = T))
    sub_event_df = new_df %>% group_by(!!sym(year.col)) %>% summarise(events = n())
    sub_catch_df = sub_catch_df %>% mutate_if(is.character, factor, levels = years)
    sub_event_df = sub_event_df %>% mutate_if(is.character, factor, levels = years)

    sub_catch_df = sub_catch_df %>% complete(!!sym(year.col), fill = list(catch = 0.0)) %>%
      pivot_wider(names_from = year.col, values_from = !year.col)
    sub_event_df = sub_event_df %>% complete(!!sym(year.col), fill = list(events = 0.0)) %>%
      pivot_wider(names_from = year.col, values_from = !year.col)

    new_record = NULL
    if(attr_ndx == 1) {
      new_record = round(catch_df - sub_catch_df, 2)
    } else if(attr_ndx == 2) {
      new_record = round(event_df - sub_event_df, 2)
    } else if(attr_ndx == 3) {
      new_record = round((sub_catch_df - catch_df) / catch_df * 100, 2)
    } else if(attr_ndx == 4) {
      new_record = round((sub_event_df - event_df) / event_df * 100, 2)
    } else if(attr_ndx == 5) {
      new_record = sub_event_df
    } else  {
      new_record = sub_catch_df
    }

    if(first_record) {

      record = new_record
      rownames(record) = rule
      return(record);
    } else {
      ## get a temp value from the first entry
      rowlabs = rownames(record)
      record =rbind(record, new_record)
      rowlabs = c(rowlabs, rule)
      if(any(duplicated(rowlabs)))
        stop("found duplicate rule labels, this isn't allowed")
      row.names(record) = rowlabs;
      return(record);
    }
  } else {
    #init calcs
    catch = sum(df[, catch.col], na.rm = T)
    events = nrow(df)
    # apply rule
    new_df = subset(df, index)
    upd_catch = sum(new_df[, catch.col], na.rm = T)
    upd_events = nrow(new_df)

    catch_percent = round((upd_catch - catch) / catch * 100, 2)
    event_percent = round((upd_events - events) / events * 100, 2)
    catch_removed = round(catch - upd_catch, 2)
    events_removed = round(events - upd_events, 2)
    new_record = NULL
    if(first_record) {
      new_record = data.frame(rule = rule, events = upd_events, catch = upd_catch, catch_removed = catch_removed, relative_catch = catch_percent, events_removed = events_removed, relative_events = event_percent)
    } else {
      ## get a temp value from the first entry
      new_record = record[1,]
      new_record$rule = rule
      new_record$catch_removed = catch_removed
      new_record$relative_catch = catch_percent
      new_record$events_removed = events_removed
      new_record$relative_events = event_percent
      new_record$events = upd_events
      new_record$catch = upd_catch
    }
    record = rbind(record, new_record)
    return(record)
  }
  return(FALSE);
}

#' apply_grooming_rule
#' @details will remove rows from a data frame based on the conditional vector index.
#' @param df the df we are manipulating
#' @param index a conditional index length(index) == nrow(df). TRUE means it stays, FALSE means it is being excluded from
#' @param keep_NA if index has NA's and this is true then they stay in the data
#' @export
#' @return a modified version of df with the FALSE values from 'index' removed
#'
apply_grooming_rule = function(df, index, keep_NA = TRUE) {
  ## business rules
  if(length(index) != nrow(df))
    stop(paste0("length of index must equal rows in dataframe"))
  if(sum(is.na(index)) >= 1) {
    cat(paste0("Found ", sum(is.na(index)), " NA's in your index "))
    if(keep_NA) {
      index[is.na(index)] = TRUE
      cat("converting them to True\n")
    } else {
      index[is.na(index)] = FALSE
      cat("converting them to False\n")
    }
  }
  if(sum(!index) == nrow(df)) {
    stop(paste0("All values in index are FALSE, this indicates an erroroneous index vector"))
  }
  new_df = subset(df, subset = index)
  return(new_df)
}
