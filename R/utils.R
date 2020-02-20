#' Get Single Row Guanrank
#' Calculates the Guanrank for a single row.
#'
#' @param row vector containing named strings for "event," "n," "time"
#' @param sf the survfit object
#' @param censored original data frame filtered for censored (no event) samples
#' @param uncensored original data frame filtered for uncensored (no event)
#' samples
#'
#' @return Guanrank (non-normalized) for the row
#'
#' @importFrom magrittr %>%
get_rank <- function(row, sf, censored, uncensored) {
  row_event <- as.logical(row["event"])
  row_time <- as.numeric(row["time"])
  row_surv <- sf$surv[sf$time == row_time]
  if (!row_event) {
    # Calculate score for a censored subject
    # Assuming that for all T > t means for all events with time greater than
    # t. Not just for all unique times (could be fewer than events)
    # sum(1 - get_unc_df$surv/row_surv) is equivalent to the 1st + 2nd sum in
    # the equation for the censored hazard rank
    geq_unc_df <- uncensored %>% dplyr::filter(time >= row_time)
    geq_unc <- ifelse(nrow(geq_unc_df) == 0, 0,
      sum(1 - geq_unc_df$surv / row_surv)
    )

    # sum(1 - 0.5*geq_cen_df$surv/row_surv) is equivalent to the 3rd sum in the
    # equation for the censored hazard rank
    geq_cen_df <- censored %>% dplyr::filter(time >= row_time)
    geq_cen <- ifelse(nrow(geq_cen_df) == 0, 0,
      sum(1 - 0.5 * geq_cen_df$surv / row_surv)
    )

    # Equivalent to the 4th sum
    lt_cen_df <- censored %>% dplyr::filter(time < row_time)
    lt_cen <- ifelse(nrow(lt_cen_df) == 0, 0, 
                     sum(0.5 * row_surv / lt_cen_df$surv))
    return(sum(geq_unc, geq_cen, lt_cen, -0.5))
  }
  else {
    # Calculate score for an uncensored subject
    # Equivalent to the 3rd sum in the hazard rank equation for uncensored
    # subjects
    gt_tot <- (uncensored %>% dplyr::filter(time > row_time) %>% nrow()) +
      (censored %>% dplyr::filter(time > row_time) %>% nrow())

    # Equivalent to the 1st sum in the hazard rank equation for uncensored
    # subjects
    eq_unc <- 0.5 * nrow(uncensored %>% dplyr::filter(time == row_time))

    # Equivalent to the 2nd and 4th sums in the hazard rank equation for
    # uncensored subjects
    leq_cen_df <- censored %>% dplyr::filter(time <= row_time)
    leq_cen <- ifelse(nrow(leq_cen_df) == 0, 0, sum(row_surv / leq_cen_df$surv))
    return(sum(eq_unc, gt_tot, leq_cen, -0.5))
  }
}

#' Apply Guanrank to Vectors or Data Frame
#'
#' @param time vector of numeric survival times (possibly censored). Will be
#' attempted to coerce to numeric via as.numeric
#' @param event vector of logicals indicating whether or not an event occurred
#' (if no event occurred, the sample is considered right-censored)
#'
#' @return vector of numeric Guanrank values, normalized to [0,1], with higher
#' value indicating higher hazard
#' @export
#' @importFrom magrittr %>%
#' @examples
#' time <- c(93, 169, 922, 176, 789, 378, 780, 47, 77, 85, 94, 91, 64, 100)
#' event <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
#'TRUE, TRUE, TRUE, TRUE)
#' vectorized_guanrank(time, event)
#'
#' # Additionally, this can be done on columns of a data frame directly
#' df <- data.frame("pfs" = time, "status" = event)
#' dplyr::mutate(df, gr = vectorized_guanrank(pfs, status))
vectorized_guanrank <- function(time, event) {
  data <- data.frame("time" = as.numeric(time), "event" = as.logical(event))

  sf <- survival::survfit(survival::Surv(time, event) ~ 1)
  surv_table <- data.frame("time" = sf$time, "surv" = sf$surv)
  censored <- data[order(data$time), ] %>% dplyr::filter(!event)
  uncensored <- data[order(data$time), ] %>% dplyr::filter(event)
  censored <- suppressMessages(censored %>% 
                                 dplyr::left_join(surv_table, copy = TRUE))
  uncensored <- suppressMessages(uncensored %>% 
                                   dplyr::left_join(surv_table))
  grs <- apply(data, 1, function(row) {
    get_rank(row, sf, censored, uncensored)
  })
  grs <- grs / max(grs)
  return(grs)
}
