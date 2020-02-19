p <- function(ta, tb, sf) {
  later_prob <- min(sf$surv[sf$time == tb], sf$surv[sf$time == ta])
  earlier_prob <- max(sf$surv[sf$time == tb], sf$surv[sf$time == ta])

  return(1 - (later_prob / earlier_prob))
}

get_rank <- function(row, sf, censored, uncensored) {
  surv_table <- data.frame('time'=sf$time, 'surv'=sf$surv)
  row_event <- as.logical(row['event'])
  row_n <- as.numeric(row['n'])
  row_time <- as.numeric(row['time'])
  row_surv <- sf$surv[sf$time == row_time]
  censored <- censored %>% left_join(surv_table, copy=TRUE)
  uncensored <- uncensored %>% left_join(surv_table)
  if(!row_event) {
    # Calculate score for a censored subject
    # Assuming that for all T > t means for all events with time greater than
    # t. Not just for all unique times (could be fewer than events)
    # sum(1 - get_unc_df$surv/row_surv) is equivalent to the 1st + 2nd sum in
    # the equation for the censored hazard rank
    geq_unc_df <- uncensored %>% dplyr::filter(time >= row_time)
    geq_unc <- ifelse(nrow(geq_unc_df)==0, 0,
                      sum(1 - geq_unc_df$surv/row_surv))

    # sum(1 - 0.5*geq_cen_df$surv/row_surv) is equivalent to the 3rd sum in the
    # equation for the censored hazard rank
    geq_cen_df <- censored %>% dplyr::filter(time >= row_time)
    geq_cen <- ifelse(nrow(geq_cen_df)==0, 0,
                      sum(1 - 0.5*geq_cen_df$surv/row_surv))

    # Equivalent to the 4th sum
    lt_cen_df <- censored %>% dplyr::filter(time < row_time)
    lt_cen <- ifelse(nrow(lt_cen_df)==0, 0, sum(0.5*row_surv/lt_cen_df$surv))
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
    eq_unc <- 0.5*nrow(uncensored %>% filter(time == row_time))

    # Equivalent to the 2nd and 4th sums in the hazard rank equation for
    # uncensored subjects
    leq_cen_df <- censored %>% dplyr::filter(time <= row_time)
    leq_cen <-  ifelse(nrow(leq_cen_df)==0, 0, sum(row_surv/leq_cen_df$surv))
    return(sum(eq_unc, gt_tot, leq_cen, -0.5))
  }
}

vectorized_guanrank <- function(id, time, event, data = NULL) {
  if(is.null(data)){
    data <- data.frame('id' = id, 'time' = time, 'event' = event)
  }
  else {
    id <- data$id
    time <- data$time
    event <- data$event
  }
  sf <- survival::survfit(survival::Surv(time, event)~1)
  censored <- data[order(data$time), ] %>% dplyr::filter(!as.logical(event))
  uncensored <- data[order(data$time), ] %>% dplyr::filter(as.logical(event))
  grs <- apply(data, 1, function(row) {
    suppressMessages(get_rank(row, sf, censored, uncensored))
    })
  grs <- grs/max(grs)
  return(grs)
}
