####get median expression level####
#' Title get_med
#'
#' @param count_norm
#' @param refDataExp
#' @param weight_table
#' @param generate_weights
#'
#' @return
#' @export
#'
#' @examples
get_med <- function(count_norm, refDataExp, weight_table, generate_weights) {

  ENSG=rownames(count_norm)

  ####calculate median for all genes####
  pickGeneDFall_tmp=count_norm %>% mutate(ENSG=ENSG) %>% left_join(select(refDataExp, chr, ENSG), by = "ENSG")
  if (generate_weights == TRUE) {
    pickGeneDFall <- pickGeneDFall_tmp %>% mutate(med = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median), var = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, var)) %>%
      select(ENSG, chr, med, var)
  } else {
    pickGeneDFall <- pickGeneDFall_tmp %>% mutate(med = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median)) %>%
      select(ENSG, med)
  }

  return(pickGeneDFall)
}
