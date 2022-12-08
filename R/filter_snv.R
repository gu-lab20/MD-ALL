####filter SNVs of interest for samples###
#' Title
#'
#' @param one_smpSNP
#' @param keepSNP
#' @param minDepth
#' @param mafRange
#'
#' @return
#' @export
#'
#' @examples
filter_snv <- function(one_smpSNP, keepSNP, minDepth, mafRange) {

  smpSNPdata.tmp= one_smpSNP %>% dplyr::select(sampleID, ID, maf, chr, start, depth, arm) %>%
    filter(data.table::inrange(maf, mafRange[1], mafRange[2]), depth > minDepth) %>% filter(chr != "Y")
  if (keepSNP[1] != FALSE) {
    smpSNPdata.tmp <- smpSNPdata.tmp %>% filter(ID %in% keepSNP)
  }
  return(smpSNPdata.tmp)
}
