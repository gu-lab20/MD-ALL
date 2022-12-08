#' get_norm_exp
#'
#' @param count_table
#' @param standard_samples
#' @param minReadCnt
#' @param samp_prop
#' @param weight_table
#' @param weight_samp_prop
#' @param generate_weights
#'
#' @return
#' @export
#'
#' @examples
get_norm_exp <- function(count_table, standard_samples, minReadCnt, samp_prop, weight_table, weight_samp_prop, generate_weights) {

  sample_name="TestSample"
  #rename count table
  names(count_table)=c("ENSG", "count")


  #check the count file format
  if (ncol(count_table) != 2 | !is.numeric(count_table[, 2])) {
    return(paste0("Incorrect count file format for the sample: ", sample_name))
  }

  # set column to join by
  # data.table::setkey(count_table, ENSG)

  #inner join diploid reference and analyzed sample
  # data.table::setkey(standard_samples, ENSG)
  # final_mat <- as.data.frame(count_table[standard_samples, nomatch = 0])
  final_mat=count_table %>% left_join(standard_samples)

  #keep genes for determining gender for later
  gender_genes = final_mat %>% filter(ENSG %in% "ENSG00000012817")

  #filter genes based on reads count; top 1-q have read count > N, filter base on weight
  keepIdx_tmp = final_mat %>%
    mutate(keep_gene = apply(.[, -1], MARGIN = 1, FUN = function(x) sum(x > minReadCnt) > (length(x) * samp_prop)), id = row_number()) %>%
    filter(keep_gene == TRUE)

  if (generate_weights == FALSE) {
    keepIdx = keepIdx_tmp %>% inner_join(weight_table, by = "ENSG") %>%
      group_by(chromosome_name) %>% mutate(weight_chr_quant = quantile(weight, 1 - weight_samp_prop)) %>%
      filter(weight >= weight_chr_quant) %>% pull(id)
  } else {
    keepIdx = keepIdx_tmp %>% pull(id)
  }

  # filter table for normalization, get rid of genes with 0 counts and keep genes for gender estimation
  count_filt <- final_mat %>% .[c(keepIdx), ] %>% .[pull(., 2) != 0, ] %>% bind_rows(gender_genes) %>% distinct(ENSG, .keep_all = TRUE)
  ENSG <- count_filt$ENSG
  count_filt <- select(count_filt, -ENSG)

  #sample Deseq normalization
  count_col <- as.data.frame(colnames(count_filt))

  dds <- DESeq2::DESeqDataSetFromMatrix(colData = count_col, countData = count_filt, design= ~ 1)
  dds_vst <- DESeq2::varianceStabilizingTransformation(dds, blind=T, fitType='local')
  count_norm <- as.data.frame(SummarizedExperiment::assay(dds_vst))


  colnames(count_norm)[1] <- sample_name
  colnames(count_norm)[2:ncol(count_norm)] <- paste0("control_", 1:(ncol(count_norm)-1))
  print(paste0("Normalization for sample: ", sample_name, " completed"))

  #Modify table for downstream analysis
  rownames(count_norm) <- ENSG

  return(count_norm)
}
