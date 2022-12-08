#' prepare_snv
#'
#' @param snv_file
#' @param centr_ref
#' @param minDepth
#' @param snv_format
#'
#' @return
#' @export
#'
#' @examples
prepare_snv <- function(snv_file, centr_ref, minDepth, snv_format = "vcf") {

  sample_name="TestSample"
  smpSNP=list()
  print(paste("Preparing file with snv information for:", sample_name))

  #prepare data from custom table
  if (snv_format == "custom") {
    #read the table
    snv_table_pre <- fread(snv_file)

    #check if all appropriate columns are present
    cols <- colnames(snv_table_pre)
    chr <- str_which(cols, "^#Chromosome$|#CHROM$|^CHR$|^chr$|^Chr$|^CHROM$|^chrom$|^Chrom$|^CHROMOSOME$|^chromosome$|^Chromosome$")
    start <- str_which(cols, "^START$|^start$|^Start$|^POS$|^pos$|^Pos$")
    depth <- str_which(cols, "^DEPTH$|^depth$|^Depth$|^DP$|^dp$|^Dp$")
    maf <- str_which(cols, "^MAF$|^maf$|^Maf$|^HET$|^het$|^Het$")
    to_keep <- as.numeric(c(chr, start, depth, maf))
    if (length(to_keep) != 4) {
      smpSNP[[sample_name]] <- "Incorrect column name in a custom snv file."
      return(smpSNP)
    }
    snv_table <- snv_table_pre[, to_keep, with = FALSE]
    data.table::setnames(snv_table, colnames(snv_table), c("chr", "start", "depth", "maf"))

    #Check some column parameters
    if (is.numeric(snv_table[, start]) == FALSE | is.numeric(snv_table[, depth]) == FALSE | is.numeric(snv_table[, maf]) == FALSE) {
      smpSNP[[sample_name]] <- "Incorrect type of a column in a custom file with snv information."
      return(smpSNP)
    }

  } else if (snv_format == "vcf") {
    snv_table <- vcf_to_snv(snv_file)
    if(is.character(snv_table)) {
      smpSNP[[sample_name]] <- snv_table
      return(smpSNP)
    }
  }

  smpSNP[[sample_name]] <- snv_table %>%
    filter(chr %in% c(c(1:22, "X"), paste0("chr", c(1:22, "X")))) %>%
    mutate(chr = sub("chr", "", chr)) %>% left_join(centr_ref, by = "chr") %>%
    mutate(chr = factor(chr, levels=c(1:22, "X")), ID=paste0(chr,"-", start), sampleID=sample_name) %>%
    mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr")))

  return(smpSNP)
}

# snv_file="tests/test.HaplotypeCaller.vcf"
