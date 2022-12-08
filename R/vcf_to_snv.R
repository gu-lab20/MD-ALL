#' vcf_to_snv
#'
#' @param vcf_file
#' @param maf_tresh
#' @param depth_tresh
#'
#' @return
#' @export
#'
#' @examples
vcf_to_snv=function(vcf_file, maf_tresh = 0.01, depth_tresh = 5, outfmt=1){
  #read the vcf files
  message("Reading in vcf file..")
  vcf_data=as.data.frame(data.table::fread(vcf_file, skip = "#CHROM", header = TRUE))

  #check whether the input is in correct format
  if (dim(vcf_data)[2] < 10) {
    vcf_final <- "Incorrect vcf file format. Incorrect number of columns"
    return(vcf_final)
  }
  if (!identical(colnames(vcf_data)[1:9], c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))) {
    vcf_final <- "Incorrect vcf file format."
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, 'INFO'], "DP=[0-9]+") != TRUE &
      str_detect(string = str_split(vcf_data[1, 'FORMAT'], pattern = ":")[[1]][3], pattern = "DP") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No depth information in the INFO column and missing information in the FORMAT column"
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, 9], "AD") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No allele depth (AD) in FORMAT column"
    return(vcf_final)
  }

  vcf_data <- vcf_data[, 1:10]
  names(vcf_data)=c("chr", "start", "ID", "ref", "var", "qual", "FILTER", "INFO", "FORMAT", "INFO2")

  # Getting depth out of the INFO column
  message("Extracting depth..")

  if (str_detect(vcf_data[1, 'INFO'], "DP=[0-9]+") == TRUE) {
    vcf_data$depth = as.numeric(gsub("DP=", "", str_extract(vcf_data$INFO, "DP=[0-9]+")))
  } else {
    DP_i=grep("DP",unlist(strsplit(vcf_data$FORMAT[1],split = ":")))
    vcf_data$depth=sapply(str_split(vcf_data$INFO2,":"),function(x){x[DP_i]})
  }

  #reference allele and alternative allele depths
  message("Extracting reference allele and alternative allele depths..")
  vcf_data$REF_ALT_num = sapply(str_split(vcf_data$INFO2, ":"), function(x) x[2])

  #extract count for alternative allele
  vcf_data$varDp = as.numeric(sapply(str_split(vcf_data$REF_ALT_num, ","), function(x) x[2]))

  #mutant allele frequency
  vcf_data$maf = vcf_data$varDp/vcf_data$depth
  message("Needed information from vcf extracted")

  #return needed columns
  if(outfmt==1){
  vcf_final <- vcf_data[,c("chr", "start", 'depth', "maf")]
  } else {
    vcf_final=vcf_data[,c("chr", "start","ref" ,'var','depth', "maf")]
  }

  message("Finished reading vcf")
  return(vcf_final)

}
# vcf_final=NULL
# vcf_file="tests/test.HaplotypeCaller.vcf"



