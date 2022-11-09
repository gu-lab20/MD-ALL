#' read_input
#'
#' @param file, String, Input file name
#' @param delimiter, String, delimiter
#' @param header, Logical, if using the first line as header
#' @param count_field_num, num, num of count column in the file
#' @return
#' @export read_input
#'
#' @examples
#' read_input(file)
read_input=function(file,delimiter="\t",header=F,count_field_num=2){
  df_out=read.table(file,sep = delimiter,stringsAsFactors = F,header=header)
  df_out=df_out[c(1,count_field_num)]
  names(df_out)=c("feature","TestSample")
  df_out
}
