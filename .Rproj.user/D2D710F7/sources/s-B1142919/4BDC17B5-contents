#' read_input
#'
#' @param file, String, Input file name
#' @param delimiter, String, delimiter
#' @param header, Logical, if using the first line as header
#'
#' @return
#' @export read_input
#'
#' @examples
#' read_input(file)
read_input=function(file,delimiter="\t",header=F){
  df_out=read.table(file,header = T,sep = delimiter,stringsAsFactors = F,header=header)
  df_out=df_out[1:2]
  names(df_df_out)=c("feature","TestSample")
}

