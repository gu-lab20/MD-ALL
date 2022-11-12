#' load_data
#'
#' @param label
#'
#' @return
#' @export
#'
#' @examples
load_data=function(label="obj_234_HTSeqCount"){
  if(label=="obj_234_HTSeqCount"){
    out=  readRDS("data-raw/obj_234_HTSeqCountMini.rds")
  }

  if(label=="obj_2042_HTSeqCount"){
    out=  readRDS("data-raw/obj_2042_HTSeqCountMini.rds")
  }

  out
}
