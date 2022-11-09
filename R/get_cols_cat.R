#' get_cols_cat
#'
#' @param value
#' @param cols_in
#' @param cols_in_default
#'
#' @return
#' @export
#'
#' @examples
get_cols_cat=function(value,cols_in,cols_in_default){
  if(!any(names(cols_in) %in% value)){
    print("Supplied color labels not inclued in values levels, use default cols instead")
    if(length(unique(value))>  length(cols_in_default)){message("Defaule color number not enough, the rest will using grey")}
    cols_out=cols_in_default[1:length(unique(value))]
    names(cols_out)=sort(unique(value))
  }

  if(any(names(cols_in) %in% value)){
    cols_out=cols_in[names(cols_in) %in% value]
    if(length(cols_out)<length(levels(as.factor(value)))){
      cols_out=cols_in[1:length(levels(as.factor(value)))]
    }
  }
  cols_out
}
