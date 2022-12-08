#' knn_pred_one
#'
#' @param indata
#' @param var_y
#' @param var_x
#' @param KNN_K
#'
#' @return
#' @export
#'
#' @examples
knn_pred_one=function(indata,var_y,var_x,KNN_K=5){
  indata=indata %>% mutate(obs=1:n())

  formula_in=formula(paste0(var_y,"~",paste0(var_x,collapse = "+")))

  df_freq=as.data.frame(table(unlist(indata[var_y]))) %>% arrange(Freq)

  if(min(df_freq$Freq)==1){
    #split data
    levels_1=as.character(df_freq$Var1[df_freq$Freq==1])
    indata_1=indata[unlist(indata[var_y]) %in% levels_1,]
    indata_rest=indata %>% filter(!obs %in% indata_1$obs)

    #get KNN model
    knnFit_rest = caret::train(formula_in, data = indata_rest, method = "knn",
                               preProcess = c("center", "scale"),
                               tuneGrid = expand.grid(k = c(KNN_K)))

    #get KNN predicted value
    indata_rest$knn_pred=predict(knnFit_rest) %>% as.vector()
    indata_1$knn_pred=predict(knnFit_rest,newdata = indata_1) %>% as.vector()
    df_knn=bind_rows(indata_rest,indata_1) %>% arrange(obs)

    if(!all(df_knn$obs==indata$obs)){stop("obs not match in original data and predicted dataset")}

    knn_pred=df_knn$knn_pred
  }

  if(min(df_freq$Freq)>1){
    knnFit = train(formula_in, data = indata, method = "knn",
                   preProcess = c("center", "scale"),
                   tuneGrid = expand.grid(k = c(KNN_K)))
    knn_pred=predict(knnFit) %>% as.vector()
  }

  knn_pred
}
