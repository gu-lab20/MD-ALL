#' gg_boxPlot
#'
#' @param df_in
#' @param var_value
#' @param var_group
#' @param var_strata
#' @param size
#' @param highlightLevel
#' @param sizeHighlight
#' @param plot_title
#' @param x_lab
#' @param y_lab
#'
#' @return
#' @export
#'
#' @examples
gg_boxPlot=function(df_in,var_value,var_group,var_strata=NULL,size=0.6,
                    highlightLevel=NULL,sizeHighlight=1.5,
                    cols=subtypeCol,
                    plot_title=NULL,x_lab=NULL,y_lab=NULL
){
  #get plot data
  df_in1=df_in[c(var_value,var_group)]
  names(df_in1)=c("var_value","var_group")

  #get labs
  if(is.null(x_lab)){x_lab=var_group}
  if(is.null(y_lab)){y_lab=var_value}
  if(is.null(plot_title)){plot_title=var_group}

  #get color
  # col_in=get_cols_cat(df_in$feature,cols_in=cols,cols_in_default)
  cols_in=cols[names(cols) %in% df_in1$var_group]

  # base boxplot Fun
  get_p=function(df){
    ggplot(df,aes(var_group,var_value,fill=var_group))+
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size=0.5) +
      scale_fill_manual(values=cols_in) +
      labs(title = plot_title, x = x_lab, y = y_lab)+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color="grey77"),
            panel.grid.minor.y = element_line(color="grey77"),

            panel.background = element_rect(color="white",fill="white"),

            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle =30, vjust = 1, hjust=1),
            legend.position = "none")
  }

  p1=get_p(df_in1)

  if(!is.null(highlightLevel)){
    df_highlight=df_in1[df_in1$var_group %in% highlightLevel,] %>% mutate(highlight=substr(highlightLevel,1,4))
    df_=df_in1[!df_in1$var_group %in% highlightLevel,] %>% mutate(highlight="Ref")

    p1=get_p(df_)+
      # geom_boxplot(data=df_highlight,aes(var_group,var_value,fill=var_group),outlier.shape = NA) +
      geom_point(data = df_highlight,aes(var_group,var_value,fill=var_group),size=sizeHighlight,color="red") +
      facet_grid(~highlight,scales = "free", space='free')
  }
  p1
}
