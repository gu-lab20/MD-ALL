#' gg_dimPlot
#'
#' @param df_in
#' @param x
#' @param y
#' @param var_col
#' @param cols
#' @param size
#' @param axis_by
#' @param highlightLevel
#' @param sizeHighlight
#' @param plot_title
#' @param x_lab
#' @param y_lab
#' @param legend_title
#'
#' @return
#' @export
#'
#' @examples
gg_dimPlot=function(df_in,x="uMAP_1",y="uMAP_2",var_col="diag",cols=subtypeCol(),
                    size=1,axis_by=5,
                    highlightLevel=NULL,sizeHighlight=2,
                    plot_title=NULL,x_lab=NULL,y_lab=NULL,legend_title="Group"){
  #get plot data
  df_in1=df_in[c(x,y,var_col)]
  names(df_in1)=c("x","y","var_col")

  #get labs
  if(is.null(x_lab)){x_lab=x}
  if(is.null(y_lab)){y_lab=y}
  if(is.null(plot_title)){plot_title=var_col}

  #get xaixs breaks
  min_x=min(df_in1$x);max_x=max(df_in1$x);x_breaks=round(seq(min_x,max_x,by=axis_by),0)
  min_y=min(df_in1$y);max_y=max(df_in1$y);y_breaks=round(seq(min_y,max_y,by=axis_by),0)

  # col_in=get_cols_cat(df_in$feature,cols_in=cols,cols_in_default)
  col_in=cols[names(cols) %in% df_in1$var_col]

  p1=ggplot(df_in1,aes(x=x,y=y,col=var_col)) +
    geom_point(data = df_in1,aes(x=x,y=y,color=var_col),size=size)+
    scale_color_manual(values = col_in) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = plot_title, x = x_lab, y = y_lab,col=legend_title)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.line = element_line(colour = "black"))+
    guides(alpha = "none")

  if(!is.null(highlightLevel)){
    df_highlight=df_in1[df_in1$var_col %in% highlightLevel,]

    p1=p1+geom_point(data=df_highlight, aes(x, y, color=var_col), size=sizeHighlight,show.legend = FALSE)+
      geom_text_repel(data=df_highlight, aes(x, y, label=var_col), fontface = "bold",show.legend = FALSE)

  }
  p1
}
