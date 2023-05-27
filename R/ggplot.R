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
    scale_color_manual(values = col_in,
                       guide = guide_legend(title.position = "top",title.hjust = 0.5) ) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = plot_title, x = x_lab, y = y_lab,col=legend_title)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black",
                                   arrow = arrow(length = unit(0.3, "lines"), type = "closed")),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.position = "bottom",
          plot.title = element_text(size=17,face="bold",hjust=0.5)
          )+
    guides(alpha = "none")

  if(!is.null(highlightLevel)){
    df_highlight=df_in1[df_in1$var_col %in% highlightLevel,]

    p1=p1+geom_point(data=df_highlight, aes(x, y),color="cyan",shape=10, size=sizeHighlight,show.legend = FALSE)+
      geom_text_repel(data=df_highlight, aes(x, y, label=var_col), color="cyan",fontface = "bold",show.legend = FALSE)

  }
  p1
}


#' gg_barplot
#'
#' @param df_plot
#' @param x
#' @param y
#' @param group.by
#' @param strata
#' @param cols
#' @param type
#' @param position
#' @param bar_width
#' @param y_lab
#' @param x_lab
#' @param legend_lab
#' @param x_tick_label_var
#' @param title
#'
#' @return
#' @export
#'
#' @examples
gg_barplot=function(df_plot,x=NULL,y=NULL,group.by=NULL,strata=NULL,cols=cols_in_default,type="identity",
                    position="stack",bar_width=0.9,
                    y_lab=NULL,x_lab=NULL,legend_lab=NULL,
                    x_tick_label_var=NULL,title=""

){

  if(is.null(y_lab)){y_lab=y};if(is.null(x_lab)){x_lab=x}

  df_plot1=df_plot[c(x,y,group.by,strata,x_tick_label_var)]
  names(df_plot1)=c("x","y",'group.by','strata',"x_tick_label_var")[!c(is.null(x),is.null(y),is.null(group.by),is.null(strata),is.null(x_tick_label_var))]

  p1=ggplot(df_plot1,aes(x=x,y=y,fill=group.by))+
    geom_bar(stat="identity",position=position,width = bar_width) +
    scale_fill_manual(values = cols) +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+xlab(x_lab) +ylab(y_lab)+labs(title=title)

  if(!is.null(strata)){
    p1=p1+facet_wrap(~strata,nrow = 1)
  }

  if(!is.null(x_tick_label_var)){
    df_x_tick_label=df_plot1 %>% select(x,x_tick_label_var) %>% distinct() %>% arrange(x)
    p1=p1+scale_x_discrete(labels=df_x_tick_label$x_tick_label_var)
  }

  if(!is.null(legend_lab)){p1=p1+labs(fill = legend_lab)}

  p1
}

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
    ggplot(df,aes(var_group,var_value))+
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size=1,aes(color=var_group)) +
      scale_color_manual(values=cols_in) +
      labs(title = plot_title, x = x_lab, y = y_lab)+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color="grey77",linetype = 4),
            panel.grid.minor.y = element_line(color="grey77",linetype = 4),
            panel.background = element_rect(color="white",fill="white"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle =30, vjust = 1, hjust=1),
            plot.title = element_text(size=17,face="bold",hjust=0.5),
            legend.position = "none")
  }

  p1=get_p(df_in1)

  if(!is.null(highlightLevel)){
    df_highlight=df_in1[df_in1$var_group %in% highlightLevel,] %>%
      mutate(highlight=substr(highlightLevel,1,4),
             var_group1=sort(unique(df_in1$var_group))[1])
    df_=df_in1[!df_in1$var_group %in% highlightLevel,] %>% mutate(highlight="Ref")

    # p1=get_p(df_)+
    #   geom_point(data = df_highlight,aes(var_group,var_value,fill=var_group),size=sizeHighlight,color="red") +
    #   facet_grid(~highlight,scales = "free", space='free')

    geo_seq <- function(a,b,n){seq(a,b,(b-a)/n)}

    y=geo_seq(min(df_in1$var_value,na.rm = T),max(df_in1$var_value,na.rm = T),3)
    label=sprintf("%.1f",y)

    df_tick=bind_rows(
      data.frame(y=y,label=label,stringsAsFactors = F),
      data.frame(y=df_highlight$var_value,label=df_highlight$var_group,stringsAsFactors = F)
    ) %>% arrange(y)

    p1=get_p(df_)+
      geom_hline(yintercept = df_highlight$var_value,color="red",linetype=4) +
      scale_y_continuous(breaks = df_tick$y,labels = df_tick$label)
  }
  p1
}

#' draw_DimPlot
#'
#' @param obj_in
#' @param group.by
#' @param cols
#' @param axis_by
#' @param reduction
#' @param highlightLevel
#' @param sizeHighlight
#' @param plot_title
#'
#' @return
#' @export
#'
#' @examples
draw_DimPlot=function(obj_in,group.by="diag",cols=subtypeCol(),axis_by=10,reduction="tsne",
                      plot_title=NULL,legend_title="Subtypes",
                      highlightLevel=NULL,sizeHighlight=2){

  # if(!reduction %in% c("tsne","umap")){stop("Reduction not defined")}

  df_plot=get_embeding_feature(obj_in = obj_in,feature = group.by,reduction = reduction)

  if(reduction=="tsne"){
    if(is.null(plot_title)){plot_title="tSNE"}
    x_lab=paste0("tSNE_1\nGene N = ", unique(df_plot$FeatureN), "; Perplexity N=", unique(df_plot$perplexityN))

    p1=gg_dimPlot(df_in=df_plot,x="tSNE_1",y="tSNE_2",var_col=group.by,cols=cols,
                  size=1,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab=x_lab,y_lab=NULL,legend_title=legend_title)  }

  if(grepl("umap",reduction)){
    if(is.null(plot_title)){plot_title="uMAP"}

    plot_title=paste0("UMAP: Gene N = ", unique(df_plot$FeatureN), "; Neighbor N=", unique(df_plot$n_neighbors))

    p1=gg_dimPlot(df_in=df_plot,x="uMAP_1",y="uMAP_2",var_col=group.by,cols=cols,
                  size=1,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab="UMAP1",y_lab="UMAP2",legend_title=legend_title)
  }
  p1
}


#' draw_BoxPlot
#'
#' @param obj_in
#' @param group.by
#' @param features
#' @param cols
#' @param assay_name_in
#' @param plot_title
#' @param x_lab
#' @param y_lab
#' @param highlightLevel
#' @param sizeHighlight
#'
#' @return
#' @export
#'
#' @examples
draw_BoxPlot=function(obj_in,group.by="diag",features=c("CRLF2"),
                      cols=subtypeCol(),
                      assay_name_in="vst",
                      plot_title="Box Plot",x_lab="Subtype",y_lab=NULL,
                      useGeneName=T,df_geneName=info_gtf_hg38,
                      highlightLevel=NULL,sizeHighlight=1.5){

  if(is.null(y_lab)){y_lab=features}

  if(!useGeneName){
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(features,group.by))
    gg_boxPlot(df_feature,var_value = features,var_group = group.by,plot_title = plot_title,x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }

  if(useGeneName){
    geneId=df_geneName$gene_id[df_geneName$gene_name %in% features][1]
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(geneId,group.by))
    gg_boxPlot(df_feature,var_value = geneId,var_group = group.by,plot_title = plot_title,x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }
}


#' gg_tilePlot
#'
#' @param df_in
#' @param x
#' @param y
#' @param var_col
#' @param cols
#' @param y_lab
#' @param x_lab
#' @param legend_lab
#' @param reverse_y
#' @param x_tick_label_var
#' @param squared
#' @param add_horizontal
#' @param add_border
#' @param border_lwd
#' @param border_linetype
#' @param border_color
#' @param x.axis
#'
#' @return
#' @export
#'
#' @examples
gg_tilePlot=function(df_in,x,y,var_col,cols=subtypeCol(),
                     title="Heatmap",
                     y_lab="",x_lab="Number of features",legend_lab="Subtype:",reverse_y=F,
                     x_tick_label_var=NULL,
                     squared=T,
                     add_border=T,border_lwd=0.3,border_linetype=1,border_color="white",
                     x.axis=T){
  #get data
  df_in1=df_in[c(x,y,var_col,x_tick_label_var)]
  names(df_in1)=c('x','y','var_col',"x_tick_label_var")[c(T,T,T,!is.null(x_tick_label_var))]

  #get lab
  if(is.null(y_lab)){y_lab=y}
  if(is.null(x_lab)){x_lab=x}

  p1=ggplot(df_in1, aes(x = x, y = y, fill = var_col)) +
    geom_tile() +
    scale_fill_manual(values = cols) +
    labs(x=x_lab,y=y_lab,title=title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 17, face = "bold",hjust = 0.5),
          legend.position = "bottom",

    )

  if(reverse_y){p1=p1 + scale_y_discrete(limits = rev(sort(unique(df_in1$y))))}

  if(!is.null(legend_lab)){p1=p1+labs(fill = legend_lab)}

  if(!is.null(x_tick_label_var)){
    df_x_tick_label=df_in1 %>% select(x,x_tick_label_var) %>% distinct() %>% arrange(x)
    p1=p1+scale_x_discrete(labels=df_x_tick_label$x_tick_label_var)
  }

  if(!x.axis){p1=p1+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())}

  if(add_border){p1=p1 +
    geom_tile(color = border_color,lwd = border_lwd,linetype = border_linetype)+
    theme(legend.key = element_rect(fill = NA, colour = NA))
  }

  if(squared){p1=p1 + coord_fixed()}


  p1
}

#' gg_heatmap
#'
#' @param df_in
#' @param x
#' @param y
#' @param var_col
#' @param y_lab
#' @param x_lab
#' @param y_tick_label_var
#' @param reverse_y
#' @param y_xais_side
#' @param title
#'
#' @return
#' @export
#'
#' @examples
gg_heatmap=function(df_in,x = "labels",y = "variable",var_col = "value",
                    y_lab=NULL,x_lab=NULL,
                    y_tick_label_var=NULL,
                    reverse_y=F,
                    y_xais_side="right",
                    title="Scores (column scaled to 0-1)"
){

  #get data
  df_in1=df_in[c(x,y,var_col,y_tick_label_var)]
  names(df_in1)=c('x','y','var_col',"y_tick_label_var")[c(T,T,T,!is.null(y_tick_label_var))]

  #get lab
  if(is.null(y_lab)){y_lab=""}
  if(is.null(x_lab)){x_lab=""}
  if(is.null(title)){title=NULL}

  p1=ggplot(df_in1, aes(x=x, y=y, fill= var_col)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("midnightblue","navyblue","deepskyblue3","lightseagreen","yellow1")) +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")+
    guides(fill = guide_colourbar(ticks = FALSE,barwidth = 15,barheight = 0.8,title.position = "top",title.hjust = 0.5,title = title))

  if(is.null(y_tick_label_var) & reverse_y){
    p1=p1 + scale_y_discrete(position=y_xais_side,limits = rev(sort(unique(df_in1$y))))}

  if((!is.null(y_tick_label_var)) & reverse_y){
    df_y_tick_label=df_in1 %>% select(y,y_tick_label_var) %>% distinct() %>% arrange(y)
    p1=p1+scale_y_discrete(position=y_xais_side,labels=rev(df_y_tick_label$y_tick_label_var),limits = rev(sort(unique(df_in1$y))))  }

  if((!is.null(y_tick_label_var)) & !reverse_y){
    df_y_tick_label=df_in1 %>% select(y,y_tick_label_var) %>% distinct() %>% arrange(y)
    p1=p1+scale_y_discrete(position=y_xais_side,labels=df_y_tick_label$y_tick_label_var)
  }

  p1
}
