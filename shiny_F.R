library(stringr)
library(dplyr)
library(DESeq2)
library(shiny)
library(Rtsne)
library(caret)
library(ggrepel)
library(shinydashboard)
library(Rphenograph)
library(reshape2)
library(tidyr)
library(tibble)
library(data.table)
library(future)

setwd("//isi-dcnl/user_data/zgu_grp/DryLab/bin/jiliu/R")
gc(rm(list=ls()))



#run vst --------------------------
# load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/for_GEPprediction_tSNE2KNN_1499.rdata")

load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/for_GEPprediction_Phenograph_1409_1batch.rdata")
for_tsne_1409_=as.data.frame(for_tsne_1409)
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/for_vst_N234.rdata")
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/gene_name.RData")
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/Diag_vst_MedianMeanValue.RData")
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/df_forboxplot.RData")

gene_list=gene_ref_$gene_name


diagCol=c()
{
  diagCol["BCL2/MYC"]="seagreen2"
  diagCol["DUX4"]='grey40'
  diagCol["ETV6-RUNX1"]="gold2"
  diagCol["HLF"]= "skyblue"
  diagCol["Hyperdiploid"]="#3E9F32"
  diagCol["iAMP21"]="lightslateblue"
  diagCol["IKZF1 N159Y"]="#CCCC33"
  diagCol["KMT2A"]="#1F78B5"
  diagCol["Low hypodiploid"]="#1E90FF"
  diagCol["MEF2D"]="#66C2A6"
  diagCol["NUTM1"]='black'
  diagCol["PAX5-ETV6"]="#808000"
  diagCol["PAX5 P80R"]="orangered"
  diagCol["PAX5alt"]="#FFA620"
  diagCol["Ph"]="magenta3"
  diagCol["TCF3-PBX1"]="darkgoldenrod4"
  diagCol["Y"]="#E6BEFF"
  diagCol["ZEB2/CEBPE"]="#D27B1C86"
  diagCol["ZNF384"]="#A8DD00"
  diagCol["_Prediction"]="red4"
  diagCol["Sample_test"]="red4"
  
}

# color -------------------
#BALL color
subtypeCol=c()
{
  subtypeCol["ETV6-RUNX1"]="gold2"
  subtypeCol["ETV6-RUNX1-like"]="pink"
  subtypeCol["KMT2A"]="#1F78B5"
  subtypeCol["Ph"]="magenta3"
  subtypeCol["DUX4"]='grey40'
  subtypeCol["TCF3-PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["BCL2/MYC"]="seagreen2"
  subtypeCol["NUTM1"]='black'
  subtypeCol["HLF"]= "skyblue"
  subtypeCol["PAX5(P80R)"]="orangered"
  subtypeCol["PAX5 P80R"]="orangered"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["LowHypo"]="#1E90FF"
  subtypeCol["Low hypodiploid"]="#1E90FF"
  subtypeCol["NearHaploid"]='blue3'
  subtypeCol["Near haploid"]='blue3'
  subtypeCol["Ph-like"]="red4"
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["PAX5-ETV6"]="#808000"
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["IKZF1(N159Y)"]="#CCCC33"
  subtypeCol["IKZF1 N159Y"]="#CCCC33"
  subtypeCol["LowHyper"]="cyan"
  subtypeCol["Bother"]='grey75'
  subtypeCol["Low hyperdiploid"]='grey75'
  subtypeCol["CRLF2(non-Ph-like)"]='grey75'
  subtypeCol["KMT2A-like"]='grey75'
  subtypeCol["ZNF384-like"]='grey75'
  subtypeCol["Other"]='grey75'
  subtypeCol["ZEB2/CEBPE"]="#D27B1C86"
  subtypeCol["Y"]="#E6BEFF"
  subtypeCol["Unknown"]="#469990"
  subtypeCol["_Prediction"]="red4"
  subtypeCol["Normal"]="pink"
}


#write_tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

#read tsv --------------
get_COH_sample=function(infile){
  COH_sample=substr(basename(infile),1,12)
}

read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}


draw_tsne_MAD=function(indata,data_diag,diag_col,N,i,axis_by,size=1){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1")) +    ylab("tSNE dimension 2")+
    theme_bw() +
    
    geom_point(data=subset(data_plot, diag != "_Prediction" ), aes(X, Y, color=diag), size=size)+
    geom_point(data=subset(data_plot, diag == "_Prediction" ), aes(X, Y, color=diag), size=3, alpha=0.75)+
    geom_text_repel(data=subset(data_plot, diag == "_Prediction" ), aes(X, Y, label=COH_sample), fontface = "bold")+
    scale_color_manual(values = diag_col) +
    
    scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = axis_by),1)) +
    
    theme(axis.title=element_text(size=10),
          #axis.text=element_blank(),
          #axis.ticks=element_blank(),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10,face="bold")
    ) +
    scale_fill_discrete(name="Experimental\nCondition")+
    guides(color = guide_legend(override.aes = list(size = 1)))
}

run_tsne=function(indata,perplexity_one,out_dims){
  indata1=t(indata)
  cat("Run tSNE: Used Feature N=", dim(indata1)[2], "; Used Sample N=", dim(indata1)[1],"; Perplexity=", perplexity_one,  "\n")
  
  set.seed(10) # Sets seed for reproducibility
  tsne_out = Rtsne(indata1, dims = out_dims, perplexity = perplexity_one, 
                   theta = 0.5, 
                   max_iter = 5000, 
                   check_duplicates = F,
                   partial_pca=T, 
                   num_threads = 4)
  
  #summary output
  df_id=data.frame(COH_sample=colnames(indata)) 
  
  df_tsne=as.data.frame(tsne_out$Y)
  names(df_tsne)=c("X","Y","Z")[1:ncol(df_tsne)]
  df_tsne$COH_sample=colnames(indata)
  
  df_tsne=df_id %>% left_join(df_tsne) %>%
    mutate(PerplexityN=perplexity_one,
           FeatureN=dim(indata1)[2]) 
  df_tsne
}

create_deseqObject_from_matrixdata=function(matrix_in,diag_in,sample_var_name,design_var){
  
  sample_name_in=names(matrix_in)[names(matrix_in) %in% unlist(diag_in[sample_var_name])]
  sample_name_in=sort(sample_name_in)
  feature_in=rownames(matrix_in)
  
  cat(paste0("Used feature N=",length(feature_in)));  cat("\n")
  cat(paste0("Used sample N=",length(sample_name_in)));  cat("\n")
  
  matrix_in1=matrix_in[,sample_name_in]
  
  row.names(matrix_in1)=feature_in
  
  diag_in1=diag_in[unlist(diag_in[sample_var_name]) %in% sample_name_in,] 
  diag_in1=diag_in1[order(unlist(diag_in1[sample_var_name])),]
  
  formula_in=formula(paste0("~",paste0(design_var,collapse = "+")))
  
  htseq_object = DESeq2::DESeqDataSetFromMatrix(countData = matrix_in1,
                                        colData = diag_in1,
                                        design = formula_in)
  cat(dim(htseq_object));cat("\n")
  htseq_object
}

run_vst=function(object_in){
  cat("Running vst\n")
  object_out=DESeq2::varianceStabilizingTransformation(object_in,blind = T)
  data_round=data.frame(assay(object_out))
  dim(data_round)
  data_round
}

runVST_FromFile=function(file){
  id="sample_test"
  # id=word(basename(file),1,sep="[.]")
  #prepare count for reference + new sample
  count_one=read.table(file,header = F,stringsAsFactors = F) 
  names(count_one)=c("feature",id)
  count_panel_with_one=count_gene_234 %>% left_join(count_one)
  row.names(count_panel_with_one)=count_panel_with_one$feature
  
  #prepare diag with panel data
  if(!id %in% diag_234$COH_sample){
    diag_in=rbind(
      diag_234 %>% mutate(source="N2239") %>% select(COH_sample,diag,source),
      data.frame(
        COH_sample=id,
        stringsAsFactors = F,
        diag="Test",
        source="for_prediction"
      )
    )
  }
  
  if(id %in% diag_234$COH_sample){
    diag_in=diag_234
  }
  
  diag_in$diag=as.factor(diag_in$diag)
  row.names(diag_in)=diag_in$COH_sample
  
  #running vst
  deseqobject=create_deseqObject_from_matrixdata(count_panel_with_one,diag_in,"COH_sample","diag")
  
  df_vst=run_vst(deseqobject)
  
  #output data
  df_vst_one=data.frame(
    feature=row.names(df_vst),
    vst=unlist(df_vst[id]),
    stringsAsFactors = F
  )
  # rownames(df_vst_one)<- gsub(id,"",rownames(df_vst_one))
  names(df_vst_one)[2]=id
  df_vst_one
}
# draw tsne--------------
draw_tsne_FromVst=function(df_vst,diag_ref,df_vst_ref,feature_in,perplexity_one){
  COH_sample=names(df_vst)[2]

  cat(paste0("Run analysis for ",COH_sample,"\n\n"))

  diag_ref_1=diag_ref[!diag_ref$COH_sample %in% COH_sample,]

  diag_in_=rbind(
    diag_ref_1 %>% mutate(source="ref") %>% select(COH_sample,diag,source),
    data.frame(
      COH_sample=COH_sample,
      diag="_Prediction",
      source="for_prediction",
      stringsAsFactors = F
    )
  )

  df_vst_ref$feature=row.names(df_vst_ref)
  df_vstAll=df_vst_ref %>% left_join(df_vst)
  row.names(df_vstAll)=df_vstAll$feature
  df_vstAll$feature=NULL


  for_tsne=df_vstAll[head(feature_in,800),]

  for_tsne_plot=run_tsne(for_tsne,perplexity_one, 2)
  for_tsne_plot1=for_tsne_plot %>% left_join(diag_in_)

  plot_out=draw_tsne_MAD(for_tsne_plot,diag_in_,diagCol,800,perplexity_one,10)
  return(plot_out)
}

tsne_df_FromVst=function(df_vst,diag_ref,df_vst_ref,feature_in,perplexity_one){
  COH_sample=names(df_vst)[2]
  
  cat(paste0("Run analysis for ",COH_sample,"\n\n"))
  
  diag_ref_1=diag_ref[!diag_ref$COH_sample %in% COH_sample,]
  
  diag_in_=rbind(
    diag_ref_1 %>% mutate(source="ref") %>% select(COH_sample,diag,source),
    data.frame(
      COH_sample=COH_sample,
      diag="_Prediction",
      source="for_prediction",
      stringsAsFactors = F
    )
  )
  
  df_vst_ref$feature=row.names(df_vst_ref)
  df_vstAll=df_vst_ref %>% left_join(df_vst)
  row.names(df_vstAll)=df_vstAll$feature
  df_vstAll$feature=NULL
  
  
  for_tsne=df_vstAll[head(feature_in,800),]
  
  for_tsne_plot=run_tsne(for_tsne,perplexity_one, 2)
  for_tsne_plot1=for_tsne_plot %>% left_join(diag_in_)
  for_tsne_plot1
}

tsne_df_fromTable=function(indata,data_diag,diag_col){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
}

#Boxplot--------------------

draw_boxplot=function(df_vst_ref, df_vst_sample, variable_type_in, gene_ref_,diag_ref){
  df_vst_ref=as.data.frame(df_vst_ref)
  df_vst_ref$feature=row.names(df_vst_ref)
  df_vst_comb=df_vst_ref %>% left_join(df_vst_sample)
  geneid=gene_ref_$gene_id[gene_ref_$gene_name==variable_type_in]
  if(!geneid %in% df_vst_comb$feature){stop("Typed gene not in top 5000 variable genes")}
  if(geneid %in% df_vst_comb$feature){
    df_vst1=as.data.frame(t(df_vst_comb[df_vst_comb$feature==geneid,]))
    names(df_vst1)[1]="value"
    df_vst1$COH_sample=row.names(df_vst1)
    df_vst2=df_vst1 %>% left_join(diag_ref) %>% filter(!grepl("EN",value))
    df_vst2$diag[df_vst2$COH_sample=="sample_test"]="_Prediction"
    df_vst2$value=as.numeric(df_vst2$value)
    table(df_vst2$diag)
    col_in=unique(subtypeCol[names(subtypeCol) %in% df_vst2$diag])
    ggplot(df_vst2,aes(x=diag,y=value))+
      geom_boxplot(aes(color=diag),coef=6,size=0.75) +
      scale_color_manual(values=col_in) +
      geom_jitter(size=0.5)+
      ylab(label = "Count") +
      theme_bw()
  
 }
}
draw_boxplot(for_tsne_1409_, x, "IRX1", gene_ref_,diag_1409)

#Neighbor prediction -----------------------------------------
get_predDiag=function(value_one){
  freq=table(value_one)
  freq=sort(freq,decreasing = T)
  list(pred_GEhalf=ifelse(any(freq > sum(freq)/2),names(freq)[freq > sum(freq)/2],"Fail_prediction"),
       pred_freq=paste0(paste0(names(freq),": ",freq," (",sprintf("%.1f",100*freq/sum(freq)),"%)"),collapse = "; ")
  )
}

get_top_neighbor=function(COH_sample,df_in,diag_in,top_gene,FeatureN_list){
  df_neighbor_featureN=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
    # print(paste0("FeatureN=",FeatureN_one))
    print("get feature df")
    input_gene=top_gene[1:FeatureN_one]
    df_in_one=df_in[row.names(df_in) %in% input_gene,]
    df_in_one_t=as.matrix(t(df_in_one))
    
    print("Get id df")
    df_diag_match=data.frame(COH_sample=row.names(df_in_one_t),stringsAsFactors = F) %>% left_join(diag_in) %>% select(COH_sample,diag)
    
    print(paste0("Run find neighbors, FeatureN=",FeatureN_one))
    neighborIdx = find_neighbors(df_in_one_t, k=80)
    neighbor_subtype = as.data.frame(t(apply(neighborIdx, 1, function(x) as.vector(df_diag_match$diag[x]))))
    neighbor_subtype_1 = neighbor_subtype[,-1] %>% mutate(COH_sample=df_diag_match$COH_sample,FeatureN=FeatureN_one)
    neighbor_subtype_2=neighbor_subtype_1[neighbor_subtype_1$COH_sample==COH_sample,c("FeatureN",paste0("V",2:21))]
    neighbor_subtype_2
  }))
  df_neighbor_featureN
}

draw_neighbor_heatmap=function(df_in,inset_x= -0.24,margin_right=11){
  
  row.names(df_in)=df_in$FeatureN
  df_in$FeatureN=NULL
  
  nrow_list=1:dim(df_in)[1]
  ncol_list=1:dim(df_in)[2]
  
  par(mar=c(5.1, 4.1, 4.1, margin_right), xpd=TRUE)
  plot(c(0, dim(df_in)[2]), c(0, dim(df_in)[1]),  type= "n", xaxt="n",yaxt="n",
       xlab = "Neighbors", ylab = "Top MAD gene number"
  )
  
  for(nrow_one in nrow_list){
    for(ncol_one in ncol_list){
      value_one=as.character(df_in[nrow_one,ncol_one])
      # print(nrow_one);print(ncol_one);print(value_one)
      col_one=subtypeCol[value_one]
      rect(ncol_one-1, nrow_one-1,ncol_one,nrow_one,  col = col_one) 
    }
  }
  
  axis(1,at=ncol_list-0.5,labels = ncol_list)
  axis(2,at=nrow_list-0.5,labels = row.names(df_in),las=2)
  
  subtypeCol=sort(subtypeCol)
  
  legend("topright", inset=c(inset_x,0), title="Subtype",
         legend=names(subtypeCol[names(subtypeCol) %in% as.character(unlist(df_in))]), 
         fill=subtypeCol[names(subtypeCol) %in% as.character(unlist(df_in))])
}

Phenograph_pred=function(df_in,diag_in,top_gene,FeatureN_one,neighbor_k,ratio_cutoff=0.5){
  # print("run phenograph")
  # 
  # print("get feature df")
  input_gene=top_gene[1:FeatureN_one]
  df_in_one=df_in[row.names(df_in) %in% input_gene,]
  df_in_one_t=as.matrix(t(df_in_one))
  
  # print("Get id df")
  df_diag_match=data.frame(COH_sample=row.names(df_in_one_t),stringsAsFactors = F) %>% left_join(diag_in) %>% select(COH_sample,diag) %>% mutate(obs=1:n())
  
  print(paste0("Run Phenograph, FeatureN=",FeatureN_one,", NeighborN=",neighbor_k))
  
  set.seed(10)
  PG_out=Rphenograph(df_in_one_t, neighbor_k)
  
  df_cluster=data.frame(
    obs=as.numeric(names(membership(PG_out[[2]]))),
    cluster=as.vector(membership(PG_out[[2]])),
    stringsAsFactors = F
  )
  
  df_cluster_detail=df_diag_match %>% left_join(df_cluster) %>% 
    group_by(diag, cluster) %>% mutate(N_diagCluster=n()) %>% 
    group_by(cluster) %>%  mutate(N_cluster=n()) %>%
    mutate(ratio = round(N_diagCluster/N_cluster, digits = 2)) %>% 
    arrange(cluster, desc(N_cluster))
  
  df_cluster_diag=df_cluster_detail %>% select(diag,cluster,ratio,N_cluster,N_diagCluster,) %>% distinct() %>% 
    filter(ratio>ratio_cutoff) %>% group_by(diag) %>% arrange(desc(N_diagCluster)) %>% 
    mutate(diag_pred=diag,
           diagCluster_index=1:n(),
           diagCluster_index_max=n(),
           diag_pred_granular=ifelse(diagCluster_index_max==1,diag,paste0(diag,"_",diagCluster_index)),
           diag_pred_freq=paste0(diag_pred,"(",N_diagCluster,";",ratio*100,"%)"),
           diagCluster_index=NULL,diagCluster_index_max=NULL,N_cluster=NULL,N_diagCluster=NULL)
  df_cluster_diag$diag=NULL
  
  
  df_diag_pred=df_cluster_detail %>% 
    left_join(df_cluster_diag %>% mutate(ratio=NULL)) %>%
    mutate(obs=NULL,FeatureN=FeatureN_one,top_neighborN=neighbor_k)
  
  df_diag_pred$diag_pred=ifelse(is.na(df_diag_pred$diag_pred),"FailPrediction",df_diag_pred$diag_pred)
  df_diag_pred
}

Phenograph_pred_list=function(df_in,diag_in,top_gene,FeatureN_list,top_neighborN_list,ratio_cutoff=0.5){
  df_Phenograph_pred=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
    bind_rows(lapply(top_neighborN_list, function(top_neighborN_one){
      Phenograph_pred(df_in,diag_in,top_gene,FeatureN_one = FeatureN_one,neighbor_k = top_neighborN_one)
    }))
  }))
  df_Phenograph_pred
}

get_Phenograph_pred_sum=function(df_pred,top_neighborN_one){
  df_pred_=df_pred %>% filter(top_neighborN==top_neighborN_one) 
  sample_list=unique(df_pred_$COH_sample)
  
  df_pred_sum=bind_rows(lapply(sample_list, function(COH_sample_one){
    df_pred_one=df_pred_ %>% filter(COH_sample==COH_sample_one) %>%
      mutate(label=paste0("FeatureN",FeatureN))
    
    df_pred_one$diag_pred_granular[df_pred_one$diag_pred_granular==""]="FailPrediction"
    df_pred_one$diag_pred_freq[df_pred_one$diag_pred_freq==""]="FailPrediction"
    
    df_pred_one_t=dcast(df_pred_one,COH_sample~label,value.var = "diag_pred_freq")
    df_pred_one_t$pred=get_predDiag(df_pred_one$diag_pred)$pred_GEhalf
    df_pred_one_t$pred_freq=get_predDiag(df_pred_one$diag_pred)$pred_freq
    
    df_pred_one_t$pred_granular=get_predDiag(df_pred_one$diag_pred_granular)$pred_GEhalf
    df_pred_one_t$pred_granular_freq=get_predDiag(df_pred_one$diag_pred_granular)$pred_freq
    
    df_pred_one_t=df_pred_one_t[c("COH_sample","pred","pred_freq",'pred_granular',"pred_granular_freq",df_pred_one$label)]
    df_pred_one_t
  }))
  df_pred_sum
}




# draw heatmap from vst data
Phenograph_pred_FromVst=function(file_vst,vst_ref,diag_in_,top_gene,FeatureN_one,neighbor_k,ratio_cutoff=0.5){
 
  #prepare diag_ref
  # df_vst_file=read_tsv(file_vst)
  df_vst_file=file_vst
  COH_sample=names(df_vst_file)[2]
  cat(paste0("Run analysis for ",COH_sample,"\n\n"))
  
  diag_in_1=diag_in_[!diag_in_$COH_sample %in% COH_sample,]
  diag_ref=rbind(
    diag_in_1 %>% mutate(source="ref") %>% select(COH_sample,diag,source),
    data.frame(
      COH_sample=COH_sample,
      diag="_Prediction",
      source="for_prediction",
      stringsAsFactors = F
    )
  )
  
  vst_ref$feature=row.names(vst_ref)
  df_vst=vst_ref %>% left_join(df_vst_file)
  row.names(df_vst)=df_vst$feature
  df_vst$feature=NULL
  

  # print("run phenograph")
  # 
  # print("get feature df")
  input_gene=top_gene[1:FeatureN_one]
  df_vst_one=df_vst[row.names(df_vst) %in% input_gene,]
  df_vst_one_t=as.matrix(t(df_vst_one))
  
  # print("Get id df")
  df_diag_match=data.frame(COH_sample=row.names(df_vst_one_t),stringsAsFactors = F) %>% left_join(diag_in_) %>% select(COH_sample,diag) %>% mutate(obs=1:n())
  
  print(paste0("Run Phenograph, FeatureN=",FeatureN_one,", NeighborN=",neighbor_k))
  
  set.seed(10)
  PG_out=Rphenograph(df_vst_one_t, neighbor_k)
  
  df_cluster=data.frame(
    obs=as.numeric(names(membership(PG_out[[2]]))),
    cluster=as.vector(membership(PG_out[[2]])),
    stringsAsFactors = F
  )
  
  df_cluster_detail=df_diag_match %>% left_join(df_cluster) %>% 
    group_by(diag, cluster) %>% mutate(N_diagCluster=n()) %>% 
    group_by(cluster) %>%  mutate(N_cluster=n()) %>%
    mutate(ratio = round(N_diagCluster/N_cluster, digits = 2)) %>% 
    arrange(cluster, desc(N_cluster))
  
  df_cluster_diag=df_cluster_detail %>% select(diag,cluster,ratio,N_cluster,N_diagCluster,) %>% distinct() %>% 
    filter(ratio>ratio_cutoff) %>% group_by(diag) %>% arrange(desc(N_diagCluster)) %>% 
    mutate(diag_pred=diag,
           diagCluster_index=1:n(),
           diagCluster_index_max=n(),
           diag_pred_granular=ifelse(diagCluster_index_max==1,diag,paste0(diag,"_",diagCluster_index)),
           diag_pred_freq=paste0(diag_pred,"(",N_diagCluster,";",ratio*100,"%)"),
           diagCluster_index=NULL,diagCluster_index_max=NULL,N_cluster=NULL,N_diagCluster=NULL)
  df_cluster_diag$diag=NULL
  
  
  df_diag_pred=df_cluster_detail %>% 
    left_join(df_cluster_diag %>% mutate(ratio=NULL)) %>%
    mutate(obs=NULL,FeatureN=FeatureN_one,top_neighborN=neighbor_k)
  
  df_diag_pred$diag_pred=ifelse(is.na(df_diag_pred$diag_pred),"FailPrediction",df_diag_pred$diag_pred)
  df_diag_pred
}

Phenograph_pred_list_FromVst=function(file_vst,vst_ref,diag_in_,top_gene,FeatureN_list,top_neighborN_list,ratio_cutoff=0.5){
  df_Phenograph_pred_=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
    bind_rows(lapply(top_neighborN_list, function(top_neighborN_one){
      Phenograph_pred_FromVst(file_vst,vst_ref,diag_in_,top_gene,FeatureN_one = FeatureN_one,neighbor_k = top_neighborN_one)
    }))
  }))
  df_Phenograph_pred_
}

# Phenograph Prediction
Phenograph_pred_out=function(file_vst,vst_ref,diag_in_,top_gene,FeatureN_list,top_neighborN_list,ratio_cutoff=0.5){
  df_Phenograph_pred=Phenograph_pred_list_FromVst(file_vst,vst_ref,diag_in_,top_gene,FeatureN_list,top_neighborN_list,ratio_cutoff=0.5)
  df_vst_file=file_vst
  COH_sample=names(df_vst_file)[2]
  df_Phenograph_pred1=df_Phenograph_pred[df_Phenograph_pred$COH_sample==COH_sample,]
  df_pred_out=bind_rows(
    get_Phenograph_pred_sum(df_Phenograph_pred1,5) %>% mutate(top_neighborN=5),
    get_Phenograph_pred_sum(df_Phenograph_pred1,10) %>% mutate(top_neighborN=10),
    get_Phenograph_pred_sum(df_Phenograph_pred1,20) %>% mutate(top_neighborN=20)
  )
  df_pred_out
}

draw_neighbor_heatmap_FromVst=function(file_vst,vst_ref,diag_in_,top_gene,FeatureN_list){
  #prepare diag_ref
  # df_vst_file=read_tsv(file_vst)
  df_vst_file=file_vst
  COH_sample=names(df_vst_file)[2]
  cat(paste0("Run analysis for ",COH_sample,"\n\n"))
  
  diag_in_1=diag_in_[!diag_in_$COH_sample %in% COH_sample,]
  diag_ref=rbind(
    diag_in_1 %>% mutate(source="ref") %>% select(COH_sample,diag,source),
    data.frame(
      COH_sample=COH_sample,
      diag="_Prediction",
      source="for_prediction",
      stringsAsFactors = F
    )
  )
  vst_ref$feature=row.names(vst_ref)
  df_vst=vst_ref %>% left_join(df_vst_file)
  row.names(df_vst)=df_vst$feature  
  df_vst$feature=NULL
  
  df_top_neighbor=get_top_neighbor(COH_sample,df_vst,diag_in_,top_gene,FeatureN_list)
  draw_neighbor_heatmap(df_top_neighbor,inset_x= -0.24,margin_right=11)
}

run_tsne_panel_basedOnTopGeneList_3D=function(indata,genelist,n_pool,i_pool){
  df_tsne_panel=bind_rows(lapply(n_pool, function(n){
    lapply(i_pool, function(i){
      feature_in=head(genelist,n)
      indata1=indata[feature_in,]
      df_tsne_temp=run_tsne(indata1,i,3)
      df_tsne_temp
    })
  }))
  df_tsne_panel
}

knn_prediction_ontSNE=function(data_in){
  perplexity_list=unique(data_in$PerplexityN)
  FeatureN_list=unique(data_in$FeatureN)
  
  pred_all=bind_rows(lapply(perplexity_list, function(perplexity_one){
    pred_=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
      for_knn_one=data_in[data_in$PerplexityN==perplexity_one & data_in$FeatureN==FeatureN_one,]
      
      df_train=for_knn_one[!for_knn_one$diag=="_Prediction",]
      df_test=for_knn_one[for_knn_one$diag=="_Prediction",]
      
      pred_one=knn_prediction_train_test_df(df_train,df_test,c("X","Y","Z"),"COH_sample",5)
      pred_one$PerplexityN=perplexity_one
      pred_one$FeatureN=FeatureN_one
      pred_one
    }))
    pred_
  }))
  pred_all
}

get_pred_summary=function(predict_out,COH_sample){
  freq=table(predict_out)
  freq=sort(freq,decreasing = T)
  
  df_knn_out=data.frame(
    Sample=COH_sample,
    KNN_pred=ifelse(any(freq > sum(freq)/2),names(freq)[freq > sum(freq)/2],"Fail_prediction"),
    KNN_pred_freq=paste0(paste0(names(freq),": ",freq," (",sprintf("%.1f",100*freq/sum(freq)),"%)"),collapse = "; "),
    stringsAsFactors = F
  )
  df_knn_out
}

knn_prediction_train_test_df=function(df_train,df_test,var_x,id_var,kNN_N){
  formula_in=formula(paste0("diag","~",paste0(var_x,collapse = "+")))
  KNN_train_fit=caret::train(formula_in, data = df_train, method = "knn",
                             preProcess = c("center", "scale"),
                             tuneGrid = expand.grid(k = c(kNN_N)))
  knn_pred = predict(KNN_train_fit,newdata = df_test) %>% as.vector()
  outdata=df_test[id_var] %>% mutate(pred_knn=knn_pred)
  outdata
}

# KNN Prediction
KNN_pred_out=function(file_vst,COH_sample,vst_ref,diag_in_,feature_in,n_pool,i_pool){
  # n_pool=seq(200,2000,100)
  # i_pool=c(20,30)
  # COH_sample=substr(basename(file_vst),1,12)
  df_vstBatch=vst_ref[unique(c(diag_in_$COH_sample,COH_sample))]
  for_tSNE=df_vstBatch[feature_in,]

  # df_vst_file=file_vst
  # COH_sample=names(df_vst_file)[2]

  cat(paste0("Run analysis for ",COH_sample,"\n\n"))

  diag_in_1=diag_in_[!diag_in_$COH_sample %in% COH_sample,]
  diag_ref=rbind(
    diag_in_1 %>% mutate(source="ref") %>% select(COH_sample,diag,source),
    data.frame(
      COH_sample=COH_sample,
      diag="_Prediction",
      source="for_prediction",
      stringsAsFactors = F
    )
  )

  for_knn=run_tsne_panel_basedOnTopGeneList_3D(for_tSNE,feature_in,n_pool,i_pool)
  df_forKNN_diag=for_knn %>% left_join(diag_ref)

  pred_all=knn_prediction_ontSNE(df_forKNN_diag)

  df_knn_out=get_pred_summary(pred_all$pred_knn,COH_sample)
  df_knn_out
}




# test------------------------

# vst_ref=vst_gene_1499
# diag_in_=diag_1449
# top_gene=top_feature_list_1499
# FeatureN_one=50
# neighbor_k=10
# FeatureN_list=seq(50,2000,900)
# top_neighborN_list=c(5,10)
# diag_in_=diag_1409
# top_gene=top_feature_list_1409

# plot_out=draw_tsne_FromVst(runVST_FromFile(file),diag_1449,vst_gene_1499,top_feature_list_1499)
# plot_out=draw_tsne_FromVst(x,diag_1409,z,top_feature_list_1409)
#   
# plot_heatmap=draw_neighbor_heatmap_FromVst(runVST_FromFile(file),for_tsne_1409_,diag_1409,top_feature_list_1409,FeatureN_list = seq(50,2000,50),top_neighborN_list = c(5,10,20))
# y=pred_out(x,for_tsne_1409_,diag_1409,top_feature_list_1409,FeatureN_list = seq(50,2000,50),top_neighborN_list = c(5,10,20),ratio_cutoff=0.5)
# tsne_table=tsne_df_FromVst(file_vst,diag_1409,for_tsne_1409_,top_feature_list_1409,20)
# data_plot=tsne_table %>% left_join(diag_1409) %>% na.omit()
# diag_col=diagCol[names(diagCol) %in% diag_1409$diag]
# plot=ggplot() +
#                 xlab(paste0("tSNE dimension 1\ntop gene num = ", 800, "; perplexity=", 20)) +    ylab("tSNE dimension 2")+
#                 theme_bw() +
# 
#                 geom_point(data=subset(data_plot, diag != "_Prediction" ), aes(X, Y, color=diag), size=1)+
#                 geom_point(data=subset(data_plot, diag == "_Prediction" ), aes(X, Y, color=diag), size=3, alpha=0.75)+
#                 geom_text_repel(data=subset(data_plot, diag == "_Prediction" ), aes(X, Y, label=COH_sample), fontface = "bold")+
#                 scale_color_manual(values = diag_col) +
# 
#                 scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = 10),1)) +
#                 scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = 10),1)) +
# 
#                 theme(axis.title=element_text(size=10),
#                       #axis.text=element_blank(),
#                       #axis.ticks=element_blank(),
#                       legend.text = element_text(size=10),
#                       legend.title = element_text(size=10,face="bold")
#                 ) +
#                 scale_fill_discrete(name="Experimental\nCondition")+
#                 guides(color = guide_legend(override.aes = list(size = 1)))
# plot                 
# 
# 
# file="test/COH000066_D1.DUX4patched.HTSeq"
# x=runVST_FromFile(file)


