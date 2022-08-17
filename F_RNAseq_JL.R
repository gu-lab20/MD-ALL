#install.packages(c("Boruta"))

# library(Rtsne)
# library(caret)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# library(ggrepel)
# library(ggpubr)
# library(DESeq2)

generate_html_stastics=function(indata,file,var_list){
  summarytools::view(summarytools::dfSummary(indata[var_list]),file,method = "panel")
}
#write_tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

#read tsv --------------
read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}

get_PU=function(dir,postfix){
  fq_1stlines=system(paste0("ls ",dir_fq_raw,"*.",postfix,"|xargs -i echo 'paste <(ls {}) <(zcat {}|head -n1)'|bash|sed 's/ /\t/g'"),intern = T)
  df_fq_raw=data.frame(
    V1=word(fq_1stlines,1,sep="\t"),
    V2=word(fq_1stlines,2,sep="\t"),
    V3=word(fq_1stlines,3,sep="\t"),
    stringsAsFactors = F
  ) %>%
    mutate(
      PU=paste0(gsub("@","",word(V2,1,sep=":")),":",word(V2,2,sep=":"),":",word(V2,3,sep=":"),".",word(V2,4,sep=":"),".",word(V3,-1,sep=":")),
      i5=word(word(PU,-1,sep="[.]"),1,sep="[+]")
    )
}


# gene list in --------------------------------
file_gene_infor="//isi-dcnl/user_data/zgu_grp/DryLab/bin/jiliu/R/reference/Homo_sapiens.GRCh38.V102.summary.txt"
dux4_genelist=readLines("//isi-dcnl/user_data/zgu_grp/DryLab/bin/jiliu/R/reference/DUX4_gene.list")

info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding" | info_gtf$gene_name %in% dux4_genelist]
geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
gene_in=codinggeneList[!codinggeneList %in% geneXYM]

rm(file_gene_infor)
rm(codinggeneList)
rm(geneXYM)

# get vst matrix --------------------------------
get_vst_matrix=function(files){
  for(file in files){
    print(match(file,files))
    id=word(basename(file),1,sep="[.]")
    count_one=read.table(file,header = T,stringsAsFactors = F) 
    if(match(file,files)==1){vst_all=count_one}
    if(match(file,files)>1){vst_all=vst_all %>% left_join(count_one)}
  }
  vst_all
}


########################################  Calculate freq, no group variable, one or multiple variables ########################################

#value_x=unlist(indata["hp"])
#var_x="hp"

cal_freq_NoGroup=function(indata,var_list,digit_n){
  bind_rows(lapply(var_list, function(var_x){
    value_x=unlist(indata[var_x])
    out1=data.frame(table(value_x))
    out2=data.frame(
      Variable=var_x,
      variable_levels=paste0(as.character(unlist(out1[,1]))),
      N_used=as.character(sum(!is.na(value_x))),
      All=
        paste0(unlist(out1[,2]),
               " (",
               sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),
               "%)"
        ),
      stringsAsFactors = F
    )
    out2$N_used[2:nrow(out2)]=NA
    bind_rows(data.frame(Variable=var_x,variable_levels=var_x,stringsAsFactors = F),out2) 
  }))
  
}


# Creat Deseq object from dataset -----------------------------------------------------------------------
#data_diag is can be a filtered sample list. 
#this function will select the collapsed samples. Default sample name should be "COH_sample"
#var_effect can be vetor of diagnosis or/and batch

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
  
  htseq_object = DESeqDataSetFromMatrix(countData = matrix_in1,
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

rlogDist=function(rlogDF, outDir, fout, threshLine=0.3){
  samples=colnames(rlogDF)
  sampleNum=length(samples);
  sampleCol=rainbow(sampleNum, alpha = 0.6)
  dens.x=c()
  dens.y=c()
  maxY=0;
  for( i in 1:sampleNum){
    sampleI=rlogDF[,i]
    density=density(sampleI, from=-5, to=20 )
    maxY=ifelse(maxY > max(density$y), maxY, max(density$y));
    dens.x=cbind(dens.x, density$x)
    dens.y=cbind(dens.y, density$y)
  }
  
  
  pdf(file = file.path(outDir,fout), width = 10, height = 12)
  par(mfrow=c(2,1), mar=c(5,5,2,2))
  plot(1,1, type='n', ylim=c(0,maxY), xlim=c(5,18), xlab = "Normalization value", ylab="Density", main = "" ,cex.lab=1.5,cex.axis=1.5)
  for( i in 1:sampleNum){
    points(dens.x[,i], dens.y[,i], type = "l", col=sampleCol[i]);
  }
  densX=rowMedians(dens.x)
  densY=rowMedians(dens.y)
  points(densX, densY, type='l', lwd = 2)
  #dist
  distDens=c()
  for( r in 1:sampleNum){
    distI=dist(rbind(densY, dens.y[,r]))
    distDens=c(distDens, as.vector(distI))
  }
  names(distDens)=samples
  write.table(x = distDens, file = file.path(outDir,"dist2Median.txt"), quote=F, row.names = T, sep = "\t")
  distDensSrt=sort(distDens, decreasing = F)
  plot(x = 1:sampleNum, y=distDensSrt, xlab = "Sample index", ylab="Distance to the median line",cex.lab=1.5,cex.axis=1.5)
  abline(h = threshLine)
  text(x = 5, y=threshLine, labels = paste("dist =",threshLine), pos = 3)
  dev.off()
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}



draw_density_plot=function(indata,dir_out,file_density,dist_cutoff){
  if(!dir.exists(dir_out)){dir.create(dir_out)}
  rlogDist(indata, dir_out,fout = file_density, threshLine = dist_cutoff)
  data_dist=read.table(paste0(dir_out,"/dist2Median.txt"),skip = 1,stringsAsFactors = F)
  names(data_dist)=c("sample_id","dist")
  data_dist
}


get_high_MAD_gene=function(indata,N_genes){
  geneMadOrd = names(sort(apply(indata,1,mad),decreasing = T))
  topMadGenes=geneMadOrd[1:N_genes]
  outdata=indata[topMadGenes,]
}

get_topMAD_features=function(df_in,N_genes){
  geneMadOrd = names(sort(apply(df_in,1,mad),decreasing = T))
  if(length(geneMadOrd)>N_genes){topMadGenes=geneMadOrd[1:N_genes]}
  if(length(geneMadOrd)<=N_genes){topMadGenes=geneMadOrd}
  topMadGenes
}


remove_high_corr_genes=function(indata,highCorrCutoff){
  set.seed(0)
  modelDF=t(indata) %>% as.data.frame()
  corrMatrix = cor(modelDF)
  highCorrIdx = findCorrelation(corrMatrix, cutoff=highCorrCutoff)
  highCorrGenes=colnames(modelDF[, highCorrIdx])
  
  modelDFnoCorr=modelDF[, -highCorrIdx]
  noHighCorrGeneDF=t(modelDFnoCorr)
  noHighCorrGeneDF
}


remove_nonCOdingXYM_genes=function(indata){
  info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
  codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding"]
  geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
  gene_in=codinggeneList[!codinggeneList %in% geneXYM]
  outdata=indata[row.names(indata) %in% gene_in,]
}

#indata=deseq_panel_vst_batchC1
#data_diag=diag_in
#sample_id="COH_sample"
#batch_var="batch3"

#batch effect correction -------------------------
correct_batch=function(indata,diag_in,sample_id){
  print(dim(indata))
  print("running batch correction")
  
  sample_matrix=names(indata)
  sample_diag=unlist(diag_in[sample_id])
  
  sample_in=sort(sample_matrix[sample_matrix %in% sample_diag])
  
  data_diag1=diag_in[unlist(diag_in[sample_id]) %in% sample_in,] %>% arrange(sample_id) %>%
    mutate(batch1=word(batch,1,sep="_"),batch2=word(batch,2,sep="_"),batch3=word(batch,3,sep="_"))
  
  indata1=indata[,unlist(data_diag1[sample_id])]
  
  modcombat = model.matrix(~1, data=data_diag1)
  
  print("1. Correct for mRNA/totalRNA")
  out_data1 = sva::ComBat(dat=as.matrix(indata1), batch=data_diag1$batch2, mod=modcombat, par.prior=TRUE) %>% as.data.frame() 
  print("2. Correct for stranded/non-stranded")
  out_data2 = sva::ComBat(dat=as.matrix(out_data1), batch=data_diag1$batch1, mod=modcombat, par.prior=TRUE) %>% as.data.frame() %>% round(., digits = 10)
  print("3. Correct for length:100bp/150bp/75bp etc.")
  out_data3 = sva::ComBat(dat=as.matrix(out_data2), batch=data_diag1$batch3, mod=modcombat, par.prior=TRUE) %>% as.data.frame() %>% round(., digits = 5)
  
  print(paste0("outdata dim: ", paste0(dim(out_data3),collapse = ",")))
  out_data3
}

correct_batch_var=function(indata,diag_in,sample_id,batch_var){
  
  print(dim(indata))
  print(paste0("running batch correction for ",batch_var))
  
  sample_matrix=colnames(data.frame(indata))
  sample_diag=unlist(diag_in[sample_id])
  sample_id_overlap=sample_matrix[sample_matrix %in% sample_diag]
  
  sample_id_in=sort(sample_id_overlap)
  
  data_diag1=diag_in[unlist(diag_in[sample_id]) %in% sample_id_in,] %>% arrange(sample_id) %>%
    mutate(batch_c=unlist(diag_in[batch_var]))
  
  indata1=indata[,unlist(data_diag1[sample_id])]
  
  modcombat = model.matrix(~1, data=data_diag1)
  
  out_data3 = sva::ComBat(dat=as.matrix(indata1), batch=data_diag1$batch_c, mod=modcombat, par.prior=TRUE) %>% 
    as.data.frame() %>% round(., digits = 8)
  
  print(paste0("outdata dim: ", paste0(dim(out_data3),collapse = ",")))
  out_data3
}

remove_lowexpression_genes=function(indata,cutoff){
  out_data=indata[apply(indata, 1, max) >= cutoff,]
}


run_boruta=function(matrix_in,diag_in,max_run){
  data_with_diagnosis=as.data.frame(t(matrix_in)) %>%
    mutate(COH_sample=colnames(matrix_in)) %>% 
    left_join(diag_in[c("COH_sample","diag")]) %>%
    mutate(COH_sample=NULL,diag=as.factor(diag)) %>%  na.omit()
  
  print(dim(data_with_diagnosis))
  
  set.seed(10)
  boruta_out = Boruta(diag~., data = data_with_diagnosis, doTrace = 2,maxRuns=max_run)
  boruta_out_fixed = TentativeRoughFix(boruta_out)
  
  importance_boruta=attStats(boruta_out)
  importance_boruta_fixed=attStats(boruta_out_fixed)
  
  data_boruta_confirmed=matrix_in[row.names(importance_boruta)[importance_boruta$decision=="Confirmed"],]
  data_boruta_confirmed_fixed=matrix_in[row.names(importance_boruta_fixed)[importance_boruta_fixed$decision=="Confirmed"],]
  
  return(list(
    fit_boruta=boruta_out,
    importance_boruta=importance_boruta,importance_boruta_fixed=importance_boruta_fixed,
    data_boruta_confirmed=data_boruta_confirmed,data_boruta_confirmed_fixed=data_boruta_confirmed_fixed))
}


# for tsne ----------------------
from_vst2tsnePanel_GEP_df=function(df_in,diag_in){
  df_vst_batchC=correct_batch(df_in,diag_in,"COH_sample","batch")
  df_out=df_vst_batchC[row.names(df_vst_batchC) %in% topMADgene,]
  df_out
}

from_vst2tsnePanel_df=function(df_in,diag_in){
  
  df_vst_batchC=correct_batch(df_in,diag_in,"COH_sample","batch")
  
  df_vst_batchC_codingNoXYM=df_vst_batchC[row.names(df_vst_batchC) %in% gene_in,]
  
  df_vst_batchC_codingNoXYM_highE=remove_lowexpression_genes(df_vst_batchC_codingNoXYM,5)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE,3000)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr=remove_high_corr_genes(df_vst_batchC_codingNoXYM_highE_MAD3000,0.75)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr,1500)
  print(dim(df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000))
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000
}


from_vstbatch2tsnePanel_df=function(df_in,gene_in){
  print("codingNoXYM")
  df_vst_batchC_codingNoXYM=df_in[row.names(df_in) %in% gene_in,]
  print("remove low exp")
  df_vst_batchC_codingNoXYM_highE=remove_lowexpression_genes(df_vst_batchC_codingNoXYM,5)
  print("top MAD10000")
  df_vst_batchC_codingNoXYM_highE_MAD10000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE,10000)
  print("remove corr")
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr=remove_high_corr_genes(df_vst_batchC_codingNoXYM_highE_MAD10000,0.75)
  
  print("top MAD5000")
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD5000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr,5000)
  
  print(paste0("FeatureNum, SampleNum: ",paste0(dim(df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD5000),collapse = ", ")))
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD5000
}


#run tsne base function -----------------------------
run_tsne_noRunPCA=function(indata,perplexity,out_dims){
  cat("Run tSNE: Used Feature N=", dim(indata)[2], "; Used Sample N=", dim(indata)[1],"; Perplexity=", perplexity,  "\n")
  
  # Sets seed for reproducibility
  set.seed(10) 
  
  # Run tSNE
  rlogTsneOut = Rtsne(indata, 
                      dims = out_dims, 
                      perplexity = perplexity, 
                      theta = 0.5, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=F, 
                      pca=F,
                      num_threads = 4)
  #summary output
  df_id=data.frame(COH_sample=row.names(indata)) 
  
  df_tsne=as.data.frame(rlogTsneOut$Y)
  names(df_tsne)=c("X","Y","Z")[1:ncol(df_tsne)]
  df_tsne$COH_sample=row.names(indata)
  
  df_tsne=df_id %>% left_join(df_tsne) %>%
    mutate(perplexity=perplexity,
           FeatureN=dim(indata)[2]) 
  df_tsne
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

run_tsne_panel_basedOnTopGeneList_2D=function(indata,genelist,n_pool,i_pool){
  df_tsne_panel=bind_rows(lapply(n_pool, function(n){
    lapply(i_pool, function(i){
      feature_in=head(genelist,n)
      indata1=indata[feature_in,]
      df_tsne_temp=run_tsne(indata1,i,2)
      df_tsne_temp
    })
  }))
  df_tsne_panel
}

#tsne based on MAD --------------------------
draw_tsne_MAD=function(indata,data_diag,diag_col,N,i,axis_by,size=1){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1\ntop gene num = ", N, "; perplexity=", i)) +    ylab("tSNE dimension 2")+
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


draw_tsne_MAD_panel=function(indata,diag_in_tsne,diagCol=subtypeCol,i_pool,file_path,file_name,col.n=6,legend_loc="bottom",size=1){
  n_list=unique(indata$FeatureN)
  n_plot=length(n_list)*length(i_pool)
  print(paste0("plot number:",n_plot))
  col.r=n_plot/col.n
  
  tsnePanelList=list()
  for (i in i_pool){
    for(N in n_list){
      run_index=(match(N,n_list)-1)*length(i_pool)+match(i,i_pool)
      print(run_index)
      df_tsne_in=indata %>% filter(FeatureN==N & PerplexityN==i)
      tsnePanelList[[run_index]]=draw_tsne_MAD(df_tsne_in,diag_in_tsne,diagCol,N,i,10,size)
    }}
  
  tsnePanel=ggarrange(plotlist=tsnePanelList,ncol=col.n,nrow=col.r,common.legend = T,legend=legend_loc)
  file_=file.path(file_path,paste0("tSNE_panel_",file_name))
  ggsave(paste0(file_,".pdf"),width = 5*col.n,height = 4.5*col.r)
  ggsave(paste0(file_,".png"),width = 5*col.n,height = 4.5*col.r,dpi = 300)
}



draw_tsne=function(indata,data_diag,diag_col=subtypeCol,N,i,axis_by=10){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1\ntop gene num = ", N, "; perplexity=", i)) +
    ylab("tSNE dimension 2")+
    theme_bw() +
    geom_point(data=data_plot, aes(X, Y, color=diag)) +
    scale_color_manual(values = diag_col) +
    scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = axis_by),1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          legend.text = element_text( size=15),
          legend.title = element_text( size=15,face="bold")
    ) 
}

draw_tsne_fromPred2D=function(file_pred2D,FeatureN_one,PerplexityN_one){
  load("../1.prepare_diag/2.data_temp/diag_2239.rdata")
  df_tsne_2D=read_tsv(file_pred2D)
  
  df_tsne_2D_1=df_tsne_2D %>% filter(FeatureN==FeatureN_one & PerplexityN==PerplexityN_one)
  diag_1367=diag_2239 %>% filter(COH_sample %in% df_tsne_2D_1$COH_sample)
  
  draw_tsne(df_tsne_2D_1,diag_1367,N = 800,i = 30)
}



from_diag_to_tsne=function(diag_in){
  load("../2.prepare_matrix/2.data_temp/vst_gene_NormalBALL_2411.rdata")
  # vstBatch_gene_=correct_batch_var(vst_gene_NormalBALL_2411,diag_in,"COH_sample","batch")
  vstBatch_gene_=correct_batch(vst_gene_NormalBALL_2411,diag_in,"COH_sample")
  for_tsne_=from_vstbatch2tsnePanel_df(vstBatch_gene_,gene_in)
  
  #test tsne
  # top800=get_topMAD_features(for_tsne_,800)
  # for_tsne_800=for_tsne_[row.names(for_tsne_) %in% top800,]
  #df_tsne=run_tsne(for_tsne_800,30,2)
  #col_in=subtypeCol[names(subtypeCol) %in% diag_in$diag]
  #draw_tsne(df_tsne,diag_in,col_in,800,30,10)
  
  #return(list(vstBatch_gene_=vstBatch_gene_,for_tsne_=for_tsne_,df_tsne=df_tsne))
  return(list(vstBatch_gene_=vstBatch_gene_,for_tsne_=for_tsne_))
}

#Draw tsne based on ggplot2 ----------------------------
draw_tsne_MAD=function(indata,data_diag,diag_col,N,i,axis_by,size=1){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1\ntop gene num = ", N, "; perplexity=", i)) +    ylab("tSNE dimension 2")+
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


draw_tsne_MAD_panel=function(indata,diag_in_tsne,diagCol=subtypeCol,i_pool,file_path,file_name,col.n=6,legend_loc="bottom",size=1){
  n_list=unique(indata$FeatureN)
  n_plot=length(n_list)*length(i_pool)
  print(paste0("plot number:",n_plot))
  col.r=n_plot/col.n
  
  tsnePanelList=list()
  for (i in i_pool){
    for(N in n_list){
      run_index=(match(N,n_list)-1)*length(i_pool)+match(i,i_pool)
      print(run_index)
      df_tsne_in=indata %>% filter(FeatureN==N & PerplexityN==i)
      tsnePanelList[[run_index]]=draw_tsne_MAD(df_tsne_in,diag_in_tsne,diagCol,N,i,10,size)
    }}
  
  tsnePanel=ggarrange(plotlist=tsnePanelList,ncol=col.n,nrow=col.r,common.legend = T,legend=legend_loc)
  file_=file.path(file_path,paste0("tSNE_panel_",file_name))
  ggsave(paste0(file_,".pdf"),width = 5*col.n,height = 4.5*col.r)
  ggsave(paste0(file_,".png"),width = 5*col.n,height = 4.5*col.r,dpi = 300)
}

draw_tsne=function(indata,data_diag,diag_col=subtypeCol,N,i,axis_by=10){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  diag_col=diag_col[names(diag_col) %in% data_diag$diag]
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1\ntop gene num = ", N, "; perplexity=", i)) +
    ylab("tSNE dimension 2")+
    theme_bw() +
    geom_point(data=data_plot, aes(X, Y, color=diag)) +
    scale_color_manual(values = diag_col) +
    scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = axis_by),1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          legend.text = element_text( size=15),
          legend.title = element_text( size=15,face="bold")
    ) 
}

draw_tsne_fromPred2D=function(file_pred2D,FeatureN_one,PerplexityN_one){
  load("../1.prepare_diag/2.data_temp/diag_2239.rdata")
  df_tsne_2D=read_tsv(file_pred2D)
  
  df_tsne_2D_1=df_tsne_2D %>% filter(FeatureN==FeatureN_one & PerplexityN==PerplexityN_one)
  diag_=diag_2239 %>% filter(COH_sample %in% df_tsne_2D_1$COH_sample)
  
  draw_tsne(df_tsne_2D_1,diag_,N = 800,i = 30)
}

from_diag_to_tsne=function(diag_in){
  load("../2.prepare_matrix/2.data_temp/vst_gene_NormalBALL_2411.rdata")
  # vstBatch_gene_=correct_batch_var(vst_gene_NormalBALL_2411,diag_in,"COH_sample","batch")
  vstBatch_gene_=correct_batch(vst_gene_NormalBALL_2411,diag_in,"COH_sample")
  for_tsne_=from_vstbatch2tsnePanel_df(vstBatch_gene_,gene_in)
  
  #test tsne
  # top800=get_topMAD_features(for_tsne_,800)
  # for_tsne_800=for_tsne_[row.names(for_tsne_) %in% top800,]
  #df_tsne=run_tsne(for_tsne_800,30,2)
  #col_in=subtypeCol[names(subtypeCol) %in% diag_in$diag]
  #draw_tsne(df_tsne,diag_in,col_in,800,30,10)
  
  #return(list(vstBatch_gene_=vstBatch_gene_,for_tsne_=for_tsne_,df_tsne=df_tsne))
  return(list(vstBatch_gene_=vstBatch_gene_,for_tsne_=for_tsne_))
}

#Draw tsne based on basic R graphics ----------------------------
#need X,Y,diag and col
draw_tsne_one=function(df_in,cex_point=0.5){
  par(mar = c(4.5,3.5,1,1))
  plot(df_in$X,df_in$Y,col=df_in$col,pch=20,
       ylab="",xlab="",cex=cex_point)
  xlab_lab=paste0("tSNE dimension 1\nTop MAD gene N = ",
                  unique(df_in$FeatureN)[1],"; Perplexity = ",unique(df_in$PerplexityN)[1])
  title(xlab=xlab_lab,line=3)
  title(ylab="tSNE dimension 2",line=2)
}

draw_tsne_legend=function(df_in,withPercentage=TRUE,ncol=4,cex_legend=1.5,position="center",header_legend=""){
  df_freq=as.data.frame(table(df_in$diag)) %>% 
    transmute(diag=as.character(Var1),
              Freq=Freq,
              percentage=paste0(diag,": ",Freq," (",round(Freq/nrow(df_in)*100,1),"%)"))
  df_legend=df_in %>% select(diag,col) %>% distinct() %>% left_join(df_freq) %>% arrange(desc(Freq))
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes = FALSE, xlab = "",ylab = "")
  
  if(withPercentage){legend(position,title = header_legend,legend = df_legend$percentage,col=df_legend$col,pch=20,bty = "n",ncol = ncol,cex = cex_legend)}
  if(!withPercentage){legend(position,title=header_legend,legend = df_legend$diag,col=df_legend$col,pch=20,bty = "n",ncol = ncol,cex = cex_legend)}
}

draw_tsne_one_withLegend=function(df_in,outfile_prefix){
  draw_plot=function(){
    layout(matrix(c(1,2),ncol=2),widths = c(1,0.38))
    par(mar = c(4.5,3.5,0.5,0))
    draw_tsne_one(df_in,cex_point=0.5)
    par(mar = c(0,0,1,0))
    draw_tsne_legend(df_in,withPercentage=F,ncol=1,cex_legend=0.6,position="topleft",header_legend = "Group")
  }
  #draw pdf
  pdf(paste0(outfile_prefix,".pdf"),width = 7.5,height = 5.5)
  draw_plot()
  dev.off()
  
  #draw png
  png(paste0(outfile_prefix,".png"),width = 2200,height = 1700,res=300)
  draw_plot()
  dev.off()
  
  draw_plot()
}

#col.n=2,height_png=2500,width_png=2000,height_pdf=10,width_pdf=8,cex_legend = 0.7
draw_tsne_panel=function(df_in,outfile_prefix,col.n=6,
                         width_pdf=20,height_pdf=16,width_png=5000,height_png=4000,
                         withPercentage=T,cex_legend=1.5,cex_point=0.5){
  PerplexityN_list=unique(df_in$PerplexityN)
  FeatureN_list=unique(df_in$FeatureN)
  
  PerplexityN_one=PerplexityN_list[1]
  FeatureN_one=FeatureN_list[1]
  
  pics_n=length(PerplexityN_list) * length(FeatureN_list)
  row.n=ceiling(pics_n/col.n)
  
  draw_plot=function(){
    layout(matrix(c(1:pics_n,rep(pics_n+1,col.n)),ncol=col.n,byrow = T),heights = c(rep(1,row.n),0.5))
    for(PerplexityN_one in PerplexityN_list){
      for(FeatureN_one in FeatureN_list){
        df_one=df_in[df_in$PerplexityN==PerplexityN_one & df_in$FeatureN==FeatureN_one,]
        draw_tsne_one(df_one)
      }}
    
    df_one_top=df_in[df_in$PerplexityN==PerplexityN_list[1] & df_in$FeatureN==FeatureN_list[1],]
    draw_tsne_legend(df_one_top,ncol = 5,withPercentage = withPercentage,cex_legend=cex_legend)
  }
  
  #draw pdf
  pdf(paste0(outfile_prefix,".pdf"),width = width_pdf,height = height_pdf)
  draw_plot()
  dev.off()
  
  #draw png
  png(paste0(outfile_prefix,".png"),width = width_png,height = height_png,res=300)
  draw_plot()
  dev.off()
}

#tsne based on boruta --------------------------

run_tsne_boruta=function(indata,N,i){
  gene_in=row.names(importance_boruta_in)[1:N]
  
  data_for_tsne=t(indata[gene_in,])
  
  
  cat("Perplexity:", i, "; top MAD gene number:", dim(data_for_tsne)[2], "\n");
  dim(data_for_tsne)
  # Run TSNE
  set.seed(10) # Sets seed for reproducibility
  rlogTsneOut = Rtsne(data_for_tsne, dims = 2, perplexity = i, 
                      theta = 0.1, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=T, 
                      num_threads = 4)
  
  tsneDF=data.frame(COH_sample=row.names(data_for_tsne)) %>% 
    mutate(X=rlogTsneOut$Y[,1], 
           Y=rlogTsneOut$Y[,2],
           perplexity=i,
           gene_n=N) 
  
  tsneDF
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


df_ball_col=data.frame(
  diag=names(subtypeCol),
  col=subtypeCol,
  obs=1:length(subtypeCol),
  stringsAsFactors = F
)

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
}

#draw_diag_barplot -----------------
draw_diag_barplot=function(diag_in,file_barplot){
  df_freq=data.frame(table(diag_in$diag)) %>% arrange(Var1) %>% arrange(Freq)
  names(df_freq)[1]="diag"
  df_freq=df_freq %>% left_join(df_ball_col) %>% 
    mutate(percetage=sprintf("%.1f",100*Freq/nrow(diag_in)),
           diag_percentage=paste0(diag," (",percetage,"%)"))
  
  max_x=max(df_freq$Freq)
  
  pdf(file_barplot,width = 7,height = 14)
  par(mar=c(8,10,8,6),oma=c(6, 4.1, 4.1, 2.1))
  plot_bar=barplot(df_freq$Freq,xlab="Frequency",xaxt='n',xlim=c(0,max_x),names.arg = df_freq$Var1,col=df_freq$col,horiz=TRUE,cex.names = 1.8)
  axis(2, at=plot_bar, labels=df_freq$diag_percentage, tick=FALSE, las=2, line=-0.5, cex.axis=1.2,font=2)
  axis(1,cex.axis=1.2,font=2)
  dev.off()
}

# generate_data_for_tsne_from_DeseqObejct_NoBatchCorrection------------------------------------
generate_data_for_tsne_from_DeseqObejct_NoBatchCorrection=function(in_object,file_gene_infor,dir_density,file_name_density){
  
  #VarianceStabilizingTransformation
  cat("Running vst\n")
  htseq_vst=DESeq2::varianceStabilizingTransformation(in_object,blind = T)
  htseq_rlog=round(data.frame(assay(htseq_vst)), digits = 2 )
  
  cat("Draw density plot\n")
  
  if(!dir.exists(dir_density)){dir.create(dir_density)}
  rlogDist(htseq_rlog, dir_density,fout = file_name_density, threshLine = 0.3)
  
  #extract coding gene
  cat("extract coding gene\n")
  info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
  codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding"]
  htseq_rlog=htseq_rlog[row.names(htseq_rlog) %in% codinggeneList,]
  
  
  #batch effect correction
  
  #modcombat = model.matrix(~1, data=data_diag)
  #htseq_correctBatch = ComBat(dat=as.matrix(htseq_rlog), batch=data_diag$batch, mod=modcombat, par.prior=TRUE) %>%
  #  as.data.frame() %>% round(., digits = 2)
  
  
  #remove genes on chr X, Y and MT
  cat("remove genes on chr X, Y and MT\n")
  geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
  htseq_noXYM=htseq_rlog[!row.names(htseq_rlog) %in% geneXYM,]
  
  #remove consistantly low-expressed genes
  htseq_highExpr=htseq_noXYM[apply(htseq_noXYM, 1, max) >= 5,]
  
  cat("running MAD\n")
  #remove low mad genes
  N=10000
  geneMadOrd = names(sort(apply(htseq_highExpr,1,mad),decreasing = T))
  topMadGenes=geneMadOrd[1:N]
  htseq_topMAD=htseq_highExpr[topMadGenes,]
  
  #remove high corr genes, find attributes that are highly corrected (ideally >0.75)
  cat("Removing correlated genes\n")
  indata=htseq_topMAD
  highCorrCutoff=0.75
  htseq_forPlot=remove_high_corr_genes(htseq_topMAD,0.75)
  dim(htseq_forPlot)
  htseq_forPlot
  
}

#PCA -----------------------------
cal_pc_basedOnTopMADFeature=function(df_in,top_feature,top_n){
  df_in1=df_in[row.names(df_in) %in% top_feature[1:top_n],]
  
  pca_fit=prcomp(df_in1)
  
  df_pc=data.frame(pca_fit$rotation)
  df_pc$NFeature=top_n
  df_pc$COH_sample=row.names(df_pc)
  df_pc
}

# KNN prediction -------------------------------------

knn_pred=function(indata,diag_in,test_level){
  
  df_=indata %>% left_join(diag_in)
  
  df_train=df_ %>% filter(!diag==test_level)
  df_test=df_ %>% filter(diag==test_level)
  
  KNN_K=5
  knnFit = train(diag ~X+Y+Z, data = df_train, method = "knn",
                 preProcess = c("center", "scale"),
                 tuneGrid = expand.grid(k = c(KNN_K)))
  knn_pred = predict(knnFit, newdata = df_test) %>% as.vector()
  outdata=df_test %>% select(COH_sample,diag) %>% mutate(diag_pred_knn=knn_pred)
  outdata
}

KNN_pred_df_with_test_level=function(indata,outcome,varX_list,id_var,KNN_K,test_level){
  df_train=indata %>% filter(!diag==test_level)
  df_test=indata %>% filter(diag==test_level)
  
  formula_in=formula(paste0(outcome,"~",paste0(varX_list,collapse = "+")))
  
  knnFit = caret::train(formula_in, data = df_train, method = "knn",
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(k = c(KNN_K)))
  knn_pred = predict(knnFit, newdata = df_test) %>% as.vector()
  outdata=df_test %>% select(COH_sample,diag) %>% mutate(diag_pred_knn=knn_pred)
  outdata
}

KNN_pred_tSNE_panel=function(indata,outcome,varX_list,id_var,KNN_K,i_pool,n_pool){
  formula_in=formula(paste0(outcome,"~",paste0(varX_list,collapse = "+")))
  df_2=bind_rows(lapply(n_pool, function(n_one){
    df_1=bind_rows(lapply(i_pool, function(i_one){
      indata1=indata[indata$FeatureN==n_one & indata$perplexity==i_one,]
      knnFit = caret::train(formula_in, data = indata1, method = "knn",
                            preProcess = c("center", "scale"),
                            tuneGrid = expand.grid(k = c(KNN_K)))
      knn_pred = predict(knnFit) %>% as.vector()
      outdata=indata1[c(id_var,outcome)] %>% mutate(diag_pred_knn=knn_pred,perplexity=i_one,FeatureN=n_one)
      outdata
    }))
    df_1
  }))
}

from_tsne2knn=function(for_tsne,diag_in,top_feature_list,featureN_list,perplexity_list,dims){
  df_pred_all=
    bind_rows(lapply(perplexity_list, function(perplexity_one){
      bind_rows(lapply(featureN_list, function(FeatureN_one){
        print(paste0("Do Feature Number:",FeatureN_one,", perplexity Number:",perplexity_one))
        
        var_x=top_feature_list[1:FeatureN_one]
        for_tsne_one=for_tsne[row.names(for_tsne) %in% var_x,]
        dim(for_tsne_one)
        
        df_tsne=run_tsne(for_tsne_one,perplexity_one,dims)
        df_tsne_diag=df_tsne %>% left_join(diag_in)
        
        var_coordinate_tsne=names(df_tsne)[names(df_tsne) %in% c("X","Y","Z")]
        formula_in=formula(paste0(outcome,"~",paste0(var_coordinate_tsne,collapse = "+")))
        KNN_train_fit=caret::train(formula_in, data = df_tsne_diag, method = "knn",
                                   preProcess = c("center", "scale"),
                                   tuneGrid = expand.grid(k = c(kNN_N)))
        
        knn_pred = predict(KNN_train_fit,newdata = df_tsne_diag) %>% as.vector()
        df_tsne_diag=df_tsne_diag[c("COH_sample",var_coordinate_tsne,"diag")]
        outdata=df_tsne_diag %>% mutate(diag_pred_knn=knn_pred)
        outdata$FeatureN=FeatureN_one
        outdata$PerplexityN=perplexity_one
        outdata
      }))
    }))
  df_pred_all
}

knn_pred_file=function(file,diag_in){
  df_tsne=read_tsv(file) %>% mutate(group=paste0(gene_n,"_",perplexity)) %>% left_join(diag_in)
  
  indata=df_tsne
  groups=unique(indata$group)
  
  for(group in groups){
    g_index=match(group,groups)
    print(g_index)
    indata_one=indata[indata$group==group,]
    predict_one=knn_pred(indata_one)
    names(predict_one)[3]=paste0("pred_",group)
    if(g_index==1){predict_out=predict_one}
    if(g_index>1){predict_out=predict_out %>% left_join(predict_one)}
  }
  predict_out$knn_pred=apply(predict_out[3:ncol(predict_out)],1,function(x){names(sort(table(x),decreasing=TRUE)[1])})
  predict_out
}

knn_pred_forDf=function(indata,diag_in){
  indata=indata %>% mutate(group=paste0(gene_n,"_",perplexity)) %>% left_join(diag_in)
  groups=unique(indata$group)
  KNN_K=5
  for(group in groups){
    g_index=match(group,groups)
    print(g_index)
    indata_one=indata[indata$group==group,]
    predict_one=knn_pred(indata_one,diag_in,"_Prediction")
    names(predict_one)[3]=paste0("pred_",group)
    if(g_index==1){predict_out=predict_one}
    if(g_index>1){predict_out=predict_out %>% left_join(predict_one)}
  }
  predict_out
}

cal_classifcation_rate=function(indata,var_y,var_predicted){
  data_freq_classcification=as.data.frame(table(indata[var_y]==unlist(indata[var_predicted]),unlist(indata[var_y])))
  data_freq_classcification1=tidyr::spread(data_freq_classcification,Var1,Freq) 
  names(data_freq_classcification1)[1:3]=c(var_y,"False","True")
  data_freq_classcification2=data_freq_classcification1%>%
    mutate(
      N_subtype=False +True,
      Right=paste0(True," (",sprintf("%.2f",100*True/N_subtype),"%)"),
      Wrong=paste0(False," (",sprintf("%.2f",100*False/N_subtype),"%)"),
      True=NULL,False=NULL
    )
  data_freq_classcification2
}

cal_err_rate=function(indata,y,y_predicted){
  x1=spread(as.data.frame(table(unlist(indata[y])==unlist(indata[y_predicted]))),Var1,Freq)
  if(ncol(x1)==2){names(x1)=c("Wrong","Right")}
  if(ncol(x1)==1){names(x1)=c("Right")}
  x1$Wrong=nrow(indata)-x1$Right
  x1$err_rate=round(x1$Wrong/nrow(indata),4)
  x1$err_freq_percentage=paste0(x1$Wrong," (",sprintf("%.2f",x1$err_rate*100),"%)")
  x1$Right=NULL
  x1
}

#PCA/tSNE+KNN prediction -------

calculate_PC_df=function(df_train,df_test,top_feature_list,used_FeatureN){
  FeatureN_one=used_FeatureN
  print(paste0("Do Feature Number:",FeatureN_one))
  
  var_x=top_feature_list[1:FeatureN_one]
  
  pc_train_fit=prcomp(df_train[var_x])
  pc_train=as.data.frame(predict(pc_train_fit))
  pc_train$COH_sample=df_train$COH_sample
  
  pc_test=as.data.frame(predict(pc_train_fit,newdata = df_test))
  pc_test$COH_sample=df_test$COH_sample
  
  pc_train=pc_train %>% left_join(diag_ref)
  pc_test=pc_test %>% left_join(diag_ref)
  return(list(pc_train=pc_train,pc_test=pc_test))
}


knn_prediction_onPC=function(pc_train,pc_test,PCN_list,kNN_N){
  df_pred1=bind_rows(lapply(PCN_list, function(PCN_one){
    print(paste0("Do PC Number:",PCN_one))
    print("Fit KNN")
    outdata=knn_prediction_train_test_df(pc_train,pc_test,paste0("PC",1:PCN_one),"COH_sample",kNN_N)
    outdata$PCN=PCN_one
    outdata
  }))
  df_pred1
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

knn_prediction_train_test_df=function(df_train,df_test,var_x,id_var,kNN_N){
  formula_in=formula(paste0("diag","~",paste0(var_x,collapse = "+")))
  KNN_train_fit=caret::train(formula_in, data = df_train, method = "knn",
                             preProcess = c("center", "scale"),
                             tuneGrid = expand.grid(k = c(kNN_N)))
  knn_pred = predict(KNN_train_fit,newdata = df_test) %>% as.vector()
  outdata=df_test[id_var] %>% mutate(pred_knn=knn_pred)
  outdata
}

get_pred_summary=function(predict_out){
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
       xlab = "Neighbors", ylab = "Top MAD gene number",
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


get_neighbor_pred=function(df_in,diag_in,top_gene,FeatureN_list,neighbor_list){
  
  df_pred_all_featureN=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
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
    neighbor_subtype_1 = neighbor_subtype[,-1]
    
    # print("Get pred")
    df_pred_all=bind_rows(lapply(neighbor_list, function(top_neighbor_one){
      print(paste0("Get prediction, top_neighbor=",top_neighbor_one))
      df_pred=bind_rows(lapply(1:nrow(neighbor_subtype_1), function(i){
        value_one=as.character(unlist(neighbor_subtype_1[i,1:top_neighbor_one]))
        freq=table(value_one)
        freq=sort(freq,decreasing = T)
        
        data.frame(
          COH_sample=df_diag_match$COH_sample[i],
          FeatureN=FeatureN_one,
          top_neighborN=top_neighbor_one,
          diag=df_diag_match$diag[i],
          pred_GEhalf=ifelse(any(freq > sum(freq)/2),names(freq)[freq > sum(freq)/2],"Fail_prediction"),
          # pred_top=names(freq[1]),
          pred_freq=paste0(paste0(names(freq),": ",freq," (",sprintf("%.1f",100*freq/sum(freq)),"%)"),collapse = "; "),
          stringsAsFactors = F
        )
      }))
      
      df_pred
    }))
    df_pred_all
  }))
  df_pred_all_featureN
}

get_errRate_df=function(df_in,pred_var="diag_pred"){
  FeatureN_list=unique(df_in$FeatureN)
  FeatureN_list=FeatureN_list[!is.na(FeatureN_list)]
  top_neighborN_list=unique(df_in$top_neighborN)
  top_neighborN_list=top_neighborN_list[!is.na(top_neighborN_list)]
  diag_list=unique(df_in$diag)
  diag_list=diag_list[!is.na(diag_list)]
  
  df_err_3=bind_rows(lapply(diag_list, function(diag_one){
    df_one=df_in %>% filter(diag==diag_one)
    df_err_2=bind_rows(lapply(top_neighborN_list, function(top_neighborN_one){
      df_err_1=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
        # print(paste0(diag_one,",",FeatureN_one,",",top_neighborN_one))
        df_one_x=df_one %>% filter(top_neighborN==top_neighborN_one & FeatureN==FeatureN_one)
        rate=table(df_one_x$diag==df_one_x[pred_var])/nrow(df_one_x)
        
        df_err=data.frame(
          diag=diag_one,
          FeatureN=FeatureN_one,
          top_neighborN=top_neighborN_one,
          err_rate=ifelse(length(rate)==2,rate[1],ifelse(names(rate)=="FALSE",1,0)),
          stringsAsFactors = F
        )
        df_err
      }))
      df_err_1
    }))
    df_err_2
  }))
  df_err_3
}

draw_errplot_phenographGroup=function(df_err,file_err_plot){
  
  FeatureN_list=unique(df_err$FeatureN)
  top_neighborN_list=unique(df_err$top_neighborN)
  diag_list=unique(df_err$diag)
  
  df_col=data.frame(top_neighborN=top_neighborN_list,
                    col=unique(subtypeCol)[1:length(top_neighborN_list)],
                    pch=1:length(top_neighborN_list),stringsAsFactors = F)
  
  df_err1=df_err %>% left_join(df_col)
  
  pdf(file_err_plot,width = 20,height = 16)
  par(mfrow=c(5,4))
  
  for(diag_one in diag_list){
    df_err_2=df_err1 %>% filter(diag==diag_one)
    
    err_min=min(df_err_2$err_rate,na.rm = T)
    
    plot(NA,xlim=c(min(df_err_2$FeatureN),max(df_err_2$FeatureN)),ylim=c(0,1),main=paste0(diag_one,": Min error rate=",round(err_min,2)),
         ylab="Error rate",xlab="Feature N")
    
    for(top_neighborN_one in top_neighborN_list){
      lines(df_err_2$FeatureN[df_err_2$top_neighborN==top_neighborN_one],
            df_err_2$err_rate[df_err_2$top_neighborN==top_neighborN_one],
            col=df_err_2$col[df_err_2$top_neighborN==top_neighborN_one])
      points(df_err_2$FeatureN[df_err_2$top_neighborN==top_neighborN_one],
             df_err_2$err_rate[df_err_2$top_neighborN==top_neighborN_one],
             col=df_err_2$col[df_err_2$top_neighborN==top_neighborN_one],
             pch=df_err_2$pch[df_err_2$top_neighborN==top_neighborN_one])
    }
  }
  
  plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt = "n",yaxt = "n",xlab="",ylab="",bty="n")
  legend(c("topleft"),legend = df_col$top_neighborN,col=df_col$col,lty=1,pch=df_col$pch,title = "top_neighborN",cex=2,bty ="n",ncol = 2)
  
  dev.off()
}

get_dfUnMatch_minimumSet=function(df_pred,df_pred_err){
  diag_list=unique(df_pred_err$diag)
  df_unmatch_minimumSet=bind_rows(lapply(diag_list, function(diag_one){
    df_pred_err_one=df_pred_err %>% filter(diag==diag_one)
    df_pred_err_one1=df_pred_err_one %>% filter(err_rate==min(df_pred_err_one$err_rate)) %>% slice_head(n=1)
    
    df_pred_one=df_pred %>% filter(diag==diag_one) %>%
      filter(FeatureN==df_pred_err_one1$FeatureN & top_neighborN==df_pred_err_one1$top_neighborN) %>%
      filter(!diag==diag_pred)
    df_pred_one
  }))
}

#Phenograph prediction -----------------------------------------
Phenograph_pred=function(df_in,diag_in,top_gene,FeatureN_one,neighbor_k=10,ratio_cutoff=0.5){
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


#prediction evaulation -----------------------

draw_errPlot_onFeatureN=function(pred_err_all,filename){
  
  pred_err_all=pred_err_all %>% arrange(FeatureN)
  
  x_min=min(pred_err_all$FeatureN)
  x_max=max(pred_err_all$FeatureN)
  
  y_min=min(pred_err_all$err_rate)
  y_max=max(pred_err_all$err_rate)
  pdf(filename,width = 30,height = 10)
  
  plot(NA,xlim=c(x_min,x_max),ylim=c(y_min,y_max),las=2,ylab="error rate",xlab="Input feature N")
  for(diag_one in key_diag){
    pred_err_one=pred_err_all[pred_err_all$diag==diag_one,]
    lines(pred_err_one$FeatureN,pred_err_one$err_rate,
          col=df_diag$col[df_diag$diag==diag_one])
    points(pred_err_one$FeatureN,pred_err_one$err_rate,
           pch=df_diag$pch[df_diag$diag==diag_one],
           col=df_diag$col[df_diag$diag==diag_one])
  }
  h_ref=seq(0.01,0.4,0.02)
  for(h_ref_one in h_ref){abline(h=h_ref_one,lty=3,col="gray")}
  axis(at=h_ref, side=2,las=2)
  axis(at=seq(0,5000,50), side=1,las=2)
  
  legend("topright",legend = df_diag$diag,col=df_diag$col,pch=df_diag$pch,ncol = 10)
  dev.off()
}

# Boruta -------------------------------------


library(Boruta)

run_boruta=function(matrix_in,diag_in,max_run){
  data_with_diagnosis=as.data.frame(t(matrix_in)) %>%
    mutate(COH_sample=colnames(matrix_in)) %>% 
    left_join(diag_in[c("COH_sample","diag")]) %>%
    mutate(COH_sample=NULL,diag=as.factor(diag)) %>%  na.omit()
  
  print(dim(data_with_diagnosis))
  
  set.seed(10)
  boruta_out = Boruta(diag~., data = data_with_diagnosis, doTrace = 2,maxRuns=max_run)
  boruta_out_fixed = TentativeRoughFix(boruta_out)
  
  importance_boruta=attStats(boruta_out)
  importance_boruta_fixed=attStats(boruta_out_fixed)
  
  data_boruta_confirmed=matrix_in[row.names(importance_boruta)[importance_boruta$decision=="Confirmed"],]
  data_boruta_confirmed_fixed=matrix_in[row.names(importance_boruta_fixed)[importance_boruta_fixed$decision=="Confirmed"],]
  
  return(list(
    fit_boruta=boruta_out,
    importance_boruta=importance_boruta,importance_boruta_fixed=importance_boruta_fixed,
    data_boruta_confirmed=data_boruta_confirmed,data_boruta_confirmed_fixed=data_boruta_confirmed_fixed))
}


extract_boruta_imp=function(in_rdata){
  load(in_rdata)
  
  in_boruta=boruta_out
  diag_one=word(basename(in_rdata),1,sep="[.]")
  
  df_imp=in_boruta$importance_boruta %>% filter(decision=="Confirmed") %>% arrange(desc(medianImp))
  df_imp$gene_id=row.names(df_imp)
  df_imp=df_imp %>% left_join(info_gtf)
  df_imp$diag=diag_one
  df_imp
}

# DE analysis -------------------------------
DE_oneGroup_vs_Others=function(in_count_matrix,in_diag,test_level,batch=c("batch2","batch3")){
  in_diag1=in_diag %>% mutate(diag1=ifelse(diag_new==test_level,diag_new,"0other"))
  in_diag1$diag1=as.factor(in_diag1$diag1)
  in_diag1$batch=as.factor(in_diag1$batch)
  in_diag1$batch2=as.factor(in_diag1$batch2)
  in_diag1$batch3=as.factor(in_diag1$batch3)
  
  deseq_object=create_deseqObject_from_matrixdata(in_count_matrix,in_diag1,"COH_sample",c("diag1",batch))
  dim(deseq_object)
  
  dds = DESeq2::DESeq(deseq_object)
  res=results(dds, contrast=c("diag1",test_level,"0other"))
  
  res_DE=as.data.frame(res@listData) %>% 
    mutate(feature=res@rownames,
           test=paste0(test_level,"_vs_","other")) %>% 
    arrange(padj) 
  names(res_DE)[7]="gene_id"
  res_DE
}


# for exon level ----------------------
to_yesNo=function(one){ifelse(any(unlist(strsplit(one,"[+]")) %in% gene_level_list),"Yes","No")}
to_feature=function(one){ifelse(any(unlist(strsplit(one,"[+]")) %in% gene_level_list),
                                unlist(strsplit(one,"[+]"))[which(unlist(strsplit(one,"[+]")) %in% gene_level_list)[1]],
                                unlist(strsplit(one,"[+]"))[1])}
to_geneName=function(one){paste0(df_gene_name$gene_name[match(unlist(strsplit(one,"[+]")),df_gene_name$feature)],collapse = "+")}


#get_refMatrix_MeanSd for celldex data-----------------------
get_refMatrix_MeanSd=function(ref_in,var){
  df_ref_matrix=as.data.frame(assays(ref_in)$logcounts)
  df_ref_phenotype=as.data.frame(colData(ref_in))
  df_ref_phenotype$sampleid=row.names(df_ref_phenotype)
  
  celltype_list=unique(unlist(df_ref_phenotype[var]))
  
  df_mean=bind_cols(lapply(celltype_list, function(celltype_one){
    sample_in=df_ref_phenotype$sampleid[unlist(df_ref_phenotype[var])==celltype_one]
    
    df_one=df_ref_matrix[sample_in]
    
    df_one$mean=apply(df_one, 1, mean,na.rm=T)
    df_one=df_one['mean']
    names(df_one)=celltype_one
    df_one
  }))
  
  df_genename=data.frame(gene_name=row.names(df_mean),stringsAsFactors = F)
  df_mean$gene_name=row.names(df_mean)
  df_mean=df_genename %>% left_join(df_mean)
  df_sd=bind_cols(lapply(celltype_list, function(celltype_one){
    sample_in=df_ref_phenotype$sampleid[unlist(df_ref_phenotype[var])==celltype_one]
    df_one=df_ref_matrix[sample_in]
    df_one$mean=apply(df_one, 1, sd,na.rm=T)
    df_one=df_one['mean']
    names(df_one)=celltype_one
    df_one
  }))
  
  df_sd$gene_name=row.names(df_sd)
  df_sd=df_genename %>% left_join(df_sd)
  
  list(df_mean=df_mean, df_sd=df_sd)
}


#matrix2tpm --------------------------------------------------
#counts:A numeric data.frame of read counts with samples (columns) and genes (rows).
matrix2tpm = function(counts) {
  len=nrow(counts)
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#get_ReferenceSampe_PhenotypeClass --------------------------------------------------------
#get ReferenceSampe and PhenotypeClass for cibersort using SummarizedExperiment object 
get_ReferenceSampe_PhenotypeClass=function(ref_in,var){
  #get TPM reference sample matrix
  df_ref_matrix=as.data.frame(assays(ref_in)$logcounts)
  df_ref_matrix=df_ref_matrix[!grepl("[/]",row.names(df_ref_matrix)),]
  df_ref_matrix=apply(df_ref_matrix, 2, function(x){2^x})
  df_ref_matrix_tpm=as.data.frame(matrix2tpm(df_ref_matrix))
  sample_n=ncol(df_ref_matrix_tpm)
  df_ref_matrix_tpm$gene_name=row.names(df_ref_matrix_tpm)
  df_ref_matrix_tpm=df_ref_matrix_tpm[c("gene_name",colnames(df_ref_matrix_tpm)[1:sample_n])]
  
  #get phenotype class
  df_id=data.frame(id=colnames(df_ref_matrix_tpm),stringsAsFactors = F) %>% filter(!id=="gene_name")
  df_ref_phenotype=as.data.frame(colData(ref_in))
  df_ref_phenotype$id=row.names(df_ref_phenotype)
  df_ref_phenotype=df_id %>% left_join(df_ref_phenotype) %>% mutate(label=1)
  
  df_ref_phenotype$cellType=gsub(" ",".",df_ref_phenotype$cellType)
  
  df_ref_phenotype1=dcast(df_ref_phenotype,cellType~id,value.var = "label")
  df_phenotypeClass=apply(df_ref_phenotype1,2,function(x){ifelse(is.na(x),2,x)})
  
  #test if id names equal
  names_phenotypeclass=colnames(df_phenotypeClass)[2:ncol(df_phenotypeClass)]
  names_refSample=colnames(df_ref_matrix_tpm)[2:(sample_n+1)]
  print(paste0("All names equal:",all(names_phenotypeclass==names_refSample)))
  #rename
  # colnames(df_ref_matrix_tpm)=c("gene_name",unlist(df_ref_phenotype[var]))
  #output
  list(ReferenceSampe=df_ref_matrix_tpm,PhenotypeClass=df_phenotypeClass)
}

#get_strandness --------------------------------------------------------
get_strandness=function(file_list){
  
  df_strand=bind_rows(lapply(file_list, function(file_one){
    id=word(file_one,2,sep="/")
    value_one=system(paste0("cat ",file_one,"|grep ++|cut -d: -f2"),intern = T)
    df_one=as.data.frame(matrix(c(id,value_one),nrow=1))
    df_one
  }))
  
  df_strand$V2=as.numeric(as.character(df_strand$V2))
  df_strand$V3=as.numeric(as.character(df_strand$V3))
  df_strand$percentage_max=apply(df_strand[c("V2","V3")],1,max)
  df_strand=df_strand %>% mutate(batch1=ifelse(percentage_max>0.7,"Stranded","Unstranded"))
  df_strand
}


#get fq length --------------------------------------------------------
get_fq_length=function(fqfile_pattern="fq_rna/*.fq.gz",delim="[.]",fqid_index=1){
  fq_length=system(paste0("ls ",fqfile_pattern,"|xargs -i echo 'paste <(echo {}) <(zcat {}|head -n2|tail -n1)'|bash"),intern = T)
  df_fq_length=data.frame(fq=word(fq_length,1,sep="\t"),
                          reads=word(fq_length,2,sep="\t"),
                          length=as.numeric(nchar(word(fq_length,2,sep="\t"))),
                          stringsAsFactors = F)
  
  df_fq_length=df_fq_length %>% mutate(
    name_in_source=basename(fq),
    fqid=word(name_in_source,fqid_index,sep=delim)
  )
  df_fq_length
}


































