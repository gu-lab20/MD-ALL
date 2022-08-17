# set -------------------------
setwd("/scratch/zuhu/project/HTB/ALL/out/COH002918_R1/")

library(stringr)
library(dplyr)
library(DESeq2)
library(ggrepel)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
gc(rm(list=ls()))


#load data and functions --------------------
source("/home/zgu_labs/bin/R_script/F_RNAseq.R")

#Arguments -------------
COH_sample="COH002914_R1"
vstBatch="/scratch/zuhu/project/HTB/ALL/out/stats/GEP_prediction/vstBatch.rdata"
outdir="GEP_prediction_PCA"

args <- commandArgs(trailingOnly = TRUE)
print(args)

COH_sample=args[1]
vstBatch=args[2]
outdir=args[3]

#load reference information --------------------
load("/home/zgu_labs/bin/R_script/BALL_GEPprediction/for_GEP_prediction_PCA_1502.rdata")
load(vstBatch)


# vst_batch=correct_batch(vst_gene_1502,diag_1502,"COH_sample")
# 
# vst_batch_codingNoXYM=vst_batch[row.names(vst_batch) %in% gene_in,]
# 
# vst_batch_codingNoXYM_highE=remove_lowexpression_genes(vst_batch_codingNoXYM,5)
# 
# vst_batch_codingNoXYM_highE_MAD10000=get_high_MAD_gene(vst_batch_codingNoXYM_highE,10000)
# 
# vst_batch_codingNoXYM_highE_MAD10000_noCorr=remove_high_corr_genes(vst_batch_codingNoXYM_highE_MAD10000,0.75)
# dim(vst_batch_codingNoXYM_highE_MAD10000_noCorr)
# vst_batch_codingNoXYM_highE_MAD10000_noCorr_MAD5000=get_high_MAD_gene(vst_batch_codingNoXYM_highE_MAD10000_noCorr,5000)
# 
# top_feature_list_1502=get_topMAD_features(vst_batch_codingNoXYM_highE_MAD10000_noCorr_MAD5000,5000)
# 
# save(diag_1502,vst_gene_1502,top_feature_list_1502,
#      file="/home/zgu_labs/bin/R_script/BALL_GEPprediction/for_GEP_prediction_PCA_1502.rdata")

diag_ref=diag_1502
top_feature_list=top_feature_list_1502


df_vstBatch_transpose=as.data.frame(t(vstBatch_forPrediction[top_feature_list,]))
df_vstBatch_transpose$COH_sample=row.names(df_vstBatch_transpose)

df_feature_train=df_vstBatch_transpose[df_vstBatch_transpose$COH_sample %in% diag_ref$COH_sample,]
df_feature_test=df_vstBatch_transpose[df_vstBatch_transpose$COH_sample %in% COH_sample,]

FeatureN_list=c(450,500,600,650,700,750)
PCN_list=c(20,21,22,23,24,25,26)

FeatureN_list=c(25,seq(40,80,20),seq(1020,2000,20))
PCN_list=c(25)


pred_all=bind_rows(lapply(FeatureN_list, function(FeatureN_one){
  list_pc=calculate_PC_df(df_feature_train,df_feature_test,top_feature_list,FeatureN_one)
  pred_=knn_prediction_onPC(list_pc$pc_train,list_pc$pc_test,PCN_list,5)
  pred_$FeatureN=FeatureN_one
  pred_
}))
cat 
df_knn_out1=get_pred_summary(pred_all$pred_knn)
df_knn_out2=get_pred_summary(pred_all$pred_knn[pred_all$PCN %in% c(23,24)])


if(!dir.exists(outdir)){dir.create(outdir)}
write_tsv(pred_all,file.path(outdir,paste0(COH_sample,".PCA.KNN.pred_detail.tsv")))

print(paste0("Output file: ",file.path(outdir,paste0(COH_sample,".PCA.KNN.pred.tsv"))))

write_tsv(df_knn_out1,file.path(outdir,paste0(COH_sample,".PCA.KNN.pred.tsv")))
#write_tsv(df_knn_out2,file.path(outdir,paste0(COH_sample,".PCA.KNN.pred_PCN23_24.tsv")))


# draw tSNE 2D ---------------------
# for_tsne_plot=run_tsne(for_pca,30)
# write_tsv(for_tsne_plot,file.path(outdir,paste0(COH_sample,".df_tsne.tsv")))
# 
# draw_tsne_MAD(for_tsne_plot,diag_in_,diagCol,1000,30,5)
# ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D.pdf")),width = 10,height=7,dpi = 300)
# ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D.png")),width = 10,height=7,dpi = 300)








