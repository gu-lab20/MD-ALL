# set -------------------------
# setwd("/scratch/zuhu/project/HTB/ALL/out_new/COH002914_R1/TRANSCRIPTOME/GEP_prediction/")
setwd("//isi-dcnl/user_data/zgu_grp/DryLab/bin/jiliu/R/")

library(dplyr)
library(Rtsne)
library(caret)
library(ggrepel)

# library(stringr)
# library(DESeq2)
gc(rm(list=ls()))

#load data and functions --------------------
source("F_RNAseq_JL.R")

#Arguments -------------
# COH_sample=file_path_sans_ext(basename("test_out/COH001000_D1.vst.tsv"))
COH_sample=tools::file_path_sans_ext("test_out//COH001000_D1.vst.tsv")

COH_sample="COH001000_D1"
vstBatch="test_out/vstBatch.rdata"
outdir="test_out/GEP_Prediction"

# args <- commandArgs(trailingOnly = TRUE)
# print(args)

# COH_sample=args[1]
# vstBatch=args[2]
# outdir=args[3]

#load reference information --------------------
# load("/home/zgu_labs/bin/R/RNAseq/GEP_prediction/for_GEPprediction_tSNE2KNN_1319.rdata")
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/for_GEPprediction_tSNE2KNN_1499.rdata")

load(vstBatch)

n_pool=seq(200,900,50)
i_pool=c(20,30)

# n_pool=c(60,seq(100,5000,50))
# i_pool=c(20,30,40,50)

# n_pool=c(4950,5000)
# i_pool=c(20)

# n_pool=c(200,250)
# i_pool=c(20,30)

diag_ref=diag_1449
feature_in=top_feature_list_1499

# generate df for tsne, use established gene list in reference cohort -------------------
df_vstBatch=vstBatch_forPrediction[unique(c(diag_ref$COH_sample,COH_sample))]

for_tSNE=df_vstBatch[feature_in,]

# run tSNE 3D ----------------------
#make diag
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

if(!dir.exists(outdir)){dir.create(outdir)}

for_knn=run_tsne_panel_basedOnTopGeneList_3D(for_tSNE,feature_in,n_pool,i_pool) 
write_tsv(for_knn,file.path(outdir,paste0(COH_sample,".df_tSNE_forKNN.tsv")))

# KNN prediction K=5 -------------------
df_forKNN_diag=for_knn %>% left_join(diag_in_)

pred_all=knn_prediction_ontSNE(df_forKNN_diag)

df_knn_out=get_pred_summary(pred_all$pred_knn)

print(paste0("outfile: ",file.path(outdir,paste0(COH_sample,".tSNE.KNN.pred_detail.tsv"))))
print(paste0("outfile: ",file.path(outdir,paste0(COH_sample,".tSNE.KNN_prediction.tsv"))))

write_tsv(pred_all,file.path(outdir,paste0(COH_sample,".tSNE.KNN.pred_detail.tsv")))
write_tsv(df_knn_out,file.path(outdir,paste0(COH_sample,".tSNE.KNN_prediction.tsv")))


# draw tSNE 2D ---------------------
for_tsne_plot=run_tsne(for_tSNE[head(feature_in,800),],30,2)
for_tsne_plot1=for_tsne_plot %>% left_join(diag_in_)
write_tsv(for_tsne_plot1,file.path(outdir,paste0(COH_sample,".df_tsne_forPlot.tsv")))

draw_tsne_MAD(for_tsne_plot,diag_in_,diagCol,800,30,10)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D.pdf")),width = 10,height=8,dpi = 300)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D.png")),width = 10,height=8,dpi = 300)







