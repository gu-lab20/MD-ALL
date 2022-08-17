# set -------------------------
# setwd("/scratch/zuhu/project/HTB/ALL/out_new/COH002914_R1/TRANSCRIPTOME/GEP_prediction/")
# setwd("/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/out_new/COH000066_D1/TRANSCRIPTOME/GEP_prediction/")
# setwd("/scratch/zuhu/project/GordanaRaca/ALL/out_new/COH002928_D1/TRANSCRIPTOME/GEP_prediction/")

library(stringr)
library(ggrepel)
library(dplyr)
library(Rphenograph)
library(reshape2)
library(Rtsne)

gc(rm(list=ls()))

#load data and functions --------------------
source("F_RNAseq_JL.R")

#Arguments -------------
COH_sample="COH001000_D1"
file_vst="test_out/COH001000_D1.vst.tsv"
# outdir="test_out/GEP_Prediction"
outdir="test_out/GEP_Prediction/test"

# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# 
# COH_sample=args[1]
# file_vst=args[2]
# outdir=args[3]


#load reference information --------------------
load("C:/Users/jiliu/OneDrive - City of Hope National Medical Center/Documents/Rotation/Gu/reference/for_GEPprediction_tSNE2KNN_1499.rdata")

diag_ref=diag_1449
top_gene=top_feature_list_1499
vstBatch_ref=vst_gene_1499

print(paste0("Run analysis for ",COH_sample))
if(!dir.exists(outdir)){dir.create(outdir)}

print("Get diag df")
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

print("Get vst values")
df_vst=read_tsv(file_vst)

vstBatch_ref$feature=row.names(vstBatch_ref)
df_for_pred=vstBatch_ref %>% left_join(df_vst)
row.names(df_for_pred)=df_for_pred$feature
df_for_pred$feature=NULL

#get neighbor heatmap ------------------------------------------------------------
print("Run phenograph")
# FeatureN_list=c(50,5000)

# FeatureN_list=seq(50,5000,50)
# top_neighborN_list=c(5,10,20)


FeatureN_list=seq(50,2000,900)
top_neighborN_list=10

df_Phenograph_pred=Phenograph_pred_list(df_for_pred,diag_in_,top_gene,FeatureN_list,top_neighborN_list,ratio_cutoff=0.5)

df_Phenograph_pred1=df_Phenograph_pred[df_Phenograph_pred$COH_sample==COH_sample,]

write_tsv(df_Phenograph_pred1,file.path(outdir,paste0(COH_sample,".PhenographCluster_pred_detail.tsv")))

df_pred_out=bind_rows(
  get_Phenograph_pred_sum(df_Phenograph_pred1,5) %>% mutate(top_neighborN=5),
  get_Phenograph_pred_sum(df_Phenograph_pred1,10) %>% mutate(top_neighborN=10),
  get_Phenograph_pred_sum(df_Phenograph_pred1,20) %>% mutate(top_neighborN=20)
)

write_tsv(df_pred_out,file.path(outdir,paste0(COH_sample,".PhenographCluster_pred.tsv")))
#get neighbor heatmap ------------------------------------------------------------
print("Get neighbor heatmap")

df_top_neighbor=get_top_neighbor(COH_sample,df_for_pred,diag_in_,top_gene,FeatureN_list)
write_tsv(df_top_neighbor,paste0(outdir,"/",COH_sample,".top_neighbor.tsv"))

pdf(paste0(outdir,"/",COH_sample,".neighbors_50-5000.pdf"),width = 10,height=12)
draw_neighbor_heatmap(df_top_neighbor,inset_x = -0.25)
dev.off()

png(paste0(outdir,"/",COH_sample,".neighbors_50-5000.png"),width = 3500,height=4800,res = 300)
draw_neighbor_heatmap(df_top_neighbor,inset_x = -0.2)
dev.off()

#run tsne ------------------------------------------------------------
print("Run tSNE")

df_tsne_perplexity10=run_tsne(df_for_pred[head(top_gene,800),],10,2)
df_tsne_perplexity20=run_tsne(df_for_pred[head(top_gene,800),],20,2)

write_tsv(bind_rows(df_tsne_perplexity10,df_tsne_perplexity20),
          file.path(outdir,paste0(COH_sample,".tsne_2D.tsv")))

draw_tsne_MAD(df_tsne_perplexity10,diag_in_,diagCol,800,10,10)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D_perplexity10.pdf")),width = 10,height=8,dpi = 300)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D_perplexity10.png")),width = 10,height=8,dpi = 300)

draw_tsne_MAD(df_tsne_perplexity20,diag_in_,diagCol,800,20,10)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D_perplexity10.pdf")),width = 10,height=8,dpi = 300)
ggsave(file.path(outdir,paste0(COH_sample,".tSNE_2D_perplexity20.png")),width = 10,height=8,dpi = 300)











