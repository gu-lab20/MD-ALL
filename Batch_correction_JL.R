# set -------------------------
# setwd("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction")

setwd("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction")

library(stringr)
library(dplyr)
# library(DESeq2)
# library(ggrepel)
gc(rm(list=ls()))

#load data and functions --------------------
source("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/F_RNAseq_JL.R")
source("functions.R")
# load("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/for_vst_N234.rdata")
# diag_1449=read_tsv("/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/R_code/19.select_GoldRef/2.data_temp/selection_6/diag_1499.tsv")
# top_feature_list_1499=readLines("/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/R_code/19.select_GoldRef/2.data_temp/selection_6/top_feature.list")
# vst_gene_1499=vst_gene_2239[diag_1449$COH_sample]
# save(diag_1449,vst_gene_1499,top_feature_list_1499,
#      file="/home/zgu_labs/bin/R/RNAseq/GEP_prediction/for_GEPprediction_tSNE2KNN_1499.rdata")
load("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/for_GEPprediction_tSNE2KNN_1499.rdata")
#Arguments -------------

vst_list=list.files("test_out/",full.names = T)

# vst_list=vst_list[grepl("HTSeq",vst_list)]

vst_list="//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/vst_list"
id_info="//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/test/df_batch.tsv"
vstBatch="//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction/test_out/vstBatch.rdata"

# args <- commandArgs(trailingOnly = TRUE)
# cat("Arguments:\n")
# print(args)
# cat("\n")

# vst_list=args[1]
# id_info=args[2]
# vstBatch=args[3]

vst_gene_in=vst_gene_1499
diag_ref=diag_1449

#get vst matrix -------------
#files_vst=readLines(vst_list)
files_vst=list.files("test_out/",full.names = T)
cat(paste0("vst file N= ",length(files_vst),"\n\n"))
df_vst=get_vst_matrix(files_vst)
cat("\n")

cat(paste0("Done reading vst files\n"))
cat(paste0("vst matrix SampleN= ",ncol(df_vst)-1,"\n\n"))

vst_gene_in$feature=row.names(vst_gene_in)
df_vst=df_vst[unique(c("feature",names(df_vst)[!names(df_vst) %in% names(vst_gene_in)]))]
vst_all=vst_gene_in %>% left_join(df_vst)
row.names(vst_all)=vst_all$feature
#get diag --------------------
id_info=read_tsv(id_info)

diag_test=id_info %>% transmute(COH_sample=COH_sample,batch=batch)
cat(paste0("Batchs in ref cohort: ",paste0(unique(diag_ref$batch),collapse = ", "),"\n"))
cat(paste0("Batchs in df for prediction: ",paste0(unique(diag_test$batch),collapse = ", "),"\n\n"))

diag_for_batch=rbind(
  diag_ref %>% select(COH_sample,batch),
  diag_test
) %>% distinct() %>% filter(COH_sample %in% names(vst_all))

#select vst
vst_all=vst_all[,diag_for_batch$COH_sample]

#run batch correction
cat(paste0("Start batch correcting\n\n"))
vstBatch_forPrediction=correct_batch(vst_all,diag_for_batch,"COH_sample")

#save batch file
if(!dir.exists(dirname(vstBatch))){dir.create(dirname(vstBatch))}
save(vstBatch_forPrediction,file=vstBatch)




