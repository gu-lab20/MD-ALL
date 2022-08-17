# set -------------------------
setwd("/scratch/zuhu/project/GordanaRaca/ALL/out/COH002956_D1/TRANSCRIPTOME/")
setwd("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction")

library(stringr)
library(dplyr)
library(DESeq2)
gc(rm(list=ls()))

#load data and functions --------------------
source("/home/zgu_labs/bin/R/RNAseq/F_RNAseq.R")

#Arguments -------------
file="/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/out_new/COH000066_D1/TRANSCRIPTOME/HTSeq_DUX4/COH000066_D1.DUX4patched.HTSeq"
outdir="/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/out_new/COH000066_D1/TRANSCRIPTOME/GEP_prediction"

args <- commandArgs(trailingOnly = TRUE)
print(args)

file=args[1]
outdir=args[2]
#load reference data --------------------
load("/home/zgu_labs/bin/R/RNAseq/GEP_prediction/for_vst_N234.rdata")

id=word(basename(file),1,sep="[.]")

#prepare count for reference + new sample-------------------
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

#running vst------------------
deseqobject=create_deseqObject_from_matrixdata(count_panel_with_one,diag_in,"COH_sample","diag")
df_vst=run_vst(deseqobject)

#output data -----------
df_vst_one=data.frame(
  feature=row.names(df_vst),
  vst=unlist(df_vst[id]),
  stringsAsFactors = F
)

names(df_vst_one)[2]=id
  
#write file ----------------
write_tsv(df_vst_one,file.path(outdir,paste0(id,".vst.tsv")))
