#set ---------------------
setwd("/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/R_code/17.phenograph")

library(dplyr)
library(stringr)
library(reshape2)
library(Rphenograph)

gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNAseq.R")

#load data
load("../1.prepare_diag/2.data_temp/diag_1552.rdata")
load("../5.tSNE_gene/2.data_temp/1552/for_tsne_1552.rdata")


Rphenograph::find_neighbors()

iris_unique <- unique(iris) # Remove duplicates
data <- as.matrix(iris_unique[,1:4])
Rphenograph_out <- Rphenograph(data, k = 45)
modularity(Rphenograph_out[[2]])
membership(Rphenograph_out[[2]])
iris_unique$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
ggplot(iris_unique, aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)) + geom_point(size = 3)+theme_bw()


df_for_tsne_1552=for_tsne_1552$for_tsne_
top_gene=get_topMAD_features(df_for_tsne_1552,800)


df_for_tsne_1552_1=df_for_tsne_1552[top_gene,]

dim(df_for_tsne_1552)

df_in=as.matrix(t(df_for_tsne_1552_1))
dim(df_in)

Rphenograph_out=Rphenograph(df_in, k = 5)

df_Rphenograph_out=data.frame(
  COH_sample=row.names(df_in),
  stringsAsFactors = F
) %>% mutate(
  Rphenograph_group=membership(Rphenograph_out[[2]])
) %>% left_join(diag_1552)

df_freq=data.frame(table(df_Rphenograph_out$Rphenograph_group,df_Rphenograph_out$diag)) %>% mutate(varname=paste0("RphenographGroup",Var1))

df_freq1=dcast(df_freq,varname~Var2,value.var = "Freq")

write_tsv(df_freq1,"df_freq_Rphenograph.tsv")


table(df_Rphenograph_out$Rphenograph_group)

ggplot(iris_unique, aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)) + geom_point(size = 3)+theme_bw()

df_tsne=for_tsne_1552$df_tsne

df_Rphenograph_out_withtsne=df_Rphenograph_out %>% left_join(df_tsne)

diag_in=df_Rphenograph_out_withtsne %>% mutate() %>% select(COH_sample)




