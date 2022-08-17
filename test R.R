library(tidyr)
library(tibble)
library(data.table)

df_vst_file=x
COH_sample=names(df_vst_file)[2]
cat(paste0("Run analysis for ",COH_sample,"\n\n"))
diag_in_=diag_1409
vst_ref=for_tsne_1409_

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


gene_ref=rename(gene_ref_,gene_id="variable")
df_MedianMeanValue_comb=left_join(gene_ref,df_MedianMeanValue)
for_box_plot=dplyr::select(df_MedianMeanValue_comb,"diag","gene_name","median","mean")
print(for_box_plot)
save(for_box_plot,file="df_forboxplot.RData")



# ------------------
df_vst=x
diag_ref=df_MedianMeanValue
df_gene_name=gene_ref_

df_vst_rename1=rename(df_vst,feature="variable")
df_vst_rename2=rename(df_vst_rename1,sample_test="median")
df_vst_diag=add_column(df_vst_rename2,diag="_Prediction")
df_vst_all=bind_rows(df_vst_diag,df_MedianMeanValue)

gene_ref=rename(gene_ref_,gene_id="variable")
df_vst_all_comb=left_join(gene_ref,df_vst_all)
for_box_plot=dplyr::select(df_vst_all_comb,"diag","gene_name","median","mean")

col_in=unique(subtypeCol[names(subtypeCol) %in% df_vst_all$diag])
variable_type_in=c("TCF3","PAX5")
df_for_boxplot=filter(for_box_plot,gene_name %in% variable_type_in)
ggplot(df_for_boxplot,aes(x=gene_name,y=median)) + 
  geom_boxplot(show.legend = F,coef = 6) + 
  geom_jitter(aes(color=diag)) + 
  scale_color_manual(values=col_in) +
  geom_text_repel(data=subset(df_for_boxplot, diag == "_Prediction" ), aes(label=COH_sample), fontface = "bold")+
  theme_bw()
