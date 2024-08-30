#' run_vst
#'
#' @param obj_in , object, object of reference
#' @param assay_name_in , string, assay name used for running vst
#' @param assay_name_out , string, assay names used for storing the vst values
#' @param var_design , string, variable in coldata used in vst running formula design. will not affecting the vst values
#'
#' @return
#' @export run_vst
#'
#' @examples
#' obj_x=run_vst(obj_in = obj_x)
run_vst=function(obj_in,assay_name_in="counts",assay_name_out="vst",var_design=NULL){
  matrix_=as.matrix(assays(obj_in$SE)[[assay_name_in]])
  df_colData=as.data.frame(as.matrix(colData(obj_in$SE)))

  if(!is.null(var_design)){formula_in=formula(paste0("~",paste0(var_design,collapse = "+")))}
  if(is.null(var_design)){formula_in=formula(paste0("~","1"))}

  obj_deseq2 = DESeq2::DESeqDataSetFromMatrix(countData = matrix_,colData = df_colData,design = formula_in)
  cat("Running vst\n")
  obj_deseq2=DESeq2::varianceStabilizingTransformation(obj_deseq2,blind = T)
  df_vst=round(as.matrix(assay(obj_deseq2)),6)

  assays(obj_in$SE)[[assay_name_out]]=df_vst
  cat("Done vst\n")
  obj_in
}


#' get_vst_values
#'
#' @param obj_in , object, object of reference
#' @param file_obj , string, file of object of reference
#' @param df_count , dataframe, count dataframe for merge
#'
#' @return
#' @export get_vst_values
#'
#' @examples
get_vst_values=function(obj_in=NULL,file_obj=NULL,df_count,out_fmt="df"){

  obj_x=obj_merge(obj_in=obj_in,file_obj = file_obj,df_in = df_count,assay_name_in = "counts")

  obj_x=run_vst(obj_in = obj_x)

  df_vst=as.data.frame(assays(obj_x$SE)[["vst"]])
  df_vst$feature=row.names(df_vst)

  df_vst_out=df_vst[c("feature",names(df_count)[2:ncol(df_count)])] %>% arrange(feature)

  if(out_fmt=="df"){out=df_vst_out}

  if(out_fmt=="obj"){out=obj_x}
  out
}



#' run_umap
#'
#' @param obj_in
#' @param n_neighbors
#' @param dims
#' @param variable_n
#' @param out_label
#' @param feature_panel
#'
#' @return
#' @export
#'
#' @examples
run_umap=function(obj_in,n_neighbors=30,dims=2,variable_n=800,out_label="umap",feature_panel="variable_genes",min_dist=0.3){
  #get features
  variable_genes_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]

  #get running data
  indata=assays(obj_in$SE[variable_genes_in,])[["vst"]]
  indata1=t(indata)
  cat("Running uMAP: Feature N=", dim(indata1)[2], "; Sample N=", dim(indata1)[1],"; n_neighbors=", n_neighbors,  "\n")

  #running umap
  set.seed(10)
  fit_umap=umap(indata1,n_neighbors=n_neighbors,random_state=123,n_components=dims,min_dist=min_dist)

  #get output
  df_umap=as.data.frame(fit_umap$layout)
  names(df_umap)=c("uMAP_1","uMAP_2","uMAP_3")[1:ncol(df_umap)]
  df_umap$n_neighbors=n_neighbors
  df_umap$FeatureN=dim(indata1)[2]
  df_umap$COH_sample=row.names(df_umap)

  #summary output
  obj_in[[out_label]]=df_umap
  obj_in[[paste0(out_label,"_fit")]]=fit_umap

  obj_in
}

#' run_PhenoGraph
#'
#' @param obj_in
#' @param feature_panel
#' @param variable_n
#' @param neighbor_k
#' @param ratio_cutoff
#' @param out_label
#'
#' @return
#' @export
#'
#' @examples
run_PhenoGraph=function(obj_in,
                        feature_panel="variable_genes",
                        variable_n = variable_n,
                        neighbor_k=neighbor_k,
                        ratio_cutoff=0.5,
                        out_label="PhenographPred"){

  variable_genes_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]
  indata=assays(obj_in$SE[variable_genes_in,])[["vst"]]

  indata1=t(indata)
  cat("Run Phenograph: Used Feature N=", dim(indata1)[2], "; Used Sample N=", dim(indata1)[1],"; Neighbor_k=", neighbor_k,  "\n")

  df_diag_match=data.frame(COH_sample=row.names(indata1),diag=obj_in$SE$diag,stringsAsFactors = F) %>% mutate(obs=1:n())

  set.seed(10)
  PG_out=Rphenograph(indata1, neighbor_k)

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
    mutate(diag_pred=ifelse(!is.na(cluster),diag,"NoClusterAsigned"),
           diagCluster_index=1:n(),
           diagCluster_index_max=n(),
           diag_pred_granular=ifelse(diagCluster_index_max==1,diag_pred,paste0(diag_pred,"_",diagCluster_index)),
           diag_pred_freq=paste0(diag_pred,"(",N_diagCluster,";",ratio*100,"%)"),
           diag_pred_granular_freq=paste0(diag_pred_granular,"(",N_diagCluster,";",ratio*100,"%)"),
    ) %>% select(diag,cluster,diag_pred,diag_pred_freq,diag_pred_granular,diag_pred_granular_freq)

  df_cluster_diag$diag=NULL

  df_diag_pred=df_cluster_detail %>%
    left_join(df_cluster_diag %>% mutate(ratio=NULL)) %>%
    mutate(obs=NULL,FeatureN=variable_n,top_neighborN=neighbor_k)

  df_diag_pred$diag_pred=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred)
  df_diag_pred$diag_pred_freq=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_freq)
  df_diag_pred$diag_pred_granular=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_granular)
  df_diag_pred$diag_pred_granular_freq=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_granular_freq)

  row.names(df_diag_pred)=df_diag_pred$COH_sample

  obj_in[[out_label]]=df_diag_pred
  obj_in
}

#' get_PhenoGraphPred
#'
#' @param obj_in
#' @param panelName
#' @param SampleLevel
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_PhenoGraphPred=function(obj_in,panelName="PhenographPred",SampleLevel="TestSample",type="value"){
  df_phenograph=obj_in[[panelName]]
  SampleLevel=gsub("[-]",".",SampleLevel)
  df_out=df_phenograph[df_phenograph$COH_sample %in% c(SampleLevel),]

  value_out=paste0("PhenoGraph Clustering Labeled Subtype: ",df_out$diag_pred,
                   " (FeatureN=",df_out$FeatureN,"; NeighborN=",df_out$top_neighborN,")")
  # if(type=="value"){out=value_out}
  # if(tolower(type) %in% c("df","dataframe")){out=df_out}
  return(list(df=df_out,value=value_out))
}

#' get_PhenoGraphPreds
#'
#' @param obj_in
#' @param feature_panel
#' @param SampleLevel
#' @param variable_n_list
#' @param neighbor_k
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
get_PhenoGraphPreds=function(obj_in,feature_panel="keyFeatures",SampleLevel ="TestSample",
                             variable_n_list=c(seq(100,1000,100),1058),neighbor_k = 10,cores=1){

  if(floor(cores)==1){
    df_out=bind_rows(lapply(variable_n_list,function(featureN){
      obj_i=run_PhenoGraph(obj_in = obj_in,feature_panel = feature_panel,variable_n = featureN,neighbor_k = neighbor_k)
      data.frame(
        featureN=featureN,
        pred=get_PhenoGraphPred(obj_i,SampleLevel =SampleLevel)$df$diag_pred,
        method="PhenoGraph",
        stringsAsFactors = F
      )
    }))
  }

  if(floor(cores)>1){
    cores=floor(cores)
    require(parallel)
    df_out=bind_rows(mclapply(variable_n_list,function(featureN){
      obj_i=run_PhenoGraph(obj_in = obj_in,feature_panel = feature_panel,variable_n = featureN,neighbor_k = neighbor_k)
      data.frame(
        featureN=featureN,
        pred=get_PhenoGraphPred(obj_i,SampleLevel =SampleLevel)$df$diag_pred,
        method="PhenoGraph",
        stringsAsFactors = F
      )
    },mc.cores = cores))
  }
  df_out
}

#' get_SVMPreds
#'
#' @param models_svm
#' @param df_in
#'
#' @return
#' @export
#'
#' @examples
get_SVMPreds=function(models_svm,df_in,id="TestSample"){
  bind_rows(lapply(1:length(models_svm),function(i){
    model_svm=models_svm[[i]]
    data.frame(
      featureN=as.numeric(gsub("FeatureN","",names(models_svm)[[i]])),
      pred=as.character(predict(model_svm,newdata = t(df_in[id]))),
      method="SVM",
      stringsAsFactors = F
    )
  }))
}

#' get_pred_label
#'
#' @param df_pred
#'
#' @return
#' @export
#'
#' @examples
get_pred_label=function(df_pred){
  (df_pred %>% mutate(n=n()) %>%  group_by(pred) %>%
     mutate(n_pred=n(),percentage=round(n_pred/n,2),label=paste0(pred,",",n_pred,"(",percentage,")")) %>%
     ungroup() %>% arrange(desc(percentage)) %>%
     select(label) %>% distinct() %>% mutate(label=paste0(label,collapse = "|")) %>%
     distinct())$label
}

#' get_pred_result
#'
#' @param df_pred
#'
#' @return
#' @export
#'
#' @examples
get_pred_result=function(df_pred,count_PAX5ETV6_as_PAX5alt=F){
  if(!count_PAX5ETV6_as_PAX5alt){
    out=(df_pred %>% mutate(n=n()) %>%  group_by(pred) %>%
           mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
           select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$out
  }

  if(count_PAX5ETV6_as_PAX5alt){
    out=(df_pred %>% mutate(pred=ifelse(pred %in% c("PAX5::ETV6","PAX5 P80R"),"PAX5alt",pred),n=n()) %>%  group_by(pred) %>%
           mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
           select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$out
  }
  out
}

#' get_GEP_pred
#'
#' @param count_matrix
#' @param featureN_PG
#'
#' @return
#' @export
#'
#' @examples
get_GEP_pred=function(count_matrix,
                      featureN_PG=c(100,200,300,400,500,600,700,800,900,1000,1058)){
  df_vsts=get_vst_values(obj_in = obj_234_HTSeq,df_count = count_matrix)

  df_listing=data.frame(id=colnames(df_vsts)[2:ncol(df_vsts)],stringsAsFactors = F) %>%
    mutate(obs=1:n())

  df_out_GEP=bind_rows(lapply(df_listing$obs,function(i){
    sample_id=df_listing$id[i]
    cat(paste0("Prediction for ",sample_id," ...\n"))

    #imputation
    df_vst_i=f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vsts[c("feature",sample_id)])

    #get obj
    obj_=obj_merge(obj_in = obj_1821,df_in = df_vst_i,assay_name_in = "vst")

    #phenograph
    cat("Running phenograph...")
    df_out_phenograph=get_PhenoGraphPreds(obj_in = obj_,feature_panel = "keyFeatures",SampleLevel = sample_id,
                                          neighbor_k = 10,
                                          # variable_n_list = c(seq(100,1000,100),1058)
                                          variable_n_list = featureN_PG
    )
    cat("Running SVM...\n")

    df_out_svm=get_SVMPreds(models_svm,df_in = df_vst_i,id = sample_id)

    df_feateure_exp=get_geneExpression(df_vst = df_vsts[c("feature",sample_id)],genes = c("CDX2","CRLF2","NUTM1"))

    data.frame(
      sample_id=sample_id,
      CDX2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"],
      CRLF2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CRLF2"],
      NUTM1=df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"],
      phenoGraph_pred=get_pred_result(df_out_phenograph),
      phenoGraph_predScore=get_pred_score(df_out_phenograph),
      phenoGraph_predLabel=get_pred_label(df_out_phenograph),

      svm_pred=get_pred_result(df_out_svm),
      svm_predScore=get_pred_score(df_out_svm),
      svm_predLabel=get_pred_label(df_out_svm),

      stringsAsFactors = F
    )
  }))
}

#' get_pred_score
#'
#' @param df_pred
#'
#' @return
#' @export
#'
#' @examples
get_pred_score=function(df_pred,count_PAX5ETV6_as_PAX5alt=F){
  if(!count_PAX5ETV6_as_PAX5alt){
    out=(df_pred %>% mutate(n=n()) %>%  group_by(pred) %>%
           mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
           select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$percentage
  }

  if(count_PAX5ETV6_as_PAX5alt){
    out=(df_pred %>% mutate(pred=ifelse(pred %in% c("PAX5::ETV6","PAX5 P80R"),"PAX5alt",pred),n=n()) %>%  group_by(pred) %>%
           mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
           select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$percentage
  }
  out
}


#' get_subtype_final
#'
#' @param id
#' @param df_feateure_exp
#' @param df_out_phenograph
#' @param df_out_svm
#' @param out_mutation
#' @param chrom_n
#' @param CNV_label
#' @param fusion_fc
#' @param fusion_c
#'
#' @return
#' @export
#'
#' @examples
get_subtype_final=function(
  id,
  df_feateure_exp=df_feateure_exp,
  df_out_phenograph=df_out_phenograph,df_out_svm=df_out_svm,
  out_mutation=out_mutation,
  chrom_n=chrom_n,CNV_label=CNV_label,
  fusion_fc=fusion_fc,fusion_c=fusion_c){

  #Gene expression
  subtype_exp=NA
  if(df_feateure_exp$Expression[df_feateure_exp$Gene=="CRLF2"]>=9){subtype_exp="CRLF2"}
  if(df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"]>=9){ subtype_exp="CDX2/UBTF"}
  if(df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"]>=9){subtype_exp="NUTM1"}
  # if(df_feateure_exp$Expression[df_feateure_exp$Gene=="HLF"]>=9){  subtype_exp="HLF"}
  # if(df_feateure_exp$Expression[df_feateure_exp$Gene=="DUX4"]>=9){ subtype_exp="DUX4"}

  #CNV
  Chr21Alteration=ifelse(grepl("21",CNV_label),"Yes","No")

  if(chrom_n>=51){
    subtype_cnv="Hyperdiploid"
  } else if (chrom_n>=47) {
    subtype_cnv="Low hyperdiploid"
  } else if (chrom_n>=39) {
    subtype_cnv=NA
  } else if (chrom_n >=31) {
    subtype_cnv="Low hypodiploid"
  } else {
    subtype_cnv="Near haploid"
  }

  #fusion
  if(nrow(fusion_fc)>0){
    fusion_fc_=fusion_fc %>% select(FusionInFile,Spanning_pairs,Spanning_unique_reads) %>% distinct()
  } else {
    fusion_fc_=data.frame()
  }

  if(nrow(fusion_c)>0){
    fusion_c_=fusion_c %>% select(FusionInFile,readsA,readsB) %>% distinct()
  } else {
    fusion_c_=data.frame()
  }

  fusion_out=bind_rows(fusion_fc[-(2:3)],fusion_c[-(2:3)]) %>% distinct()
  subtype_fusion=sort(unique(fusion_out$PossibleSubtype))

  #GEP
  pg_=get_pred_result(df_out_phenograph)
  score_phenograph=get_pred_score(df_out_phenograph)
  subtype_phenograph_label=get_pred_label(df_out_phenograph)
  gp_all=df_out_phenograph$pred

  svm_=get_pred_result(df_out_svm)
  score_svm=get_pred_score(df_out_svm)
  subtype_svm_label=get_pred_label(df_out_svm)
  svm_all=df_out_svm$pred

  gep_all=sort(unique(c(gp_all,svm_all)))

  #GEP for PAX5alt
  pg1=get_pred_result(df_out_phenograph,count_PAX5ETV6_as_PAX5alt = T)
  score_phenograph1=get_pred_score(df_out_phenograph,count_PAX5ETV6_as_PAX5alt = T)

  svm1=get_pred_result(df_out_svm,count_PAX5ETV6_as_PAX5alt = T)
  score_svm1=get_pred_score(df_out_svm,count_PAX5ETV6_as_PAX5alt = T)

  if(svm_==pg_){subtype_GEP=svm_} else {subtype_GEP=sort(c(svm_,pg_))}
  if(svm1==pg1){subtype_GEP1=svm1} else {subtype_GEP1=sort(c(svm1,pg1))}

  #summarise ----
  like_subtypes=c("Ph-like","KMT2A-like","ZNF384-like","ETV6::RUNX1-like","PAX5alt")
  names(like_subtypes)=c("Ph","KMT2A","ZNF384","ETV6::RUNX1","PAX5alt")

  subtype_final=NA;confidence=NA;stage=NA

  if(is.na(subtype_final)){
    if(any(subtype_exp %in% c("CDX2/UBTF","HLF","NUTM1"))){
      subtype_exp_=subtype_fusion[subtype_fusion %in% c("CDX2/UBTF","HLF","NUTM1")]
      if(any(subtype_exp_ %in% subtype_GEP)){
        subtype_final=subtype_exp_[subtype_exp_ %in% subtype_GEP]
        subtype_final=subtype_final[1]
        stage="HighGeneExp"
      }
    }
  }

  if(is.na(subtype_final)){
    if(any(subtype_fusion %in% c('ETV6::RUNX1','KMT2A','TCF3::PBX1','MEF2D','ZNF384','PAX5::ETV6','BCL2/MYC','DUX4',"NUTM1","CDX2/UBTF","HLF","Ph","ZEB2/CEBP"))){
      subtype_fusion_=subtype_fusion[subtype_fusion %in% c('ETV6::RUNX1','KMT2A','TCF3::PBX1','MEF2D','ZNF384','PAX5::ETV6','BCL2/MYC','DUX4',"NUTM1","CDX2/UBTF","HLF","Ph","ZEB2/CEBP")]
      if(any(subtype_fusion_ %in% subtype_GEP)){
        subtype_final=subtype_fusion_[subtype_fusion_ %in% subtype_GEP]
        subtype_final=subtype_final[1]
        stage="Fusion"
      }
    }
  }

  if(is.na(subtype_final)){if(subtype_cnv %in% "Near haploid"     & ("Hyperdiploid"   %in% subtype_GEP | "Low hypodiploid" %in% subtype_GEP)){subtype_final="Near haploid";stage="Aneuploid"}}
  if(is.na(subtype_final)){if(subtype_cnv %in% "Hyperdiploid"     & "Hyperdiploid"    %in% subtype_GEP){subtype_final="Hyperdiploid";    stage="Aneuploid"}}
  if(is.na(subtype_final)){if(subtype_cnv %in% "Low hyperdiploid" & "Hyperdiploid"    %in% subtype_GEP){subtype_final="Low hyperdiploid";stage="Aneuploid";confidence="Low"}}
  if(is.na(subtype_final)){if(subtype_cnv %in% "Low hypodiploid"  & "Low hypodiploid" %in% subtype_GEP){subtype_final="Low hypodiploid"; stage="Aneuploid"}}
  if(is.na(subtype_final)){if(Chr21Alteration %in% "Yes"          & "iAMP21"          %in% subtype_GEP){subtype_final="iAMP21";          stage="Aneuploid"}}

  if(is.na(subtype_final)){if("IKZF1 N159Y" %in% out_mutation$subtype_mutation & "IKZF1 N159Y" %in% subtype_GEP){subtype_final="IKZF1 N159Y";stage="Mutation"}}
  if(is.na(subtype_final)){if("PAX5 P80R" %in% out_mutation$subtype_mutation & "PAX5 P80R" %in% subtype_GEP){subtype_final="PAX5 P80R";stage="Mutation"}}
  if(is.na(subtype_final)){if("ZEB2/CEBP" %in% out_mutation$subtype_mutation & "ZEB2/CEBP" %in% subtype_GEP){subtype_final="ZEB2/CEBP";stage="Mutation"}}

  if(is.na(subtype_final)){if(!"ETV6::RUNX1" %in% subtype_fusion &  "ETV6::RUNX1" %in% subtype_GEP){subtype_final="ETV6::RUNX1-like";stage="like"}}
  if(is.na(subtype_final)){if(!"KMT2A" %in% subtype_fusion &  "KMT2A" %in% subtype_GEP){subtype_final="KMT2A-like";stage="like"}}
  if(is.na(subtype_final)){if(!"ZNF384" %in% subtype_fusion &  "ZNF384" %in% subtype_GEP){subtype_final="ZNF384-like";stage="like"}}

  if(is.na(subtype_final)){if(!"Ph" %in% subtype_fusion & "Ph" %in% subtype_GEP & "Ph-like" %in% subtype_fusion){subtype_final="Ph-like";stage="Ph-like"}}
  if(is.na(subtype_final)){if(!"Ph" %in% subtype_fusion & "Ph" %in% subtype_GEP & "CRLF2" %in% subtype_exp){subtype_final="Ph-like";stage="Ph-like"}}
  if(is.na(subtype_final)){if(!"Ph" %in% subtype_fusion & svm_=="Ph" & pg_=="Ph"){subtype_final="Ph-like";stage="Ph-like"}}

  if(is.na(subtype_final)){if("PAX5alt" %in% subtype_GEP & "PAX5alt" %in% subtype_fusion){subtype_final="PAX5alt";stage="PAX5alt"}}
  if(is.na(subtype_final)){if("PAX5alt" %in% subtype_GEP & grepl("PAX5",out_mutation$out_text_BALLmutation)){subtype_final="PAX5alt";stage="PAX5alt"}}
  if(is.na(subtype_final)){if(svm_=="PAX5alt" & pg_=="PAX5alt"){subtype_final="PAX5alt";stage="PAX5alt"}}

  if(is.na(subtype_final)){if(subtype_exp %in% "CRLF2" & (!"Ph" %in% subtype_GEP) & "CRLF2(non-Ph-like)" %in% subtype_fusion){subtype_final="CRLF2(non-Ph-like)";confidence="Low";stage="CRLF2"}}

  if(is.na(subtype_final)){
    x=(data.frame(g=c(pg_,svm_),s=c(score_phenograph,score_svm)) %>% arrange(desc(s)) %>% slice_head(n=1) %>% filter(s>0.5))$g
    if(length(x)>=1){
    if(x %in% names(like_subtypes)){x=like_subtypes[names(like_subtypes)==x]}
    subtype_final=x;stage="GEPalone"
    }
  }

  if(is.na(subtype_final)){subtype_final="Other";confidence="Low";stage="Other"}

  #Modify subtype_final
  if("iAMP21" %in% subtype_GEP & length(subtype_GEP)==2){
    if(((pg_=="iAMP21" & score_phenograph > 0.9) | (svm_=="iAMP21" & score_svm > 0.9)) & "Yes" %in% Chr21Alteration){
      if(subtype_final!="iAMP21"){stage="Modify_iAMP21"}
      subtype_final="iAMP21"
    }
  }

  if("Ph" %in% subtype_GEP & length(subtype_GEP)==2){
    if(((pg_=="Ph" & score_phenograph > 0.9) | (svm_=="Ph" & score_svm > 0.9)) & ("CRLF2" %in% subtype_exp | "Ph-like" %in% subtype_fusion)){
      if(subtype_final!="Ph-like"){stage="Modify_Ph-like"}
      subtype_final="Ph-like";
    }
  }

  if("NUTM1" %in% subtype_GEP & length(subtype_GEP)==2){
    if(((pg_=="NUTM1" & score_phenograph > 0.9) | (svm_=="NUTM1" & score_svm > 0.9)) & ("NUTM1" %in% subtype_exp | "NUTM1" %in% subtype_fusion)){
      if(subtype_final!="NUTM1"){stage="Modify_NUTM1"}
      subtype_final="NUTM1"
    }
  }

  if(svm_=="Low hypodiploid" & pg_=="Low hypodiploid" & grepl("TP53",out_mutation$out_text_BALLmutation)){
    if(subtype_final!="Low hypodiploid"){stage="Modify_LowHypodiploid"}
    subtype_final="Low hypodiploid"
  }

  if("Ph" %in% subtype_fusion  & "Ph" %in% gep_all){
    if(subtype_final!="Ph"){stage="Modify_Ph"}
    subtype_final="Ph";stage="Modify_Ph"
  }

  #Modify confidence
  define_confidence1=function(subtype_final,pg_,svm_,score_svm,score_phenograph,subtype_fusion){
    confidence=NA

    # Feature gene expression subtypes
    if(subtype_final %in% c("CDX2/UBTF","NUTM1","HLF")){
      if(((pg_==subtype_final & score_phenograph > 0.9) | (svm_==subtype_final & score_svm > 0.9)) & (subtype_final %in% subtype_exp)){confidence="High"}
      if(pg_==subtype_final & pg_==svm_ & subtype_final  %in%  subtype_exp){confidence="High"}
    }

    # Fusion subtypes
    if(subtype_final %in% c("NUTM1",'ETV6::RUNX1','KMT2A','TCF3::PBX1','MEF2D','ZNF384','PAX5::ETV6','BCL2/MYC','DUX4',"NUTM1","HLF","Ph","ZEB2/CEBP","HLF","NUTM1","CDX2/UBTF")){
      if(pg_==subtype_final & pg_==svm_ & subtype_final %in% subtype_fusion){confidence="High"}
      if(((pg_==subtype_final & score_phenograph > 0.9) | (svm_==subtype_final & score_svm > 0.9)) & subtype_final %in% subtype_fusion){confidence="High"}
    }

    #Aneuploid change subtypes
    Aneuploid_subtypes=c("Near haploid","Hyperdiploid","Low hypodiploid")
    Aneuploid_GEPs=list(c("Hyperdiploid","Low hypodiploid"),"Hyperdiploid","Low hypodiploid")
    if(subtype_final %in% Aneuploid_subtypes){
      subtype_final1=Aneuploid_GEPs[[match(subtype_final,Aneuploid_subtypes)]]
      if(pg_ %in% subtype_final1 & svm_ %in% subtype_final1 & subtype_final %in% subtype_cnv){confidence="High"}
      if(((pg_ %in% subtype_final1 & score_phenograph > 0.9) | (svm_ %in% subtype_final1 & score_svm > 0.9)) & subtype_final %in% subtype_cnv){confidence="High"}
    }

    #iAMP21
    if(subtype_final == "iAMP21"){
      if(pg_==subtype_final & score_phenograph > 0.9 & svm_==subtype_final & score_svm > 0.9 & Chr21Alteration == "Yes"){confidence="High"}
    }

    # Mutation subtypes
    if(subtype_final %in% c("PAX5 P80R","IKZF1 N159Y","ZEB2/CEBP")){
      if(pg_==subtype_final & pg_==svm_ & subtype_final  %in%  out_mutation$subtype_mutation){confidence="High"}
      if(((pg_==subtype_final & score_phenograph > 0.9) | (svm_==subtype_final & score_svm > 0.9)) & (subtype_final %in% out_mutation$subtype_mutation)){confidence="High"}
    }

    # -like subytpes
    like_subtypes=c("Ph-like","KMT2A-like","ZNF384-like","ETV6::RUNX1-like","PAX5alt")
    names(like_subtypes)=c("Ph","KMT2A","ZNF384","ETV6::RUNX1","PAX5alt")

    if(subtype_final %in% like_subtypes){
      subtype_final1=names(like_subtypes)[like_subtypes==subtype_final]
      if(pg_==subtype_final1 & pg_==svm_ & ((pg_==subtype_final1 & score_phenograph > 0.9) | (svm_==subtype_final1 & score_svm > 0.9))){confidence="High"}
    }

    # DUX4
    if(subtype_final %in% c("DUX4")){
      if(subtype_final==pg_ & pg_==svm_ & score_svm>0.9 & score_phenograph>0.9){confidence="High"}
    }

    # Ph-like
    if(subtype_final %in% c("Ph-like")){
      if("Ph"==pg_ & pg_==svm_ & "Ph-like" %in%  subtype_fusion){confidence="High"}
      if("Ph"==pg_ & pg_==svm_ & "CRLF2" %in%  subtype_exp){confidence="High"}
    }

    #PAX5alt
    if(subtype_final %in% c("PAX5alt")){
      if(subtype_final==pg_ & pg_==svm_ & "PAX5" %in% subtype_fusion){confidence="High"}
      if(subtype_final==pg_ & pg_==svm_ & grepl("PAX5",out_mutation$out_text_BALLmutation)){confidence="High"}
    }

    confidence
  }

  if(is.na(confidence)){
    confidence=define_confidence1(subtype_final,pg_,svm_,score_svm,score_phenograph,subtype_fusion)
  }

  if(is.na(confidence)){
    confidence="Low"
  }

  #Modify PAX5alt
  if((!(subtype_final=='PAX5::ETV6' | subtype_final=="PAX5 P80R")) & ((!is.na(confidence) & confidence=="Low")|is.na(confidence))){
    if("PAX5alt" %in%  pg1 & "PAX5alt" %in% svm1){subtype_final='PAX5alt';stage="PAX5merge"}
    if("PAX5alt" %in%  pg1 & score_phenograph1>0.9 & "PAX5alt" %in% svm1 & score_svm1>0.9){subtype_final='PAX5alt';stage="PAX5merge";confidence="High"}
  }

  #get output table ----
  df_sum=
    data.frame(
      sample_id=id,
      CDX2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"],
      CRLF2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CRLF2"],
      # HLF=df_feateure_exp$Expression[df_feateure_exp$Gene=="HLF"],
      NUTM1=df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"],
      # DUX4=df_feateure_exp$Expression[df_feateure_exp$Gene=="DUX4"],

      fusion_fc=  ifelse(
        nrow(fusion_fc_)>=1,
        paste0(paste0(fusion_fc_$FusionInFile,"(",
                      paste0(fusion_fc_$Spanning_pairs,";",fusion_fc_$Spanning_unique_reads),")"),collapse = ","),
        NA),
      fusion_c=  ifelse(
        nrow(fusion_c_)>=1,
        paste0(paste0(fusion_c_$FusionInFile,"(",
                      paste0(fusion_c_$readsA,";",fusion_c_$readsB),")"),collapse = ","),NA),
      Mutation_BALL=gsub("\n","|",out_mutation$out_text_BALLmutation),
      Mutation_Sub_def=out_mutation$out_text_SubtypeDefiningMutation,

      RNAseqCNV_ChromN=chrom_n,
      RNAseqCNV_label=gsub("\n",";",CNV_label),
      Chr21Alteration=Chr21Alteration,

      subtype_phenograph=pg_,
      score_phenograph=score_phenograph,
      subtype_phenograph_label=subtype_phenograph_label,

      subtype_svm=svm_,
      score_svm=score_svm,
      subtype_svm_label=subtype_svm_label,

      subtype_exp=subtype_exp,
      subtype_mutation=out_mutation$subtype_mutation,
      subtype_cnv=subtype_cnv,

      subtype_fusion=paste0(subtype_fusion,collapse = ","),

      subtype_GEP=paste0(subtype_GEP,collapse = ","),

      subtype_final_possible=subtype_final,
      confidence=confidence,
      stage=stage,

      stringsAsFactors = F
    )
  df_sum
}

#' run_one_sample
#'
#' @param sample_id
#' @param file_count
#' @param file_vcf
#' @param file_fusioncatcher
#' @param file_cicero
#' @param featureN_PG
#' @param minReadCnt
#' @param minDepth
#' @param mafmin
#' @param mafmax
#'
#' @return
#' @export
#'
#' @examples
run_one_sample=function(sample_id="",file_count,file_vcf,file_fusioncatcher="",file_cicero="",featureN_PG=c(100),minReadCnt=3,minDepth=20,mafmin=0.1,mafmax=0.85){
  df_count=read_input(file_count,delimiter = "\t",header = F)
  df_vst=get_vst_values(obj_in = obj_234_HTSeq,df_count = df_count)
  #imputation
  df_vst_i=f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vst)

  #get feature gene expression
  df_feateure_exp=get_geneExpression(df_vst = df_vst,genes = c("CDX2","CRLF2","NUTM1"))

  # Add testing sample to reference dataset for subtype prediction
  obj_=obj_merge(obj_in = obj_1821,df_in = df_vst_i,assay_name_in = "vst")

  #Get phenograph prediction
  df_out_phenograph=get_PhenoGraphPreds(obj_in = obj_,feature_panel = "keyFeatures",SampleLevel = "TestSample",
                                        neighbor_k = 10,
                                        # variable_n_list = c(seq(100,1000,100),1058)
                                        variable_n_list = c(100)
  )

  #Get SVM prediction
  df_out_svm=get_SVMPreds(models_svm,df_in = df_vst_i)

  df_pred=bind_rows(df_out_phenograph,df_out_svm) %>% mutate(N=sprintf("%04d",featureN))

  # Get RNAseqCNV
  RNAseqCNV_out=run_RNAseqCNV(df_count = df_count,snv_file = file_vcf,minReadCnt = minReadCnt,minDepth = minDepth,mafRange = c(mafmin,mafmax))
  CNV_label=paste0(RNAseqCNV_out$df_cnv_out$gender,";\n",RNAseqCNV_out$df_cnv_out$chrom_n,",",RNAseqCNV_out$df_cnv_out$alterations)
  chrom_n=RNAseqCNV_out$df_cnv_out$chrom_n

  # Get subtype defining mutation
  out_mutation=get_BALL_mutation(file_vcf)

  # Get fusions
  fusion_fc=get_BALL_fusion(file_fusioncatcher,type = "fc")

  fusion_c=get_BALL_fusion(file_cicero,type = "c")

  # Get summary
  df_sum=get_subtype_final(
    id="TestSample",
    df_feateure_exp = df_feateure_exp,
    df_out_phenograph = df_out_phenograph,df_out_svm = df_out_svm,
    out_mutation = out_mutation,
    chrom_n = chrom_n,CNV_label = CNV_label,
    fusion_fc = fusion_fc,fusion_c = fusion_c)
  df_sum
}

#' run_multiple_samples
#'
#' @param file_listing
#' @param featureN_PG
#' @param minReadCnt
#' @param minDepth
#' @param mafmin
#' @param mafmax
#'
#' @return
#' @export
#'
#' @examples
run_multiple_samples=function(file_listing,featureN_PG=c(100),minReadCnt=3,minDepth=20,mafmin=0.1,mafmax=0.85){
  df_listing=as.data.frame(vroom::vroom(file_listing,progress = FALSE,show_col_types=F))

  df_listing=df_listing %>% mutate(obs=1:n())

  files_count=df_listing$count

  for (i in df_listing$obs){
    file_count=files_count[i]
    df_count=read_input(file_count)
    names(df_count)=c("feature",df_listing$id[i])
    if (i==1){df_counts=df_count}
    if (i >1){df_counts=df_counts %>% left_join(df_count)}
  }

  df_vsts=get_vst_values(obj_in = obj_234_HTSeq,df_count = df_counts)

  # df_vsts=df_vsts %>% sample_frac(0.99)

  df_vsts_i=list()
  df_outs_phenograph=list()
  df_outs_svm=list()
  RNAseqCNV_outs=list()
  i=1
  for(i in df_listing$obs){
    sample_id=df_listing$id[i]

    #imputation
    df_vst_i=f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vsts[c("feature",sample_id)])

    df_vsts_i[[sample_id]]=df_vst_i

    #get obj
    obj_=obj_merge(obj_in = obj_1821,df_in = df_vst_i,assay_name_in = "vst")

    #phenograph
    df_outs_phenograph[[sample_id]]=get_PhenoGraphPreds(obj_in = obj_,feature_panel = "keyFeatures",SampleLevel = sample_id,
                                                        neighbor_k = 10,
                                                        # variable_n_list = c(seq(100,1000,100),1058)
                                                        variable_n_list = featureN_PG
    )

    #svm
    df_outs_svm[[sample_id]]=get_SVMPreds(models_svm,df_in = df_vst_i,id = sample_id)

    #RNAseqCNV
    RNAseqCNV_outs[[sample_id]]=run_RNAseqCNV(df_count = read_input(files_count[i]),snv_file = df_listing$vcf[i],
                                              genome_version = "hg38",
                                              minReadCnt = minReadCnt,
                                              minDepth = minDepth,
                                              mafRange = c(mafmin,mafmax))
  }

  df_sums=bind_rows(lapply(df_listing$obs, function(i){
    sample_id=df_listing$id[i]

    #gene exp
    df_feateure_exp=get_geneExpression(df_vst = df_vsts[c("feature",sample_id)],genes = c("CDX2","CRLF2","NUTM1"))

    #prediction
    df_out_phenograph=df_outs_phenograph[[sample_id]]
    df_out_svm=df_outs_svm[[sample_id]]

    #RNAseqCNV
    RNAseqCNV_out=RNAseqCNV_outs[[sample_id]]
    CNV_label=paste0(RNAseqCNV_out$df_cnv_out$gender,";\n",RNAseqCNV_out$df_cnv_out$chrom_n,",",RNAseqCNV_out$df_cnv_out$alterations)
    chrom_n=RNAseqCNV_out$df_cnv_out$chrom_n

    #mutation
    out_mutation=get_BALL_mutation(df_listing$vcf[i])

    #fusioncatcher
    fusion_fc=get_BALL_fusion(df_listing$fusioncatcher[i],type = "fc")

    #cicero
    fusion_c=get_BALL_fusion(df_listing$cicero[i],type = "c")

    #sum
    df_sum=get_subtype_final(
      id=sample_id,
      df_feateure_exp = df_feateure_exp,
      df_out_phenograph = df_out_phenograph,df_out_svm = df_out_svm,
      out_mutation = out_mutation,
      chrom_n = chrom_n,CNV_label = CNV_label,
      fusion_fc = fusion_fc,fusion_c = fusion_c)
    df_sum
  }))

  list(df_counts=df_counts,df_vsts=df_vsts,df_vsts_i=df_vsts_i,df_outs_phenograph=df_outs_phenograph,df_outs_svm=df_outs_svm,RNAseqCNV_outs=RNAseqCNV_outs,df_sums=df_sums)
}


#' run_countMatrix
#'
#' @param file_countMatrix
#' @param featureN_PG
#'
#' @return
#' @export
#'
#' @examples
run_countMatrix=function(file_countMatrix,featureN_PG=c(1000,1058)){
  df_counts=as.data.frame(vroom::vroom(file_countMatrix))
  names(df_counts)[1]="feature"

  ids=names(df_counts)[-1]
  df_vsts=get_vst_values(obj_in = obj_234_HTSeq, df_count = df_counts)


  df_vsts_i=list()
  df_outs_phenograph=list()
  df_outs_svm=list()
  for(i in 1:length(ids)){
    sample_id=ids[i]

    #imputation
    df_vst_i=f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vsts[c("feature",sample_id)])

    df_vsts_i[[sample_id]]=df_vst_i

    #get obj
    obj_=obj_merge(obj_in = obj_1821,df_in = df_vst_i,assay_name_in = "vst")

    #phenograph
    df_outs_phenograph[[sample_id]]=get_PhenoGraphPreds(obj_in = obj_,feature_panel = "keyFeatures",SampleLevel = sample_id,
                                                        neighbor_k = 10,
                                                        # variable_n_list = c(seq(100,1000,100),1058)
                                                        variable_n_list = featureN_PG
    )

    #svm
    df_outs_svm[[sample_id]]=get_SVMPreds(models_svm,df_in = df_vst_i,id = sample_id)
  }

  df_sums=bind_rows(lapply(1:length(ids), function(i){
    sample_id=ids[i]

    #gene exp
    df_feateure_exp=get_geneExpression(df_vst = df_vsts[c("feature",sample_id)],genes = c("CDX2","CRLF2","NUTM1"))

    #prediction
    df_out_phenograph=df_outs_phenograph[[sample_id]]
    df_out_svm=df_outs_svm[[sample_id]]

    #get GEP
    pg_=get_pred_result(df_out_phenograph)
    score_phenograph=get_pred_score(df_out_phenograph)
    subtype_phenograph_label=get_pred_label(df_out_phenograph)
    gp_all=df_out_phenograph$pred

    svm_=get_pred_result(df_out_svm)
    score_svm=get_pred_score(df_out_svm)
    subtype_svm_label=get_pred_label(df_out_svm)
    svm_all=df_out_svm$pred

    #sum
    df_sum=
      data.frame(
        sample_id=sample_id,
        CDX2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"],
        CRLF2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CRLF2"],
        NUTM1=df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"],

        subtype_phenograph=pg_,
        score_phenograph=score_phenograph,
        subtype_phenograph_label=subtype_phenograph_label,

        subtype_svm=svm_,
        score_svm=score_svm,
        subtype_svm_label=subtype_svm_label
      )
  }))

  list(df_counts=df_counts,df_vsts=df_vsts,df_vsts_i=df_vsts_i,df_outs_phenograph=df_outs_phenograph,df_outs_svm=df_outs_svm,df_sums=df_sums)
}


#' run_multiple_umap
#'
#' @param outs_multiple
#' @param n_neighbors
#' @param variable_n
#'
#' @return
#' @export
#'
#' @examples
run_multiple_umap=function(outs_multiple,n_neighbors=10,variable_n=1058){
  sample_ids=outs_multiple$df_sums$sample_id

  print(paste0("UMAP for: ",paste0(sample_ids,collapse = ",")))

  umaps=list()

  for(sample_id in sample_ids){
    message(paste0("Running UMAP for ",sample_id,"\n"))
    obj_=obj_merge(obj_in = obj_1821,df_in = outs_multiple$df_vsts_i[[sample_id]],assay_name_in = "vst")
    obj_=run_umap(obj_in = obj_,out_label = "umap",n_neighbors = n_neighbors,variable_n = variable_n,feature_panel = "keyFeatures")
    umaps[[sample_id]]=draw_DimPlot(obj_,group.by = "diag_raw",reduction = "umap",highlightLevel = "TestSample")+
      theme(legend.position = "bottom")
  }

  umaps

}


