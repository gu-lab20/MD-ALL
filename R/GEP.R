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
get_pred_result=function(df_pred){
  (df_pred %>% mutate(n=n()) %>%  group_by(pred) %>%
     mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
     select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$out
}

#' get_pred_score
#'
#' @param df_pred
#'
#' @return
#' @export
#'
#' @examples
get_pred_score=function(df_pred){
  (df_pred %>% mutate(n=n()) %>%  group_by(pred) %>%
     mutate(n_pred=n(),percentage=round(n_pred/n,2),out=ifelse(percentage>0.5,pred,"Unclassified")) %>% ungroup() %>%
     select(pred,percentage,out) %>% distinct() %>% arrange(desc(percentage)) %>% slice_head(n=1))$percentage
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
  if(df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"]>=9){subtype_exp="CDX2/UBTF"}
  if(df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"]>=9.5){subtype_exp="NUTM1"}


  #CNV
  Chr21Alteration=ifelse(grepl("21",CNV_label),"Yes","No")

  if(chrom_n>=51){
    subtype_cnv="Hyperdiploid"
  } else if (chrom_n>=47) {
    subtype_cnv="Low Hyperdiploid"
  } else if (chrom_n>=39) {
    subtype_cnv=NA
  } else if (chrom_n >=31) {
    subtype_cnv="Low Hypodiploid"
  } else {
    subtype_cnv="Near Haploid"
  }

  #fusion
  fusion_fc_top=fusion_fc %>% slice_head(n=2)
  fusion_c_top=fusion_c %>% slice_head(n=2)
  fusion_out=bind_rows(fusion_fc[-(2:3)],fusion_c[-(2:3)]) %>% distinct()


  #GEP
  gep_=NA;subtype_GEP_def=NA
  pg_=get_pred_result(df_out_phenograph)
  score_phenograph=get_pred_score(df_out_phenograph)
  subtype_phenograph_label=get_pred_label(df_out_phenograph)
  gp_all=df_out_phenograph$pred



  svm_=get_pred_result(df_out_svm)
  score_svm=get_pred_score(df_out_svm)
  subtype_svm_label=get_pred_label(df_out_svm)
  svm_all=df_out_svm$pred


  gep_all=sort(unique(c(gp_all,svm_all)))



  # pg_="y";svm_="x"
  if(pg_==svm_){gep_=svm_}
  if(!pg_==svm_){
    if(any(c(pg_,svm_) %in% fusion_out$PossibleSubtype)){
      gep_=c(pg_,svm_)[c(pg_,svm_) %in% fusion_out$PossibleSubtype]
    } else {
      if(all(c(svm_,pg_)=="Unclassified")){gep_="Unclassified"}
      else if ("Unclassified" %in% c(svm_,pg_)){
        gep_=c(svm_,pg_)[-match("Unclassified",c(svm_,pg_))]
      } else {
        gep_=c(svm_,pg_)
      }
    }
  }

  subtype_GEP=gep_

  if("ETV6::RUNX1" %in% subtype_GEP){gep_all=unique(c(gep_all,"ETV6::RUNX1-like"))}
  if("Ph" %in% subtype_GEP){gep_all=unique(c(gep_all,"Ph-like"))}

  if("ETV6::RUNX1" %in% subtype_GEP){subtype_GEP=unique(c(subtype_GEP,"ETV6::RUNX1-like"))}
  if("Ph" %in% subtype_GEP){subtype_GEP=unique(c(subtype_GEP,"Ph-like"))}

  #fusion sum
  fusion_fc_=NA;fusion_c_=NA;subtype_fusion_GEP=NA;subtype_fusion_alone=NA;df_fusion_fc=data.frame();df_fusion_c=data.frame()
  if(any(gep_all %in% c(fusion_fc$PossibleSubtype,fusion_c$PossibleSubtype))){

    if(any(gep_all %in% fusion_fc$PossibleSubtype)){
      df_fusion_fc=fusion_fc %>% filter(PossibleSubtype %in% gep_all) %>%
        arrange(desc(Spanning_unique_reads)) %>% slice_head(n=2) %>% arrange(FusionInFile)
      fusion_fc_=df_fusion_fc$FusionInFile
    }

    if(any(gep_all %in% fusion_c$PossibleSubtype)){
      df_fusion_c=fusion_c %>% filter(PossibleSubtype %in% gep_all) %>%
        arrange(desc(readsB)) %>% slice_head(n=2) %>% arrange(FusionInFile)
      fusion_c_=df_fusion_c$FusionInFile
    }

    df_fusion_=bind_rows(df_fusion_fc,df_fusion_c) %>%
      select(FusionInFile,PossibleSubtype) %>% distinct() %>% arrange(FusionInFile)

    subtype_fusion_GEP=sort(unique(df_fusion_$PossibleSubtype))

  }

  if(!any(gep_all %in% c(fusion_fc$PossibleSubtype,fusion_c$PossibleSubtype))){

    fusion_def_subtypes=c("BCL2/MYC","MEF2D","CRLF2","ETV6::RUNX1","ETV6::RUNX1-like")

    df_fusion_out=fusion_out %>%
      group_by(method) %>% slice_head(n=2) %>% arrange(FusionInFile)

    fusion_fc_=df_fusion_out$FusionInFile[df_fusion_out$method=="fc"]
    fusion_c_=df_fusion_out$FusionInFile[df_fusion_out$method=="c"]

    subtype_fusion_alone=unique(df_fusion_out$PossibleSubtype)
  }

  #summarise ----
  subtype_final=NA;confidence="High"

  subtype_final=ifelse(any(subtype_GEP %in% subtype_exp),intersect(subtype_GEP,subtype_exp),NA)
  if(is.na(subtype_final)){if("Hyperdiploid" %in% gep_all & "Hyperdiploid" %in% subtype_cnv){subtype_final="Hyperdiploid"}}

  if(is.na(subtype_final)){subtype_final=ifelse(any(subtype_GEP %in% subtype_fusion_GEP),intersect(subtype_GEP,subtype_fusion_GEP),NA)}

  if(is.na(subtype_final)){if("iAMP21" %in% gep_all & Chr21Alteration=="Yes"){subtype_final="iAMP21"}}
  if(is.na(subtype_final)){if( "Near Haploid" %in% subtype_cnv){subtype_final="Near Haploid"}}
  if(is.na(subtype_final)){if( "Low hypodiploid" %in% subtype_cnv){subtype_final="Low hypodiploid"}}

  if(is.na(subtype_final)){subtype_final=ifelse(any(gep_all %in% subtype_fusion_GEP),intersect(gep_all,subtype_fusion_GEP),NA)}
  if(is.na(subtype_final)){subtype_final=ifelse(any(subtype_GEP %in% subtype_fusion_alone),intersect(subtype_GEP,subtype_fusion_alone),NA)}
  if(is.na(subtype_final)){subtype_final=ifelse(any(subtype_GEP %in% subtype_cnv),intersect(subtype_GEP,subtype_cnv),NA)}
  if(is.na(subtype_final)){subtype_final=ifelse(any(subtype_GEP %in% out_mutation$subtype_mutation),intersect(subtype_GEP,out_mutation$subtype_mutation),NA)}


  if(is.na(subtype_final)){
    if(svm_=="DUX4" & pg_=="DUX4" & score_phenograph==1 & score_svm==1){subtype_final="DUX4"}
  }

  if(is.na(subtype_final)){
    if(any(c('ETV6::RUNX1',"RUNX1::ETV6") %in% c(fusion_fc_,fusion_c_))){subtype_final="ETV6::RUNX1"}
  }

  if(is.na(subtype_final)){
    if(('ETV6::RUNX1' %in% gp_all & 'ETV6::RUNX1' %in% svm_all) & (pg_=='ETV6::RUNX1' | svm_=='ETV6::RUNX1')){
      subtype_final="ETV6::RUNX1-like"
    }
  }

  if(is.na(subtype_final)){
    if ("IKZF1 N159Y" %in% gp_all & "IKZF1 N159Y" %in% svm_all &
        (pg_=="IKZF1 N159Y" | svm_=="IKZF1 N159Y")){subtype_final="IKZF1 N159Y"}
  }

  if(is.na(subtype_final)){if(pg_=="KMT2A" & svm_=="KMT2A"){subtype_final="KMT2A-like"}}

  if(is.na(subtype_final)){if(pg_=="Low hypodiploid" & svm_=="Low hypodiploid"){subtype_final="Low hypodiploid"}}


  if(is.na(subtype_final)){if(pg_=="PAX5alt" & svm_=="PAX5alt"){subtype_final="PAX5alt"}}
  if(is.na(subtype_final)){if((pg_=="PAX5alt" | svm_=="PAX5alt") & grepl("PAX5",out_mutation$out_text_BALLmutation)){subtype_final="PAX5alt"}}

  if(is.na(subtype_final)){if(pg_=="Ph" & svm_=="Ph"){subtype_final="Ph-like"}}
  if(is.na(subtype_final)){if("Ph" %in% gp_all & "Ph" %in% svm_all & "Ph" %in% subtype_GEP){subtype_final="Ph-like"}}

  if(is.na(subtype_final)){if(pg_=="ZNF384" & svm_=="ZNF384"){subtype_final="ZNF384-like"}}

  if(is.na(subtype_final)){
    if(("CRLF2(non-Ph-like)" %in% fusion_fc_top$PossibleSubtype | "CRLF2(non-Ph-like)" %in% fusion_c_top$PossibleSubtype ) &
       "CRLF2" %in% subtype_exp &
       ((!"Ph" %in% pg_)|(!"Ph" %in% svm_))){subtype_final="CRLF2(non-Ph-like)"}
  }

  #Overwrite
  if(pg_=="Hyperdiploid" & "Hyperdiploid" %in% svm_all & "Low Hyperdiploid" %in% subtype_cnv){subtype_final="Low Hyperdiploid"}

  if(any(c('BCL2::IGH','IGH::BCL2',"IGH::MYC","MYC::IGH","IGH::CASC11","FGFR3::IGH","IGH::FGFR3") %in% c(fusion_c_top$FusionInFile,fusion_fc_top$FusionInFile)) &
     min(c(score_phenograph,score_svm))<0.8){
    subtype_final="BCL2/MYC"
    fusion_fc_=fusion_fc$FusionInFile[fusion_fc$PossibleSubtype=="BCL2/MYC"]
    fusion_c_=fusion_c$FusionInFile[fusion_c$PossibleSubtype=="BCL2/MYC"]}
  if(any(c("ETV6::RUNX1","RUNX1::ETV6") %in% c(fusion_c_top$FusionInFile,fusion_fc_top$FusionInFile))){subtype_final="ETV6::RUNX1"}

  if(("PAX5alt" %in% pg_ & "PAX5alt" %in% svm_) & grepl("PAX5",out_mutation$out_text_BALLmutation)){subtype_final="PAX5alt"}
  if(("PAX5alt" %in% pg_ & "PAX5alt" %in% svm_) & "iAMP21" %in% subtype_final){subtype_final="PAX5alt"}

  if((any(c("PAX5alt","Ph") %in% pg_) & any(c("PAX5alt","Ph") %in% svm_)) & "iAMP21" %in% subtype_final){subtype_final=paste0(gep_,collapse = "|");confidence="Medium"}

  if(("Ph" %in% pg_ & "Ph" %in% svm_) & "Ph-like" %in% subtype_fusion_GEP){subtype_final="Ph-like"}
  if(any(c("BCR::ABL1","ABL1::BCR") %in% c(fusion_fc_top$FusionInFile,fusion_c_top$FusionInFile))){subtype_final="Ph"}

  if("Near Haploid" %in% subtype_cnv & any(c("Hyperdiploid","Low hypodiploid") %in% c(pg_,svm_))){subtype_final="Near Haploid"}

  if("NUTM1" %in% subtype_fusion_GEP & "NUTM1" %in% subtype_exp){subtype_final="NUTM1"}
  if("TCF3::PBX1" %in% gp_all & "TCF3::PBX1" %in% svm_all & "TCF3::PBX1" %in% subtype_fusion_GEP){subtype_final="TCF3::PBX1"}
  if("ZNF384" %in% gp_all & "ZNF384" %in% svm_all & "ZNF384" %in% subtype_fusion_GEP){subtype_final="ZNF384"}

  if(is.na(subtype_final)){if(length(subtype_fusion_alone)==1 & "Ph-like" %in% subtype_fusion_alone){subtype_final="Ph-like"}}

  if(is.na(subtype_final)){subtype_final=paste0(gep_,collapse = "|");confidence="Medium"}

  # if("Low Hyperdiploid" %in% pg_)

  #get output table ----
  df_sum=
    data.frame(
      sample_id=id,
      CDX2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CDX2"],
      CRLF2=df_feateure_exp$Expression[df_feateure_exp$Gene=="CRLF2"],
      # PAX5=df_feateure_exp$Expression[df_feateure_exp$Gene=="PAX5"],
      # HLF=df_feateure_exp$Expression[df_feateure_exp$Gene=="HLF"],
      NUTM1=df_feateure_exp$Expression[df_feateure_exp$Gene=="NUTM1"],
      # MEGF10=df_feateure_exp$Expression[df_feateure_exp$Gene=="MEGF10"],
      # DUX4=df_feateure_exp$Expression[df_feateure_exp$Gene=="DUX4"],

      fusion_fc=paste0(fusion_fc_,collapse = ","),
      fusion_c=paste0(fusion_c_,collapse = ","),

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

      subtype_fusion_GEP=paste0(subtype_fusion_GEP,collapse = ","),
      subtype_fusion_alone=paste0(subtype_fusion_alone,collapse = ","),

      subtype_GEP=paste0(gep_,collapse = "|"),

      subtype_final=subtype_final,
      confidence=confidence,

      stringsAsFactors = F
    )
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
  df_listing=as.data.frame(vroom::vroom(file_listing))

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
    RNAseqCNV_outs[[sample_id]]=run_RNAseqCNV(df_count = read_input(files_count[i]),snv_file = df_listing$VCF[i],
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
    out_mutation=get_BALL_mutation(df_listing$VCF[i])

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
run_countMatrix=function(file_countMatrix,featureN_PG=c(100,1058)){
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


