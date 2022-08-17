#args <- commandArgs(trailingOnly = TRUE)
#print(args)

#file_tsv=args[1]

#write_tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

#read tsv --------------
read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}

########################################  Calculate freq, with group variable, one variable ########################################

#cal_freq_one(indata,group,var_cat_one,"freq1")

#group="all"

#############calculation for one group variable
cal_freq_one=function(indata,group,var_cat_one,freq_need="freq1"){
  value_group=unlist(indata[group])
  value_cat_one=unlist(indata[var_cat_one])
  
  factor_group=factor(value_group)
  factor_cat=factor(value_cat_one)
  
  levels(factor_group)[length(levels(factor_group))]
  levels(factor_cat)[length(levels(factor_cat))]
  
  freq_raw=table(value_cat_one,value_group)
  freq_sum_col=apply(freq_raw,2,sum)
  
  #0
  freq1=data.frame(Variable=rep(var_cat_one,length(levels(factor_cat))+1),
                   variable_levels=c(var_cat_one,row.names(freq_raw)),
                   index=c(NA,paste0(var_cat_one,"_",row.names(freq_raw))),
                   stringsAsFactors = F)
  freq1$N_used=NA;freq1$N_used[2]=sum(freq_raw)
  
  #1
  if (freq_need=="freq1"){
    percentage_raw1=freq_raw
    for (i in 1:length(levels(factor_cat))){
      percentage_raw1[i,]=percentage_raw1[i,]/freq_sum_col
    }
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw1[]=cal(percentage_raw1)
    percentage_raw1[]=paste0(freq_raw," (",percentage_raw1,")")
    percentage_raw1=as.data.frame(percentage_raw1,stringsAsFactors = F)
    percentage_raw1=dcast(percentage_raw1,value_cat_one~value_group,value.var = 'Freq')
    
    names(percentage_raw1)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw1)
    
  }
  #2
  if (freq_need=="freq2"){
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw2=freq_raw/sum(freq_raw)
    percentage_raw2[]=cal(percentage_raw2)
    percentage_raw2[]=paste0(freq_raw," (",percentage_raw2,")")
    percentage_raw2=as.data.frame(percentage_raw2,stringsAsFactors = F)
    percentage_raw2=dcast(percentage_raw2,value_cat_one~value_group,value.var = 'Freq')
    names(percentage_raw2)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw2)
  }
  freq_out[]=lapply(freq_out, as.character)
  freq_out
}


#freq_base1=cal_freq_set_cat(data_merge,'all',var_cat_base)
cal_freq_set_cat=function(indata,group,var_cat_all){
  freq_out_set=
    bind_rows(lapply(var_cat_all,function(var_cat_one){
      cal_freq_one(indata,group,var_cat_one,"freq1")
    }))
  freq_out_set
}

########################################  Calculate chisq ########################################
#indata=data_1
#var_group="V3"
#var_freq="V1"
cal_chisq=function(var_freq_list,var_group,indata){
  freq_group  = factor(unlist(indata[var_group]))
  
  bind_rows(lapply(var_freq_list, function(var_freq){
    print(var_freq)
    freq_var    = factor(unlist(indata[var_freq]))
    
    n_min=min(table(freq_group,freq_var))
    
    x1=data.frame(table(freq_group,freq_var))
    x1$value=paste0("value_",x1$freq_var)
    x2=tidyr::spread(x1[-2], value, Freq)
    
    if(any(apply(x2[2:ncol(x2)],1,sum)==0) |  any(apply(x2[2:ncol(x2)],2,sum)==0)){
      chisq_p=NA
    } else{
      
      if(n_min>4){
        chisq_p=chisq.test(freq_var,freq_group)$p.value}
      
      if(n_min<=4){
        chisq_p=fisher.test(table(freq_group,freq_var),simulate.p.value = T,alternative = "two.sided")$p.value
      }
    }
    data.frame(
      Variable=var_freq,
      P_raw=signif(chisq_p,digits = 6),
      P=convert_p(chisq_p),
      stringsAsFactors = F)
  }))
  
}


######################################## Merge freq and chisq, one or multiple x var ########################################

#indata=indata
#group="source1"
#var_cat_one="points_sum"

#cal_freq_chisq_one(indata,"source1","points_sum","freq1")


cal_freq_chisq=function(indata,group,var_cat_list,freq_type){
  bind_rows(lapply(var_cat_list, function(var_cat_one){
    out_freq=cal_freq_one(indata,group,var_cat_one,freq_type)
    out_chisq=cal_chisq(var_cat_one,group,indata)
    out_freq=out_freq %>% left_join(out_chisq)
    n=1:nrow(out_freq)
    out_freq$P=ifelse(n==2,out_freq$P,NA)
    out_freq
  }))
}


########################################  Convert p values ########################################
convert_p=function(P_raw){
  P=ifelse(P_raw<0.001,"<0.001",
           ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                  ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                         ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                )))))
  P
}

#Generate intake --------------------

#get PU--------
get_PU=function(dir,postfix="fq.gz",read_delim="_",read_index=5,read_redundancy="R",fqid_delim="_",fqid_index,fqid_redundancy){
  fq_1stlines=system(paste0("ls ",dir,"*.",postfix,"|xargs -i echo 'paste <(ls {}) <(zcat {}|head -n1)'|bash|sed 's/ /\t/g'"),intern = T)
  df_fq_raw=data.frame(
    file=word(fq_1stlines,1,sep="\t"),
    V2=word(fq_1stlines,2,sep="\t"),
    V3=word(fq_1stlines,3,sep="\t"),
    stringsAsFactors = F
  ) %>%
    mutate(
      PU=paste0(gsub("@","",word(V2,1,sep=":")),":",word(V2,2,sep=":"),":",word(V2,3,sep=":"),".",word(V2,4,sep=":"),".",word(V3,-1,sep=":")),
      i7=word(word(PU,-1,sep="[.]"),1,sep="[+]"),
      name_in_source=basename(file),
      read=gsub(read_redundancy,"",word(name_in_source,read_index,sep=read_delim)),
      fqid=gsub(fqid_redundancy,"",word(name_in_source,fqid_index,sep=fqid_delim))
    )
}

# the input variables in fq included: file,name_in_source,sample_id,order,read,PL,ID,PU,LB,target_name,project_name,subproject,tissue,COH_ID,COH_sample,SM
generate_intake=function(fq1){
  data_read_1=fq1 %>% filter(read==1) %>% group_by(COH_sample,SM) %>% mutate(x=1,TmpId=cumsum(x),x=NULL)
  data_read_2=fq1 %>% filter(read==2) %>% group_by(COH_sample,SM) %>% mutate(x=1,TmpId=cumsum(x),x=NULL)
  
  cohid_multiple=names(table(data_read_1$COH_sample))[table(data_read_1$COH_sample)>=2]
  
  fq3=rbind(data_read_1,data_read_2) %>% arrange(sample_id,SM,TmpId) %>%
    mutate(ID=paste0('B',sprintf("%03d",TmpId)))
  
  intake_exome=fq3 %>% select(file,name_in_source,sample_id,read,ID,SM,PL,PU,LB,target_name,project_name,subproject,tissue,COH_ID,COH_sample)
  intake_exome
}

get_fq_link=function(df_intake,dir_old,dir_new,outdir=""){
  df_ln=df_intake %>% select(name_in_source,read,COH_sample) %>%
    mutate(fq_new=paste0(COH_sample,".R",read,".fq.gz")) %>%
    mutate(fq_name_old=paste0(dir_old,name_in_source),
           fq_name_new=paste0(dir_new,fq_new)) %>%
    mutate(code=paste0("ln -s ",fq_name_old," ",fq_name_new))
  write.table(df_ln["code"],paste0(outdir,"ln.code"),row.names = F,col.names = F,quote = F)
  df_ln
}

#Get HTSeq count matrix --------------------

get_htseq_matrix=function(files_htseq){
  for (file_htseq in files_htseq){
    id=word(basename(file_htseq),1,sep="[.]")
    df=read.table(file_htseq,header = F,sep = "\t") %>%
      filter(grepl("ENSG",V1))
    names(df)=c("feature",id)
    print(match(file_htseq,files_htseq))
    if (match(file_htseq,files_htseq)==1){count_all=df}
    if (match(file_htseq,files_htseq)>1){count_all=count_all %>% left_join(df)}
  }
  count_all=count_all %>% arrange(feature)
  count_all
}

get_Featurecounts_matrix=function(files_featurecounts){
  for (file_featurecounts in files_featurecounts){
    df=read.table(file_featurecounts,header = T,sep = "\t")
    df=df[c(1,ncol(df))]
    
    names(df)=c("feature",word(basename(file_featurecounts),1,sep="[.]"))
    print(match(file_featurecounts,files_featurecounts))
    if (match(file_featurecounts,files_featurecounts)==1){count_all_featurecounts=df}
    if (match(file_featurecounts,files_featurecounts)>1){count_all_featurecounts=count_all_featurecounts %>% left_join(df)}
  }
  count_all_featurecounts=count_all_featurecounts %>% arrange(feature)
  count_all_featurecounts
}

get_count_matrix=function(files){
  for(file in files){
    print(match(file,files))
    id=word(basename(file),1,sep="[.]")
    count_one=read.table(file,header = F,stringsAsFactors = F) 
    names(count_one)=c("feature",id)
    if(match(file,files)==1){count_all=count_one}
    if(match(file,files)>1){count_all=count_all %>% left_join(count_one)}
  }
  count_all
}

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

get_junction_matrix=function(files){
  for(file in files){
    print(match(file,files))
    id=word(basename(file),1,sep="[.]")
    count_one=read_tsv(file)
    names(count_one)=c("feature",id)
    if(match(file,files)==1){out_all=count_one}
    if(match(file,files)>1){out_all=out_all %>% left_join(count_one)}
  }
  out_all
}

#serach vcf file to extract AF with position information, one record --------------------------------------
#chr="chr1"
#start=10489025
#end=10499925
#vcf_file="/home/zuhu/snp142Common.vcf.gz"
#bcftools="/home/zuhu/miniconda3/bin/bcftools"

search_vcf_onerecord=function(chr,start,end,vcf_file,bcftools){
  para_f="'%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/AF\\n'"
  search_region=paste0(chr,":",start,"-",end)
  
  command=paste0(bcftools," query -f ",para_f," ",vcf_file," -r ",search_region," 2> /dev/null")
  search_out=system(command,intern=T)
  
  data_search_out=data.frame(raw_out=search_out,stringsAsFactors = F)
  
  data_search_out=tidyr::separate(data_search_out, raw_out, sep = "\t",into = c("chr","pos","ref","alt","qual","filter","AF"))

  data_search_out=tidyr::separate_rows(data_search_out,alt,sep=",") %>% filter(!alt=="")
  
  data_search_out = data_search_out %>% mutate(start=as.numeric(pos),
                                               length_ref=nchar(ref),
                                               length_alt=nchar(alt),
                                               end=start+ifelse(length_ref>=length_alt,length_ref,length_alt)-1,
                                               length_ref=NULL,length_alt=NULL,pos=NULL)
  
  if(nrow(data_search_out)==0){data_search_out=data.frame(chr=NA,ref=NA,alt=NA,qual=NA,filter=NA,AF=NA,start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  
  names(data_search_out)=paste0(names(data_search_out),"_database")
  data_search_out$chr=chr
  data_search_out$start=start
  data_search_out$end=end
  data_search_out
}

#serach vcf file to extract AF with position information, using data set, with chr, start and end --------------------------------------

search_vcf_dataset=function(indata,vcf_file,bcftools){
  search_out=bind_rows(lapply(1:nrow(indata), function(i){
    data_one=indata[i,]
    
    chr=data_one$chr
    start=data_one$start
    end=data_one$end
    
    search_vcf_onerecord(chr,start,end,vcf_file,bcftools)
  }))
  
  search_out1=indata %>% left_join(search_out) 
  search_out1
}

#serach bb file to extract AF with position information, one record --------------------------------------

search_bb_one=function(chr,start,end,bb_in,bigBedToBed){
  start=as.numeric(start)-1
  command=paste0(bigBedToBed," -chrom=",chr," -start=",start," -end=",end," ",bb_in," stdout")
  search_out=system(command,intern=T)
  data_search_out=data.frame(search_out=search_out,stringsAsFactors = F)
  
  get_max_maf=function(x){
    if(x ==""){x1=NA}
    if(x!=""){x1=max(as.numeric(unlist(strsplit(x,","))))}
    x1
  }
  
  data_search_out=data_search_out %>% 
    mutate(chr=stringr::word(search_out,1,sep="\t"),
           pos=as.numeric(stringr::word(search_out,2,sep="\t"))+1,
           ref=stringr::word(search_out,5,sep="\t"),
           alt=stringr::word(search_out,7,sep="\t"),
           V10=stringr::word(search_out,10,sep="\t"),
    )
  
  data_search_out=tidyr::separate_rows(data_search_out,alt,sep=",") %>% filter(!alt=="")
  data_search_out$AF=unlist(lapply(data_search_out$V10, get_max_maf)) 
  data_search_out = data_search_out %>% mutate(start=pos,
                                               length_ref=nchar(ref),
                                               length_alt=nchar(alt),
                                               end=pos-1+ifelse(length_ref>=length_alt,length_ref,length_alt),
                                               length_ref=NULL,length_alt=NULL,V10=NULL,search_out=NULL,pos=NULL)
  
  if(nrow(data_search_out)==0){data_search_out=data.frame(chr=chr,ref="",alt="",AF=NA,start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  
  names(data_search_out)=paste0(names(data_search_out),"_database")
  data_search_out$chr=chr
  data_search_out$start=start+1
  data_search_out$end=end
  data_search_out
  
  
  
}

#serach bb file to extract AF with position information, using data set, with chr, start and end --------------------------------------
search_bb_dataset=function(indata,bb_in,bigBedToBed){
  search_out=bind_rows(lapply(1:nrow(indata), function(i){
    data_one=indata[i,]
    
    chr=data_one$chr
    start=data_one$start
    end=data_one$end
    
    search_bb_one(chr,start,end,bb_in,bigBedToBed)
  }))
  
  search_out1=indata %>% left_join(search_out) 
  search_out1
}


#possion test --------------------------------
#fit_test=rateratio.test::rateratio.test(c(data_one$Group1,data_one$Group2),c(data_one$n_group1,data_one$n_group2))

#fisher exact test -------------------------------

# scale to [min_target, max_target] --------------------------
scale_target=function(value,min_target,max_target,min_value,max_value){
  ifelse(value==-Inf | value== Inf,NA,
         ((value-min_value)/(max_value-min_value))*(max_target-min_target)+min_target
  )
}





















