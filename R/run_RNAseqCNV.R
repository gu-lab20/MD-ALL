#' run_RNAseqCNV
#'
#' @param df_count
#' @param mafRange
#' @param genome_version
#' @param snv_format
#' @param model_gend
#' @param model_dip
#' @param model_alter
#' @param model_alter_noSNV
#' @param dpRatioChromEdge
#' @param generate_weights
#' @param minReadCnt
#' @param samp_prop
#' @param weight_samp_prop
#' @param minDepth
#' @param weight_tab
#'
#' @return
#' @export
#'
#' @examples
run_RNAseqCNV=function(df_count,
                       snv_file,
                       mafRange = c(0.1, 0.9),
                       genome_version="hg38",
                       snv_format="vcf",
                       model_gend = model_gender,
                       model_dip = model_dipl,
                       model_alter = model_alt,
                       model_alter_noSNV = model_noSNV,
                       generate_weights = FALSE,
                       minReadCnt=3,
                       samp_prop = 0.8,
                       weight_samp_prop = 1,
                       minDepth = 20,
                       weight_tab=weight_table){

  #get ref dataset ----
  if (genome_version == "hg38") {
    referData = gene_annot_hg38
    keptSNP = dbSNP_hg38
    par_region = pseudoautosomal_regions_hg38
    centr_ref = centromeres_hg38
    diploid_standard = diploid_standard_hg38
  } else if (genome_version == "hg19") {
    referData = gene_annot_hg19
    keptSNP = dbSNP_hg19
    par_region = pseudoautosomal_regions_hg19
    centr_ref = centromeres_hg19
    diploid_standard = diploid_standard_hg19
  }

  #create standard samples ----
  cols_ind <- c(1:20, ncol(diploid_standard))

  standard_samples=as.data.frame(diploid_standard)[,cols_ind]

  #Create estimation table
  est_table <- data.frame(sample = character(),
                          gender = factor(levels = c("female", "male")),
                          chrom_n = integer(),
                          alterations = character(), stringsAsFactors = FALSE)

  #calculate normalized count values with DESeq2 normalization method for batch of samples from the input
  count_norm <- get_norm_exp(count_table = df_count,
                             standard_samples = standard_samples,
                             minReadCnt = minReadCnt,
                             samp_prop = samp_prop,
                             weight_table = weight_tab,
                             weight_samp_prop = weight_samp_prop,
                             generate_weights)


  #calculate median gene expression across diploid reference and analyzed sample
  pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = referData, generate_weights = generate_weights)

  #load SNP data
  smpSNP <- prepare_snv(snv_file = snv_file, centr_ref = centr_ref, snv_format = snv_format, minDepth = minDepth)

  #filter SNP data base on dpSNP database
  smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keptSNP, minDepth = minDepth, mafRange = mafRange)

  #analyze chromosome-level metrics (out-dated)
  smpSNPdata <- calc_chrom_lvl(smpSNPdata.tmp)
  #arm-level metrics
  smpSNPdata_a_2 <- calc_arm(smpSNPdata.tmp)

  #select sample
  sample_name="TestSample"
  # count_norm_samp <- count_norm %>% select(!!quo(tidyselect::all_of(sample_name))) %>% mutate(ENSG = rownames(.))
  count_norm_samp <- count_norm %>% select(sample_name) %>% mutate(ENSG = rownames(.))

  #join reference data and weight data
  count_ns <- count_transform(count_ns = count_norm_samp, pickGeneDFall, refDataExp = referData, weight_table = weight_tab)

  #remove PAR regions
  count_ns <- remove_par(count_ns = count_ns, par_reg = par_region)

  #Calculate metrics for chromosome arms
  feat_tab <- get_arm_metr(count_ns = count_ns, smpSNPdata = smpSNPdata_a_2, sample_name = sample_names, centr_ref = centr_ref)

  #Export the per-arm median of log2 fold change of expression
  log2fold_arm_s <- select(feat_tab, chr, arm, arm_med) %>% mutate(arm_med = round(arm_med, 3))
  colnames(log2fold_arm_s)[3] <- sample_name
  log2fold_arm=log2fold_arm_s

  #estimate gender
  count_ns_gend <- count_norm_samp %>% filter(ENSG %in% "ENSG00000012817") %>%  select(ENSG, !!quo(sample_name)) %>% tidyr::spread(key = ENSG, value = !!quo(sample_name))
  gender <- ifelse(predict(model_gend, newdata = count_ns_gend, type = "response") > 0.5, "male", "female")

  #preprocess data for karyotype estimation and diploid level adjustement
  # model diploid level
  feat_tab$chr_status <- randomForest:::predict.randomForest(model_dip, feat_tab, type = "class")
  #exclude non-informative regions  and
  #if the model was not able to call changes
  #(mainly due problematic density graphs on chromosome X) change the value to unknown
  feat_tab_dipl <- feat_tab %>%
    filter(arm != "p" | !chr %in% c(13, 14, 15, 21)) %>% mutate(chr_status = ifelse(is.na(chr_status), "unknown", as.character(chr_status))) %>%
    metr_dipl()

  #model alteration on chromosome arms an in case of problematic SNV graph, use model without this information included
  print(paste0("Estimating chromosome arm CNV",": ", sample_name))
  feat_tab_alt <- feat_tab_dipl %>%
    filter(chr_status != "unknown") %>%
    mutate(alteration = as.character(randomForest:::predict.randomForest(model_alter, ., type = "class")),
           alteration_prob = apply(randomForest:::predict.randomForest(model_alter, ., type = "prob"), 1, max))

  if (any(feat_tab_dipl$chr_status == "unknown")) {
    feat_tab_alt <- feat_tab_dipl %>%
      filter(chr_status == "unknown") %>%
      mutate(alteration = as.character(randomForest:::predict.randomForest(model_alter_noSNV, ., type = "class")),
             alteration_prob = apply(randomForest:::predict.randomForest(model_alter_noSNV, ., type = "prob"), 1, max)) %>%
      bind_rows(feat_tab_alt)
  }

  feat_tab_alt <- colour_code(feat_tab_alt, conf_tresh = 0.85) %>% group_by(chr) %>%
    mutate(alteration = as.character(alteration), chr_alt = as.character(ifelse(length(unique(alteration)) == 1, unique(alteration), "ab")))

  #estimate karyotype
  kar_list <- gen_kar_list(feat_tab_alt = feat_tab_alt, sample_name = sample_name, gender = gender)
  est_table <- rbind(est_table, kar_list)

  df_cnv_out=cbind(est_table , status = "not checked", comments = "none")

  return(list(df_cnv_out=df_cnv_out,feat_tab_alt=feat_tab_alt,smpSNPdata=smpSNPdata,count_ns=count_ns))
}
