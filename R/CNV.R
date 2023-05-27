####Calculate peak statistics for whole chromosome####
#' calc_chrom_lvl
#'
#' @param smpSNPdata.tmp
#'
#' @return
#' @export
#'
#' @examples
calc_chrom_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>% filter(snvOrd<=1000) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
    ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata)
}

####Calculate peak statistics for arms separately#s###
#' Title calc_arm_lvl
#'
#' @param smpSNPdata.tmp
#'
#' @return
#' @export
#'
#' @examples
calc_arm_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata_a=smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
    ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata_a)
}

####find max on y axis from density vector of allele frequency####
#' Title densityMaxY
#'
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
densityMaxY <- function(vec){
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    return(max(d$y))
  }
}

####find peak in a vector####
#' Title findPeak
#'
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
findPeak <- function(vec){
  vec=vec[between(vec, 0.1, 0.9)]
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    maxVec=max(d$y)
    maxPos=d$x[d$y == maxVec]
    # peak=d$x[which(d$y==max(d$y[which(diff(sign(diff(d$y) ))==-2)]))]
    return(maxPos[1])
  }
}

####find peak distance####
#' Title find_peak_dist
#'
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
find_peak_dist <- function(vec) {
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    #Choose the two highest peaks from the MAF density curve
    peaks_max = d$x[which(d$y %in% sort(d$y[which(diff(sign(diff(d$y) ))==-2)], decreasing = TRUE)[1:2])]
    #making sure the distance is measured between symetric peaks
    #symmetry treshold set at 1.2
    if (length(na.omit(peaks_max)) > 1 & data.table::inrange(max(peaks_max - 0.5), 0.42 - min(peaks_max), 0.58 - min(peaks_max))) {
      dist_peak <- abs(peaks_max[1] - peaks_max[2])
      return(dist_peak)
    } else {
      return(0)
    }
  }
}

#### calculate chromosomal statistics ####
# calc arm extended
#' Title calc_arm
#'
#' @param smpSNPdata.tmp
#'
#' @return
#' @export
#'
#' @examples
calc_arm <- function(smpSNPdata.tmp) {

  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peak_m_dist = abs(peak - 0.5), y_0.5 = find_y_0.5(maf), peakdist = find_peak_dist(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red')) %>%
    ungroup() %>% filter(chr %in% c(1:22, "X")) %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata)
}

######height of density curve on 0.5 on x axis
#' Title find_y_0.5
#'
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
find_y_0.5 <- function(vec) {
  vec=vec[between(vec, 0.1, 0.9)]
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    peak_0.5=d$y[which.min(abs(d$x - 0.5))]
    return(peak_0.5)
  }
}

####normalize normalized counts (against median of expression for each gene) and join weight values####
#beta-needs cleaning
#' Title count_transform
#'
#' @param count_ns
#' @param pickGeneDFall
#' @param refDataExp
#' @param weight_table
#'
#' @return
#' @export
#'
#' @examples
count_transform <- function(count_ns, pickGeneDFall, refDataExp, weight_table) {
  count_ns_tmp = count_ns %>% left_join(select(pickGeneDFall, ENSG, med), by = "ENSG") %>%
    mutate(count_nor_med=log2(.[, 1] / med) ) %>% filter(med != 0)
  sENSGinfor=refDataExp[match(count_ns_tmp$ENSG, refDataExp$ENSG), ] %>% select(chr, end, start)

  #keeping only the genes which have weights calculated for geom_poit and boxplot
  count_ns = cbind(sENSGinfor, count_ns_tmp) %>% inner_join(weight_table, by = "ENSG")
}

####filter out par regions####
#' Title remove_par
#'
#' @param count_ns
#' @param par_reg
#'
#' @return
#' @export
#'
#' @examples
remove_par <- function(count_ns, par_reg) {
  #### get rid of PAR regions

  parX = filter(par_reg, chr == "X")
  parX <- as.data.frame(parX)
  count_ns = count_ns %>% filter(chr != "X" | ! data.table::inrange(start, parX[1, 1], parX[1, 2]) & ! data.table::inrange(start, parX[2, 1], parX[2, 2]) &
                                   ! data.table::inrange(end, parX[1, 1], parX[1, 2]) & !  data.table::inrange(end, parX[2, 1], parX[2, 2]))
  parY = filter(par_reg, chr == "Y")
  parY <- as.data.frame(parY)
  count_ns = count_ns %>% filter(chr != "Y" | ! data.table::inrange(start, parY[1, 1], parY[1, 2]) & ! data.table::inrange(start, parY[2, 1], parY[2, 2]) &
                                   ! data.table::inrange(end, parY[1, 1], parY[1, 2]) & !  data.table::inrange(end, parY[2, 1], parY[2, 2]))
}

####get arm metrics####
#' Title get_arm_metr
#'
#' @param count_ns
#' @param smpSNPdata
#' @param sample_name
#' @param centr_ref
#'
#' @return
#' @export
#'
#' @examples
get_arm_metr <- function(count_ns, smpSNPdata, sample_name, centr_ref) {

  #calculate weighted median for every chromosome and use only 1:22
  summ_arm <- count_ns %>% filter(!is.infinite(count_nor_med)) %>% filter(chr %in% c(1:22, "X")) %>% left_join(centr_ref, by = "chr") %>% mutate(chr = factor(chr, levels = c(1:22, "X"))) %>%
    mutate(arm = ifelse(end < cstart, "p", ifelse(end > cend, "q", "centr"))) %>% filter(arm != "centr") %>% group_by(chr, arm) %>%

    # get rid of 21 p arm
    filter(!(chr == 21 & arm == "p")) %>%

    mutate(arm_med = weighted_quantile(x = count_nor_med,w = weight, probs = 0.5, na.rm = TRUE),
           up_quart = weighted_quantile(x = count_nor_med, w = weight, probs = 0.75, na.rm = TRUE),
           low_quart = weighted_quantile(x = count_nor_med, w = weight, probs = 0.25, na.rm = TRUE)) %>%

    ungroup() %>% distinct(chr, arm, arm_med, up_quart, low_quart) %>%
    left_join(distinct(.data = smpSNPdata, chr, arm, peakdist, peak_m_dist, peak_max, y_0.5), by = c("chr", "arm")) %>% mutate(chr = factor(chr, levels = c(1:22, "X")), sd = sd(arm_med), mean_all = mean(arm_med)) %>% ungroup() %>%

    mutate(sds_median = (arm_med - mean_all)/sd, sds_025 = (low_quart - mean_all)/sd, sds_075 = (up_quart-mean_all)/sd,
           n_02_04 = sum(data.table::inrange(peakdist, 0.2, 0.4) & peak_m_dist > 0.08), n_04 = sum(data.table::inrange(peakdist, 0.4, 0.9) & peak_m_dist > 0.08)) %>%

    arrange(chr) %>%

    select(chr, arm, arm_med, up_quart, low_quart, peak_max, peak_m_dist, peakdist, y_0.5, sd, sds_median, sds_025, sds_075, n_02_04, n_04)
  return(summ_arm)

}

#### Function for calculating weighted quantiles ####
## The code is taken from package spatstat, version 1.63.3. The function is no longer part of spatstat package
#' Title weighted_quantile
#'
#' @param x vector of numerical values
#' @param w vector of numerical weights, need to be larger than one
#' @param probs the quantile to be calculated
#' @param na.rm logical; remove NA values
#'
#' @return
#' @export
#'
#' @examples
weighted_quantile <- function (x, w, probs = seq(0, 1, 0.25), na.rm = TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if (anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  stopifnot(all(w >= 0))
  if (all(w == 0))
    stop("All weights are zero", call. = FALSE)
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  result <- numeric(length(probs))
  for (i in seq_along(result)) {
    p <- probs[i]
    lefties <- which(Fx <= p)
    if (length(lefties) == 0) {
      result[i] <- x[1]
    }
    else {
      left <- max(lefties)
      result[i] <- x[left]
      if (Fx[left] < p && left < length(x)) {
        right <- left + 1
        y <- x[left] + (x[right] - x[left]) * (p - Fx[left])/(Fx[right] -
                                                                Fx[left])
        if (is.finite(y))
          result[i] <- y
      }
    }
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE),
                          "%")
  return(result)
}

#### calculate statistics with diploid knowldedge ####
#' Title metr_dipl
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
metr_dipl <- function(data) {

  if (sum(data$chr_status == "dipl") > 3) {
    sd_chr = "dipl"
  } else {
    sd_chr = "no_dipl"
  }

  sd_dipl <- data %>% filter(chr_status == sd_chr) %>% mutate(sd_dipl = sd(arm_med)) %>% pull(sd_dipl) %>% unique()
  mean_dipl <- data %>% filter(chr_status == "dipl") %>% mutate(mean_dipl = mean(arm_med)) %>% pull(mean_dipl) %>% unique()

  data_mod <- data %>% mutate(sd_dipl = sd_dipl, mean_dipl = mean_dipl, sds_median_dipl = (arm_med - mean_dipl)/sd_dipl, sds_025_dipl = (low_quart - mean_dipl)/sd_dipl, sds_075_dipl = (up_quart-mean_dipl)/sd_dipl)

  return(data_mod)
}

### Create colour coding for estimation values ####
#' Title colour_code
#'
#' @param data
#' @param conf_tresh
#'
#' @return
#' @export
#'
#' @examples
colour_code <- function(data, conf_tresh) {
  data_col <- data %>% mutate(colour_arm = factor(ifelse(alteration_prob < conf_tresh, "low", "high"), levels = c("low", "high"))) %>% group_by(chr) %>% mutate(min_prob = min(alteration_prob)) %>% ungroup() %>% mutate(colour_chr = factor(ifelse(min_prob < conf_tresh, "low", "high"), levels = c("low", "high"), ordered = TRUE)) %>%
    select(-min_prob)
  return(data_col)
}

### Generate report for sample ###
#' Title gen_kar_list
#'
#' @param feat_tab_alt
#' @param sample_name
#' @param gender
#'
#' @return
#' @export
#'
#' @examples
gen_kar_list <- function(feat_tab_alt, sample_name, gender) {

  ##extract only the chromosomes which have only one call and have high quality
  tab_mod <- feat_tab_alt %>% select(chr, arm, alteration, colour_arm) %>% group_by(chr) %>% mutate(alt_types = length(unique(alteration)), qual_types = length(unique(colour_arm))) %>% ungroup()

  # whole chromosome changes
  alt_chr_whole <- tab_mod %>% filter(alt_types == 1 & qual_types == 1 & (colour_arm != "high" | alteration != 0)) %>%
    mutate(alt_sign = ifelse(alteration %in% c(1, 2), "+", ifelse(alteration == -1, "-", "")), qual_sign = ifelse(colour_arm == "high", "", "?"), alt_str = paste0(qual_sign, chr, alt_sign)) %>%
    distinct(chr, alt_str, alteration)
  double_chr_whole <- alt_chr_whole %>% filter(alteration == 2)
  alt_chr_whole_fin <- alt_chr_whole %>% bind_rows(double_chr_whole) %>% arrange(chr)

  # arm only changes
  alt_arm <- tab_mod %>% filter((alt_types == 2 | qual_types == 2) & (colour_arm != "high" | alteration != 0)) %>%
    mutate(alt_sign = ifelse(alteration %in% c(1, 2), "+", ifelse(alteration == -1, "-", "")), qual_sign = ifelse(colour_arm == "high", "", "?"), alt_str = paste0(qual_sign, chr, arm, alt_sign)) %>%
    distinct(chr, alt_str, alteration, arm)
  double_arm <- alt_arm %>% filter(alteration == 2)
  alt_arm_fin <- alt_arm %>% bind_rows(double_arm) %>% arrange(chr, arm) %>% select(-arm)

  #fuse the tables
  alterations <- bind_rows(alt_chr_whole_fin, alt_arm_fin) %>% arrange(chr) %>% pull(alt_str) %>% paste0(collapse = ", ")

  #chromosome number
  chrom_diff <- tab_mod %>% filter(alt_types == 1, qual_types == 1, colour_arm == "high") %>% distinct(chr, alteration) %>% pull(alteration) %>% as.numeric(.) %>% sum(.)
  chrom_n = 46 + chrom_diff

  #fill in empty alteration vector
  if (alterations == "") {
    alterations <- "none"
  }


  kar_table <- data.frame(sample = sample_name,
                          gender = factor(gender, levels = c("female", "male")),
                          chrom_n = chrom_n,
                          alterations = alterations, stringsAsFactors = FALSE)

  return(kar_table)
}

####adjust for diploid level based on diploid chromosomes####
#' Title adjust_dipl
#'
#' @param feat_tab_alt
#' @param count_ns
#'
#' @return
#' @export
#'
#' @examples
adjust_dipl <- function(feat_tab_alt, count_ns) {

  if (sum(feat_tab_alt$chr_status == "dipl", na.rm = TRUE) < 1) {
    print("Unable to adjust expression level")
    return(count_ns)
  }

  baseline_shift = median(feat_tab_alt$arm_med[feat_tab_alt$chr_status == "dipl"], na.rm = TRUE)
  count_ns$count_nor_med <- count_ns$count_nor_med - baseline_shift
  return(count_ns)
}

####Calculate weighted boxplot values####
#' Title get_box_wdt
#'
#' @param count_ns
#' @param scaleCols
#'
#' @return
#' @export
#'
#' @examples
get_box_wdt <- function(count_ns, scaleCols) {
  box_wdt <- count_ns %>% filter(chr %in% c(1:22, "X"))  %>% group_by(chr) %>% mutate(med_weig = weighted_quantile(x = count_nor_med, w = weight, probs = 0.5, na.rm = TRUE), low = weighted_quantile(x = count_nor_med, w = weight, probs = 0.25),
                                                                                      high = weighted_quantile(x = count_nor_med, w = weight, probs = 0.75), IQR = abs(high - low), max = high + IQR*1.5, min = low - IQR*1.5) %>%
    distinct(chr, .keep_all = TRUE)

  colours <- c()
  for(i in 1:nrow(box_wdt)) {
    colours[i] <- scaleCols$colour[which.min(abs(box_wdt$med_weig[i] - scaleCols$med_expr))]
  }

  box_wdt <- box_wdt %>% ungroup() %>% mutate(medianCol = colours) %>% select(chr, med_weig, low, high, min, max, medianCol) %>%
    mutate(chr = factor(x = chr, levels = c(1:22, "X")), pos = 0.5)

  return(box_wdt)
}

#### change ylim accordingly values in graph
#' Title adjust_ylim
#'
#' @param box_wdt
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
adjust_ylim <- function(box_wdt, ylim) {
  box_wdt_noy <- filter(box_wdt, !chr %in% "Y")
  if (any(box_wdt_noy$max > ylim[2])) {
    ylim[2] <- max(box_wdt_noy$max)*1.1
  }

  if (any(box_wdt_noy$min < ylim[1])) {
    ylim[1] <- min(box_wdt_noy$min)*1.1
  }
  return(ylim)
}

####Apply limit to datapoints and normalize point position####
#' Title prep_expr
#'
#' @param count_ns
#' @param dpRatioChrEdge
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
prep_expr <- function(count_ns, dpRatioChrEdge, ylim) {
  count_ns_final= count_ns %>% select(chr, end, count_nor_med, weight) %>%
    bind_rows(dpRatioChrEdge) %>% filter(chr %in% c(1:22, "X"), between(count_nor_med, ylim[1], ylim[2]) ) %>%
    mutate(chr=factor(chr, levels = c(1:22, "X"))) %>%
    arrange(chr, end) %>% group_by(chr) %>%
    mutate(normPos=scales::rescale(end, from = range(end, na.rm = TRUE)))
  return(count_ns_final)
}

####plot expression boxplot and point plot####
#' Title plot_exp
#'
#' @param count_ns_final
#' @param box_wdt
#' @param sample_name
#' @param ylim
#' @param estimate
#' @param feat_tab_alt
#' @param gender
#'
#' @return
#' @export
#'
#' @examples
plot_exp <- function(count_ns_final, box_wdt, sample_name, ylim, estimate, feat_tab_alt, gender,size1=3,size2=3) {
  # ggplot() + ylim(ylim) + ylab("log2 fold change of expression") + scale_fill_identity()+
  #   geom_point(data = count_ns_final, aes(x = normPos, y = count_nor_med),size = size1,alpha = 0.32, show.legend = FALSE)+
  #   #scale_size(range = c(2, 6)) +
  #   #scale_alpha(range = c(0.22, 0.4)) +
  #   geom_boxplot(data = box_wdt,
  #                aes(ymin = min, lower = low, middle = med_weig, upper = high, ymax = max, fill=medianCol, x = pos),
  #                alpha=0.75, outlier.colour = NA, stat = "identity", show.legend = FALSE) +
  #   geom_hline(yintercept = 0, colour = "red") +
  #
  #   labs(title = paste0(sample_name),
  #        subtitle = paste0("estimated gender: ", gender)) +
  #   geom_point(data = data.frame(x = c(0.5, 0.5),
  #                                y = c(ylim[2], ylim[2]),
  #                                point_col = c("low", "high"),
  #                                chr = factor(c(1, 1), levels = c(1:22, "X"))),
  #              mapping = aes(x = x, y = y, color = point_col), shape = 0, size = size2, stroke = 2) +
  #
  #   geom_label(data = distinct(feat_tab_alt, chr, colour_chr, chr_alt),
  #              aes(x = 0.5, y = ylim[2], color = colour_chr, label = chr_alt),
  #              label.size = 2,
  #              show.legend = FALSE) +
  #   scale_color_manual(limits = c("low", "high"), values=c("orangered", "black")) +
  #
  #   guides(color = guide_legend(title = "Quality"))+
  #   facet_grid(.~chr) +
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5, size = 18),
  #         plot.subtitle = element_text(hjust = 0.5, size = 12),
  #         axis.title.x=element_blank(),
  #         axis.text.x = element_blank(),
  #         axis.title.y = element_text(size = 15),
  #         axis.text.y = element_text(size = 12),
  #         axis.ticks = element_blank(),
  #         legend.justification = "top",
  #         legend.text = element_text(size = 15),
  #         legend.title = element_text(size = 17)) +
  #   guides(size = "none")

  gp_expr <- ggplot() + ylim(ylim) + ylab("log2 fold change of expression") +
    scale_fill_identity()+
    geom_point(data = count_ns_final, aes(x = normPos, y = count_nor_med), size = size1,alpha = 0.32, show.legend = FALSE)+
    #scale_size(range = c(2, 6)) +
    #scale_alpha(range = c(0.22, 0.4)) +
    geom_boxplot(data = box_wdt,
                 aes(ymin = min, lower = low, middle = med_weig, upper = high, ymax = max, fill=medianCol, x = pos),
                 alpha=0.75, outlier.colour = NA, stat = "identity", show.legend = FALSE) +
    geom_hline(yintercept = 0, colour = "red") +

    labs(title = paste0(sample_name),
         subtitle = paste0("estimated gender: ", gender))

  if (estimate == TRUE) {
    gp_expr <- gp_expr +

      geom_point(data = data.frame(x = c(0.5, 0.5),
                                   y = c(ylim[2], ylim[2]),
                                   point_col = c("low", "high"),
                                   chr = factor(c(1, 1), levels = c(1:22, "X"))),
                 mapping = aes(x = x, y = y, color = point_col), shape = 0, size = size2, stroke = 2) +

      geom_label(data = distinct(feat_tab_alt, chr, colour_chr, chr_alt),
                 aes(x = 0.5, y = ylim[2], color = colour_chr, label = chr_alt),
                 label.size = 2,
                 show.legend = FALSE) +
      scale_color_manual(limits = c("low", "high"), values=c("orangered", "black")) +

      guides(color = guide_legend(title = "Quality"))
  }

  gp_expr <- gp_expr +
    facet_grid(.~chr) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 12),
          axis.ticks = element_blank(),
          legend.justification = "top",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    guides(size = "none")
  gp_expr
}

####plot snv density plots####
#' Title plot_snv
#'
#' @param smpSNPdata
#' @param sample_name
#' @param estimate
#'
#' @return
#' @export
#'
#' @examples
plot_snv <- function(smpSNPdata, sample_name, estimate) {

  missedChr=c(1:22, "X")[table(smpSNPdata$chr) < 15]
  if(length(missedChr) > 0){
    tmpSNPdata=data.frame(sampleID = sample_name, ID=paste0(missedChr, "-1"), maf=0.5, chr=factor(missedChr, levels = c(1:22, "X")), start=1,
                          depth=100, snvOrd=1, snvNum=1, peak_max=0, peak=0, peakCol="red", stringsAsFactors = F)
    smpSNPdata = bind_rows(smpSNPdata, tmpSNPdata)
  }
  if(nrow(smpSNPdata)<500){
    gp.maf=ggplot()+annotate("text", x = 1, y = 1, label = "Low coverage")+theme_void()
  }else{
    snvNumDensityMaxY=smpSNPdata %>% select(chr, snvNum, peak_max, peakdist) %>% unique()
    yAxisMax=snvNumDensityMaxY %>% filter(snvNum > 100) %>% .$peak_max %>% max()
    snvNumDF = snvNumDensityMaxY %>% mutate(x=0.5, y=yAxisMax*1.05)
    peakdist_dat = snvNumDensityMaxY %>% mutate(x = 0.5, y = yAxisMax*1.1, label = round(peakdist, 3))
    gp.maf=ggplot(data=smpSNPdata) + xlab("Minor allele frequency") + ylab("Density") +
      geom_density(aes(maf, color=peakCol), show.legend = FALSE) +
      geom_vline(xintercept = c(1/3, 0.5, 2/3), alpha = 0.4, size = 0.5)+
      geom_label(data = peakdist_dat, aes(x, y, label = label), fill = "white", color = "black", vjust=0)+
      scale_color_identity(guide = guide_legend(override.aes = list(color = "white")))+
      scale_x_continuous(breaks = round(c(1/3, 2/3), 3), labels = c("1/3", "2/3"), minor_breaks = NULL, limits = c(0,1)) +
      scale_y_continuous(breaks = c(seq(from = 1, to = floor(yAxisMax)), yAxisMax*1.15), labels = c(seq(from = 1, to = floor(yAxisMax)), "peak dist."), limits = c(0, yAxisMax*1.25)) +
      facet_grid(.~chr, scales="free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 40, vjust = 0.5),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.title.x = element_text(size = 17),
            axis.title.y = element_text(size = 15),
            axis.text.y = element_text(size = 12),
            plot.margin = unit(c(0,1,1,1), "lines")
      )
    if (estimate == TRUE) {
      #need to add identical legend as for the expression in order for the graphs to align correctly with ggarrange
      gp.maf <- gp.maf +
        geom_point(data = data.frame(x = c(0.5, 0.5), y = c(0, 0), point_col = c("white", "white"), chr = factor(c(1, 1), levels = c(1:22, "X"))), mapping = aes(x = x, y = y, color = point_col), shape = 20, size = 5, alpha = 0) +
        guides(color = guide_legend(
          title = "Quality"
        )) +
        theme(
          legend.text = element_text(color = "white", size = 15),
          legend.title = element_text(color = "white", size = 17),
          legend.key = element_blank(),
          legend.justification = "top"
        )
    }
  }
  return(gp.maf)
}


#' vcf_to_snv
#'
#' @param vcf_file
#' @param maf_tresh
#' @param depth_tresh
#'
#' @return
#' @export
#'
#' @examples
vcf_to_snv=function(vcf_file, maf_tresh = 0.01, depth_tresh = 5, outfmt=1){
  #read the vcf files
  message("Reading in vcf file..")
  vcf_data=as.data.frame(data.table::fread(vcf_file, skip = "#CHROM", header = TRUE))

  #check whether the input is in correct format
  if (dim(vcf_data)[2] < 10) {
    vcf_final <- "Incorrect vcf file format. Incorrect number of columns"
    return(vcf_final)
  }
  if (!identical(colnames(vcf_data)[1:9], c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))) {
    vcf_final <- "Incorrect vcf file format."
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, 'INFO'], "DP=[0-9]+") != TRUE &
      str_detect(string = str_split(vcf_data[1, 'FORMAT'], pattern = ":")[[1]][3], pattern = "DP") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No depth information in the INFO column and missing information in the FORMAT column"
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, 9], "AD") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No allele depth (AD) in FORMAT column"
    return(vcf_final)
  }

  vcf_data <- vcf_data[, 1:10]
  names(vcf_data)=c("chr", "start", "ID", "ref", "var", "qual", "FILTER", "INFO", "FORMAT", "INFO2")

  # Getting depth out of the INFO column
  message("Extracting depth..")

  if (str_detect(vcf_data[1, 'INFO'], "DP=[0-9]+") == TRUE) {
    vcf_data$depth = as.numeric(gsub("DP=", "", str_extract(vcf_data$INFO, "DP=[0-9]+")))
  } else {
    DP_i=grep("DP",unlist(strsplit(vcf_data$FORMAT[1],split = ":")))
    vcf_data$depth=sapply(str_split(vcf_data$INFO2,":"),function(x){x[DP_i]})
  }

  #reference allele and alternative allele depths
  message("Extracting reference allele and alternative allele depths..")
  vcf_data$REF_ALT_num = sapply(str_split(vcf_data$INFO2, ":"), function(x) x[2])

  #extract count for alternative allele
  vcf_data$varDp = as.numeric(sapply(str_split(vcf_data$REF_ALT_num, ","), function(x) x[2]))

  #mutant allele frequency
  vcf_data$maf = vcf_data$varDp/vcf_data$depth
  message("Needed information from vcf extracted")

  #return needed columns
  if(outfmt==1){
    vcf_final <- vcf_data[,c("chr", "start", 'depth', "maf")]
  } else {
    vcf_final=vcf_data[,c("chr", "start","ref" ,'var','depth', "maf")]
  }

  message("Finished reading vcf")
  return(vcf_final)

}

#' get_norm_exp
#'
#' @param count_table
#' @param standard_samples
#' @param minReadCnt
#' @param samp_prop
#' @param weight_table
#' @param weight_samp_prop
#' @param generate_weights
#'
#' @return
#' @export
#'
#' @examples
get_norm_exp <- function(count_table, standard_samples, minReadCnt, samp_prop, weight_table, weight_samp_prop, generate_weights) {

  sample_name="TestSample"
  #rename count table
  names(count_table)=c("ENSG", "count")


  #check the count file format
  if (ncol(count_table) != 2 | !is.numeric(count_table[, 2])) {
    return(paste0("Incorrect count file format for the sample: ", sample_name))
  }

  # set column to join by
  # data.table::setkey(count_table, ENSG)

  #inner join diploid reference and analyzed sample
  # data.table::setkey(standard_samples, ENSG)
  # final_mat <- as.data.frame(count_table[standard_samples, nomatch = 0])
  final_mat=count_table %>% left_join(standard_samples)

  #keep genes for determining gender for later
  gender_genes = final_mat %>% filter(ENSG %in% "ENSG00000012817")

  #filter genes based on reads count; top 1-q have read count > N, filter base on weight
  keepIdx_tmp = final_mat %>%
    mutate(keep_gene = apply(.[, -1], MARGIN = 1, FUN = function(x) sum(x > minReadCnt) > (length(x) * samp_prop)), id = row_number()) %>%
    filter(keep_gene == TRUE)

  if (generate_weights == FALSE) {
    keepIdx = keepIdx_tmp %>% inner_join(weight_table, by = "ENSG") %>%
      group_by(chromosome_name) %>% mutate(weight_chr_quant = quantile(weight, 1 - weight_samp_prop)) %>%
      filter(weight >= weight_chr_quant) %>% pull(id)
  } else {
    keepIdx = keepIdx_tmp %>% pull(id)
  }

  # filter table for normalization, get rid of genes with 0 counts and keep genes for gender estimation
  count_filt <- final_mat %>% .[c(keepIdx), ] %>% .[pull(., 2) != 0, ] %>% bind_rows(gender_genes) %>% distinct(ENSG, .keep_all = TRUE)
  ENSG <- count_filt$ENSG
  count_filt <- select(count_filt, -ENSG)

  #sample Deseq normalization
  count_col <- as.data.frame(colnames(count_filt))

  dds <- DESeq2::DESeqDataSetFromMatrix(colData = count_col, countData = count_filt, design= ~ 1)
  dds_vst <- DESeq2::varianceStabilizingTransformation(dds, blind=T, fitType='local')
  count_norm <- as.data.frame(SummarizedExperiment::assay(dds_vst))


  colnames(count_norm)[1] <- sample_name
  colnames(count_norm)[2:ncol(count_norm)] <- paste0("control_", 1:(ncol(count_norm)-1))
  print(paste0("Normalization for sample: ", sample_name, " completed"))

  #Modify table for downstream analysis
  rownames(count_norm) <- ENSG

  return(count_norm)
}


####get median expression level
#' Title get_med
#'
#' @param count_norm
#' @param refDataExp
#' @param weight_table
#' @param generate_weights
#'
#' @return
#' @export
#'
#' @examples
get_med <- function(count_norm, refDataExp, weight_table, generate_weights) {

  ENSG=rownames(count_norm)

  ####calculate median for all genes
  pickGeneDFall_tmp=count_norm %>% mutate(ENSG=ENSG) %>% left_join(select(refDataExp, chr, ENSG), by = "ENSG")
  if (generate_weights == TRUE) {
    pickGeneDFall <- pickGeneDFall_tmp %>% mutate(med = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median), var = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, var)) %>%
      select(ENSG, chr, med, var)
  } else {
    pickGeneDFall <- pickGeneDFall_tmp %>% mutate(med = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median)) %>%
      select(ENSG, med)
  }

  return(pickGeneDFall)
}

####filter SNVs of interest for samples###
#' Title
#'
#' @param one_smpSNP
#' @param keepSNP
#' @param minDepth
#' @param mafRange
#'
#' @return
#' @export
#'
#' @examples
filter_snv <- function(one_smpSNP, keepSNP, minDepth, mafRange) {

  smpSNPdata.tmp= one_smpSNP %>% dplyr::select(sampleID, ID, maf, chr, start, depth, arm) %>%
    filter(data.table::inrange(maf, mafRange[1], mafRange[2]), depth > minDepth) %>% filter(chr != "Y")
  if (keepSNP[1] != FALSE) {
    smpSNPdata.tmp <- smpSNPdata.tmp %>% filter(ID %in% keepSNP)
  }
  return(smpSNPdata.tmp)
}

#' prepare_snv
#'
#' @param snv_file
#' @param centr_ref
#' @param minDepth
#' @param snv_format
#'
#' @return
#' @export
#'
#' @examples
prepare_snv <- function(snv_file, centr_ref, minDepth, snv_format = "vcf") {

  sample_name="TestSample"
  smpSNP=list()
  print(paste("Preparing file with snv information for:", sample_name))

  #prepare data from custom table
  if (snv_format == "custom") {
    #read the table
    snv_table_pre <- fread(snv_file)

    #check if all appropriate columns are present
    cols <- colnames(snv_table_pre)
    chr <- str_which(cols, "^#Chromosome$|#CHROM$|^CHR$|^chr$|^Chr$|^CHROM$|^chrom$|^Chrom$|^CHROMOSOME$|^chromosome$|^Chromosome$")
    start <- str_which(cols, "^START$|^start$|^Start$|^POS$|^pos$|^Pos$")
    depth <- str_which(cols, "^DEPTH$|^depth$|^Depth$|^DP$|^dp$|^Dp$")
    maf <- str_which(cols, "^MAF$|^maf$|^Maf$|^HET$|^het$|^Het$")
    to_keep <- as.numeric(c(chr, start, depth, maf))
    if (length(to_keep) != 4) {
      smpSNP[[sample_name]] <- "Incorrect column name in a custom snv file."
      return(smpSNP)
    }
    snv_table <- snv_table_pre[, to_keep, with = FALSE]
    data.table::setnames(snv_table, colnames(snv_table), c("chr", "start", "depth", "maf"))

    #Check some column parameters
    if (is.numeric(snv_table[, start]) == FALSE | is.numeric(snv_table[, depth]) == FALSE | is.numeric(snv_table[, maf]) == FALSE) {
      smpSNP[[sample_name]] <- "Incorrect type of a column in a custom file with snv information."
      return(smpSNP)
    }

  } else if (snv_format == "vcf") {
    snv_table <- vcf_to_snv(snv_file)
    if(is.character(snv_table)) {
      smpSNP[[sample_name]] <- snv_table
      return(smpSNP)
    }
  }

  smpSNP[[sample_name]] <- snv_table %>%
    filter(chr %in% c(c(1:22, "X"), paste0("chr", c(1:22, "X")))) %>%
    mutate(chr = sub("chr", "", chr)) %>% left_join(centr_ref, by = "chr") %>%
    mutate(chr = factor(chr, levels=c(1:22, "X")), ID=paste0(chr,"-", start), sampleID=sample_name) %>%
    mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr")))

  return(smpSNP)
}

# snv_file="tests/test.HaplotypeCaller.vcf"


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

  #get ref dataset
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

  #create standard samples
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

#' Weight calculation
#'
#' Get weights for CNA analysis
#'
#' The function provides an option to create customized weight matrix based on specific data (cancer) type. The weights are based on the correlation of gene expression with CNV and expression variance across the cohort. To calculate gene expression-CNV correlation
#' reliable set of CNVs on arm level is required.  For the function to work properly, each chromosome arm should have
#' minimum two CNV states in the provided CNV data. If this is not the case the parameter fill_cor should be set as TRUE to artificially fill out the correlation coefficients.
#'
#' @param CNV_data path to a table with known CNVs on arm-level. The information is used to generate CNV-gene expression level correlation coeffiecient, that are
#' subsequently used for gene weight calculation. The more samples are provided, the more robust the calculated correlation coeffiecient will be.
#' There are two optional input formats. Either in wide format, where each column represents one chromosome arm e.g. columns (1p,1q,2p,2q,3p,3q..) and
#' one additional column with sample name. The values in chromosome arm column are the CNV states on that chromosome arm, if the CNV state is unknown the value should be set to NA. Also the CNV_data file
#' can be in "long" format. Here the table is required to have three columns: sample: sample name, arm: chromosome arm (format: "1p", "19q", "Xp") and CNV: (CNV state, e.g. +1, -1, +2, 0)
#' @param standard_samp character vector with standard sample names, which will be used for normalization, and adjustment of diploid level. These samples should not
#' carry any CNVs. Providing standard samples is important for datasets, where high numbers of CNVs (e.g. 10 gained chromosomes or 10 deleted chromosomes) of one type are present is some samples.
#'  The sample names must be included in the first column of metadata table. If no standard samples samples are input, the data will be normalized as a batch and no diploid level adjustment will
#'  be performed.
#' @param metadata path to a metadata table with three columns. First colum: sample names, second column: file names of count files, third column: file names of snv files. There should be no header.
#'  More information is included in the package README file.
#' @param CNV_data_format string, in which format is the CNV_data provided, can be either "wide" or "long"
#' @param fill_cor if TRUE (default), genes for which correlation coefficients cannot be estimated and and those for which adjusted p-value by Benjamini-Hochberg method is higher than 0.1, the correlation coefficient value will be set
#' to 0.1. This was determined as the mean gene expression-CNV correlation coefficient from qcute lymphoblastic leukemia data on all genes. This parameter may help in creating weights for datasets, for which the coverage by CNVs
#' of genome is low and thus correlation coeffiecients may be constructed only for small proportion of genes. If FALSE, genes without determined correlation coefficients and those for which adjusted p-value is higher than 0.1 will be excluded.
#' @param output_file path to a file, where gene weight matrix will be saved.
get_weights <- function(CNV_data, count_dir, standard_samp = NULL, metadata, CNV_data_format, fill_cor = TRUE) {

  if(!CNV_data_format %in% c("long", "wide")) {
    stop("CNV_data_format has to be set to either long or wide format")
  }


  CNV_table <- fread(CNV_data)

  if(CNV_data_format == "wide") {
    CNV_table <- pivot_longer(data = CNV_table, names_to = "arm", values_to = "CNV", cols = c(paste0(c(1:22, "X"), "q"), paste0(c(1:12, 16:20, "X"), "p")))
  }


  #function for geometric
  gm_mean = function(a){prod(a)^(1/length(a))}

  #read metadata table
  metadata <- fread(metadata, header = FALSE)

  #get file paths to count files
  metadata$count_path <- file.path(count_dir, pull(metadata, 2))

  if (any(duplicated(metadata[, 2]))) {
    stop("Duplicated samples")
  }

  if (!is.null(standard_samples)) {

    #Check whether the samples are present in the sample table
    if (all(standard_samp %in% pull(metadata, 1)) == FALSE) {
      stop("The input standard samples are not in metadata table.")
    }
    #Create standard sample table
    standard_table <- create_standard(standard_samples = standard_samp, sample_table = metadata)
  }

  if (!is.null(standard_samp)) {
    for (i in 1:nrow(metadata)) {
      if (i == 1) {
        count_norm <- get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = FALSE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% .[, !duplicated(colnames(.))] %>% select(1) %>%
          mutate(ENSG = rownames(.))
      } else {
        count_norm <- inner_join(count_norm, get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = FALSE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% .[, !duplicated(colnames(.))] %>%
                                   select(1) %>%  mutate(ENSG = rownames(.)))
      }
    }

  } else {
    count_norm <- get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = TRUE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% mutate(ENSG = row.names(.))
  }

  #annotate the normalized counts with chromosome and chromosome arm
  count_norm_an <-  count_norm %>% left_join(as.data.frame(refDataExp), by = "ENSG") %>% left_join(mutate(centr_ref, chr = as.character(chr)), by = "chr") %>% mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr"))) %>%
    mutate(arm = paste0(chr, arm)) %>% select(-start, -end, -cstart, -cend, -chr)

  if (!is.null(standard_samp)) {

    #calculate pseudo reference in diploid samples for second normalization
    count_dipl_norm <-  count_norm_an[, which(colnames(count_norm_an) %in% c(standard_samp, "arm", "ENSG"))]
    pseudo_ref_norm <- apply(X = select(count_dipl_norm, -arm, - ENSG), MARGIN = 1, function(x) gm_mean(x))
    ref_dipl = cbind(select(count_dipl_norm, arm, ENSG), pseudo_ref = as.numeric(pseudo_ref_norm)) %>% filter(pseudo_ref != 0)

    #second normalization, especially for samples with high numbers of CNA
    count_norm_an_cor <- count_norm_an
    non_dipl <- pull(metadata[!V1 %in% standard_samp], V1)
    for (i in non_dipl) {
      sample_vst <- count_norm_an[, c(non_dipl, "ENSG", "arm")]
      colnames(sample_vst)[1] <- "count"
      size_fac = sample_vst %>% right_join(ref_dipl, by = c("ENSG", "arm")) %>% inner_join(filter(CNV_table, sample == i), by = c("arm")) %>% filter(CNV == 0) %>% mutate(size_facs = count/pseudo_ref) %>% summarise(median(size_facs))
      count_norm_an_cor[, i] <- count_norm_an_cor[, i] / size_fac
    }
    count_norm_an <- count_norm_an_cor
  }

  #annotate the data with alteration information
  count_final <- gather(count_norm_an, key = "sample", value = "vst", which(!colnames(count_norm_an) %in% c("arm", "ENSG"))) %>% inner_join(select(CNV_table, sample, arm, CNV), by = c("sample", "arm"))

  #get per-gene CNA expression correlation
  genes <- unique(count_final$ENSG)

  for (i in 1:length(genes)) {
    gene_data <- count_final %>% filter(ENSG == genes[i])
    suppressWarnings(cor <- cor.test(x = gene_data$vst, y = as.numeric(gene_data$CNV)))

    if (i != 1) {
      gene_cor <- rbind(gene_cor, data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value))
    } else {
      gene_cor <- data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value, stringsAsFactors = FALSE)
    }
  }

  # calculate adjusted p values, filter reliable genes and optionally fill out the rest of genes with base correleation coefficient of 0.1
  gene_cor_padj <- gene_cor %>% filter(!is.na(pearson_r) & !is.na(p)) %>% mutate(p_adj = p.adjust(p, method = "hochberg"))
  gene_cor_filt <- gene_cor_padj %>% filter(p_adj < 0.1)

  if (nrow(gene_cor_filt) < 1000 & fill_cor == FALSE) {
    stop("There is not enough data to generate robust correlation coefficients for most of the genes. Either expand the data or set fill_cor to TRUE to artificially assign a baseline correlation coefficient.")
  }

  if (fill_cor == TRUE) {
    gene_fill <- gene_cor %>% filter(is.na(pearson_r) | is.na(p)) %>% mutate(cor = 0.1)
    gene_fill <- gene_fill %>% bind_rows(gene_cor_padj %>% filter(p_adj  > 0.1) %>% mutate(cor = 0.1) %>% select(-p_adj))
  }

  #get per-gene variance and depth
  var <- apply(select(count_dipl_norm, -arm, -ENSG), 1, var)
  depth <- apply(select(count_dipl_norm, -arm, -ENSG), 1, mean)
  var_depth_table <- data.frame(ENSG = count_dipl_norm$ENSG, var = var, depth = depth, arm = count_dipl_norm$arm, stringsAsFactors = FALSE)

  #calculate weights
  weight_table <- var_depth_table %>% right_join(gene_cor, by = "ENSG") %>% mutate(weight = scales::rescale(pearson_r^5, to = c(1,100))*scales::rescale(1/var, to = c(1,100))) %>% mutate(chromosome_name = sub("p|q", "", arm)) %>%
    arrange(desc(weight)) %>% dplyr::select(ENSG, weight, chromosome_name)

  write.table(weight_table, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}


#' Title get_RNAseqCNV_plot
#'
#' @param RNAseqCNV_out
#' @param adjust
#' @param estimate_lab
#' @param scale_cols
#'
#' @return
#' @export
#'
#' @examples
get_RNAseqCNV_plot=function(RNAseqCNV_out,
                            sample_name="TestSample",
                            adjust=TRUE,
                            estimate_lab = TRUE,
                            scale_cols = scaleCols,
                            dpRatioChromEdge = dpRatioChrEdge,
                            size1=3,size2=3){
  #get parameters
  gender=as.character(RNAseqCNV_out$df_cnv_out$gender)[1]
  count_ns=RNAseqCNV_out$count_ns
  feat_tab_alt=RNAseqCNV_out$feat_tab_alt
  smpSNPdata=RNAseqCNV_out$smpSNPdata

  #adjust the gene expression according to the estimation of which chromosomes are diploid
  if (adjust == TRUE) {
    count_ns <-  adjust_dipl(feat_tab_alt, count_ns)
  }

  #calculate box plots
  box_wdt <- get_box_wdt(count_ns = count_ns, scaleCols = scale_cols)

  #adjust y axis limits
  ylim <- adjust_ylim(box_wdt = box_wdt, ylim = c(-0.4, 0.4))

  count_ns_final <- prep_expr(count_ns = count_ns, dpRatioChrEdge = dpRatioChromEdge, ylim = ylim)

  # filter low weighted genes for clearer visualization
  #count_ns_final <- filter_expr(count_ns_final = count_ns_final, cutoff = 0.6)


  gg_exp <- plot_exp(count_ns_final = count_ns_final, box_wdt = box_wdt, sample_name = sample_name, ylim = ylim,
                     estimate = estimate_lab, feat_tab_alt = feat_tab_alt, gender = gender,size1 = size1,size2 = size2)

  gg_snv <- plot_snv(smpSNPdata, sample_name = sample_name, estimate = estimate_lab)

  # fig <- arrange_plots(gg_exp = gg_exp, gg_snv = gg_snv)

  fig = plot_grid(gg_exp,gg_snv, ncol = 1,rel_heights = c(0.75,0.25),align = "v")

  fig
}


