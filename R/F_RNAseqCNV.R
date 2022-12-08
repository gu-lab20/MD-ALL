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
plot_exp <- function(count_ns_final, box_wdt, sample_name, ylim, estimate, feat_tab_alt, gender) {
  gp_expr <- ggplot() + ylim(ylim) + ylab("log2 fold change of expression") +
    scale_fill_identity()+
    geom_point(data = count_ns_final, aes(x = normPos, y = count_nor_med, size = 0.2), alpha = 0.32, show.legend = FALSE)+
    #scale_size(range = c(2, 6)) +
    #scale_alpha(range = c(0.22, 0.4)) +
    geom_boxplot(data = box_wdt, aes(ymin = min, lower = low, middle = med_weig, upper = high, ymax = max, fill=medianCol, x = pos), alpha=0.75, outlier.colour = NA, stat = "identity", show.legend = FALSE)+
    geom_hline(yintercept = 0, colour = "red")+
    labs(title = paste0(sample_name),
         subtitle = paste0("estimated gender: ", gender))
  if (estimate == TRUE) {
    gp_expr <- gp_expr +
      geom_point(data = data.frame(x = c(0.5, 0.5), y = c(ylim[2], ylim[2]), point_col = c("low", "high"), chr = factor(c(1, 1), levels = c(1:22, "X"))),
                 mapping = aes(x = x, y = y, color = point_col), shape = 0, size = 4, stroke = 2) +
      geom_label(data = distinct(feat_tab_alt, chr, colour_chr, chr_alt), aes(x = 0.5, y = ylim[2], color = colour_chr, label = chr_alt), label.size = 2, show.legend = FALSE) +
      scale_color_manual(limits = c("low", "high"), values=c("orangered", "black")) +
      guides(color = guide_legend(
        title = "Quality"
      ))
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

####arrange expression and snv graph####
#' Title arrange_plots
#'
#' @param gg_exp
#' @param gg_snv
#'
#' @return
#' @export
#'
#' @examples
arrange_plots <- function(gg_exp, gg_snv) {
  fig <- ggarrange(plotlist =list(expr=gg_exp, maf=gg_snv)
                   , ncol = 1, nrow = 2, heights = c(3, 1), align='v')
  return(fig)
}





