xytitle <- 8
textsize <- 8
anglenum <- 0
vj <- 0
hj <- 0.5
legendsize <- 0.4
titlesize <- 9
theme4_backline <- theme(axis.title.x = element_text(size=xytitle,face="bold",colour="black"),
                         axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj,colour="black" ),
                         axis.title.y = element_text(size=xytitle,face="bold",colour="black"),
                         axis.text.y = element_text( size = textsize,colour="black"),
                         plot.title = element_text(size = titlesize),
                         plot.subtitle = element_text(size = textsize),
                         panel.background = element_rect(fill=NA, colour=NA),
                         legend.background = element_rect(fill=NA, colour=NA),
                         legend.title = element_text( size = textsize),
                         legend.text = element_text( size = textsize),
                         legend.key.size =  unit(0.5, "cm"),
                         panel.grid.major = element_line(color = "#EDEDED", size = 0.5),
                         panel.grid.minor = element_line(color = "#EDEDED", size = 0.35),
                         panel.border = element_rect(colour = "black", fill = NA, size = 0.35),
                         strip.background = element_blank(),
                         axis.line = element_blank( ),
                         axis.ticks = element_line(colour = "black", size = 0.25),
                         axis.ticks.length.x=unit(0.08, "cm"),
                         axis.ticks.length.y=unit(0.08, "cm")
)
theme2_nobackplus <- theme(
  axis.title.x = element_text(size=xytitle,face="bold",colour="black"),
  axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj,colour="black" ),
  axis.title.y = element_text(size=xytitle,face="bold",colour="black"),
  axis.text.y = element_text( size = textsize,colour="black"),
  plot.title = element_text(size = titlesize),
  plot.subtitle = element_text(size = textsize),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_rect(fill = NA, colour = NA),
  legend.background = element_rect(fill = NA, colour = NA),
  legend.title = element_text( size = textsize),
  legend.text = element_text( size = textsize),
  legend.key.size =  unit(0.5, "cm"),
  axis.line.x = element_line(
    colour = "darkgrey",
    size = 0.15,
    linetype = "dashed"
  ),
  axis.line.y = element_line(
    colour = "black",
    size = 0.15
  ),
  axis.ticks = element_line(colour = "black", size = 0.25),
  axis.ticks.length.x = unit(0.08, "cm"),
  axis.ticks.length.y = unit(0.08, "cm")
)
theme2_nobacknoerror <- theme(
  axis.title.x = element_text(size=xytitle,face="bold",colour="black"),
  axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj,colour="black" ),
  axis.title.y = element_text(size=xytitle,face="bold",colour="black"),
  axis.text.y = element_text( size = textsize,colour="black"),
  plot.title = element_text(size = titlesize),
  plot.subtitle = element_text(size = textsize),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_rect(fill = NA, colour = NA),
  legend.background = element_rect(fill = NA, colour = NA),
  legend.title = element_text( size = textsize),
  legend.text = element_text( size = textsize),
  legend.key.size =  unit(0.5, "cm"),
  axis.line.x = element_line(
    colour = "black",
    size = 0.25
  ),
  axis.line.y = element_line(
    colour = "black",
    size = 0.25
  ),
  axis.ticks = element_line(colour = "black", size = 0.25),
  axis.ticks.length.x = unit(0.08, "cm"),
  axis.ticks.length.y = unit(0.08, "cm")
)
theme1_noback <- theme(axis.title.x = element_text(size=xytitle,face="bold",colour="black"),
                       axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj,colour="black" ),
                       axis.title.y = element_text(size=xytitle,face="bold",colour="black"),
                       axis.text.y = element_text( size = textsize,colour="black"),
                       plot.title = element_text(size = titlesize),
                       plot.subtitle = element_text(size = textsize),
                       panel.grid.major = element_blank( ),
                       panel.grid.minor = element_blank( ),
                       panel.border = element_blank( ),
                       panel.background = element_rect(fill=NA, colour=NA),
                       legend.background = element_rect(fill=NA, colour=NA),
                       legend.title = element_text( size = textsize),
                       legend.text = element_text( size = textsize),
                       legend.key.size =  unit(0.5, "cm"),
                       axis.line.x = element_blank(),
                       axis.line.y = element_line(colour="black", size=0.25),
                       axis.ticks = element_line(colour = "black", size = 0.25),
                       axis.ticks.length.x=unit(0.08, "cm"),
                       axis.ticks.length.y=unit(0.08, "cm")
)
themelegendright <- theme(legend.position = "right",
                          legend.title = element_text(size = textsize),
                          legend.key.size = unit(legendsize, "cm"),
                          legend.text = element_text(size = textsize),
                          panel.background = element_rect(fill=NA, colour=NA))
themelegendleft <- theme(legend.position = "left",
                         legend.title = element_text(size = textsize),
                         legend.key.size = unit(legendsize, "cm"),
                         legend.text = element_text(size = textsize,angle=90),
                         panel.background = element_rect(fill=NA, colour=NA))
themelegendtop <- theme(legend.position = "top",
                        legend.title = element_text(size = textsize),
                        legend.key.size = unit(legendsize, "cm"),
                        legend.text = element_text(size = textsize),
                        panel.background = element_rect(fill=NA, colour=NA))
themelegendnone <- theme(legend.position = "none",
                         legend.title = element_text(size = textsize),
                         legend.key.size = unit(legendsize, "cm"),
                         legend.text = element_text(size = textsize),
                         panel.background = element_rect(fill=NA, colour=NA))
themetext30 <- theme(axis.text.x = element_text( size = textsize, angle = 30,
                                                 vjust = 0.9,hjust= 0.9 )
)
themetext90 <- theme(axis.text.x = element_text( size = textsize, angle = 90,hjust = 1, vjust = 0.5 )
)
themey90 <- theme( axis.text.y = element_text( size = textsize, angle = 90,vjust = 0,hjust= 0.5))
themeynone <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank())
open_xlsx <- function(file,i){
  a <- as.data.frame(read_excel(file,sheet=i))
  return(a)
}
jaccard_index <- function(x, y) {
  intersection <- length(intersect(x, y))
  union <- length(base::union(x, y))
  if (union == 0) {
    return(NA)
  } else {
    return(intersection / union)
  }
}
indi_prot <- function(data,cate,metadata){

  if(cate=="ms"){
    d1 <- melt(as.matrix(data))
    colnames(d1) <- c("X1","OriginalID","value")
    d1 <- left_join(d1,metadata,by="OriginalID")
    d1re <- d1 %>%
      dplyr::group_by(SubjectID, X1) %>%
      dplyr::summarise(
        zero_count = sum(value == 0, na.rm = TRUE),
        total_count = n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(percentage=(total_count-zero_count)/total_count,
                    number=(total_count-zero_count))
    freq <- d1re %>% dplyr::group_by(SubjectID, number) %>%
      dplyr::summarise(
        count = n(),
        .groups = "drop"
      ) %>% dplyr::group_by(SubjectID) %>%
      dplyr::mutate(percentage = count/ sum(count))%>%
      dplyr::group_by(SubjectID) %>%
      complete(number = c(0, 1, 2, 3), fill = list(count = 0, percentage = 0)) %>%
      dplyr::ungroup()
  }else {
    data <- left_join(data,metadata,by=c("SampleID"="OriginalID"))
    d1re <- data %>%
      dplyr::group_by(SubjectID,Assay) %>%
      dplyr::summarise(
        NPXlow = sum(NPX < LOD, na.rm = TRUE),
        total_count = n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(percentage=(total_count-NPXlow)/total_count,
                    number=(total_count-NPXlow))
    freq <- d1re %>% dplyr::group_by(SubjectID, number) %>%
      dplyr::summarise(
        count = n(),
        .groups = "drop"
      ) %>% dplyr::group_by(SubjectID) %>%
      dplyr::mutate(percentage = count/ sum(count))%>%
      dplyr::group_by(SubjectID) %>%
      complete(number = c(0, 1, 2, 3), fill = list(count = 0, percentage = 0)) %>%
      dplyr::ungroup()
  }
  return(freq)
}
plot_indipro <- function(freq1_2,tit){
  colorindi <- c("#1F77B4FF", "#AEC7E8FF","#FF7F0EFF", "#FFBB78FF","#2CA02CFF", "#98DF8AFF", "#D62728FF","#FF9896FF", "#9467BDFF", "#C5B0D5FF",
                 "#8C564BFF", "#C49C94FF", "#E377C2FF", "#F7B6D2FF" )
  colorindi <- paletteer_d("ggsci::nrc_npg")
  donor_breaks <- if (is.factor(freq1_2$SubjectID)) {
    levels(freq1_2$SubjectID)
  } else {
    unique(as.character(freq1_2$SubjectID))
  }
  donor_breaks <- donor_breaks[donor_breaks %in% as.character(freq1_2$SubjectID)]
  ggplot(freq1_2, aes(x = number, y = percentage,color = SubjectID)) +
    geom_xspline(spline_shape = 0.15,size=3,alpha=1)+
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(limits = c(0, 3), expand = c(0, 0), breaks = seq(0, 3, by = 1)) +
    scale_color_manual(
      values = colorindi,
      breaks = donor_breaks,
      labels = paste0("Donor", seq_along(donor_breaks)),
      name = "DonorID"
    ) +
    labs(
      title = tit,
      x = "Number of samples",
      y = "Percentage of detected proteins"
    ) + theme2_nobacknoerror+themey90 + themelegendnone
}
impute_perseus_log2 <- function(x, shift = 1.8, width = 0.3, seed = NULL,
                                min_non_na = 3) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  na_mask <- is.na(x) | !is.finite(x)
  if (!is.null(seed)) set.seed(seed)
  params <- matrix(NA_real_, nrow = ncol(x), ncol = 2,
                   dimnames = list(colnames(x), c("mu", "sd")))
  for (j in seq_len(ncol(x))) {
    col <- x[, j]
    obs <- col[is.finite(col)]
    if (length(obs) < min_non_na) next
    mu <- mean(obs, na.rm = TRUE)
    sdv <- stats::sd(obs, na.rm = TRUE)
    params[j, ] <- c(mu, sdv)
    miss_idx <- which(!is.finite(col))
    if (length(miss_idx) > 0 && is.finite(sdv) && sdv > 0) {
      imp_vals <- stats::rnorm(length(miss_idx),
                               mean = mu - shift * sdv,
                               sd   = width * sdv)
      col[miss_idx] <- imp_vals
      x[, j] <- col
    }
  }
  attr(x, "perseus_params_by_col") <- params
  attr(x, "imputation") <- list(
    method = "Perseus_downshift_normal_on_log2",
    shift = shift, width = width,
    min_non_na = min_non_na
  )
  attr(x, "imputed_mask") <- na_mask
  return(x)
}
explore_cv <- function(npx, na.rm = TRUE) {
  sqrt(exp((log(2) * sd(npx, na.rm = na.rm))^2) - 1)
}
cv_detection_keep <- function(prop, isna, keepfilter = TRUE) {
  use_detection_filter <- (isTRUE(isna) || identical(isna, "NA")) && isTRUE(keepfilter)
  if (use_detection_filter) !is.na(prop) & prop >= 0.6 else rep(TRUE, length(prop))
}
olinkcv <- function(cate,pro,isna=TRUE,cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "3072"){
    prot <- pro
    o3072cv <- lod3072 %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- o3072cv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
      ) %>% dplyr::group_by(Assay, SubjectID) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    o3072cv_kept <- o3072cv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay, SubjectID),by = c("Assay", "SubjectID"))
    if(isna == TRUE){
      o3072cv_kept <- o3072cv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      o3072cv_kept <- o3072cv_kept
       }
    o3072cv <- o3072cv_kept %>%  dplyr::group_by(SubjectID,Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-Explore", cate1=cate1input,.groups = "drop") %>%
      dplyr::filter(!is.na(CV))
    colnames(o3072cv)  <- c("SubjectID", "protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,o3072cv)
    cvlow <- Inf
    o3072cv1 <-  o3072cv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(SubjectID, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-Explore") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(o3072cv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,o3072cv1)
  }else{
    prot <- pro
    ohtcv <- lodht %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- ohtcv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
    ) %>%
      dplyr::group_by(Assay, SubjectID) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    ohtcv_kept <- ohtcv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay, SubjectID),by = c("Assay", "SubjectID"))
    if(isna == TRUE){
      ohtcv_kept <- ohtcv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      ohtcv_kept <- ohtcv_kept
    }
    ohtcv <- ohtcv_kept %>%  dplyr::group_by(SubjectID,Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-HT",cate1=cate1input ,.groups = "drop") %>%
      dplyr::filter(!is.na(CV))
    colnames(ohtcv) <- c("SubjectID", "protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,ohtcv)
    cvlow <- Inf
    ohtcv1 <-  ohtcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(SubjectID, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-HT") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(ohtcv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,ohtcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
mscv <- function(dflistna,cate,pro,isna="0",cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "ot480"){
    prot <- pro
    otcv <- as.data.frame(dflistna[["MS-OT480"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- otcv %>%
      dplyr::group_by(protein, SubjectID) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
    if (isna == "0") {
      otcv_kept <- otcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein, SubjectID),by = c("protein", "SubjectID"))
      otcv_kept <- otcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else if (isna == "NA") {
      otcv_kept <- otcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein, SubjectID),by = c("protein", "SubjectID"))
    }else if (isna == "impute"){
      df <- as.data.frame(dflistna[["MS-OT480"]])
      df <- as.matrix(df[pro,])
      colmap <- data.frame(
        col = colnames(df),
        replicate = as.integer(str_replace(colnames(df), "_.*", "")),
        SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
        stringsAsFactors = FALSE)
      subjects <- unique(colmap$SubjectID)
      res_list <- list()
      for (s in subjects) {
        idx <- which(colmap$SubjectID == s)
        subX <- df[, idx, drop = FALSE]
        thr <- ceiling(ncol(subX) * 0.6)
        keep <- rowSums(!is.na(subX)) >= thr
        x_filtered <- subX[keep, , drop = FALSE]
        x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
        df_imp <- as.data.frame(x_imp, check.names = FALSE)
        df_imp$protein <- rownames(df_imp)
        long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
        res_list[[s]] <- long
      }
      res_all <- bind_rows(res_list)
      res_all <- left_join(res_all,metadata,by=c("variable"="OriginalID"))
      otcv_kept <- res_all
    }
    otcv <- otcv_kept %>% dplyr::group_by(SubjectID,protein) %>% dplyr::summarise(CV = explore_cv(value),cate="MS-OT480", cate1=cate1input, .groups = "drop") %>%
      dplyr::filter(!is.nan(CV))
    bindcv2 <- rbind(bindcv2,otcv)
    cvlow <- Inf
    otcv1 <-  otcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(SubjectID, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-OT480") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(otcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,otcv1)
  }else{
    prot <- pro
    astcv <- as.data.frame(dflistna[["MS-Astral"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- astcv %>%
      dplyr::group_by(protein, SubjectID) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
    if (isna == "0") {
      astcv_kept <- astcv %>% dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(protein, SubjectID),
        by = c("protein", "SubjectID"))
      astcv_kept <- astcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else if (isna == "NA") {
          astcv_kept <- astcv %>%
            dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein, SubjectID),by = c("protein", "SubjectID"))
    }else if (isna == "impute"){
          df <- as.data.frame(dflistna[["MS-Astral"]])
          df <- as.matrix(df[pro,])
          colmap <- data.frame(
            col = colnames(df),
            replicate = as.integer(str_replace(colnames(df), "_.*", "")),
            SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
            stringsAsFactors = FALSE)
          subjects <- unique(colmap$SubjectID)
          res_list <- list()
          for (s in subjects) {
            idx <- which(colmap$SubjectID == s)
            subX <- df[, idx, drop = FALSE]
            thr <- ceiling(ncol(subX) * 0.6)
            keep <- rowSums(!is.na(subX)) >= thr
            x_filtered <- subX[keep, , drop = FALSE]
            x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
            df_imp <- as.data.frame(x_imp, check.names = FALSE)
            df_imp$protein <- rownames(df_imp)
            long <- data.table::melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
            res_list[[s]] <- long
          }
          res_all <- bind_rows(res_list)
          res_all <- left_join(res_all,metadata,by=c("variable"="OriginalID"))
          astcv_kept <- res_all
        }
    astcv <- astcv_kept %>% dplyr::group_by(SubjectID,protein) %>% dplyr::summarise(CV = explore_cv(value),cate="MS-Astral", cate1=cate1input, .groups = "drop") %>%
      dplyr::filter(!is.nan(CV))
    bindcv2 <- rbind(bindcv2,astcv)
    cvlow <- Inf
    astcv1 <-  astcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(SubjectID, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-Astral") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(astcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,astcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
somacvf <- function(pro,isna = TRUE,cate1input,keepfilter=TRUE){
  somalod <- readRDS("data/processed/SomaScan_eLOD.rds")
  somalod <- somalod[,c(1,11)]
  bindcv2 <- NULL
  bindcv <- NULL
  overlap5soma <- dplyr::filter(somer_protein,EntrezGeneSymbol %in% pro)$SeqId
  somacv <- as.data.frame((soma)) %>% rownames_to_column(var = "protein") %>% dplyr::filter(protein %in% overlap5soma) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  somacv <- left_join(somacv,somalod,by=c("protein"="SeqId"))
  keep_tbl <- somacv %>% dplyr::mutate(
    valid    = !is.na(value) & !is.na(eLod),
    detected = valid & (value > eLod)
  ) %>%
    dplyr::group_by(protein, SubjectID) %>%
    dplyr::summarise(
      n_total  = sum(valid),
      n_detect = sum(detected),
      prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
      keep     = cv_detection_keep(prop, isna, keepfilter),
      .groups  = "drop")
  somacv_kept <- somacv %>%
    dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(protein, SubjectID),by = c("protein", "SubjectID"))
  if(isna == TRUE){
    somacv_kept <- somacv_kept %>% dplyr::mutate(value = if_else(value < eLod, NA_real_, value))
  }else{
    somacv_kept <- somacv_kept
  }
  somacv <- somacv_kept %>% dplyr::group_by(SubjectID,protein) %>% dplyr::summarise( CV=explore_cv(log2(value)),cate="Somascan",
                                                                                cate1 = cate1input, .groups = "drop")
  bindcv2 <- rbind(bindcv2,somacv)
  cvlow <- Inf
  somacv1 <-  somacv %>%
    dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(SubjectID, CV_rounded) %>%
    dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Somascan") %>%
    dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(somacv$protein)),cate1 = cate1input)
  bindcv <- rbind(bindcv,somacv1)
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
olinkcvall <- function(cate,pro,isna=TRUE,cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "3072"){
    prot <- pro
    o3072cv <- lod3072 %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- o3072cv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
    ) %>%
      dplyr::group_by(Assay) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    o3072cv_kept <- o3072cv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay),by = c("Assay"))
    if(isna == TRUE){
      o3072cv_kept <- o3072cv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      o3072cv_kept <- o3072cv_kept
    }
    o3072cv <- o3072cv_kept %>%  dplyr::group_by(Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-Explore", cate1=cate1input,.groups = "drop")
    colnames(o3072cv)  <- c("protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,o3072cv)
    cvlow <- Inf
    o3072cv1 <-  o3072cv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by( CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-Explore") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(o3072cv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,o3072cv1)
  }else{
    prot <- pro
    ohtcv <- lodht %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- ohtcv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
    ) %>%
      dplyr::group_by(Assay) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    ohtcv_kept <- ohtcv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay),by = c("Assay"))
    if(isna == TRUE){
      ohtcv_kept <- ohtcv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      ohtcv_kept <- ohtcv_kept
    }
    ohtcv <- ohtcv_kept %>%  dplyr::group_by(Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-HT",cate1=cate1input ,.groups = "drop")
    colnames(ohtcv) <- c( "protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,ohtcv)
    cvlow <- Inf
    ohtcv1 <-  ohtcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by( CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-HT") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(ohtcv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,ohtcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
mscvall <- function(dflistna,cate,pro,isna="0",cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "ot480"){
    prot <- pro
    otcv <- as.data.frame(dflistna[["MS-OT480"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- otcv %>%
      dplyr::group_by(protein) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
    if (isna == "0") {
      otcv_kept <- otcv %>% dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%
                                                dplyr::select(protein),by = c("protein"))
      otcv_kept <- otcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else if (isna == "NA") {
      otcv_kept <- otcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein),by = c("protein"))
    }else if (isna == "impute"){
      df <- as.data.frame(dflistna[["MS-OT480"]])
      df <- as.matrix(df[pro,])
      colmap <- data.frame(
        col = colnames(df),
        replicate = as.integer(str_replace(colnames(df), "_.*", "")),
        SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
        stringsAsFactors = FALSE)
        thr <- ceiling(ncol(df) * 0.6)
        keep <- rowSums(!is.na(df)) >= thr
        x_filtered <- df[keep, , drop = FALSE]
        x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
        df_imp <- as.data.frame(x_imp, check.names = FALSE)
        df_imp$protein <- rownames(df_imp)
        long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
      res_all <- left_join(long,metadata,by=c("variable"="OriginalID"))
      otcv_kept <- res_all
    }
    otcv <- otcv_kept %>% dplyr::group_by(protein) %>%
      dplyr::summarise(CV = explore_cv(value),cate="MS-OT480", cate1=cate1input, .groups = "drop")
    bindcv2 <- rbind(bindcv2,otcv)
    cvlow <- Inf
    otcv1 <-  otcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by( CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-OT480") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(otcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,otcv1)
  }else{
    prot <- pro
    astcv <- as.data.frame(dflistna[["MS-Astral"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- astcv %>%
      dplyr::group_by(protein) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
    if (isna == "0") {
      astcv_kept <- astcv %>% dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%
                                                dplyr::select(protein),by = c("protein"))
      astcv_kept <- astcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else if (isna == "NA") {
      astcv_kept <- astcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein),by = c("protein"))
    }else if (isna == "impute"){
      df <- as.data.frame(dflistna[["MS-Astral"]])
      df <- as.matrix(df[pro,])
      colmap <- data.frame(
        col = colnames(df),
        replicate = as.integer(str_replace(colnames(df), "_.*", "")),
        SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
        stringsAsFactors = FALSE)
      thr <- ceiling(ncol(df) * 0.6)
      keep <- rowSums(!is.na(df)) >= thr
      x_filtered <- df[keep, , drop = FALSE]
      x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
      df_imp <- as.data.frame(x_imp, check.names = FALSE)
      df_imp$protein <- rownames(df_imp)
      long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
      res_all <- left_join(long,metadata,by=c("variable"="OriginalID"))
      astcv_kept <- res_all
    }
    astcv <- astcv_kept %>% dplyr::group_by(protein) %>% dplyr::summarise(CV = explore_cv(value),cate="MS-Astral", cate1=cate1input, .groups = "drop")
    bindcv2 <- rbind(bindcv2,astcv)
    cvlow <- Inf
    astcv1 <-  astcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by( CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-Astral") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(astcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,astcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
somacvfall <- function(pro,isna = TRUE,cate1input,keepfilter=TRUE){
  somalod <- readRDS("data/processed/SomaScan_eLOD.rds")
  somalod <- somalod[,c(1,11)]
  bindcv2 <- NULL
  bindcv <- NULL
  overlap5soma <- dplyr::filter(somer_protein,EntrezGeneSymbol %in% pro)$SeqId
  somacv <- as.data.frame((soma)) %>% rownames_to_column(var = "protein") %>% dplyr::filter(protein %in% overlap5soma) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  somacv <- left_join(somacv,somalod,by=c("protein"="SeqId"))
  keep_tbl <- somacv %>% dplyr::mutate(
    valid    = !is.na(value) & !is.na(eLod),
    detected = valid & (value > eLod)
  ) %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(
      n_total  = sum(valid),
      n_detect = sum(detected),
      prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
      keep     = cv_detection_keep(prop, isna, keepfilter),
      .groups  = "drop")
  somacv_kept <- somacv %>%
    dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(protein),by = c("protein"))
  if(isna == TRUE){
    somacv_kept <- somacv_kept %>% dplyr::mutate(value = if_else(value < eLod, NA_real_, value))
  }else{
    somacv_kept <- somacv_kept
  }
  somacv <- somacv_kept %>% dplyr::group_by(protein) %>% dplyr::summarise( CV=explore_cv(log2(value)),cate="Somascan",
                                                                                     cate1 = cate1input, .groups = "drop")
  bindcv2 <- rbind(bindcv2,somacv)
  cvlow <- Inf
  somacv1 <-  somacv %>%
    dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by( CV_rounded) %>%
    dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Somascan") %>%
    dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(somacv$protein)),cate1 = cate1input)
  bindcv <- rbind(bindcv,somacv1)
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
olinkcvindi <- function(cate,pro,isna=TRUE,cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "3072"){
    prot <- pro
    o3072cv <- lod3072 %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- o3072cv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
    ) %>% dplyr::group_by(Assay) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    o3072cv_kept <- o3072cv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay),by = c("Assay"))
    if(isna == TRUE){
      o3072cv_kept <- o3072cv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      o3072cv_kept <- o3072cv_kept
    }
    o3072cv <- o3072cv_kept %>%  dplyr::group_by(Time,Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-Explore", cate1=cate1input,.groups = "drop") %>%
      dplyr::filter(!is.na(CV))
    colnames(o3072cv)  <- c("Time", "protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,o3072cv)
    cvlow <- Inf
    o3072cv1 <-  o3072cv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(Time, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-Explore") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(o3072cv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,o3072cv1)
  }else{
    prot <- pro
    ohtcv <- lodht %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
    keep_tbl <- ohtcv %>% dplyr::mutate(
      valid    = !is.na(NPX) & !is.na(LOD),
      detected = valid & (NPX > LOD)
    ) %>%
      dplyr::group_by(Assay) %>%
      dplyr::summarise(
        n_total  = sum(valid),
        n_detect = sum(detected),
        prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
        keep     = cv_detection_keep(prop, isna, keepfilter),
        .groups  = "drop")
    ohtcv_kept <- ohtcv %>%
      dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(Assay),by = c("Assay"))
    if(isna == TRUE){
      ohtcv_kept <- ohtcv_kept %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
    }else{
      ohtcv_kept <- ohtcv_kept
    }
    ohtcv <- ohtcv_kept %>%  dplyr::group_by(Time,Assay) %>% dplyr::summarise(CV = explore_cv(NPX),cate="Olink-HT",cate1=cate1input ,.groups = "drop") %>%
      dplyr::filter(!is.na(CV))
    colnames(ohtcv) <- c("Time", "protein","CV","cate","cate1")
    bindcv2 <- rbind(bindcv2,ohtcv)
    cvlow <- Inf
    ohtcv1 <-  ohtcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(Time, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Olink-HT") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(ohtcv$protein)),cate1=cate1input )
    bindcv <- rbind(bindcv,ohtcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
mscvindi <- function(dflistna,cate,pro,isna="0",cate1input,keepfilter=TRUE){
  bindcv2 <- NULL
  bindcv <- NULL
  if(cate == "ot480"){
    prot <- pro
    otcv <- as.data.frame(dflistna[["MS-OT480"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- otcv %>%
      dplyr::group_by(protein) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
    if (isna == "0") {
      otcv_kept <- otcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein),by = c("protein"))
      otcv_kept <- otcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else if (isna == "NA") {
      otcv_kept <- otcv %>%
        dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein),by = c("protein"))
    }else if (isna == "impute"){
      df <- as.data.frame(dflistna[["MS-OT480"]])
      df <- as.matrix(df[pro,])
      colmap <- data.frame(
        col = colnames(df),
        replicate = as.integer(str_replace(colnames(df), "_.*", "")),
        SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
        stringsAsFactors = FALSE)
      thr <- ceiling(ncol(df) * 0.6)
      keep <- rowSums(!is.na(df)) >= thr
      x_filtered <- df[keep, , drop = FALSE]
      x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
      df_imp <- as.data.frame(x_imp, check.names = FALSE)
      df_imp$protein <- rownames(df_imp)
      long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
      res_all <- left_join(long,metadata,by=c("variable"="OriginalID"))
      otcv_kept <- res_all
    }
    otcv <- otcv_kept %>% dplyr::group_by(Time,protein) %>% dplyr::summarise(CV = explore_cv(value),cate="MS-OT480", cate1=cate1input, .groups = "drop") %>%
      dplyr::filter(!is.nan(CV))
    bindcv2 <- rbind(bindcv2,otcv)
    cvlow <- Inf
    otcv1 <-  otcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(Time, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-OT480") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(otcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,otcv1)
  }else{
    prot <- pro
    astcv <- as.data.frame(dflistna[["MS-Astral"]]) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
      melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
    keep_tbl <- astcv %>%
      dplyr::group_by(protein) %>%
      dplyr::summarise(
        n_total  = n(),
        n_non_na = sum(!is.na(value)),
        prop     = if_else(n_total > 0, n_non_na / n_total, NA_real_),
        keep     = (!is.na(prop) & prop >= 0.6) | !keepfilter,
        .groups  = "drop"
      )
        if (isna == "0") {
          astcv_kept <- astcv %>% dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%
                                                      dplyr::select(protein),by = c("protein"))
          astcv_kept <- astcv_kept %>%dplyr::mutate(value = ifelse(is.na(value), 0, value))
        } else if (isna == "NA") {
          astcv_kept <- astcv %>%
            dplyr::inner_join(keep_tbl %>%dplyr::filter(keep) %>%dplyr::select(protein),by = c("protein"))
        }else if (isna == "impute"){
          df <- as.data.frame(dflistna[["MS-Astral"]])
          df <- as.matrix(df[pro,])
          colmap <- data.frame(
            col = colnames(df),
            replicate = as.integer(str_replace(colnames(df), "_.*", "")),
            SubjectID = as.character(str_replace(colnames(df), ".*_", "")),
            stringsAsFactors = FALSE)
          thr <- ceiling(ncol(df) * 0.6)
          keep <- rowSums(!is.na(df)) >= thr
          x_filtered <- df[keep, , drop = FALSE]
          x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
          df_imp <- as.data.frame(x_imp, check.names = FALSE)
          df_imp$protein <- rownames(df_imp)
          long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
          res_all <- left_join(long,metadata,by=c("variable"="OriginalID"))
          astcv_kept <- res_all
        }
    astcv <- astcv_kept %>% dplyr::group_by(Time,protein) %>% dplyr::summarise(CV = explore_cv(value),cate="MS-Astral", cate1=cate1input, .groups = "drop") %>%
      dplyr::filter(!is.nan(CV))
    bindcv2 <- rbind(bindcv2,astcv)
    cvlow <- Inf
    astcv1 <-  astcv %>%
      dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(Time, CV_rounded) %>%
      dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="MS-Astral") %>%
      dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(astcv$protein)),cate1=cate1input  )
    bindcv <- rbind(bindcv,astcv1)
  }
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
somacvindi <- function(pro,isna = TRUE,cate1input,keepfilter=TRUE){
  somalod <- readRDS("data/processed/SomaScan_eLOD.rds")
  somalod <- somalod[,c(1,11)]
  bindcv2 <- NULL
  bindcv <- NULL
  overlap5soma <- dplyr::filter(somer_protein,EntrezGeneSymbol %in% pro)$SeqId
  somacv <- as.data.frame((soma)) %>% rownames_to_column(var = "protein") %>% dplyr::filter(protein %in% overlap5soma) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  somacv <- left_join(somacv,somalod,by=c("protein"="SeqId"))
  keep_tbl <- somacv %>% dplyr::mutate(
    valid    = !is.na(value) & !is.na(eLod),
    detected = valid & (value > eLod)
  ) %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(
      n_total  = sum(valid),
      n_detect = sum(detected),
      prop     = if_else(n_total > 0, n_detect / n_total, NA_real_),
      keep     = cv_detection_keep(prop, isna, keepfilter),
      .groups  = "drop")
  somacv_kept <- somacv %>%
    dplyr::inner_join(keep_tbl %>% dplyr::filter(keep) %>% dplyr::select(protein),by = c("protein"))
  if(isna == TRUE){
    somacv_kept <- somacv_kept %>% dplyr::mutate(value = if_else(value < eLod, NA_real_, value))
  }else{
    somacv_kept <- somacv_kept
  }
  somacv <- somacv_kept %>% dplyr::group_by(Time,protein) %>% dplyr::summarise( CV=explore_cv(log2(value)),cate="Somascan",
                                                                                     cate1 = cate1input, .groups = "drop")
  bindcv2 <- rbind(bindcv2,somacv)
  cvlow <- Inf
  somacv1 <-  somacv %>%
    dplyr::mutate(CV_rounded = round(abs(CV), 1)) %>%   dplyr::group_by(Time, CV_rounded) %>%
    dplyr::summarise(count = n(), .groups = "drop")  %>% dplyr::mutate(cate="Somascan") %>%
    dplyr::filter(CV_rounded <= cvlow) %>%  dplyr::mutate(percent = count/length(unique(somacv$protein)),cate1 = cate1input)
  bindcv <- rbind(bindcv,somacv1)
  cvresult <- list()
  cvresult[["vol"]] <- bindcv2
  cvresult[["bar"]] <- bindcv
  return(cvresult)
}
summarise_soma_cv_to_protein <- function(cv_df, strata_cols = character()) {
  required_cv_cols <- c("protein", "CV")
  missing_cv_cols <- setdiff(required_cv_cols, colnames(cv_df))
  if (length(missing_cv_cols) > 0) {
    stop("SomaScan CV data are missing required columns: ",
         paste(missing_cv_cols, collapse = ", "))
  }
  missing_strata_cols <- setdiff(strata_cols, colnames(cv_df))
  if (length(missing_strata_cols) > 0) {
    stop("SomaScan CV data are missing requested strata columns: ",
         paste(missing_strata_cols, collapse = ", "))
  }
  if (!exists("somer_protein", inherits = TRUE)) {
    stop("Object 'somer_protein' is required for SeqId-to-protein mapping.")
  }
  soma_map_source <- get("somer_protein", inherits = TRUE)
  required_map_cols <- c("SeqId", "EntrezGeneSymbol")
  missing_map_cols <- setdiff(required_map_cols, colnames(soma_map_source))
  if (length(missing_map_cols) > 0) {
    stop("Object 'somer_protein' is missing columns: ",
         paste(missing_map_cols, collapse = ", "))
  }
  soma_map <- soma_map_source %>%
    dplyr::transmute(
      somer = as.character(SeqId),
      protein_name = trimws(as.character(EntrezGeneSymbol))
    ) %>%
    dplyr::filter(
      !is.na(somer), nzchar(somer),
      !is.na(protein_name), nzchar(protein_name),
      !toupper(protein_name) %in% c("NA", "N/A", "NULL")
    ) %>%
    dplyr::distinct()
  ambiguous_somers <- soma_map %>%
    dplyr::count(somer, name = "n_protein_names") %>%
    dplyr::filter(n_protein_names > 1) %>%
    dplyr::pull(somer)
  if (length(ambiguous_somers) > 0) {
    message(
      "Excluded ", length(ambiguous_somers),
      " SomaScan SeqId(s) with ambiguous protein-name mapping."
    )
    soma_map <- soma_map %>%
      dplyr::filter(!somer %in% ambiguous_somers)
  }
  input_somers <- unique(as.character(cv_df$protein))
  unmapped_somers <- setdiff(input_somers, soma_map$somer)
  if (length(unmapped_somers) > 0) {
    message(
      "Excluded ", length(unmapped_somers),
      " SomaScan SeqId(s) without a unique EntrezGeneSymbol mapping."
    )
  }
  cv_mapped <- cv_df %>%
    dplyr::mutate(
      somer = as.character(protein),
      CV = suppressWarnings(as.numeric(CV))
    ) %>%
    dplyr::select(-protein) %>%
    dplyr::inner_join(soma_map, by = "somer") %>%
    dplyr::filter(is.finite(CV))
  preserve_cols <- intersect(c("cate", "cate1"), colnames(cv_mapped))
  protein_group_cols <- unique(c(strata_cols, preserve_cols, "protein_name"))
  somer_group_cols <- unique(c(protein_group_cols, "somer"))
  somer_cv <- cv_mapped %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(somer_group_cols))) %>%
    dplyr::summarise(somer_CV = mean(CV), .groups = "drop")
  protein_cv <- somer_cv %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(protein_group_cols))) %>%
    dplyr::summarise(
      CV = mean(somer_CV),
      n_somers = dplyr::n_distinct(somer),
      cv_sd = ifelse(n_somers > 1, stats::sd(somer_CV), NA_real_),
      cv_min = min(somer_CV),
      cv_max = max(somer_CV),
      somer_ids = paste(sort(unique(somer)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(protein = protein_name)
  protein_cv
}
safe_cortest <- function(x, y, method = "spearman",	use = "pairwise.complete.obs") {
  idx <- is.finite(x) & is.finite(y)
  if (sum(idx) >= 5) {
    return(cor.test(x[idx], y[idx], method = method,	use = "pairwise.complete.obs"))
  } else {
    return(NA_real_)
  }
}
mscorrdiffindi3time <- function(pro = prolist[["MS-OT480"]],
                              cate = "MS-OT480",
                              data = dflistna[["MS-OT480"]],
                              isna = "0"
                    ){
  prot <- pro
  data1 <- as.data.frame(data) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  data2 <- data1 %>%
    dplyr::group_by(protein, SubjectID) %>%
    dplyr::summarise(n_total  = n(),n_non_na = sum(!is.na(value)),prop= if_else(n_total > 0, n_non_na / n_total, NA_real_),
                     keep= !is.na(prop) & prop >= 0.6,  .groups  = "drop" )
  if (isna == "0") {
    data3 <- data1%>%dplyr::inner_join(data2%>%dplyr::filter(keep)%>%dplyr::select(protein,SubjectID),by = c("protein","SubjectID"))
    data3 <- data3%>%dplyr::mutate(value = ifelse(is.na(value),0,value))
  } else if (isna == "NA") {
    data3 <- data1%>%dplyr::inner_join(data2 %>%dplyr::filter(keep)%>%dplyr::select(protein,SubjectID),by = c("protein","SubjectID"))
  } else if (isna == "impute"){
    df <- as.data.frame(data)
    df <- as.matrix(df[pro,])
    colmap <- data.frame(col = colnames(df),replicate = as.integer(str_replace(colnames(df), "_.*", "")),
                         SubjectID = as.character(str_replace(colnames(df), ".*_", "")),stringsAsFactors = FALSE)
    subjects <- unique(colmap$SubjectID)
    res_list <- list()
    for (s in subjects) {
      idx <- which(colmap$SubjectID == s)
      subX <- df[, idx, drop = FALSE]
      thr <- ceiling(ncol(subX) * 0.6)
      keep <- rowSums(!is.na(subX)) >= thr
      x_filtered <- subX[keep, , drop = FALSE]
      x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
      df_imp <- as.data.frame(x_imp, check.names = FALSE)
      df_imp$protein <- rownames(df_imp)
      res_list[[s]] <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
    }
    res_all <- bind_rows(res_list)
    res_all <- left_join(res_all,metadata,by=c("variable"="OriginalID"))
    data3 <- res_all
  }
  df_wide <- data3 %>%
    dplyr::select( protein, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)%>% column_to_rownames(var="protein")%>%
    as.matrix()
  R <- cor(df_wide, use = "pairwise.complete.obs", method = "spearman")
  R[lower.tri(R, diag = TRUE)] <- NA
  out <- as.data.frame(as.table(R), stringsAsFactors = FALSE) %>%
    tidyr::drop_na(Freq) %>%
    mutate(
      cate = cate,
      Var1_prefix = sub("_.*", "", Var1),
      Var2_prefix = sub("_.*", "", Var2),
      Pair        = paste(pmin(Var1_prefix, Var2_prefix),"vs",pmax(Var1_prefix, Var2_prefix)),
      S1 = sub(".*_", "", Var1),
      S2 = sub(".*_", "", Var2)
    )
  return(out)
}
olinkcorrdiffindi3time <- function(pro = prolist[["Olink-Explore"]],cate = "Olink-Explore",
                                   data = dflistna[["Olink-Explore"]],isna = "NA",
                                   keepfilter = TRUE
                 ){
  if(cate == "Olink-Explore"){
    prot <- pro
    data1 <- lod3072 %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
  }else{
    prot <- pro
    data1 <- lodht %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
  }
  data2 <- data1 %>% dplyr::mutate(valid = !is.na(NPX) & !is.na(LOD), detected = valid & (NPX > LOD)) %>%
    dplyr::group_by(Assay, SubjectID) %>%
    dplyr::summarise( n_total  = sum(valid),n_detect = sum(detected),
                      prop  = if_else(n_total > 0, n_detect / n_total, NA_real_),
                      keep = cv_detection_keep(prop, isna, keepfilter),  .groups  = "drop")
  data3 <- data1 %>% dplyr::inner_join(data2 %>% dplyr::filter(keep) %>% dplyr::select(Assay, SubjectID),by = c("Assay", "SubjectID"))
  if(isna == "NA"){
    data3 <- data3 %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
  }else{
    data3 <- data3
  }
  df_wide <- data3 %>% dplyr::select( Assay, SampleID, NPX) %>%
    pivot_wider(names_from = SampleID, values_from = NPX)%>% column_to_rownames(var="Assay")%>%
    as.matrix()
  R <- cor(df_wide, use = "pairwise.complete.obs", method = "spearman")
  R[lower.tri(R, diag = TRUE)] <- NA
  out <- as.data.frame(as.table(R), stringsAsFactors = FALSE) %>%
    tidyr::drop_na(Freq) %>%
    mutate(cate = cate,Var1_prefix = sub("_.*", "", Var1),  Var2_prefix = sub("_.*", "", Var2),
            Pair = paste(pmin(Var1_prefix, Var2_prefix),"vs",pmax(Var1_prefix, Var2_prefix)),
            S1 = sub(".*_", "", Var1), S2 = sub(".*_", "", Var2) )
  return(out)
}
somacorrdiffindi3time <- function(pro = prolist[["Somascan"]],cate = "Somascan",
                                  data = dflistna[["Somascan"]],isna = "NA",
                                  keepfilter = TRUE
){
  somalod <- readRDS("data/processed/SomaScan_eLOD.rds")
  somalod <- somalod[,c(1,11)]
  overlap5soma <- dplyr::filter(somer_protein,EntrezGeneSymbol %in% pro)$SeqId
  data1 <- as.data.frame((soma)) %>% rownames_to_column(var = "protein") %>% dplyr::filter(protein %in% overlap5soma) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  data1 <- left_join(data1,somalod,by=c("protein"="SeqId"))
  data2 <- data1 %>% dplyr::mutate( valid = !is.na(value) & !is.na(eLod), detected = valid & (value > eLod) ) %>%
    dplyr::group_by(protein, SubjectID) %>%
    dplyr::summarise(n_total  = sum(valid),n_detect = sum(detected),prop = if_else(n_total > 0, n_detect / n_total, NA_real_),
                     keep = cv_detection_keep(prop, isna, keepfilter),.groups  = "drop")
  data3 <- data1 %>% dplyr::inner_join(data2 %>% dplyr::filter(keep) %>% dplyr::select(protein, SubjectID),by = c("protein", "SubjectID"))
  if(isna == "NA"){
    data3 <- data3 %>% dplyr::mutate(value = if_else(value < eLod, NA_real_, value))
  }else{
    data3 <- data3
  }
  data3$value <- log2(data3$value)
  df_wide <- data3 %>% dplyr::select( protein, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)%>% column_to_rownames(var="protein")%>%
    as.matrix()
  R <- cor((df_wide), use = "pairwise.complete.obs", method = "spearman")
  R[lower.tri(R, diag = TRUE)] <- NA
  out <- as.data.frame(as.table(R), stringsAsFactors = FALSE) %>%
    tidyr::drop_na(Freq) %>%
    mutate( cate = cate,Var1_prefix = sub("_.*", "", Var1),  Var2_prefix = sub("_.*", "", Var2),
            Pair = paste(pmin(Var1_prefix, Var2_prefix),"vs",pmax(Var1_prefix, Var2_prefix)),
            S1 = sub(".*_", "", Var1), S2 = sub(".*_", "", Var2) )
  return(out)
}
msmatrixpro <- function( pro = prolist[["MS-OT480"]],
                         cate = "MS-OT480",
                         data = dflistna[["MS-OT480"]],
                         isna = "NA" ){
  prot <- pro
  data1 <- as.data.frame(data) %>% rownames_to_column(var = "protein")  %>%  dplyr::filter(protein %in% prot) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  data2 <- data1 %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(n_total  = n(),n_non_na = sum(!is.na(value)),prop= if_else(n_total > 0, n_non_na / n_total, NA_real_),
                     keep= !is.na(prop) & prop >= 0.6,  .groups  = "drop" )
  if (isna == "0") {
    data3 <- data1%>%dplyr::inner_join(data2%>%dplyr::filter(keep)%>%dplyr::select(protein),by = c("protein"))
    data3 <- data3%>%dplyr::mutate(value = ifelse(is.na(value),0,value))
  } else if (isna == "NA") {
    data3 <- data1%>%dplyr::inner_join(data2 %>%dplyr::filter(keep)%>%dplyr::select(protein),by = c("protein"))
  } else if (isna == "impute"){
    df <- as.data.frame(data)
    df <- as.matrix(df[pro,])
    colmap <- data.frame(col = colnames(df),replicate = as.integer(str_replace(colnames(df), "_.*", "")),
                         SubjectID = as.character(str_replace(colnames(df), ".*_", "")),stringsAsFactors = FALSE)
    thr <- ceiling(ncol(df) * 0.6)
    keep <- rowSums(!is.na(df)) >= thr
    x_filtered <- df[keep, , drop = FALSE]
    x_imp <- impute_perseus_log2(x_filtered,seed = 2025)
    df_imp <- as.data.frame(x_imp, check.names = FALSE)
    df_imp$protein <- rownames(df_imp)
    long <- melt(df_imp, id.vars = "protein",variable.name = "variable",value.name = "value")
    res_all <- left_join(long,metadata,by=c("variable"="OriginalID"))
    data3 <- res_all
  }
  df_wide <- data3 %>%
    dplyr::select( protein, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)%>% column_to_rownames(var="protein")%>%
    as.matrix()
  df_wide <- df_wide[,metadata$OriginalID]
  return(df_wide)
}
olinkmatrixpro <- function(pro = prolist[["Olink-Explore"]],cate = "Olink-Explore",
                                 data = dflistna[["Olink-Explore"]],isna = "NA",
                                 keepfilter = TRUE
){
  prot <- pro
  if(cate == "Olink-Explore"){
    data1 <- lod3072 %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
  }else{
    data1 <- lodht %>%  dplyr::filter(Assay %in% prot) %>% left_join(metadata,by=c("SampleID"="OriginalID"))
  }
  data2 <- data1 %>% dplyr::mutate(valid = !is.na(NPX) & !is.na(LOD), detected = valid & (NPX > LOD)) %>%
    dplyr::group_by( Assay) %>%
    dplyr::summarise( n_total  = sum(valid),n_detect = sum(detected),
                      prop  = if_else(n_total > 0, n_detect / n_total, NA_real_),
                      keep = cv_detection_keep(prop, isna, keepfilter),  .groups  = "drop")
  data3 <- data1 %>% dplyr::inner_join(data2 %>% dplyr::filter(keep) %>% dplyr::select( Assay),by = c("Assay"))
  if(isna == "NA"){
    data3 <- data3 %>% dplyr::mutate(NPX = if_else(NPX < LOD, NA_real_, NPX))
  }else{
    data3 <- data3
  }
  df_wide <- data3 %>%
    dplyr::select( Assay, SampleID, NPX) %>%
    pivot_wider(names_from = SampleID, values_from = NPX)%>% column_to_rownames(var="Assay")%>%
    as.matrix()
  df_wide <- df_wide[,metadata$OriginalID]
  return(df_wide)
}
somamatrixpro <- function(pro = prolist[["Somascan"]],cate = "Somascan",
                        data = dflistna[["Somascan"]],isna = "NA",
                        keepfilter = TRUE
){
  somalod <- readRDS("data/processed/SomaScan_eLOD.rds")
  somalod <- somalod[,c(1,11)]
  overlap5soma <- dplyr::filter(somer_protein,EntrezGeneSymbol %in% pro)$SeqId
  data1 <- as.data.frame((soma)) %>% rownames_to_column(var = "protein") %>% dplyr::filter(protein %in% overlap5soma) %>%
    melt() %>% left_join(metadata,by=c("variable"="OriginalID"))
  data1 <- left_join(data1,somalod,by=c("protein"="SeqId"))
  data2 <- data1 %>% dplyr::mutate( valid = !is.na(value) & !is.na(eLod), detected = valid & (value > eLod) ) %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(n_total  = sum(valid),n_detect = sum(detected),prop = if_else(n_total > 0, n_detect / n_total, NA_real_),
                     keep = cv_detection_keep(prop, isna, keepfilter),.groups  = "drop")
  data3 <- data1 %>% dplyr::inner_join(data2 %>% dplyr::filter(keep) %>% dplyr::select(protein),by = c("protein"))
  if(isna == "NA"){
    data3 <- data3 %>% dplyr::mutate(value = if_else(value < eLod, NA_real_, value))
  }else{
    data3 <- data3
  }
  data3$value <- log2(data3$value)
  df_wide <- data3 %>%
    dplyr::select( protein, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)%>% column_to_rownames(var="protein")%>%
    as.matrix()
  df_wide <- df_wide[,metadata$OriginalID]
  return(df_wide)
}
procorr <- function(input){
  bindresult <- NULL
  for (i in 1:3) {
    df1_orig <- input[[i]]
    for (j in (i+1):4) {
      df2_orig <- input[[j]]
      overlapprotein <- intersect(rownames(df1_orig), rownames(df2_orig))
      df1_sub <- df1_orig[overlapprotein, metadata$OriginalID, drop = FALSE]
      df2_sub <- df2_orig[overlapprotein, metadata$OriginalID, drop = FALSE]
      results <- data.frame(
        Protein = overlapprotein,
        Somer = NA,
        Correlation = NA,
        p_value = NA,
        cate = paste0(names(input)[i], " vs ", names(input)[j]),
        stringsAsFactors = FALSE
      )
      for(a in seq_along(overlapprotein)){
        prot <- overlapprotein[a]
        x <- as.numeric(df1_sub[prot, ])
        y <- as.numeric(df2_sub[prot, ])
        test <- safe_cortest(x, y, method = "spearman", use = "pairwise.complete.obs")
        bad <- is.null(test) || inherits(test, "try-error") ||
          !is.list(test) ||
          any(!is.finite(c(
            suppressWarnings(as.numeric(test$estimate)),
            suppressWarnings(as.numeric(test$p.value))
          )))
        if(isTRUE(bad)){
          results$Correlation[a] <- NA
          results$p_value[a] <- NA
        } else {
          results$Correlation[a] <- test$estimate
          results$p_value[a] <- test$p.value
        }
      }
      bindresult <- rbind(bindresult, results)
    }
  }
  bindresult2 <- NULL
  for (i in 1:4) {
    df1_orig <- input[[i]]
    df2_orig <- input[[5]]
    df1_sub <- df1_orig[, metadata$OriginalID, drop = FALSE]
    df2_sub <- df2_orig[, metadata$OriginalID, drop = FALSE]
    protsoma <- somer_protein %>%
      dplyr::filter(SeqId %in% rownames(df2_sub)) %>%
      dplyr::pull(EntrezGeneSymbol) %>%
      unique()
    overlapprotein <- intersect(rownames(df1_sub), protsoma)
    somaprotein <- somer_protein %>% dplyr::filter(EntrezGeneSymbol %in% overlapprotein)
    somaprotein <- somaprotein[, c("SeqId", "EntrezGeneSymbol")]
    results <- data.frame(
      Protein = somaprotein$EntrezGeneSymbol,
      Somer = somaprotein$SeqId,
      Correlation = NA,
      p_value = NA,
      cate = paste0(names(input)[i], " vs ", names(input)[5]),
      stringsAsFactors = FALSE
    )
    for(a in seq_along(somaprotein$EntrezGeneSymbol)){
      somer <- somaprotein$SeqId[a]
      prot <- unique(dplyr::filter(somer_protein, SeqId %in% somer)$EntrezGeneSymbol)
      x <- as.numeric(df1_sub[prot, ])
      y <- as.numeric(df2_sub[somer, ])
      test <- safe_cortest(x, y, method = "spearman", use = "pairwise.complete.obs")
      bad <- is.null(test) || inherits(test, "try-error") ||
        !is.list(test) ||
        any(!is.finite(c(
          suppressWarnings(as.numeric(test$estimate)),
          suppressWarnings(as.numeric(test$p.value))
        )))
      if(isTRUE(bad)){
        results$Correlation[a] <- NA
        results$p_value[a] <- NA
      } else {
        results$Correlation[a] <- test$estimate
        results$p_value[a] <- test$p.value
      }
    }
    bindresult2 <- rbind(bindresult2, results)
  }
  bindresult2 <- rbind(bindresult2, bindresult)
  bindresultsample <- bindresult2
  return(bindresultsample)
}
mutl_anova_test <- function(df_ancova,dependent,independent){
  fit <- aov(as.formula(paste(dependent, "~", paste(independent, collapse = "*"))),
             data=df_ancova)
  result2 <- TukeyHSD(fit)
  result2 <- do.call(rbind, TukeyHSD(fit))
  result2 <- rownames_to_column(as.data.frame(result2))
  result <- data.frame(do.call("rbind", strsplit(as.character(result2$rowname), "[:-]")))
  result <- cbind(result,result2)
  result$mark <- ifelse(result$`p adj` < 0.05, "*", NA)
  return(result)
}
tissuelistfunc <- function(inputdir, backprotein, outputdir){
  hpaannbind <- inputdir
  hpaannbind <- hpaannbind %>% dplyr::filter(protein1 %in% backprotein)
  data <- dplyr::filter(hpaannbind,RNA.tissue.specificity=="Tissue enriched")
  a <- as.data.frame(data[,"RNA.tissue.specific.nTPM"])
  a<- separate(a, colnames(a)[1],
               into =c("RNA_Tissue1","RNA_Tissue2","RNA_Tissue3","RNA_Tissue4","RNA_Tissue5","RNA_Tissue6","RNA_Tissue7"),
               sep = ";")
  a[a == ""] <- NA
  c <- NULL
  b <- as.data.frame(na.omit(a[,"RNA_Tissue1"]))
  b <- separate(b, colnames(b)[1],
                into =c("RNA_Tissue1.1","RNA_Tissue1.2"),
                sep = ": ")
  c <- c(c,b[,1])
  c <- as.data.frame(c)
  c <- cbind(c,data[,"protein1"])
  colnames(c) <- c("tissue","gene")
  result <- spread(c, gene, gene)
  tissuegene <- list()
  for (i in result$tissue) {
    e <- dplyr::filter(result, tissue == i)
    e <- column_to_rownames(e,var="tissue")
    e <- t(e)
    e <- na.omit(e)
    tissuegene[[i]] <- e[,1]
  }
  return(tissuegene)
}
sclistfunc <- function(inputdir, backprotein, outputdir){
  hpaannbind <- inputdir
  hpaannbind <- hpaannbind %>% dplyr::filter(protein1 %in% backprotein)
  a <- dplyr::filter(hpaannbind,RNA.single.cell.type.specificity=="Cell type enriched")
  a <- dplyr::filter(a,
                     RNA.single.cell.type.distribution == "Detected in all" |
                       RNA.single.cell.type.distribution == "Detected in many" |
                       RNA.single.cell.type.distribution == "Detected in some" |
                       RNA.single.cell.type.distribution == "Detected in single")
  a1 <- as.data.frame(a[,c("protein1","RNA.single.cell.type.specific.nTPM")])
  a1 <- separate(a1, colnames(a1)[2],
                 into =c("RNA_Tissue1","RNA_Tissue2","RNA_Tissue3","RNA_Tissue4","RNA_Tissue5",
                         "RNA_Tissue6","RNA_Tissue7","RNA_Tissue8","RNA_Tissue9","RNA_Tissue10"),
                 sep = ";")
  a1[a1 == ""] <- NA
  a2 <- as.data.frame(lapply(a1, function(x) sapply(strsplit(as.character(x), ":"), "[", 1)))
  a2 <- a2[!duplicated(a2$protein1),]
  rownames(a2) <- a2$protein1
  a2 <- a2[,-1]
  c <- melt(as.matrix(a2))
  c <- na.omit(c)
  colnames(c) <- c("X1","X2","value")
  singlecellgene <- list()
  for (i in names(table(c$value))) {
    b <- dplyr::filter(c,value==i)
    b <- as.character( b$X1)
    singlecellgene[[i]] <- b
  }
  return(singlecellgene)
}
secretomefunc <- function(backprotein){
  hpa_proteinatlas_file <- file.path("data", "reference", "HPA", "proteinatlas.tsv")
  HPA_cite <- read.delim(hpa_proteinatlas_file,
                         sep = "\t",
                         header = TRUE,
                         fill = TRUE)
  HPA_cite <- HPA_cite[,c(1,3,5,19,16,9,59,58,70,71,57,20,21,22,23)]
  HPA_cite$HPAgene <- HPA_cite$Gene
  HPA_sec <- dplyr::filter(HPA_cite,Secretome.location != "")
  HPA_sec <- HPA_sec[,c(1,7,8)]
  HPA_sec_ordered <- HPA_sec %>%
    dplyr::arrange(Secretome.location, Secretome.function) %>%
    dplyr::mutate(Sequence = row_number())
  HPA_sec_ordered[HPA_sec_ordered==""] <- "No annotated function"
  HPA_sec_ordered <- HPA_sec_ordered %>% dplyr::filter(Gene %in% backprotein)
  secregenelo <- list()
  secregenefun <- list()
  for (i in unique(HPA_sec_ordered$Secretome.function)) {
    ge <-  dplyr::filter(HPA_sec_ordered,Secretome.function==i) $ Gene
    secregenefun[[i]] <- as.character(ge)
  }
  for (i in unique(HPA_sec_ordered$Secretome.location)) {
    ge <-  dplyr::filter(HPA_sec_ordered,Secretome.location==i) $ Gene
    secregenelo[[i]] <- as.character(ge)
  }
  list(location = secregenelo, function_sets = secregenefun)
}
calc_fold_enrichment <- function(query_list,ref_list,background,adjust_p = TRUE,min_overlap = 0) {
  if (!is.list(query_list) || is.null(names(query_list))) {
    stop("query_list must be a named list.")
  }
  if (!is.list(ref_list) || is.null(names(ref_list))) {
    stop("ref_list must be a named list.")
  }
  background <- unique(as.character(background))
  background <- background[!is.na(background) & background != ""]
  N <- length(background)
  if (N == 0) {
    stop("background is empty; fold enrichment cannot be calculated.")
  }
  query_list <- lapply(query_list, function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & x != ""]
    intersect(x, background)
  })
  ref_list <- lapply(ref_list, function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & x != ""]
    intersect(x, background)
  })
  res <- list()
  idx <- 1
  for (q_name in names(query_list)) {
    q_genes <- query_list[[q_name]]
    n <- length(q_genes)
    for (r_name in names(ref_list)) {
      r_genes <- ref_list[[r_name]]
      K <- length(r_genes)
      k <- length(intersect(q_genes, r_genes))
      expected <- if (N > 0) n * K / N else NA_real_
      fold_enrichment <- if (!is.na(expected) && expected > 0) {
        k / expected
      } else {
        NA_real_
      }
      gene_ratio <- if (n > 0) k / n else NA_real_
      bg_ratio <- if (N > 0) K / N else NA_real_
      p_value <- if (n > 0 && K > 0 && N > 0) {
        phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
      } else {
        NA_real_
      }
      res[[idx]] <- data.frame(
        query_set = q_name,
        ref_set = r_name,
        query_size = n,
        ref_size = K,
        background_size = N,
        overlap = k,
        expected_overlap = expected,
        gene_ratio = gene_ratio,
        bg_ratio = bg_ratio,
        fold_enrichment = fold_enrichment,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  res <- do.call(rbind, res)
  if (adjust_p) {
    res$padj <- p.adjust(res$p_value, method = "BH")
  }
  if (!is.null(min_overlap) && min_overlap > 0) {
    res <- res[res$overlap >= min_overlap, , drop = FALSE]
  }
  res <- res[order(res$query_set, -res$fold_enrichment, res$p_value), , drop = FALSE]
  rownames(res) <- NULL
  return(res)
}
tsnebind <- function(mo,per,til ){

  set.seed(100)
  inputdata <- t(mo)
  tSNE_fit <-  Rtsne(inputdata,
                     perplexity = per,
                     check_duplicates = FALSE)
  tSNE_df <- tSNE_fit$Y %>%
    as.data.frame()
  tSNE_df$name <- rownames(inputdata)
  tSNE_df <- left_join(tSNE_df,metadata,by=c("name"="OriginalID"))
  color <- paletteer_d("ggthemes::Classic_20", n = 14)
  color <- c("#D4C96E", "#807961","#8DACB6", "#7D839C","#CA9791",
             "#55794D", "#EFF3F9","#924741", "#E0DC96", "#DA9F58",
             "#AC9965", "#C49C94FF", "#E377C2FF", "#F7B6D2FF" )
  point_sizes <- c(
    Day0 = 1.3,
    Day3 = 2.5,
    Day5 = 4
  )
  tSNE_df %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(fill = SubjectID,
                   shape = Sex,
                   size  = Time),
               color="#525A73",
               stroke = 0.5,
               alpha  = 1) +
    scale_fill_manual(values = color,
                      name   = "Subject",
                      guide  = guide_legend(override.aes = list(
                        shape = 21,
                        color = "#525A73",
                        size  = 3
                      ))
    ) +
    scale_shape_manual(values = c(F = 24, M = 21),
                       name   = "Sex") +
    scale_size_manual(
      values = point_sizes,
      name   = "Visit Time"
    ) +
    labs(x     = "tSNE1",
         y     = "tSNE2",
         title = til) +
    theme4_backline +
    themelegendnone
}
overlapheatmap <-  function(inpuli,cate)    {
  n <- length(inpuli)
  overlap_matrix <- matrix(0, nrow = n, ncol = n)
  jaccard_matrix <- matrix(0, nrow = n, ncol = n)
  lens <- sapply(inpuli, length)
  new_names <- paste0(names(inpuli), " (", lens, ")")
  rownames(overlap_matrix) <- new_names
  colnames(overlap_matrix) <- new_names
  rownames(jaccard_matrix) <- new_names
  colnames(jaccard_matrix) <- new_names
  for (i in 1:n) {
    for (j in 1:n) {
      overlap_count <- length(intersect(inpuli[[i]], inpuli[[j]]))
      overlap_matrix[i, j] <- overlap_count
      jaccard_matrix[i, j] <- jaccard_index(inpuli[[i]], inpuli[[j]])
    }
  }
  df_overlap <- melt(overlap_matrix, varnames = c("List1", "List2"), value.name = "Overlap")
  df_jaccard <- melt(jaccard_matrix, varnames = c("List1", "List2"), value.name = "Jaccard")
  df <- merge(df_overlap, df_jaccard, by = c("List1", "List2"))
  ggplot(df, aes(x = List1, y = List2)) +
    geom_tile(aes( fill = Jaccard), color = "white", size = 0.7) +
    geom_text(aes(label = Overlap), size = 2.5, color = "black", vjust = 0.5, hjust = 0.5) +
    scale_size_continuous(range = c(7, 19)) +
    scale_fill_gradient2(name = "Jaccard index",
                         low = "#FBF7E4", mid = "#78BEB4", high = "#A396AD", midpoint = 0.5,  na.value = "grey", guide = "colorbar"
    ) +
    labs(title = cate, x = "", y = "")+themetext90
}
safe_spearman <- function(x, y, min_n = 3) {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  n <- sum(keep)
  if (n < min_n) return(list(rho = NA_real_, p = NA_real_, n = n))
  xk <- x[keep]; yk <- y[keep]
  if (sd(xk) == 0 || sd(yk) == 0) return(list(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(xk, yk, method = "spearman",	use = "pairwise.complete.obs"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = n)
}
secretome_coverage <- function(reference, protein_sets) {
  coverage <- vapply(reference, function(genes) {
    vapply(protein_sets, function(proteins) length(intersect(proteins, genes)) / length(genes), numeric(1))
  }, numeric(length(protein_sets)))
  rownames(coverage) <- paste0(names(protein_sets), " (", lengths(protein_sets), ")")
  colnames(coverage) <- paste0(names(reference), " (", lengths(reference), ")")
  coverage_df <- reshape2::melt(coverage)
  names(coverage_df) <- c("X1", "X2", "value")
  coverage_df$value[coverage_df$value == 0] <- NA_real_
  list(heatmap = ggplot(coverage_df, aes(X1, X2)))
}
