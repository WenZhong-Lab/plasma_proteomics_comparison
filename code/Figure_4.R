# Figure 4: complete analysis and manuscript figures
# Run from the repository root.
source(file.path("code", "setup.R"), encoding = "UTF-8")
wdresult <- paste0(file.path(repo_dir, "result", "Fig4"), "/")
wdplot <- paste0(file.path(repo_dir, "plot", "Fig4"), "/")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

# Within-platform variation
# Run Figure_2.R first to produce the current detection-classification artifact.
# 03. Within-platform variation -----------------------------------------------
plot_prefix <- "03-"

load_registered_assays(wddata)
ot480 <- (diaori1)
astral <- (astori1)


rownames(olink3072)[rownames(olink3072) == "NTproBNP"] <- "NT-proBNP"
rownames(olink3072)[rownames(olink3072) == "WARS"] <- "WARS1"

detection_sets_file <- file.path(repo_dir, "result", "Fig2",
  "01-01_detection_classification_sets.rds"
)
detection_sets <- readRDS(detection_sets_file)
prolist <- detection_sets$detected
individual_prolist <- detection_sets$individual_detected
broad_prolist <- detection_sets$broad_detected

dflistna <- read_processed_matrices(wddata)
names(dflistna)[names(dflistna) == "Olink-3072"] <- "Olink-Explore"
overlap5 <- Reduce(intersect, list(prolist[[1]], prolist[[2]], prolist[[3]], prolist[[4]], prolist[[5]]))

metadata <- read_processed_metadata(wddata)
metadata <- metadata[1:30, ]
metadata$somaid <- paste0("S", metadata$OriginalID)
metadata$Time <- dplyr::recode(metadata$Time, "visit1" = "Day0", "visit2" = "Day3", "visit3" = "Day5")
metadata$Time <- factor(metadata$Time, levels = c("Day0", "Day3", "Day5"))
metadata$SubjectID <- factor(metadata$SubjectID, levels = unique(metadata$SubjectID))

format_p_stars <- function(p) dplyr::case_when( is.na(p) ~ NA_character_, p < 0.001 ~ "***", p < 0.01 ~ "**",
  p < 0.05 ~ "*", TRUE ~ ""
)

soma_cv_result_to_protein <- function(soma_cv_result, strata_cols, plot_cols) {
  somer_vol <- soma_cv_result[["vol"]]
  protein_audit <- summarise_soma_cv_to_protein( cv_df = somer_vol, strata_cols = strata_cols
  )

  missing_plot_cols <- setdiff(plot_cols, colnames(protein_audit))
  list( vol = protein_audit %>% dplyr::select(dplyr::all_of(plot_cols)), protein_audit = protein_audit,
    somer_vol = somer_vol, somer_bar = soma_cv_result[["bar"]]
  )
}
# 01. Calculate within-individual CV ------------------------------------------

full_iv_msot <- mscv(dflistna, cate = "ot480", pro = prolist[["MS-OT480"]], isna = "NA", cate1input = "Detected_MS-OT480")
full_iv_as <- mscv(dflistna, cate = "as", pro = prolist[["MS-Astral"]], isna = "NA", cate1input = "Detected_MS-Astral")
full_iv_o3072 <- olinkcv( cate = "3072", pro = prolist[["Olink-Explore"]], isna = FALSE,
  cate1input = "Detected_Olink-Explore"
)
full_iv_oht <- olinkcv(cate = "ht", pro = prolist[["Olink-HT"]], isna = FALSE, cate1input = "Detected_Olink-HT")
full_iv_soma <- somacvf(pro = prolist[["Somascan"]], isna = FALSE, cate1input = "Detected_Somascan") %>%
  soma_cv_result_to_protein( strata_cols = c("SubjectID", "cate", "cate1"),
  plot_cols = c("SubjectID", "protein", "CV", "cate", "cate1")
)

full_iv_cv_df <- rbind(full_iv_msot[[1]], full_iv_as[[1]], full_iv_o3072[[1]], full_iv_oht[[1]], full_iv_soma[[1]])
full_iv_cv_df$cate <- factor(full_iv_cv_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT",
  "Somascan"
))
full_iv_cv_df$name <- "Within-individual CV"

# 02. Calculate total CV ------------------------------------------------------

full_tot_msot <- mscvall(dflistna, cate = "ot480", pro = prolist[["MS-OT480"]], isna = "NA",
  cate1input = "Detected_MS-OT480"
)
full_tot_as <- mscvall(dflistna, cate = "as", pro = prolist[["MS-Astral"]], isna = "NA", cate1input = "Detected_MS-Astral")
full_tot_o3072 <- olinkcvall( cate = "3072", pro = prolist[["Olink-Explore"]], isna = FALSE,
  cate1input = "Detected_Olink-Explore"
)
full_tot_oht <- olinkcvall(cate = "ht", pro = prolist[["Olink-HT"]], isna = FALSE, cate1input = "Detected_Olink-HT")
full_tot_soma <- somacvfall(pro = prolist[["Somascan"]], isna = FALSE, cate1input = "Detected_Somascan") %>%
  soma_cv_result_to_protein( strata_cols = c("cate", "cate1"), plot_cols = c("protein", "CV", "cate", "cate1")
)

full_tot_cv_df <- rbind( full_tot_msot[[1]], full_tot_as[[1]], full_tot_o3072[[1]], full_tot_oht[[1]], full_tot_soma[[1]]
)
full_tot_cv_df$cate <- factor(full_tot_cv_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore",
  "Olink-HT", "Somascan"
))
full_tot_cv_df$name <- "Total CV"

# 03. Calculate inter-individual CV -------------------------------------------

full_ii_msot <- mscvindi(dflistna, cate = "ot480", pro = prolist[["MS-OT480"]], isna = "NA",
  cate1input = "Detected_MS-OT480"
)
full_ii_as <- mscvindi(dflistna, cate = "as", pro = prolist[["MS-Astral"]], isna = "NA", cate1input = "Detected_MS-Astral")
full_ii_o3072 <- olinkcvindi( cate = "3072", pro = prolist[["Olink-Explore"]], isna = FALSE,
  cate1input = "Detected_Olink-Explore"
)
full_ii_oht <- olinkcvindi(cate = "ht", pro = prolist[["Olink-HT"]], isna = FALSE, cate1input = "Detected_Olink-HT")
full_ii_soma <- somacvindi(pro = prolist[["Somascan"]], isna = FALSE, cate1input = "Detected_Somascan") %>%
  soma_cv_result_to_protein( strata_cols = c("Time", "cate", "cate1"),
  plot_cols = c("Time", "protein", "CV", "cate", "cate1")
)

full_ii_cv_df <- rbind(full_ii_msot[[1]], full_ii_as[[1]], full_ii_o3072[[1]], full_ii_oht[[1]], full_ii_soma[[1]])
full_ii_cv_df$cate <- factor(full_ii_cv_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT",
  "Somascan"
))
full_ii_cv_df$name <- "Inter-individual CV"

full_soma_cv_protein_audit <- dplyr::bind_rows( dplyr::mutate(full_iv_soma$protein_audit, name = "Within-individual CV"),
  dplyr::mutate(full_tot_soma$protein_audit, name = "Total CV"),
  dplyr::mutate(full_ii_soma$protein_audit, name = "Inter-individual CV")
)

# 04. Calculate CV for the five-platform shared detected proteins -------------

ov_pro <- Reduce(intersect, list(prolist[[1]], prolist[[2]], prolist[[4]], prolist[[3]], prolist[[5]]))
cv_col <- c( "#4E79A7FF", "#F28E2BFF", "#F1CE63FF", "#59A14FFF", "#499894FF", "#86BCB6FF", "#B07AA1FF",
  "#9D7660FF", "#D7B5A6FF", "#79706EFF"
)
# 04.1 Within-individual CV ---------------------------------------------------

ov_iv_msot <- mscv(dflistna, cate = "ot480", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-OT480")
ov_iv_as <- mscv(dflistna, cate = "as", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-Astral")
ov_iv_oe <- olinkcv(cate = "3072", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-Explore")
ov_iv_oht <- olinkcv(cate = "ht", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-HT")
ov_iv_soma <- somacvf(pro = ov_pro, isna = FALSE, cate1input = "Detected_Somascan") %>%
  soma_cv_result_to_protein( strata_cols = c("SubjectID", "cate", "cate1"),
  plot_cols = c("SubjectID", "protein", "CV", "cate", "cate1")
)

ov_iv_df <- rbind(ov_iv_msot[[1]], ov_iv_as[[1]], ov_iv_oe[[1]], ov_iv_oht[[1]], ov_iv_soma[[1]])
ov_iv_df$cate <- factor(ov_iv_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan"
))
ov_iv_df$name <- "Within-individual CV"

ov_iv_aov_df <- ov_iv_df
ov_iv_aov_df$cate <- gsub("-", "_", ov_iv_aov_df$cate)
ov_iv_ano <- mutl_anova_test(ov_iv_aov_df, c("CV"), c("cate", "SubjectID"))
ov_iv_ano <- ov_iv_ano %>% dplyr::mutate( `p adj` = suppressWarnings(as.numeric(`p adj`)), mark = case_when(
  `p adj` <= 0.001 ~ "***", `p adj` <= 0.01 ~ "**", `p adj` <= 0.05 ~ "*", TRUE ~ NA_character_
)
)
ov_iv_ano <- ov_iv_ano[-c(1:55), ] %>% dplyr::filter(X1 == X3)

# 04.2 Total CV ---------------------------------------------------------------

ov_tot_msot <- mscvall(dflistna, cate = "ot480", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-OT480")
ov_tot_as <- mscvall(dflistna, cate = "as", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-Astral")
ov_tot_oe <- olinkcvall(cate = "3072", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-Explore")
ov_tot_oht <- olinkcvall(cate = "ht", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-HT")
ov_tot_soma <- somacvfall(pro = ov_pro, isna = FALSE, cate1input = "Detected_Somascan") %>%
  soma_cv_result_to_protein( strata_cols = c("cate", "cate1"), plot_cols = c("protein", "CV", "cate", "cate1")
)

ov_tot_df <- rbind(ov_tot_msot[[1]], ov_tot_as[[1]], ov_tot_oe[[1]], ov_tot_oht[[1]], ov_tot_soma[[1]])
ov_tot_df$cate <- factor(ov_tot_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan"
))
ov_tot_df$name <- "Total CV"

# 04.3 Inter-individual CV ----------------------------------------------------

ov_ii_msot <- mscvindi(dflistna, cate = "ot480", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-OT480")
ov_ii_as <- mscvindi(dflistna, cate = "as", pro = ov_pro, isna = "NA", cate1input = "Detected_MS-Astral")
ov_ii_oe <- olinkcvindi(cate = "3072", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-Explore")
ov_ii_oht <- olinkcvindi(cate = "ht", pro = ov_pro, isna = FALSE, cate1input = "Detected_Olink-HT")
ov_ii_soma <- somacvindi(pro = ov_pro, isna = FALSE, cate1input = "Detected_Somascan") %>% soma_cv_result_to_protein(
  strata_cols = c("Time", "cate", "cate1"), plot_cols = c("Time", "protein", "CV", "cate", "cate1")
)

ov_ii_df <- rbind(ov_ii_msot[[1]], ov_ii_as[[1]], ov_ii_oe[[1]], ov_ii_oht[[1]], ov_ii_soma[[1]])
ov_ii_df$cate <- factor(ov_ii_df$cate, levels = c( "MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan"
))
ov_ii_df$Time <- factor(ov_ii_df$Time, levels = c("Day0", "Day3", "Day5"))
ov_ii_df$name <- "Inter-individual CV"

overlap_soma_cv_protein_audit <- dplyr::bind_rows( dplyr::mutate(ov_iv_soma$protein_audit, name = "Within-individual CV"),
  dplyr::mutate(ov_tot_soma$protein_audit, name = "Total CV"),
  dplyr::mutate(ov_ii_soma$protein_audit, name = "Inter-individual CV")
)

# 06. Prepare CV comparison tables --------------------------------------------

cvbind_plot_w_cm <- 33
cvbind_plot_h_cm <- 10
cvbind_fill_cols <- c( "Total CV" = "#86BCB6FF", "Within-individual CV" = "#B07AA1FF", "Inter-individual CV" = "#D7B5A6FF"
)

detected_cv_df <- dplyr::bind_rows(full_ii_cv_df, full_tot_cv_df, full_iv_cv_df) %>%
  dplyr::mutate(panel = "Detected proteins")
detected_overlap_cv_df <- dplyr::bind_rows(ov_ii_df, ov_tot_df, ov_iv_df) %>%
  dplyr::mutate(panel = "Detected protein overlap")

belowlod_oe <- setdiff(rownames(olink3072), prolist[["Olink-Explore"]])
belowlod_oht <- setdiff(rownames(olinkht), prolist[["Olink-HT"]])

bl_iv_oe <- olinkcv( cate = "3072", pro = belowlod_oe, isna = FALSE, keepfilter = FALSE,
  cate1input = "BelowLOD_Olink-Explore"
)[[1]]
bl_iv_oht <- olinkcv( cate = "ht", pro = belowlod_oht, isna = FALSE, keepfilter = FALSE, cate1input = "BelowLOD_Olink-HT"
)[[1]]
bl_iv_df <- rbind(bl_iv_oe, bl_iv_oht)
bl_iv_df$name <- "Within-individual CV"

bl_tot_oe <- olinkcvall( cate = "3072", pro = belowlod_oe, isna = FALSE, keepfilter = FALSE,
  cate1input = "BelowLOD_Olink-Explore"
)[[1]]
bl_tot_oht <- olinkcvall( cate = "ht", pro = belowlod_oht, isna = FALSE, keepfilter = FALSE,
  cate1input = "BelowLOD_Olink-HT"
)[[1]]
bl_tot_df <- rbind(bl_tot_oe, bl_tot_oht)
bl_tot_df$name <- "Total CV"

bl_ii_oe <- olinkcvindi( cate = "3072", pro = belowlod_oe, isna = FALSE, keepfilter = FALSE,
  cate1input = "BelowLOD_Olink-Explore"
)[[1]]
bl_ii_oht <- olinkcvindi( cate = "ht", pro = belowlod_oht, isna = FALSE, keepfilter = FALSE,
  cate1input = "BelowLOD_Olink-HT"
)[[1]]
bl_ii_df <- rbind(bl_ii_oe, bl_ii_oht)
bl_ii_df$name <- "Inter-individual CV"

undetected_cv_df <- dplyr::bind_rows(bl_ii_df, bl_tot_df, bl_iv_df) %>%
  dplyr::mutate(panel = "Undetected predefined protein")
# 07. Calculate CV for additional detection sets ------------------------------

build_cv_bundle <- function(protein_input, set_label) {
  pro_for <- function(platform) if (is.list(protein_input)) {
    protein_input[[platform]]
  } else {
    protein_input
  }
  iv_msot <- mscv(dflistna, cate = "ot480", pro = pro_for("MS-OT480"), isna = "NA",
    cate1input = paste0(set_label, "_MS-OT480")
  )
  iv_as <- mscv(dflistna, cate = "as", pro = pro_for("MS-Astral"), isna = "NA", cate1input = paste0(set_label, "_MS-Astral")
  )
  iv_oe <- olinkcv( cate = "3072", pro = pro_for("Olink-Explore"), isna = FALSE,
    cate1input = paste0(set_label, "_Olink-Explore")
  )
  iv_oht <- olinkcv(cate = "ht", pro = pro_for("Olink-HT"), isna = FALSE, cate1input = paste0( set_label, "_Olink-HT"
  ))
  iv_soma <- somacvf(pro = pro_for("Somascan"), isna = FALSE, cate1input = paste0(set_label, "_Somascan")) %>%
    soma_cv_result_to_protein( strata_cols = c("SubjectID", "cate", "cate1"),
    plot_cols = c("SubjectID", "protein", "CV", "cate", "cate1")
  )

  tot_msot <- mscvall(dflistna, cate = "ot480", pro = pro_for("MS-OT480"), isna = "NA",
    cate1input = paste0(set_label, "_MS-OT480")
  )
  tot_as <- mscvall(dflistna, cate = "as", pro = pro_for("MS-Astral"), isna = "NA",
    cate1input = paste0(set_label, "_MS-Astral")
  )
  tot_oe <- olinkcvall( cate = "3072", pro = pro_for("Olink-Explore"), isna = FALSE,
    cate1input = paste0(set_label, "_Olink-Explore")
  )
  tot_oht <- olinkcvall( cate = "ht", pro = pro_for("Olink-HT"), isna = FALSE, cate1input = paste0(set_label, "_Olink-HT")
  )
  tot_soma <- somacvfall(pro = pro_for("Somascan"), isna = FALSE, cate1input = paste0( set_label, "_Somascan"
  )) %>% soma_cv_result_to_protein( strata_cols = c("cate", "cate1"), plot_cols = c("protein", "CV", "cate", "cate1")
  )

  ii_msot <- mscvindi(dflistna, cate = "ot480", pro = pro_for("MS-OT480"), isna = "NA",
    cate1input = paste0(set_label, "_MS-OT480")
  )
  ii_as <- mscvindi(dflistna, cate = "as", pro = pro_for("MS-Astral"), isna = "NA",
    cate1input = paste0(set_label, "_MS-Astral")
  )
  ii_oe <- olinkcvindi( cate = "3072", pro = pro_for("Olink-Explore"), isna = FALSE,
    cate1input = paste0(set_label, "_Olink-Explore")
  )
  ii_oht <- olinkcvindi( cate = "ht", pro = pro_for("Olink-HT"), isna = FALSE, cate1input = paste0(set_label, "_Olink-HT")
  )
  ii_soma <- somacvindi(pro = pro_for("Somascan"), isna = FALSE, cate1input = paste0( set_label, "_Somascan"
  )) %>% soma_cv_result_to_protein( strata_cols = c("Time", "cate", "cate1"),
    plot_cols = c("Time", "protein", "CV", "cate", "cate1")
  )

  iv <- dplyr::bind_rows( iv_msot[[1]], iv_as[[1]], iv_oe[[1]], iv_oht[[1]], iv_soma[[1]]
  ) %>% dplyr::mutate(name = "Within-individual CV")
  tot <- dplyr::bind_rows( tot_msot[[1]], tot_as[[1]], tot_oe[[1]], tot_oht[[1]], tot_soma[[1]]
  ) %>% dplyr::mutate(name = "Total CV")
  ii <- dplyr::bind_rows( ii_msot[[1]], ii_as[[1]], ii_oe[[1]], ii_oht[[1]], ii_soma[[1]]
  ) %>% dplyr::mutate(name = "Inter-individual CV")
  data <- dplyr::bind_rows(ii, tot, iv)
  data$cate <- factor(data$cate, levels = c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan"))
  data$name <- factor(data$name, levels = c("Total CV", "Within-individual CV", "Inter-individual CV"))
  audit <- dplyr::bind_rows( dplyr::mutate(iv_soma$protein_audit, name = "Within-individual CV"),
    dplyr::mutate(tot_soma$protein_audit, name = "Total CV"),
    dplyr::mutate(ii_soma$protein_audit, name = "Inter-individual CV")
  )
  list(data = data, audit = audit)
}
individual_detected_overlap <- Reduce(intersect, individual_prolist)
broad_detected_overlap <- Reduce(intersect, broad_prolist)

individual_full_cv <- build_cv_bundle(individual_prolist, "IndividualDetected")
broad_full_cv <- build_cv_bundle(broad_prolist, "BroadlyDetected")
individual_overlap_cv <- build_cv_bundle(individual_detected_overlap, "IndividualDetectedOverlap")
broad_overlap_cv <- build_cv_bundle(broad_detected_overlap, "BroadlyDetectedOverlap")

# 03-02a. CV distributions across platform-specific and shared protein sets ---

cv_comparison_tables <- list( platform = dplyr::bind_rows(detected_cv_df,
  individual_full_cv$data %>% dplyr::mutate(panel = "Individual detected proteins"),
  broad_full_cv$data %>% dplyr::mutate(panel = "Broad detected proteins"), undetected_cv_df),
  overlap = dplyr::bind_rows(detected_overlap_cv_df,
  individual_overlap_cv$data %>% dplyr::mutate(panel = "Individual detected protein overlap"),
  broad_overlap_cv$data %>% dplyr::mutate(panel = "Broad detected protein overlap"))
)
cv_comparison_panel_levels <- list(
  platform = c("Detected proteins", "Individual detected proteins", "Broad detected proteins",
  "Undetected predefined protein"),
  overlap = c("Detected protein overlap", "Individual detected protein overlap", "Broad detected protein overlap")
)
cv_comparison_files <- c(platform = "F4a_cv_combined_by_platform.pdf", overlap = "S3b_detection_overlap_cv_by_platform.pdf")
cv_comparison_plots <- list()
for (comparison_name in names(cv_comparison_tables)) {
  plot_df <- cv_comparison_tables[[comparison_name]] %>% dplyr::mutate(
    panel = factor(panel, levels = cv_comparison_panel_levels[[comparison_name]]),
    cate = factor(cate, levels = c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")),
    name = factor(name, levels = names(cvbind_fill_cols)))
  median_df <- plot_df %>% dplyr::group_by(cate, panel, name) %>%
    dplyr::summarise(median_label = sprintf("%.2f", median(abs(CV), na.rm = TRUE) * 100), .groups = "drop")

  violin_df <- plot_df %>% dplyr::filter(CV * 100 <= 200)
  p <- ggplot(plot_df, aes(x = panel, y = CV * 100, fill = name)) +
    geom_violin(data = violin_df, trim = FALSE, colour = "transparent", alpha = 1, width = 0.9,
    position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.12, size = 0.25, colour = "white", outlier.shape = NA, alpha = 0.5,
    position = position_dodge(width = 0.8, preserve = "single")) +
    geom_text(data = median_df, aes(x = panel, y = 180, label = median_label, group = name),
    position = position_dodge2(width = 0.8, preserve = "single"), colour = "black", size = 2.8,
    vjust = 0, angle = 45, hjust = 0) + coord_cartesian(ylim = c(0, 200), clip = "off") +
    scale_fill_manual(values = cvbind_fill_cols, breaks = names(cvbind_fill_cols), name = "") +
    facet_wrap(~cate, nrow = 1, scales = "free_x") +
    labs(x = "", y = "CV (%)", fill = "") + theme2_nobacknoerror + themetext30 + themelegendtop +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), legend.text = element_text(size = 10), strip.text = element_text(size = 10))
  cv_comparison_plots[[comparison_name]] <- p
  ggsave(paste0(wdplot, cv_comparison_files[[comparison_name]]), pdf_text_grob(p),
    width = cvbind_plot_w_cm, height = cvbind_plot_h_cm, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
}
p_cv_combined_by_platform <- cv_comparison_plots$platform
p_detection_overlap_cv_by_platform <- cv_comparison_plots$overlap

# 03-02. Compare within- and inter-individual CV ------------------------------

cv_scatter_platform_order <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
cv_scatter_platform_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral", "Olink-Explore" = "Olink-Explore",
  "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan"
)
cv_scatter_category_levels <- c( "Undetected predefined protein", "Detected proteins", "Individual detected proteins",
  "Broadly detected proteins"
)
cv_scatter_colour <- c( "Undetected predefined protein" = "#D3D3D3", "Detected proteins" = "#475272",
  "Individual detected proteins" = "#AEC7E8FF", "Broadly detected proteins" = "#D28770"
)
cv_scatter_category_labels <- c( "Undetected predefined protein" = "Undetected predefined protein",
  "Detected proteins" = "Detected only", "Individual detected proteins" = "Individual detected only",
  "Broadly detected proteins" = "Broadly detected"
)

cv_scatter_membership_list <- vector("list", length(cv_scatter_platform_order))
for (platform_index in seq_along(cv_scatter_platform_order)) {
  platform <- cv_scatter_platform_order[[platform_index]]
  detected_set <- prolist[[platform]]
  individual_set <- individual_prolist[[platform]]
  broad_set <- broad_prolist[[platform]]
  detected_membership_df <- data.frame( cate = platform, protein = detected_set,
    protein_detection_category = dplyr::case_when( detected_set %in% broad_set ~ "Broadly detected proteins",
    detected_set %in% individual_set ~ "Individual detected proteins", TRUE ~ "Detected proteins"
  ), stringsAsFactors = FALSE
  )

  undetected_set <- switch(platform, "Olink-Explore" = belowlod_oe, "Olink-HT" = belowlod_oht, character()
  )
  undetected_membership_df <- if (length(undetected_set) > 0) {
    data.frame( cate = platform, protein = undetected_set, protein_detection_category = "Undetected predefined protein",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame( cate = character(), protein = character(), protein_detection_category = character(), stringsAsFactors = FALSE
    )
  }
  cv_scatter_membership_list[[platform_index]] <- dplyr::bind_rows( detected_membership_df, undetected_membership_df
  )
}
cv_scatter_membership_df <- dplyr::bind_rows(cv_scatter_membership_list) %>%
  dplyr::mutate( cate = factor(cate, levels = cv_scatter_platform_order), protein_detection_category = factor(
  protein_detection_category, levels = cv_scatter_category_levels
)
)

summarise_cv_for_scatter <- function(cv_df, group_col, output_col) {
  cv_df %>% dplyr::filter(!is.na(protein), is.finite(CV)) %>%
    dplyr::mutate(cate = as.character(cate), cv_percent = abs(CV) * 100) %>% dplyr::group_by(cate, protein) %>%
    dplyr::summarise( "{output_col}" := stats::median(cv_percent, na.rm = TRUE), "{group_col}" := dplyr::n(),
    .groups = "drop"
  )
}
cv_scatter_within_df <- summarise_cv_for_scatter( dplyr::bind_rows(full_iv_cv_df, bl_iv_df), group_col = "n_individuals",
  output_col = "within_individual_cv"
)
cv_scatter_inter_df <- summarise_cv_for_scatter( dplyr::bind_rows(full_ii_cv_df, bl_ii_df), group_col = "n_sampling_times",
  output_col = "inter_individual_cv"
)
cv_scatter_plot_df <- cv_scatter_membership_df %>% dplyr::left_join(cv_scatter_within_df, by = c("cate", "protein")) %>%
  dplyr::left_join(cv_scatter_inter_df, by = c("cate", "protein")) %>% dplyr::filter( is.finite(within_individual_cv),
  is.finite(inter_individual_cv)
)

cv_scatter_plots <- vector("list", length(cv_scatter_platform_order))
for (platform_index in seq_along(cv_scatter_platform_order)) {
  platform <- cv_scatter_platform_order[[platform_index]]
  show_legend <- identical(platform, "Olink-Explore")
  platform_df <- dplyr::filter(cv_scatter_plot_df, as.character(cate) == platform)
  point_size <- if (platform %in% c("Olink-Explore", "Olink-HT", "Somascan")) {
    0.55
  } else {
    0.7
  }
  is_olink_platform <- platform %in% c("Olink-Explore", "Olink-HT")
  if (is_olink_platform) {
    set.seed(if (identical(platform, "Olink-Explore")) {
      103L
    } else {
      104L
    })
    point_df <- platform_df[sample.int(nrow(platform_df)), , drop = FALSE]
  } else {
    point_df <- platform_df
  }
  colour_scale <- scale_colour_manual( values = cv_scatter_colour, breaks = cv_scatter_category_levels,
    labels = cv_scatter_category_labels, name = "Protein detection category\nMutually exclusive categories"
  )

  p_density <- ggplot( platform_df, aes(x = within_individual_cv, colour = protein_detection_category)
  ) + geom_density(linewidth = 0.5, adjust = 0.9, na.rm = TRUE, show.legend = FALSE) +
    coord_cartesian(xlim = c(0, 200), expand = FALSE) + colour_scale + labs(
    title = unname(cv_scatter_platform_labels[[platform]]), x = NULL, y = "Density"
  ) + theme2_nobacknoerror + theme( plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 6),
    plot.margin = margin(3, 3, 0, 3, unit = "mm")
  )

  point_layers <- if (is_olink_platform) {
    list( geom_point( data = point_df, aes(alpha = protein_detection_category), size = point_size
    ), scale_alpha_manual( values = c( "Undetected predefined protein" = 1, "Detected proteins" = 1,
      "Individual detected proteins" = 1, "Broadly detected proteins" = 0.35
    ), guide = "none"
    )
    )
  } else {
    list( geom_point( data = dplyr::filter(point_df, protein_detection_category == "Broadly detected proteins"),
      size = point_size, alpha = 0.35
    ), geom_point( data = dplyr::filter(point_df, protein_detection_category != "Broadly detected proteins"),
      size = point_size, alpha = 1
    )
    )
  }
  p_scatter <- ggplot( platform_df, aes( x = within_individual_cv, y = inter_individual_cv,
    colour = protein_detection_category
  )
  ) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.3, linetype = "dashed") + point_layers +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) +
    colour_scale + guides(colour = guide_legend(ncol = 1)) + labs( x = "Within-individual CV (%)",
    y = "Inter-individual CV (%)"
  ) +
    theme4_backline + theme( axis.title = element_text(size = 8), axis.text = element_text(size = 7),
    legend.position = if (show_legend) {
    "bottom"
  } else {
    "none"
  }, legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.height = grid::unit(3.5, "mm"),
    plot.margin = margin(0, 3, 3, 3, unit = "mm")
  )
  cv_scatter_plots[[platform_index]] <- p_density / p_scatter + patchwork::plot_layout(heights = c(1.1, 4))
}
p_cv_within_vs_inter_by_platform <- patchwork::wrap_plots(cv_scatter_plots, nrow = 1)

ggsave( paste0(wdplot, "F4c_within_vs_inter_individual_cv_scatter_by_platform.pdf"), pdf_text_grob(p_cv_within_vs_inter_by_platform),
  width = 37, height = 12, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)
# 03-09. Compare detection categories within each platform --------------------

cv_category_levels <- c( "Detected proteins", "Individual detected proteins", "Broadly detected proteins",
  "Undetected predefined protein"
)
cv_category_colours <- c( "Detected proteins" = "#475272", "Individual detected proteins" = "#AEC7E8FF",
  "Broadly detected proteins" = "#D28770", "Undetected predefined protein" = "#D3D3D3"
)
cv_category_platform_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral", "Olink-Explore" = "Olink-Explore",
  "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan"
)

cv_category_membership_df <- dplyr::bind_rows(lapply(names(prolist), function(platform) {
  detected_set <- prolist[[platform]]
  individual_set <- individual_prolist[[platform]]
  broad_set <- broad_prolist[[platform]]
  dplyr::bind_rows( data.frame( cate = platform, protein = detected_set,
    detection_category = "Detected proteins", stringsAsFactors = FALSE
  ), data.frame( cate = platform, protein = individual_set,
    detection_category = "Individual detected proteins", stringsAsFactors = FALSE
  ), data.frame( cate = platform, protein = broad_set,
    detection_category = "Broadly detected proteins", stringsAsFactors = FALSE
  )
  )
}))

build_cv_category_data <- function(detected_cv_df, undetected_cv_df, cv_name) {
  detected_category_df <- detected_cv_df %>% dplyr::mutate(cate = as.character(cate)) %>% dplyr::inner_join(
    cv_category_membership_df, by = c("cate", "protein"), relationship = "many-to-many"
  )
  undetected_category_df <- undetected_cv_df %>% dplyr::mutate( cate = as.character(cate),
    detection_category = "Undetected predefined protein"
  )
  dplyr::bind_rows(detected_category_df, undetected_category_df) %>% dplyr::mutate( cv_type = cv_name,
    cate = dplyr::recode(cate, !!!cv_category_platform_labels),
    detection_category = factor(detection_category, levels = cv_category_levels)
  )
}
cv_category_plot_df <- dplyr::bind_rows( build_cv_category_data(full_tot_cv_df, bl_tot_df, "Total CV"),
  build_cv_category_data(full_iv_cv_df, bl_iv_df, "Within-individual CV"),
  build_cv_category_data(full_ii_cv_df, bl_ii_df, "Inter-individual CV")
) %>% dplyr::filter(!is.na(CV), is.finite(CV)) %>% dplyr::mutate(
  cate = factor(cate, levels = unname(cv_category_platform_labels)),
  cv_type = factor(cv_type, levels = c("Total CV", "Within-individual CV", "Inter-individual CV"))
)

cv_category_pairwise_tests <- cv_category_plot_df %>%
  dplyr::group_by(cv_type, cate) %>%
  dplyr::group_modify(function(.x, .y) {
  present_categories <- cv_category_levels[cv_category_levels %in% as.character(unique(.x$detection_category))]
  if (length(present_categories) < 2) {
    return(data.frame())
  }
  pair_list <- combn(present_categories, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pair_list, function(pair) {
    pair_df <- .x %>%
      dplyr::filter(detection_category %in% pair) %>% dplyr::mutate(detection_category = droplevels(detection_category))
    raw_p <- tryCatch( stats::wilcox.test(CV ~ detection_category, data = pair_df, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    data.frame( category_1 = pair[1], category_2 = pair[2], n_category_1 = sum(pair_df$detection_category == pair[1]),
      n_category_2 = sum(pair_df$detection_category == pair[2]),
      median_cv_category_1 = stats::median(pair_df$CV[pair_df$detection_category == pair[1]]) * 100,
      median_cv_category_2 = stats::median(pair_df$CV[pair_df$detection_category == pair[2]]) * 100, raw_p = raw_p,
      test = "Wilcoxon rank-sum test (unadjusted)", stringsAsFactors = FALSE
    )
  }))
}) %>%
  dplyr::ungroup() %>% dplyr::mutate( fdr_bh = stats::p.adjust(raw_p, method = "BH"),
  fdr_family = "All 63 pairwise comparisons shown in 03-09"
)

plot_cv_category_metric <- function(metric_name, significance_column, y_cap = 200) {
  metric_df <- cv_category_plot_df %>% dplyr::filter(cv_type == metric_name)
  display_df <- metric_df %>% dplyr::filter(CV * 100 <= y_cap)
  test_df <- cv_category_pairwise_tests %>% dplyr::filter(cv_type == metric_name) %>% dplyr::group_by(cate) %>%
    dplyr::mutate( comparison_index = dplyr::row_number(), x1 = match(category_1, cv_category_levels),
    x2 = match(category_2, cv_category_levels), xmid = (x1 + x2) / 2, y = y_cap * (0.58 + 0.06 * comparison_index),
    significant = !is.na(.data[[significance_column]]) & .data[[significance_column]] < 0.05,
    significance_label = ifelse(significant, "*", "")
  ) %>% dplyr::ungroup()
  tick_height <- y_cap * 0.015

  p <- ggplot(display_df, aes(x = detection_category, y = CV * 100, fill = detection_category)) +
    geom_violin(trim = FALSE, colour = "transparent", width = 0.85, na.rm = TRUE) +
    geom_boxplot(width = 0.13, size = 0.25, colour = "white", outlier.shape = NA, alpha = 0.55, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, y_cap), clip = "off") +
    scale_fill_manual(values = cv_category_colours, drop = FALSE, name = "") + scale_x_discrete(drop = FALSE) +
    labs( x = "Protein detection set", y = "CV (%)", title = metric_name
  ) + facet_wrap(~cate, nrow = 1) + theme2_nobacknoerror + themetext30 + themelegendtop + theme(
    axis.text.x = element_text(size = 7, angle = 25, hjust = 1), plot.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10)
  )

  if (nrow(test_df) > 0) {
    p <- p + geom_segment( data = test_df, inherit.aes = FALSE, linewidth = 0.15, aes(x = x1, xend = x2, y = y, yend = y)
    ) + geom_segment( data = test_df, inherit.aes = FALSE, linewidth = 0.15,
      aes(x = x1, xend = x1, y = y - tick_height, yend = y)
    ) + geom_segment( data = test_df, inherit.aes = FALSE, linewidth = 0.15,
      aes(x = x2, xend = x2, y = y - tick_height, yend = y)
    ) +
      geom_text( data = dplyr::filter(test_df, significant), inherit.aes = FALSE,
      size = 4, vjust = 0, aes(x = xmid, y = y - tick_height * 0.60, label = significance_label)
    )
  }

}
p_cv_category_comparison_raw_p <- plot_cv_category_metric("Total CV", "raw_p") /
  plot_cv_category_metric("Within-individual CV", "raw_p") /
  plot_cv_category_metric("Inter-individual CV", "raw_p") +
  patchwork::plot_annotation( caption = paste(
  "All available CV observations are retained (within-individual: protein × individual;",
  "inter-individual: protein × sampling time). Detection sets are nested; pairwise Wilcoxon",
  "rank-sum tests are descriptive; * indicates raw P < 0.05 without multiple-testing correction."
)
)

ggsave( paste0(wdplot, "S3a_within_platform_cv_by_detection_category.pdf"), pdf_text_grob(p_cv_category_comparison_raw_p), width = 42,
  height = 30, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)
plot_prefix <- "04.2-"
platform_order <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
target_platforms <- platform_order[1:4]
platform_display_labels <- setNames(platform_order, platform_order)
platform_display_labels[["Somascan"]] <- "SomaScan"

scatter_display_limits <- c(0, 200)
scatter_height_cm <- 12
category_levels <- c("Detected proteins", "Individual detected proteins", "Broadly detected proteins")
category_colours <- c("Detected proteins" = "#475272", "Individual detected proteins" = "#AEC7E8FF",
  "Broadly detected proteins" = "#D28770")
category_display_labels <- c("Detected proteins" = "Detected only",
  "Individual detected proteins" = "Individual detected only", "Broadly detected proteins" = "Broad detected")

cv_payload <- list(full_iv_cv_df = full_iv_cv_df, full_ii_cv_df = full_ii_cv_df,
  full_iv_soma = full_iv_soma, full_ii_soma = full_ii_soma, prolist = prolist,
  individual_prolist = individual_prolist, broad_prolist = broad_prolist)
detected_sets <- cv_payload$prolist
individual_sets <- cv_payload$individual_prolist
broad_sets <- cv_payload$broad_prolist
for (platform in platform_order) {
  for (set_list in list(detected_sets, individual_sets, broad_sets)) {
    values <- set_list[[platform]]
  }
}
for (audit_name in c("full_iv_soma", "full_ii_soma")) {
  audit <- cv_payload[[audit_name]]$protein_audit
}
# 03. Summarise within- and inter-individual CV -------------------------------

summarise_existing_protein_cv <- function(cv_df, grouping_column, value_column, count_column) {
  required_cols <- c("cate", "protein", "CV", grouping_column)
  finite_df <- cv_df %>%
    dplyr::filter(!is.na(protein), is.finite(CV)) %>%
    dplyr::mutate(cate = as.character(cate), protein = as.character(protein))
  finite_df %>% dplyr::mutate(cv_percent = abs(CV) * 100) %>% dplyr::group_by(cate, protein) %>%
    dplyr::summarise("{value_column}" := stats::median(cv_percent), "{count_column}" := dplyr::n(), .groups = "drop")
}
within_cv_df <- summarise_existing_protein_cv(cv_payload$full_iv_cv_df, "SubjectID", "within_individual_cv", "n_individuals")
inter_cv_df <- summarise_existing_protein_cv(cv_payload$full_ii_cv_df, "Time", "inter_individual_cv", "n_sampling_times")
detected_membership_df <- dplyr::bind_rows(lapply(platform_order, function(platform) {
  protein_keys <- detected_sets[[platform]]
  data.frame(cate = rep(platform, length(protein_keys)), protein = protein_keys,
    protein_detection_category = dplyr::case_when( protein_keys %in% broad_sets[[platform]] ~ category_levels[3],
    protein_keys %in% individual_sets[[platform]] ~ category_levels[2], TRUE ~ category_levels[1]
  ))
}))
protein_cv_df <- detected_membership_df %>%
  dplyr::left_join(within_cv_df, by = c("cate", "protein"), relationship = "one-to-one") %>%
  dplyr::left_join(inter_cv_df, by = c("cate", "protein"), relationship = "one-to-one") %>%
  dplyr::mutate( comparable_cv = is.finite(within_individual_cv) & is.finite(inter_individual_cv),
  delta_cv = inter_individual_cv - within_individual_cv
)

soma_selection_audit_df <- protein_cv_df %>% dplyr::filter(cate == "Somascan") %>%
  dplyr::mutate(selected = comparable_cv & inter_individual_cv < within_individual_cv, selection_status = dplyr::case_when(
  !comparable_cv ~ "Incomplete CV", selected ~ "Selected: inter < within",
  inter_individual_cv == within_individual_cv ~ "Not selected: equal CV", TRUE ~ "Not selected: inter > within"
))
soma_selected_df <- soma_selection_audit_df %>% dplyr::filter(selected) %>% dplyr::transmute(
  protein, soma_within_individual_cv = within_individual_cv,
  soma_inter_individual_cv = inter_individual_cv, soma_delta_cv = delta_cv,
  soma_n_individuals = n_individuals, soma_n_sampling_times = n_sampling_times,
  soma_detection_category = protein_detection_category
)
comparison_df <- dplyr::bind_rows(lapply(target_platforms, function(platform) {
  target_cv <- protein_cv_df %>% dplyr::filter(cate == platform) %>% dplyr::select(-cate)
  soma_selected_df %>%
    dplyr::mutate(cate = platform, shared_detected = protein %in% detected_sets[[platform]]) %>%
    dplyr::left_join(target_cv, by = "protein", relationship = "one-to-one") %>% dplyr::mutate(
    comparable_cv = shared_detected & is.finite(within_individual_cv) & is.finite(inter_individual_cv),
    cv_direction = dplyr::case_when( !shared_detected ~ "Not shared Detected", !comparable_cv ~ "Incomplete CV",
    inter_individual_cv < within_individual_cv ~ "Inter < within",
    inter_individual_cv > within_individual_cv ~ "Inter > within", TRUE ~ "Equal CV"
  )
  )
}))
comparison_summary_df <- dplyr::bind_rows(lapply(target_platforms, function(platform) {
  d <- dplyr::filter(comparison_df, cate == platform)
  usable <- dplyr::filter(d, comparable_cv)
  n_usable <- nrow(usable)
  data.frame( platform = unname(platform_display_labels[[platform]]), n_soma_selected = nrow(soma_selected_df),
    n_shared_detected = sum(d$shared_detected), n_comparable = n_usable, n_not_shared_detected = sum(!d$shared_detected),
    n_shared_incomplete_cv = sum(d$shared_detected & !d$comparable_cv),
    n_inter_lower = sum(usable$cv_direction == "Inter < within"),
    n_inter_higher = sum(usable$cv_direction == "Inter > within"), n_equal_cv = sum(usable$cv_direction == "Equal CV"),
    pct_inter_lower = if (n_usable) {
    100 * sum(usable$cv_direction == "Inter < within") / n_usable
  } else {
    NA_real_
  },
    pct_inter_higher = if (n_usable) {
    100 * sum(usable$cv_direction == "Inter > within") / n_usable
  } else {
    NA_real_
  },
    pct_equal_cv = if (n_usable) {
    100 * sum(usable$cv_direction == "Equal CV") / n_usable
  } else {
    NA_real_
  },
    median_within_individual_cv = if (n_usable) {
    stats::median(usable$within_individual_cv)
  } else {
    NA_real_
  },
    median_inter_individual_cv = if (n_usable) {
    stats::median(usable$inter_individual_cv)
  } else {
    NA_real_
  },
    median_delta_cv = if (n_usable) {
    stats::median(usable$delta_cv)
  } else {
    NA_real_
  }, n_outside_display = sum(usable$within_individual_cv > scatter_display_limits[2] |
    usable$inter_individual_cv > scatter_display_limits[2])
  )
}))

scatter_df <- comparison_df %>% dplyr::filter(comparable_cv) %>%
  dplyr::mutate(cate = factor(cate, levels = target_platforms),
  protein_detection_category = factor(protein_detection_category, levels = category_levels))

scatter_width_cm <- 37
plot_platforms <- platform_order
soma_reference_scatter_df <- protein_cv_df %>% dplyr::filter(cate == "Somascan", protein %in% soma_selected_df$protein)
plot_scatter_df <- dplyr::bind_rows( scatter_df %>% dplyr::mutate(cate = as.character(cate)), soma_reference_scatter_df
) %>% dplyr::mutate(cate = factor(cate, levels = plot_platforms),
  protein_detection_category = factor(protein_detection_category, levels = category_levels))
# 04. Compare the SomaScan-selected proteins ----------------------------------

plot_selected_cv_scatter <- function(platform, show_legend = FALSE) {
  d <- dplyr::filter(plot_scatter_df, as.character(cate) == platform)
  target_display <- unname(platform_display_labels[[platform]])
  subtitle_text <- if (identical(platform, "Somascan")) {
    paste0("Selected N = ", nrow(soma_selected_df))
  } else {
    counts <- dplyr::filter(comparison_summary_df, .data$platform == .env$target_display)
    paste0("Shared N = ", counts$n_shared_detected, "; usable N = ", counts$n_comparable)
  }
  colour_scale <- function() scale_colour_manual(
    values = category_colours, breaks = category_levels, labels = category_display_labels,
    name = "Protein detection category\nMutually exclusive categories", drop = FALSE
  )

  density_df <- d %>% dplyr::group_by(protein_detection_category) %>% dplyr::filter(dplyr::n() >= 2L) %>% dplyr::ungroup()
  p_density <- ggplot(density_df, aes(x = within_individual_cv, colour = protein_detection_category))
  if (nrow(density_df)) {
    p_density <- p_density + stat_density(geom = "line", position = "identity", linewidth = 0.5, adjust = 0.9,
      na.rm = TRUE, show.legend = FALSE)
  } else {
    p_density <- p_density + geom_blank(data = data.frame(x = scatter_display_limits, y = c(0, 1)),
      inherit.aes = FALSE, aes(x = x, y = y)) + annotate("text", x = mean(scatter_display_limits), y = 0.5,
      label = "Density unavailable\n<2 proteins per category", size = 1.8, colour = "grey40")
  }
  p_density <- p_density + coord_cartesian(xlim = scatter_display_limits, expand = FALSE) + colour_scale() +
    labs(title = unname(platform_display_labels[[platform]]), subtitle = subtitle_text,
    x = NULL, y = "Density") + theme2_nobacknoerror + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 6),
    plot.margin = margin(3, 3, 0, 3, unit = "mm"))
  point_size <- if (platform %in% c("Olink-Explore", "Olink-HT", "Somascan")) {
    0.55
  } else {
    0.7
  }
  p_scatter <- ggplot(d, aes(x = within_individual_cv, y = inter_individual_cv, colour = protein_detection_category)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.3, linetype = "dashed")
  if (platform %in% c("Olink-Explore", "Olink-HT") && nrow(d)) {
    set.seed(if (identical(platform, "Olink-Explore")) {
      103L
    } else {
      104L
    })
    point_df <- d[sample.int(nrow(d)), , drop = FALSE]
    p_scatter <- p_scatter + geom_point(data = point_df, aes(alpha = protein_detection_category), size = point_size) +
      scale_alpha_manual(values = c("Detected proteins" = 1, "Individual detected proteins" = 1,
      "Broadly detected proteins" = 0.35), guide = "none")
  } else {
    p_scatter <- p_scatter +
      geom_point(data = dplyr::filter(d, protein_detection_category == category_levels[3]),
      size = point_size, alpha = 0.35) +
      geom_point(data = dplyr::filter(d, protein_detection_category != category_levels[3]), size = point_size, alpha = 1)
  }
  p_scatter <- p_scatter +
    coord_cartesian(xlim = scatter_display_limits, ylim = scatter_display_limits) + colour_scale() +
    guides(colour = guide_legend(ncol = 1)) +
    labs(x = "Within-individual CV (%)", y = "Inter-individual CV (%)") + theme4_backline +
    theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
    legend.position = if (show_legend) {
    "bottom"
  } else {
    "none"
  }, legend.title = element_text(size = 7), legend.text = element_text(size = 7),
    legend.key.height = grid::unit(3.5, "mm"), plot.margin = margin(0, 3, 3, 3, unit = "mm"))
  if (!nrow(d)) {
    p_scatter <- p_scatter + annotate("text", x = mean(scatter_display_limits), y = mean(scatter_display_limits),
      label = "No comparable shared proteins", size = 2.5, colour = "grey40")
  }
  p_density / p_scatter + patchwork::plot_layout(heights = c(1.1, 4))
}
legend_platform <- if (any(as.character(plot_scatter_df$cate) == "Olink-Explore")) {
  "Olink-Explore"
} else {
  present <- plot_platforms[plot_platforms %in% as.character(plot_scatter_df$cate)]
  if (length(present)) {
    present[1]
  } else {
    NA_character_
  }
}
cv_scatter_plots <- setNames(lapply(plot_platforms, function(platform) {
  plot_selected_cv_scatter(platform, show_legend = identical(platform, legend_platform))
}), plot_platforms)

scatter_output_file <- paste0(wdplot, plot_prefix, "05_somascan_selected_cross_platform_cv_scatter.pdf")

# 08. Compare proteins selected on the other platforms ------------------------

reciprocal_source_platforms <- platform_order[1:4]
reciprocal_target_platforms <- platform_order
reciprocal_required_cols <- c("cate", "protein", "protein_detection_category", "within_individual_cv", "inter_individual_cv",
  "n_individuals", "n_sampling_times", "comparable_cv", "delta_cv")
reciprocal_selection_audit_df <- protein_cv_df %>%
  dplyr::filter(cate %in% reciprocal_source_platforms) %>% dplyr::rename(source_platform = cate) %>% dplyr::mutate(
  selected = is.finite(within_individual_cv) & is.finite(inter_individual_cv) &
  within_individual_cv > inter_individual_cv, selection_status = dplyr::case_when(
  !is.finite(within_individual_cv) | !is.finite(inter_individual_cv) ~ "Incomplete CV",
  selected ~ "Selected: within > inter", within_individual_cv == inter_individual_cv ~ "Not selected: equal CV",
  TRUE ~ "Not selected: within < inter"
)
)
reciprocal_selected_df <- reciprocal_selection_audit_df %>%
  dplyr::filter(selected) %>% dplyr::transmute( source_platform, protein, source_within_individual_cv = within_individual_cv,
  source_inter_individual_cv = inter_individual_cv, source_delta_cv = delta_cv, source_n_individuals = n_individuals,
  source_n_sampling_times = n_sampling_times, source_detection_category = protein_detection_category
)

reciprocal_comparison_df <- dplyr::bind_rows(lapply(reciprocal_source_platforms, function(source) {
  selected_df <- dplyr::filter(reciprocal_selected_df, source_platform == .env$source)
  dplyr::bind_rows(lapply(reciprocal_target_platforms, function(target) {
    target_cv <- protein_cv_df %>%
      dplyr::filter(cate == .env$target) %>% dplyr::select(-cate)
    selected_df %>% dplyr::mutate( target_platform = target, is_source_reference = source == target,
      shared_detected = protein %in% detected_sets[[target]]
    ) %>% dplyr::left_join(target_cv, by = "protein", relationship = "one-to-one") %>% dplyr::mutate(
      comparable_cv = shared_detected & is.finite(within_individual_cv) &
      is.finite(inter_individual_cv), cv_direction = dplyr::case_when( !shared_detected ~ "Not shared Detected",
      !comparable_cv ~ "Incomplete CV", within_individual_cv > inter_individual_cv ~ "Within > inter",
      within_individual_cv < inter_individual_cv ~ "Within < inter", TRUE ~ "Equal CV"
    )
    )
  }))
}))
reciprocal_summary_df <- dplyr::bind_rows(lapply(reciprocal_source_platforms, function(source) {
  n_selected <- sum(reciprocal_selected_df$source_platform == source)
  dplyr::bind_rows(lapply(reciprocal_target_platforms, function(target) {
    d <- dplyr::filter(reciprocal_comparison_df, source_platform == .env$source, target_platform == .env$target)
    usable <- dplyr::filter(d, comparable_cv)
    n_usable <- nrow(usable)
    data.frame( source_platform = source, target_platform = target,
      is_source_reference = source == target, n_source_selected = n_selected,
      n_shared_detected = sum(d$shared_detected), n_comparable = n_usable, n_not_shared_detected = sum(!d$shared_detected),
      n_shared_incomplete_cv = sum(d$shared_detected & !d$comparable_cv),
      n_within_higher = sum(usable$cv_direction == "Within > inter"),
      n_within_lower = sum(usable$cv_direction == "Within < inter"), n_equal_cv = sum(usable$cv_direction == "Equal CV"),
      pct_within_higher = if (n_usable) {
      100 * sum(usable$cv_direction == "Within > inter") / n_usable
    } else {
      NA_real_
    },
      pct_within_lower = if (n_usable) {
      100 * sum(usable$cv_direction == "Within < inter") / n_usable
    } else {
      NA_real_
    },
      pct_equal_cv = if (n_usable) {
      100 * sum(usable$cv_direction == "Equal CV") / n_usable
    } else {
      NA_real_
    },
      median_within_individual_cv = if (n_usable) {
      stats::median(usable$within_individual_cv)
    } else {
      NA_real_
    },
      median_inter_individual_cv = if (n_usable) {
      stats::median(usable$inter_individual_cv)
    } else {
      NA_real_
    },
      median_delta_cv = if (n_usable) {
      stats::median(usable$delta_cv)
    } else {
      NA_real_
    }, n_outside_display = sum(usable$within_individual_cv > scatter_display_limits[2] |
      usable$inter_individual_cv > scatter_display_limits[2])
    )
  }))
}))

reciprocal_reference_summary_df <- dplyr::filter(reciprocal_summary_df, is_source_reference)
reciprocal_scatter_df <- reciprocal_comparison_df %>%
  dplyr::filter(comparable_cv) %>% dplyr::mutate(
  source_platform = factor(source_platform, levels = reciprocal_source_platforms),
  target_platform = factor(target_platform, levels = reciprocal_target_platforms),
  protein_detection_category = factor(protein_detection_category, levels = category_levels)
)

reciprocal_scatter_width_cm <- 37
reciprocal_scatter_height_cm <- 48
plot_reciprocal_cv_scatter <- function(source, target) {
  d <- dplyr::filter(reciprocal_scatter_df, as.character(source_platform) == .env$source,
    as.character(target_platform) == .env$target)
  counts <- dplyr::filter(reciprocal_summary_df, source_platform == .env$source, target_platform == .env$target)
  subtitle_text <- if (source == target) {
    paste0("Source reference; selected N = ", counts$n_source_selected)
  } else {
    paste0("Shared N = ", counts$n_shared_detected, "; usable N = ", counts$n_comparable)
  }
  colour_scale <- function() scale_colour_manual(
    values = category_colours, breaks = category_levels, labels = category_display_labels,
    name = "Protein detection category (displayed platform)", drop = FALSE
  )
  density_df <- d %>% dplyr::group_by(protein_detection_category) %>% dplyr::filter(dplyr::n() >= 2L) %>% dplyr::ungroup()
  p_density <- ggplot(density_df, aes(x = within_individual_cv, colour = protein_detection_category))
  if (nrow(density_df)) {
    p_density <- p_density + stat_density(geom = "line", position = "identity", linewidth = 0.5, adjust = 0.9,
      na.rm = TRUE, show.legend = FALSE)
  } else {
    p_density <- p_density + geom_blank(data = data.frame(x = scatter_display_limits, y = c(0, 1)),
      inherit.aes = FALSE, aes(x = x, y = y)) + annotate("text", x = mean(scatter_display_limits), y = 0.5,
      label = "Density unavailable\n<2 proteins per category", size = 1.8, colour = "grey40")
  }
  p_density <- p_density + coord_cartesian(xlim = scatter_display_limits, expand = FALSE) + colour_scale() +
    labs(title = unname(platform_display_labels[[target]]), subtitle = subtitle_text,
    x = NULL, y = "Density") + theme2_nobacknoerror + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 6),
    plot.margin = margin(3, 3, 0, 3, unit = "mm"))
  point_size <- if (target %in% c("Olink-Explore", "Olink-HT", "Somascan")) {
    0.55
  } else {
    0.7
  }
  p_scatter <- ggplot(d, aes(x = within_individual_cv, y = inter_individual_cv, colour = protein_detection_category)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.3, linetype = "dashed")
  if (target %in% c("Olink-Explore", "Olink-HT") && nrow(d)) {
    set.seed(if (identical(target, "Olink-Explore")) {
      103L
    } else {
      104L
    })
    point_df <- d[sample.int(nrow(d)), , drop = FALSE]
    p_scatter <- p_scatter + geom_point(data = point_df, aes(alpha = protein_detection_category), size = point_size) +
      scale_alpha_manual(values = c("Detected proteins" = 1, "Individual detected proteins" = 1,
      "Broadly detected proteins" = 0.35), guide = "none")
  } else {
    p_scatter <- p_scatter +
      geom_point(data = dplyr::filter(d, protein_detection_category == category_levels[3]),
      size = point_size, alpha = 0.35) +
      geom_point(data = dplyr::filter(d, protein_detection_category != category_levels[3]), size = point_size, alpha = 1)
  }
  p_scatter <- p_scatter +
    coord_cartesian(xlim = scatter_display_limits, ylim = scatter_display_limits) + colour_scale() +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "Within-individual CV (%)", y = "Inter-individual CV (%)") + theme4_backline +
    theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
    legend.position = "bottom", legend.title = element_text(size = 7),
    legend.text = element_text(size = 7), legend.key.height = grid::unit(3.5, "mm"),
    plot.margin = margin(0, 3, 3, 3, unit = "mm"))
  if (!nrow(d)) {
    empty_label <- if (!counts$n_source_selected) {
      "No source-selected proteins"
    } else {
      "No comparable shared proteins"
    }
    p_scatter <- p_scatter + annotate("text", x = mean(scatter_display_limits), y = mean(scatter_display_limits),
      label = empty_label, size = 2.5, colour = "grey40")
  }
  p_density / p_scatter + patchwork::plot_layout(heights = c(1.1, 4))
}
reciprocal_plot_rows <- list()
for (source in reciprocal_source_platforms) {
  n_selected <- sum(reciprocal_selected_df$source_platform == source)
  row_label <- paste0("Selected in ", platform_display_labels[[source]],
    ": within-individual CV > inter-individual CV; N = ", n_selected)

  p_source_label <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = row_label, hjust = 0, size = 3.2, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(2, 3, 0, 3, unit = "mm"))
  panel_list <- lapply(reciprocal_target_platforms, function(target) {
    plot_reciprocal_cv_scatter(source, target)
  })
  reciprocal_plot_rows[[paste0(source, "_label")]] <- p_source_label
  reciprocal_plot_rows[[source]] <- patchwork::wrap_plots(panel_list, nrow = 1)
}

reciprocal_plot_file <- paste0(wdplot, plot_prefix, "10_reciprocal_selected_cross_platform_cv_scatter.pdf")
# 13. Five-source, five-platform combined figure ------------------------------

five_source_plot_rows <- reciprocal_plot_rows
soma_row_label <- paste0("Selected in SomaScan: within-individual CV > inter-individual CV; N = ", nrow(soma_selected_df))
five_source_plot_rows[["Somascan_label"]] <- ggplot() +
  annotate("text", x = 0, y = 0.5, label = soma_row_label, hjust = 0, size = 3.2, fontface = "bold") +
  xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(2, 3, 0, 3, unit = "mm"))
five_source_plot_rows[["Somascan"]] <- patchwork::wrap_plots(lapply(cv_scatter_plots[platform_order], function(p) {
  p & guides(colour = "none", alpha = "none")
}), nrow = 1)
p_five_source_cross_platform_cv <- patchwork::wrap_plots(five_source_plot_rows,
  ncol = 1, heights = rep(c(0.4, 10.4), length(platform_order)), guides = "collect") &
  theme(legend.position = "bottom")

# 13. Save the five-source CV comparison --------------------------------------

five_source_plot_file <- paste0(wdplot, "S4_five_source_cross_platform_cv_scatter.pdf")
ggsave(five_source_plot_file, pdf_text_grob(p_five_source_cross_platform_cv), width = 30,height = 44,
  units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# Sample-wise correlations ---------------------------------------------------
plot_prefix <- "04-"
format_p_stars <- function(p) dplyr::case_when( is.na(p) ~ NA_character_, p < 0.001 ~ "***", p < 0.01 ~ "**",
  p < 0.05 ~ "*", TRUE ~ ""
)

options(stringsAsFactors = FALSE)

platform_levels <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
violin_fill3 <- c("#86BCB6FF", "#B07AA1FF", "#D7B5A6FF")
platform_fill5 <- c("#D9B381", "#B5252E", "#264653", "#AEC7E8FF", "#69B190FF")
names(platform_fill5) <- platform_levels
pair_fill3 <- c("#9D7660FF", "#E1C367FF", "#C97231FF")

overlap461 <- Reduce(intersect, list(prolist[[1]], prolist[[2]], prolist[[3]], prolist[[4]], prolist[[5]]))
individual_detected_overlap <- Reduce(intersect, individual_prolist)
broad_detected_overlap <- Reduce(intersect, broad_prolist)
# 01-03. Sample-wise correlations ---------------------------------------------

samplewise_catename <- "Sample-wise correlation"

samplewise_full_res_list <- list( mscorrdiffindi3time(
  pro = prolist[["MS-OT480"]], cate = "MS-OT480", data = dflistna[["MS-OT480"]], isna = "NA"
), mscorrdiffindi3time( pro = prolist[["MS-Astral"]], cate = "MS-Astral", data = dflistna[["MS-Astral"]], isna = "NA"
), olinkcorrdiffindi3time( pro = prolist[["Olink-Explore"]], cate = "Olink-Explore",
  data = dflistna[["Olink-Explore"]], isna = FALSE
), olinkcorrdiffindi3time( pro = prolist[["Olink-HT"]], cate = "Olink-HT", data = dflistna[["Olink-HT"]], isna = FALSE
), somacorrdiffindi3time( pro = prolist[["Somascan"]], cate = "Somascan", data = dflistna[["Somascan"]], isna = FALSE
)
)

samplewise_full_master_df <- do.call(rbind, samplewise_full_res_list) %>% dplyr::mutate(
  Pair_new = gsub(" ", "", paste0("v", gsub(" vs ", "_vs_v", Pair))),
  cate = factor(cate, levels = platform_levels), subject = paste0(S1, "_", S2)
)

samplewise_target_cols <- c("cate", "subject", "Pair_new", "Freq")

samplewise_full_df_total <- samplewise_full_master_df %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Total correlation")
samplewise_full_df_interindi <- samplewise_full_master_df %>% dplyr::filter(Var1_prefix == Var2_prefix) %>%
  dplyr::select(all_of(samplewise_target_cols)) %>% dplyr::mutate(cate2 = "Inter-individual correlation")
samplewise_full_df_intervisit <- samplewise_full_master_df %>%
  dplyr::filter(S1 == S2) %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Within-individual cross-visit correlation")
samplewise_full_df_bind <- dplyr::bind_rows( samplewise_full_df_total, samplewise_full_df_interindi,
  samplewise_full_df_intervisit
)
samplewise_full_df_bind$cate2 <- factor(samplewise_full_df_bind$cate2, levels = c( "Total correlation",
  "Within-individual cross-visit correlation", "Inter-individual correlation"
))

samplewise_full_median <- samplewise_full_df_bind %>% dplyr::group_by( cate2, cate
) %>% dplyr::summarise(medianva = sprintf("%.2f", median(Freq)), .groups = "drop")
samplewise_full_median$cate2 <- factor(samplewise_full_median$cate2, levels = levels(samplewise_full_df_bind$cate2)
)

samplewise_ov_catename <- "Sample-wise correlation"

samplewise_ov_res_list <- list(
  mscorrdiffindi3time(pro = overlap461, cate = "MS-OT480", data = dflistna[["MS-OT480"]], isna = "NA"),
  mscorrdiffindi3time(pro = overlap461, cate = "MS-Astral", data = dflistna[["MS-Astral"]], isna = "NA"),
  olinkcorrdiffindi3time( pro = overlap461, cate = "Olink-Explore", data = dflistna[["Olink-Explore"]], isna = FALSE
), olinkcorrdiffindi3time(pro = overlap461, cate = "Olink-HT", data = dflistna[["Olink-HT"]], isna = FALSE),
  somacorrdiffindi3time(pro = overlap461, cate = "Somascan", data = dflistna[["Somascan"]], isna = FALSE)
)

samplewise_ov_master_df <- do.call(rbind, samplewise_ov_res_list) %>% dplyr::mutate(
  Pair_new = gsub(" ", "", paste0("v", gsub(" vs ", "_vs_v", Pair))), cate = factor(cate, levels = platform_levels),
  subject = paste0(S1, "_", S2)
)

samplewise_ov_df_total <- samplewise_ov_master_df %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Total correlation")
samplewise_ov_df_interindi <- samplewise_ov_master_df %>%
  dplyr::filter(S1 != S2 & Var1_prefix == Var2_prefix) %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Inter-individual correlation")
samplewise_ov_df_intervisit <- samplewise_ov_master_df %>%
  dplyr::filter(S1 == S2) %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Within-individual cross-visit correlation")
samplewise_ov_df_bind <- dplyr::bind_rows( samplewise_ov_df_total, samplewise_ov_df_interindi, samplewise_ov_df_intervisit
)
samplewise_ov_df_bind$cate2 <- factor(samplewise_ov_df_bind$cate2, levels = c( "Total correlation",
  "Within-individual cross-visit correlation", "Inter-individual correlation"
))

samplewise_ov_median <- samplewise_ov_df_bind %>%
  dplyr::group_by( cate2, cate
) %>% dplyr::summarise(medianva = sprintf("%.2f", median(Freq)), .groups = "drop")
samplewise_ov_median$cate2 <- factor(samplewise_ov_median$cate2, levels = levels(samplewise_ov_df_bind$cate2))

belowlod_oe <- setdiff(rownames(olink3072), prolist[["Olink-Explore"]])
belowlod_oht <- setdiff(rownames(olinkht), prolist[["Olink-HT"]])

plot_w_cm <- 33
plot_h_cm <- 10
corrbind_fill_cols <- c( "Total correlation" = "#86BCB6FF",
  "Within-individual cross-visit correlation" = "#B07AA1FF", "Inter-individual correlation" = "#D7B5A6FF"
)

samplewise_bl_res_list <- list( olinkcorrdiffindi3time(
  pro = belowlod_oe, cate = "Olink-Explore", data = dflistna[["Olink-Explore"]], isna = FALSE, keepfilter = FALSE
), olinkcorrdiffindi3time( pro = belowlod_oht, cate = "Olink-HT", data = dflistna[["Olink-HT"]],
  isna = FALSE, keepfilter = FALSE
)
)

samplewise_bl_master_df <- do.call(rbind, samplewise_bl_res_list) %>% dplyr::mutate(
  Pair_new = gsub(" ", "", paste0("v", gsub(" vs ", "_vs_v", Pair))),
  cate = factor(cate, levels = platform_levels), subject = paste0(S1, "_", S2)
)
samplewise_bl_df_total <- samplewise_bl_master_df %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Total correlation")
samplewise_bl_df_interindi <- samplewise_bl_master_df %>% dplyr::filter(Var1_prefix == Var2_prefix) %>%
  dplyr::select(all_of(samplewise_target_cols)) %>% dplyr::mutate(cate2 = "Inter-individual correlation")
samplewise_bl_df_intervisit <- samplewise_bl_master_df %>%
  dplyr::filter(S1 == S2) %>% dplyr::select(all_of(samplewise_target_cols)) %>%
  dplyr::mutate(cate2 = "Within-individual cross-visit correlation")
samplewise_bl_df_bind <- dplyr::bind_rows( samplewise_bl_df_total, samplewise_bl_df_interindi, samplewise_bl_df_intervisit
)
samplewise_bl_df_bind$cate2 <- factor(samplewise_bl_df_bind$cate2, levels = c( "Total correlation",
  "Within-individual cross-visit correlation", "Inter-individual correlation"
))

build_samplewise_corr_bundle <- function(protein_input) {
  protein_for <- function(platform) if (is.list(protein_input)) {
    protein_input[[platform]]
  } else {
    protein_input
  }
  res_list <- list( mscorrdiffindi3time( pro = protein_for("MS-OT480"), cate = "MS-OT480", data = dflistna[["MS-OT480"]],
    isna = "NA"
  ), mscorrdiffindi3time( pro = protein_for("MS-Astral"), cate = "MS-Astral", data = dflistna[["MS-Astral"]], isna = "NA"
  ), olinkcorrdiffindi3time( pro = protein_for("Olink-Explore"), cate = "Olink-Explore",
    data = dflistna[["Olink-Explore"]], isna = FALSE
  ), olinkcorrdiffindi3time( pro = protein_for("Olink-HT"), cate = "Olink-HT", data = dflistna[["Olink-HT"]], isna = FALSE
  ), somacorrdiffindi3time( pro = protein_for("Somascan"), cate = "Somascan", data = dflistna[["Somascan"]], isna = FALSE
  )
  )
  master_df <- do.call(rbind, res_list) %>%
    dplyr::mutate( Pair_new = gsub(" ", "", paste0("v", gsub(" vs ", "_vs_v", Pair))),
    cate = factor(cate, levels = platform_levels), subject = paste0(S1, "_", S2)
  )
  total_df <- master_df %>% dplyr::select(all_of(samplewise_target_cols)) %>% dplyr::mutate(cate2 = "Total correlation")
  interindividual_df <- master_df %>% dplyr::filter(Var1_prefix == Var2_prefix) %>%
    dplyr::select(all_of(samplewise_target_cols)) %>%
    dplyr::mutate(cate2 = "Inter-individual correlation")
  withinindividual_df <- master_df %>% dplyr::filter(S1 == S2) %>% dplyr::select(all_of(samplewise_target_cols)) %>%
    dplyr::mutate(cate2 = "Within-individual cross-visit correlation")
  dplyr::bind_rows(total_df, interindividual_df, withinindividual_df) %>% dplyr::mutate(cate2 = factor(cate2, levels = c(
    "Total correlation", "Within-individual cross-visit correlation", "Inter-individual correlation"
  )))
}
make_samplewise_corr_plot <- function(plot_df, output_stub, save_output = TRUE) {
  median_df <- plot_df %>%
    dplyr::group_by(cate, panel, cate2) %>%
    dplyr::summarise(mediancorr = sprintf("%.2f", median(Freq, na.rm = TRUE)), .groups = "drop")
  p <- ggplot(plot_df, aes(x = panel, y = Freq, fill = cate2)) +
    geom_violin( trim = FALSE, color = "transparent", alpha = 1, width = 0.9, position = position_dodge(width = 0.8)
  ) + geom_boxplot( data = plot_df, aes(x = panel, y = Freq, fill = cate2),
    position = position_dodge(width = 0.8, preserve = "single"), width = 0.12, size = 0.25, colour = "white",
    outlier.shape = NA, alpha = 0.5
  ) + geom_text( data = median_df, aes(x = panel, y = 1, label = mediancorr, group = cate2),
    position = position_dodge2(width = 0.8, preserve = "single"), color = "black", size = 2.8, vjust = 0,
    angle = 45, hjust = 0
  ) + coord_cartesian(ylim = c(-0.3, 1.2), clip = "off") + scale_fill_manual( values = corrbind_fill_cols, labels = c(
    "Total correlation" = "Total correlation",
    "Within-individual cross-visit correlation" = "Within-individual correlation",
    "Inter-individual correlation" = "Inter-individual correlation"
  ), name = ""
  ) + labs(x = "", y = "Spearman correlation coefficient", fill = "") + facet_wrap(~cate, nrow = 1, scales = "free_x") +
    theme2_nobacknoerror + themetext30 + themelegendtop + theme(
    axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), legend.text = element_text(size = 10), strip.text = element_text(size = 10)
  )
  output_filename <- if (grepl("^[FS][0-9]", output_stub)) {
    paste0(output_stub, ".pdf")
  } else {
    paste0(plot_prefix, output_stub, ".pdf")
  }

  ggsave(paste0(wdplot, output_filename), pdf_text_grob(p), width = plot_w_cm, height = plot_h_cm, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

}

samplewise_individual_full_df <- build_samplewise_corr_bundle(individual_prolist)

samplewise_individual_overlap_df <- build_samplewise_corr_bundle(individual_detected_overlap)

samplewise_broad_full_df <- build_samplewise_corr_bundle(broad_prolist)

samplewise_broad_overlap_df <- build_samplewise_corr_bundle(broad_detected_overlap)

samplewise_corrbind_plot_df <- dplyr::bind_rows( samplewise_full_df_bind %>% dplyr::mutate(panel = "Detected proteins"),
  samplewise_individual_full_df %>% dplyr::mutate(panel = "Individual detected proteins"),
  samplewise_broad_full_df %>% dplyr::mutate(panel = "Broadly detected proteins"),
  samplewise_bl_df_bind %>% dplyr::mutate(panel = "Undetected predefined protein")
) %>% dplyr::mutate( panel = factor(panel, levels = c( "Detected proteins", "Individual detected proteins",
  "Broadly detected proteins", "Undetected predefined protein"
)), cate = factor(cate, levels = platform_levels)
)
p_samplewise_corrbind_bytech <- make_samplewise_corr_plot( samplewise_corrbind_plot_df, "F4b_samplewise_corr_bytech"
)

samplewise_overlap_sets_plot_df <- dplyr::bind_rows(
  samplewise_ov_df_bind %>% dplyr::mutate(panel = "Detected protein overlap"),
  samplewise_individual_overlap_df %>% dplyr::mutate(panel = "Individual detected protein overlap"),
  samplewise_broad_overlap_df %>% dplyr::mutate(panel = "Broadly detected protein overlap")
) %>% dplyr::mutate( panel = factor(panel, levels = c( "Detected protein overlap", "Individual detected protein overlap",
  "Broadly detected protein overlap"
)), cate = factor(cate, levels = platform_levels)
)
p_samplewise_overlap_sets_corr_bytech <- make_samplewise_corr_plot( samplewise_overlap_sets_plot_df,
  "S3c_detection_overlap_samplewise_corr_bytech"
)

# CV inputs for Olink/SomaScan comparisons -----------------------------------
  soma_cv_result_to_protein <- function(soma_cv_result, strata_cols, plot_cols) {
    somer_vol <- soma_cv_result[["vol"]]
    protein_audit <- summarise_soma_cv_to_protein( cv_df = somer_vol, strata_cols = strata_cols
    )
    missing_plot_cols <- setdiff(plot_cols, colnames(protein_audit))
    list( vol = protein_audit %>% dplyr::select(dplyr::all_of(plot_cols)), protein_audit = protein_audit,
      somer_vol = somer_vol, somer_bar = soma_cv_result[["bar"]]
    )
  }
  bls_soma_ref <- prolist[["Somascan"]]

  bls_group_levels <- c( "Detected proteins, SomaScan-shared", "Detected proteins, Olink-only",
    "Undetected predefined protein, SomaScan-shared", "Undetected predefined protein, Olink-only",
    "SomaScan, detected proteins", "SomaScan, undetected predefined protein", "SomaScan-only"
  )

  bls_col <- c( "Detected proteins, SomaScan-shared" = "#475272", "Detected proteins, Olink-only" = "#68A39D",
    "Undetected predefined protein, SomaScan-shared" = "#ADC482", "Undetected predefined protein, Olink-only" = "#D9B381",
    "SomaScan, detected proteins" = "#C27877", "SomaScan, undetected predefined protein" = "#7D576F",
    "SomaScan-only" = "#B89E8A"
  )

  bls_olink_grps <- bls_group_levels[1:4]
  bls_soma_grps <- bls_group_levels[5:7]

  bls_olink_pairs <- list( c("Detected proteins, SomaScan-shared", "Undetected predefined protein, SomaScan-shared"),
    c("Detected proteins, Olink-only", "Undetected predefined protein, Olink-only")
  )
  bls_soma_pairs <- list( c("SomaScan, detected proteins", "SomaScan, undetected predefined protein"),
    c("SomaScan, detected proteins", "SomaScan-only"), c("SomaScan, undetected predefined protein", "SomaScan-only")
  )
  bls_soma_lbl <- c("Detected proteins", "Undetected predefined protein", "SomaScan-only")

  bls_olink_ylbl <- c( "Detected proteins, SomaScan-shared" = "Detected proteins\n(SomaScan-shared)",
    "Detected proteins, Olink-only" = "Detected proteins\n(Olink-only)",
    "Undetected predefined protein, SomaScan-shared" = "Undetected predefined protein\n(SomaScan-shared)",
    "Undetected predefined protein, Olink-only" = "Undetected predefined protein\n(Olink-only)"
  )
  bls_soma_ylbl <- c( "SomaScan, detected proteins" = "Detected proteins",
    "SomaScan, undetected predefined protein" = "Undetected predefined protein", "SomaScan-only" = "SomaScan-only"
  )

  olink_cv_three <- function(olink_cate, pro_set, group_label) {
    d_iv <- olinkcv( cate = olink_cate, pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    )[[1]]
    d_iv$name <- "Within-individual CV"
    d_tot <- olinkcvall( cate = olink_cate, pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    )[[1]]
    d_tot$name <- "Total CV"
    d_ii <- olinkcvindi( cate = olink_cate, pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    )[[1]]
    d_ii$name <- "Inter-individual CV"
    out <- dplyr::bind_rows(d_iv, d_tot, d_ii)
    out$group <- group_label
    out
  }
  soma_cv_three <- function(pro_set, group_label) {
    d_iv <- somacvf( pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    ) %>% soma_cv_result_to_protein( strata_cols = c("SubjectID", "cate", "cate1"),
      plot_cols = c("SubjectID", "protein", "CV", "cate", "cate1")
    )
    d_iv <- d_iv$protein_audit
    d_iv$name <- "Within-individual CV"

    d_tot <- somacvfall( pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    ) %>% soma_cv_result_to_protein( strata_cols = c("cate", "cate1"), plot_cols = c("protein", "CV", "cate", "cate1")
    )
    d_tot <- d_tot$protein_audit
    d_tot$name <- "Total CV"

    d_ii <- somacvindi( pro = pro_set, isna = FALSE, keepfilter = FALSE, cate1input = group_label
    ) %>% soma_cv_result_to_protein( strata_cols = c("Time", "cate", "cate1"),
      plot_cols = c("Time", "protein", "CV", "cate", "cate1")
    )
    d_ii <- d_ii$protein_audit
    d_ii$name <- "Inter-individual CV"

    out <- dplyr::bind_rows(d_iv, d_tot, d_ii)
    out$group <- group_label
    out
  }
  bls_build_cv <- function(olink_cate, olink_full_rownames, olink_robust, facet_label) {
    belowlod <- setdiff(olink_full_rownames, olink_robust)
    over_ov <- intersect(olink_robust, bls_soma_ref)
    over_only <- setdiff(olink_robust, bls_soma_ref)
    below_ov <- intersect(belowlod, bls_soma_ref)
    below_only <- setdiff(belowlod, bls_soma_ref)
    soma_only <- setdiff(bls_soma_ref, union(over_ov, below_ov))
    res <- dplyr::bind_rows( olink_cv_three(olink_cate, over_ov, bls_group_levels[1]),
      olink_cv_three(olink_cate, over_only, bls_group_levels[2]), olink_cv_three(olink_cate, below_ov, bls_group_levels[3]),
      olink_cv_three(olink_cate, below_only, bls_group_levels[4]), soma_cv_three(over_ov, bls_group_levels[5]),
      soma_cv_three(below_ov, bls_group_levels[6]), soma_cv_three(soma_only, bls_group_levels[7])
    )
    res$facet <- facet_label
    res
  }

  bls_plot_df <- dplyr::bind_rows(
    bls_build_cv("3072", rownames(olink3072), prolist[["Olink-Explore"]], "Olink-Explore vs Soma"),
    bls_build_cv("ht", rownames(olinkht), prolist[["Olink-HT"]], "Olink-HT vs Soma")
  )

  bls_plot_df$name <- factor(bls_plot_df$name, levels = c( "Total CV", "Within-individual CV", "Inter-individual CV"
  ))
  bls_plot_df$group <- factor(bls_plot_df$group, levels = bls_group_levels)
  bls_plot_df$facet <- factor(bls_plot_df$facet, levels = c("Olink-Explore vs Soma", "Olink-HT vs Soma"))

  combined_reproducibility_soma_cv_protein_audit <- bls_plot_df %>% dplyr::filter(cate == "Somascan") %>%
    dplyr::select( facet, group, name, SubjectID, Time, protein, CV, n_somers, cv_sd, cv_min, cv_max, somer_ids, cate, cate1
  )

# 16. Plot within-individual CV for Olink and SomaScan ------------------------

cv_only_required <- c("bls_plot_df", "prolist", "olink3072", "olinkht")

cv_only_cols <- c("name", "cate", "facet", "SubjectID", "protein", "CV")
cv_only_soma_detected <- unique(prolist[["Somascan"]])
cv_only_olink_detected <- lapply(prolist[c("Olink-Explore", "Olink-HT")], unique)
cv_only_olink_universe <- list( "Olink-Explore" = rownames(olink3072), "Olink-HT" = rownames(olinkht)
)
cv_only_olink_undetected <- lapply(names(cv_only_olink_detected), function(platform) {
  setdiff(cv_only_olink_universe[[platform]], cv_only_olink_detected[[platform]])
})
names(cv_only_olink_undetected) <- names(cv_only_olink_detected)

cv_only_source <- bls_plot_df %>%
  dplyr::filter(name == "Within-individual CV") %>% dplyr::filter(
  (cate == "Olink-Explore" & facet == "Olink-Explore vs Soma") |
  (cate == "Olink-HT" & facet == "Olink-HT vs Soma") |
  (cate == "Somascan" & facet == "Olink-Explore vs Soma")
) %>% dplyr::transmute(platform = as.character(cate), SubjectID, protein, cv_pct = CV * 100)
cv_only_olink_labels <- c( "All detected\nproteins", "SomaScan-shared\ndetected proteins",
  "All undetected\nproteins", "SomaScan-shared\nundetected proteins"
)
cv_only_olink_sets <- lapply(names(cv_only_olink_detected), function(platform) {
  setNames(list( cv_only_olink_detected[[platform]], intersect(cv_only_olink_detected[[platform]], cv_only_soma_detected),
    cv_only_olink_undetected[[platform]], intersect(cv_only_olink_undetected[[platform]], cv_only_soma_detected)
  ), cv_only_olink_labels)
})
names(cv_only_olink_sets) <- names(cv_only_olink_detected)
cv_only_soma_sets <- list( "All SomaScan\ndetected proteins" = cv_only_soma_detected,
  "Olink-Explore\nshared detected\nproteins" = intersect( cv_only_soma_detected, cv_only_olink_detected[["Olink-Explore"]]
), "Olink-Explore\nshared undetected\nproteins" = intersect( cv_only_soma_detected,
  cv_only_olink_undetected[["Olink-Explore"]]
), "Olink-HT\nshared detected\nproteins" = intersect( cv_only_soma_detected, cv_only_olink_detected[["Olink-HT"]]
), "Olink-HT\nshared undetected\nproteins" = intersect( cv_only_soma_detected, cv_only_olink_undetected[["Olink-HT"]]
)
)

cv_only_group_data <- function(platform_key, protein_sets) {
  source_df <- dplyr::filter(cv_only_source, platform == platform_key)
  dplyr::bind_rows(lapply(names(protein_sets), function(group_label) {
    source_df %>%
      dplyr::filter(protein %in% protein_sets[[group_label]]) %>% dplyr::mutate(group = group_label)
  }))
}
cv_only_olink_df <- dplyr::bind_rows(lapply(names(cv_only_olink_sets), function(platform) {
  cv_only_group_data(platform, cv_only_olink_sets[[platform]])
})) %>%
  dplyr::mutate( platform = factor(platform, levels = c("Olink-Explore", "Olink-HT")),
  group = factor(group, levels = cv_only_olink_labels)
)
cv_only_soma_df <- cv_only_group_data("Somascan", cv_only_soma_sets) %>%
  dplyr::mutate(group = factor(group, levels = names(cv_only_soma_sets)))
cv_only_summary <- function(df) {
  df %>% dplyr::group_by(platform, group, .drop = FALSE) %>% dplyr::summarise(
    n_proteins_with_cv = dplyr::n_distinct(protein[is.finite(cv_pct)]), n_finite_cv = sum(is.finite(cv_pct)),
    n_nonfinite_cv = sum(!is.finite(cv_pct)),
    median_cv_pct = if (any(is.finite(cv_pct))) {
    median(cv_pct[is.finite(cv_pct)])
  } else {
    NA_real_
  }, .groups = "drop"
  ) %>%
    dplyr::mutate(median_label = ifelse( is.finite(median_cv_pct), sprintf("%.2f", floor(median_cv_pct * 100) / 100), "N/A"
  ))
}
cv_only_olink_summary_df <- cv_only_summary(cv_only_olink_df)
cv_only_soma_summary_df <- cv_only_summary(cv_only_soma_df)

cv_only_olink_pairs <- combn(cv_only_olink_labels, 2, simplify = FALSE)
cv_only_wilcoxon <- function(df, pair_list) {
  dplyr::bind_rows(lapply(unique(as.character(df$platform)), function(platform_key) {
    platform_df <- dplyr::filter(df, as.character(platform) == platform_key)
    dplyr::bind_rows(lapply(seq_along(pair_list), function(pair_idx) {
      pair <- pair_list[[pair_idx]]
      first_df <- dplyr::filter(platform_df, as.character(group) == pair[1], is.finite(cv_pct))
      second_df <- dplyr::filter(platform_df, as.character(group) == pair[2], is.finite(cv_pct))
      test_error <- ""
      test <- if (nrow(first_df) && nrow(second_df)) {
        tryCatch( stats::wilcox.test(first_df$cv_pct, second_df$cv_pct,
          alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE
        ),
          error = function(e) {
          test_error <<- conditionMessage(e)
          NULL
        })
      } else {
        NULL
      }
      raw_p <- if (is.null(test)) {
        NA_real_
      } else {
        unname(test$p.value)
      }
      status <- if (!nrow(first_df) || !nrow(second_df)) {
        "empty_group"
      } else if (nzchar(test_error)) {
        "test_error"
      } else if (!is.finite(raw_p)) {
        "not_estimable"
      } else {
        "ok"
      }
      data.frame( platform = platform_key, pair_idx = pair_idx, group1 = pair[1], group2 = pair[2],
        n_cv_group1 = nrow(first_df), n_cv_group2 = nrow(second_df),
        n_proteins_group1 = dplyr::n_distinct(first_df$protein), n_proteins_group2 = dplyr::n_distinct(second_df$protein),
        n_overlapping_proteins = length(intersect(first_df$protein, second_df$protein)),
        n_overlapping_cv = nrow(dplyr::inner_join( dplyr::select(first_df, SubjectID, protein),
        dplyr::select(second_df, SubjectID, protein), by = c("SubjectID", "protein")
      )),
        wilcoxon_w = if (is.null(test)) {
        NA_real_
      } else {
        unname(test$statistic)
      }, raw_p = raw_p, p_adjust_method = "none", status = status, test_error = test_error, stringsAsFactors = FALSE
      )
    }))
  })) %>% dplyr::mutate(star = dplyr::case_when( !is.finite(raw_p) ~ "", raw_p < 0.001 ~ "***",
    raw_p < 0.01 ~ "**", raw_p < 0.05 ~ "*", TRUE ~ ""
  ))
}
cv_only_olink_tests_df <- cv_only_wilcoxon(cv_only_olink_df, cv_only_olink_pairs)

cv_only_density_polygon <- function(df, ylim, width = 0.9) {
  group_levels <- levels(df$group)
  density_list <- lapply(split(df, list(df$platform, df$group), drop = TRUE), function(group_df) {
    values <- group_df$cv_pct[is.finite(group_df$cv_pct)]
    if (length(values) < 2L) {
      return(NULL)
    }
    bandwidth <- stats::bw.nrd0(values)
    full_grid <- seq(min(values) - 3 * bandwidth, max(values) + 3 * bandwidth, length.out = 512L)
    visible_grid <- seq(ylim[1], ylim[2], length.out = 2048L)
    y_grid <- sort(unique(c(full_grid, visible_grid)))
    density_values <- numeric(length(y_grid))

    for (first in seq.int(1L, length(y_grid), by = 64L)) {
      idx <- seq.int(first, min(first + 63L, length(y_grid)))
      density_values[idx] <- rowMeans(stats::dnorm(outer(y_grid[idx], values, "-"), sd = bandwidth))
    }
    data.frame( platform = group_df$platform[1], group = group_df$group[1], y = y_grid, density = density_values
    )
  })
  density_df <- dplyr::bind_rows(density_list)
  density_df <- density_df %>% dplyr::group_by(platform) %>%
    dplyr::mutate(half_width = density / max(density) * width / 2) %>% dplyr::ungroup()
  dplyr::bind_rows(lapply( split(density_df, list(density_df$platform, density_df$group), drop = TRUE),
    function(group_df) {
    center <- match(as.character(group_df$group[1]), group_levels)
    dplyr::bind_rows( dplyr::mutate(group_df, x = center - half_width),
      dplyr::mutate(group_df[nrow(group_df):1, ], x = center + half_width)
    )
  })) %>%
    dplyr::mutate(group = factor(group, levels = group_levels))
}
cv_only_make_violin <- function(df, summary_df, fill_cols, facet_platform = FALSE, ylim = c(0, 200), tests_df = NULL,
  log2_display = FALSE, y_breaks = NULL) {
  density_polygon_df <- cv_only_density_polygon(df, ylim = ylim)
  if (log2_display) {
    density_polygon_df <- density_polygon_df %>% dplyr::filter(y >= ylim[1], y <= ylim[2])
  }
  p <- ggplot(dplyr::filter(df, is.finite(cv_pct)), aes(x = group, y = cv_pct, fill = group)) +
    geom_polygon( data = density_polygon_df, inherit.aes = FALSE,
    aes(x = x, y = y, fill = group, group = interaction(platform, group)), colour = NA, alpha = 1
  ) +
    geom_boxplot(width = 0.12, linewidth = 0.25, colour = "white", outlier.shape = NA, alpha = 0.5) + geom_text(
    data = summary_df, inherit.aes = FALSE, aes(x = group, y = Inf, label = median_label), size = 2.2, vjust = 1.3
  ) +
    scale_fill_manual(values = fill_cols, guide = "none") + scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = if (is.null(y_breaks)) {
    seq(ylim[1], ylim[2], length.out = 5)
  } else {
    y_breaks
  }) +
    labs(x = NULL, y = "Within-individual CV (%)") + theme2_nobacknoerror + themelegendnone + theme(
    text = element_text(size = 8), axis.text = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5), axis.title = element_text(size = 8, face = "plain"),
    strip.text = element_text(size = 9, face = "bold"), strip.background = element_blank(), panel.grid = element_blank(),
    plot.caption = element_text(size = 6, hjust = 0), plot.margin = unit(c(1, 1, 1, 1), "mm"),
    panel.spacing = unit(0.8, "cm")
  )
  if (log2_display) {
    log2_cv_display <- scales::trans_new( "log2_cv_display", transform = function(x) log2(pmax(x, .Machine$double.xmin)),
      inverse = function(x) 2^x, domain = c(0, Inf)
    )
    p <- p + coord_trans(y = log2_cv_display, ylim = ylim, expand = FALSE, clip = "on")
  } else {
    p <- p + coord_cartesian(ylim = ylim, expand = FALSE, clip = "on")
  }
  if (facet_platform) {
    p <- p + facet_wrap(~platform, nrow = 1, drop = FALSE)
  }
  if (!is.null(tests_df)) {
    sig_df <- tests_df %>% dplyr::filter(status == "ok", star != "") %>% dplyr::arrange(platform, pair_idx) %>%
      dplyr::group_by(platform) %>%
      dplyr::mutate( x1 = match(group1, levels(df$group)), x2 = match(group2, levels(df$group)), xmid = (x1 + x2) / 2,
      y_line = ylim[2] - diff(ylim) * 0.09 - (dplyr::n() - dplyr::row_number()) * diff(ylim) * 0.045,
      y_star = y_line - diff(ylim) * 0.004
    ) %>% dplyr::ungroup()
    if (nrow(sig_df)) {
      p <- p +
        geom_segment( data = sig_df, inherit.aes = FALSE, colour = "black", linewidth = 0.3,
        aes(x = x1, xend = x2, y = y_line, yend = y_line)
      ) +
        geom_text( data = sig_df, inherit.aes = FALSE, colour = "black", size = 3,
        hjust = 0.5, vjust = 0, aes(x = xmid, y = y_star, label = star)
      )
    }
  }
  p
}
cv_only_olink_fill <- setNames(c("#475272", "#68A39D", "#ADC482", "#D9B381"), cv_only_olink_labels)
cv_only_soma_fill <- c( "All SomaScan\ndetected proteins" = "#69B190FF",
  "Olink-Explore\nshared detected\nproteins" = "#264653", "Olink-Explore\nshared undetected\nproteins" = "#A8B5BA",
  "Olink-HT\nshared detected\nproteins" = "#AEC7E8FF", "Olink-HT\nshared undetected\nproteins" = "#DFE9F6"
)
p_cv_only_olink <- cv_only_make_violin( cv_only_olink_df, cv_only_olink_summary_df, cv_only_olink_fill,
  facet_platform = TRUE, tests_df = cv_only_olink_tests_df
)

ggplot2::ggsave(paste0(wdplot, "S5a_olink_within_individual_cv_violin.pdf"), pdf_text_grob(p_cv_only_olink),
  width = 20, height = 7, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)

# 18. Plot SomaScan CV on a log2 axis starting at 2% --------------------------

cv_soma_from_2_ylim <- c(2, 100)
cv_soma_from_2_breaks <- c(2, 4, 8, 16, 32, 64, 100)
cv_soma_from_2_pairs <- combn(names(cv_only_soma_sets), 2, simplify = FALSE)
cv_soma_from_2_tests_df <- cv_only_wilcoxon(cv_only_soma_df, cv_soma_from_2_pairs)
cv_only_soma_fill <- c( "All SomaScan\ndetected proteins" = "#69B190FF",
  "Olink-Explore\nshared detected\nproteins" = "#264653", "Olink-Explore\nshared undetected\nproteins" = "#A8B5BA",
  "Olink-HT\nshared detected\nproteins" = "#AEC7E8FF", "Olink-HT\nshared undetected\nproteins" = "#DFE9F6"
)
p_cv_only_somascan_from_2 <- cv_only_make_violin( cv_only_soma_df, cv_only_soma_summary_df, cv_only_soma_fill,
  ylim = cv_soma_from_2_ylim, log2_display = TRUE, y_breaks = cv_soma_from_2_breaks
)

cv_soma_from_2_plot_file <- paste0(wdplot, "S5b_somascan_within_individual_cv_violin_log2_from_2.pdf")
cv_soma_from_2_timestamp <- format(Sys.time(), "%Y%m%dT%H%M%S")

ggplot2::ggsave(cv_soma_from_2_plot_file, pdf_text_grob(p_cv_only_somascan_from_2), width = 15, height = 7, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)
cv_soma_from_2_tests_export_df <- cv_soma_from_2_tests_df %>%
  dplyr::mutate( platform = dplyr::recode(platform, "Somascan" = "SomaScan"),
  group1 = gsub("\n", " ", group1, fixed = TRUE), group2 = gsub("\n", " ", group2, fixed = TRUE),
  method = "Two-sided unpaired Wilcoxon rank-sum; exact=FALSE; correct=TRUE", observation_unit = "individual x protein",
  independence_limitation = "Repeated individuals and overlapping protein sets; not accounted for"
)

# Correlations between SOMAmer targets ---------------------------------------
plot_prefix <- "05.6-"
protein_set_key <- "detected"
col_platform_internal <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
  "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF", "Overlap" = "#D1C1E1FF"
)
platform_display_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral",
  "Olink-Explore" = "Olink-Explore", "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan", "Overlap" = "Overlap"
)

# 01. Load matched detection and SOMAmer inputs -------------------------------

detection_sets_file <- file.path(repo_dir, "result", "Fig2",
  "01-01_detection_classification_sets.rds"
)

detection_sets <- readRDS(detection_sets_file)
detected_sets <- detection_sets[[protein_set_key]]
dflistna <- read_processed_matrices(wddata)

if ("Olink-3072" %in% names(dflistna) && !"Olink-Explore" %in% names(dflistna)) {
  names(dflistna)[names(dflistna) == "Olink-3072"] <- "Olink-Explore"
}
load_registered_assays(wddata)

required_platforms <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
external_platforms <- setdiff(required_platforms, "Somascan")
for (platform in required_platforms) {
  platform_matrix <- dflistna[[platform]]
}
soma_map_df <- somer_protein %>%
  dplyr::transmute( somamer_id = as.character(.data$SeqId), protein = as.character(.data$EntrezGeneSymbol)
) %>% dplyr::filter(!is.na(.data$somamer_id), .data$somamer_id != "", !is.na(.data$protein), .data$protein != "") %>%
  dplyr::distinct()

ambiguous_somamer_ids <- soma_map_df %>% dplyr::count(.data$somamer_id, name = "n_proteins") %>%
  dplyr::filter(.data$n_proteins > 1)
soma_map_df <- soma_map_df %>% dplyr::distinct(.data$somamer_id, .keep_all = TRUE)

soma_matrix <- dflistna[["Somascan"]]
soma_detected_proteins <- intersect( as.character(detected_sets[["Somascan"]]), unique(soma_map_df$protein)
)
soma_target_df <- soma_map_df %>%
  dplyr::filter(.data$protein %in% soma_detected_proteins, .data$somamer_id %in% rownames(soma_matrix))
make_correlation_record <- function(x, y, base_fields) {
  correlation_result <- safe_spearman(x, y, min_n = 3)
  analysis_status <- if (correlation_result$n < 3) {
    "insufficient_paired_samples"
  } else if (!is.finite(correlation_result$rho)) {
    "constant_or_nonfinite_measurement"
  } else {
    "ok"
  }
  dplyr::bind_cols( tibble::as_tibble(base_fields), tibble::tibble(
    spearman_correlation_coefficient = correlation_result$rho, raw_p_value = correlation_result$p,
    paired_sample_count = correlation_result$n, analysis_status = analysis_status
  )
  )
}
same_protein_records <- list()
same_protein_index <- 0L
soma_targets_by_protein <- split(soma_target_df$somamer_id, soma_target_df$protein)

for (protein_name in names(soma_targets_by_protein)) {
  somamer_ids <- sort(unique(soma_targets_by_protein[[protein_name]]))
  if (length(somamer_ids) < 2) {
    next
  }
  somamer_pairs <- utils::combn(somamer_ids, 2, simplify = FALSE)
  for (somamer_pair in somamer_pairs) {
    shared_samples <- colnames(soma_matrix)
    same_protein_index <- same_protein_index + 1L
    same_protein_records[[same_protein_index]] <- make_correlation_record( soma_matrix[somamer_pair[1], shared_samples],
      soma_matrix[somamer_pair[2], shared_samples], list( protein = protein_name, somamer_a = somamer_pair[1],
      somamer_b = somamer_pair[2], statistical_unit = "SOMAmer pair within the same protein"
    )
    )
  }
}
same_protein_correlation_df <- dplyr::bind_rows(same_protein_records)
same_protein_summary_df <- same_protein_correlation_df %>% dplyr::group_by(.data$protein) %>% dplyr::summarise(
  n_somamers = length(unique(c(.data$somamer_a, .data$somamer_b))), n_somamer_pairs = dplyr::n(),
  n_valid_correlations = sum(.data$analysis_status == "ok"),
  median_spearman_correlation_coefficient = if (any(.data$analysis_status == "ok")) {
  stats::median(.data$spearman_correlation_coefficient[.data$analysis_status == "ok"])
} else {
  NA_real_
},
  min_spearman_correlation_coefficient = if (any(.data$analysis_status == "ok")) {
  min(.data$spearman_correlation_coefficient[.data$analysis_status == "ok"])
} else {
  NA_real_
},
  max_spearman_correlation_coefficient = if (any(.data$analysis_status == "ok")) {
  max(.data$spearman_correlation_coefficient[.data$analysis_status == "ok"])
} else {
  NA_real_
}, .groups = "drop"
)

same_protein_plot_df <- dplyr::filter(same_protein_correlation_df, analysis_status == "ok")
# 03. Plot correlations between targets of the same protein -------------------

p_same_protein_target_correlation <- ggplot2::ggplot( same_protein_plot_df,
  ggplot2::aes(x = "SomaScan target pairs", y = .data$spearman_correlation_coefficient)
) + ggplot2::geom_violin( fill = unname(col_platform_internal["Somascan"]), colour = "transparent",
  alpha = 1, width = 0.8, trim = FALSE
) +
  ggplot2::geom_boxplot( width = 0.22, size = 0.2, colour = "white", outlier.shape = NA, alpha = 0.6
) +

  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") + ggplot2::coord_cartesian(ylim = c(-1, 1)) +
  ggplot2::labs( x = NULL, y = "Spearman correlation coefficient", subtitle = paste0( "Median = ",
  sprintf("%.2f", stats::median(same_protein_plot_df$spearman_correlation_coefficient))
)
) + theme2_nobacknoerror + ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5)) + themelegendnone

ggsave( paste0(wdplot, "F4d_somascan_same_protein_target_correlation_distribution.pdf"),
  pdf_text_grob(p_same_protein_target_correlation), width = 6, height = 7, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)

saveRDS(list(cv = list(detected = detected_cv_df, undetected = undetected_cv_df, individual = individual_full_cv, broad = broad_full_cv, detected_overlap = detected_overlap_cv_df, individual_overlap = individual_overlap_cv, broad_overlap = broad_overlap_cv),
  cv_comparison = protein_cv_df,
  sample_correlation = samplewise_corrbind_plot_df,
  overlap_sample_correlation = samplewise_overlap_sets_plot_df,
  within_individual_cv = list(olink = cv_only_olink_df, somascan = cv_only_soma_df),
  cv_tests = list(olink = cv_only_olink_tests_df, somascan = cv_soma_from_2_tests_df),
  somamer_pair_correlation = same_protein_correlation_df),
  paste0(wdresult, "variation_results.rds"))

