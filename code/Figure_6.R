# Figure 6: complete analysis and manuscript figures
# Run from the repository root.
source(file.path("code", "setup.R"), encoding = "UTF-8")
wdresult <- paste0(file.path(repo_dir, "result", "Fig6"), "/")
wdplot <- paste0(file.path(repo_dir, "plot", "Fig6"), "/")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

is_08_1_wrapper_run <- FALSE
scope_key <- "detected"

options(stringsAsFactors = FALSE)

library(ggbreak)
plot_prefix <- "08-"
in08 <- wdresult

db06_display_effects <- c("Individual", "Sex")
db06_display_pvals <- c("indi_pval", "Sex_pval")

load_registered_assays(wddata)
ot480 <- (diaori1)
astral <- (astori1)


rownames(olink3072)[rownames(olink3072) == "NTproBNP"] <- "NT-proBNP"
rownames(olink3072)[rownames(olink3072) == "WARS"] <- "WARS1"

detection_sets_file <- file.path(repo_dir, "result", "Fig2", "01-01_detection_classification_sets.rds")
detection_sets <- readRDS(detection_sets_file)
prolist <- detection_sets[[scope_key]]
scope_label <- switch( scope_key, detected = "Detected proteins", individual_detected = "Individual detected proteins",
  broad_detected = "Broadly detected proteins"
)
dflistna <- read_processed_matrices(wddata)
names(dflistna)[names(dflistna) == "Olink-3072"] <- "Olink-Explore"
overlap5 <- Reduce(intersect, unname(prolist))

metadata <- read_processed_metadata(wddata)
metadata <- metadata[1:30,]
metadata$somaid <- paste0("S",metadata$OriginalID)
metadata$Time <- dplyr::recode(metadata$Time, "visit1" = "Day0", "visit2" = "Day3", "visit3" = "Day5")
metadata$Time <- factor(metadata$Time, levels = c("Day0", "Day3", "Day5"))
metadata$SubjectID <- factor(metadata$SubjectID ,levels= unique(metadata$SubjectID ))

# 01. Prepare imputed expression matrices -------------------------------------

c8_imputed_matrix_list <- list()

c8_imputed_matrix_list[["MS-OT480"]] <- msmatrixpro(pro = prolist[["MS-OT480"]], cate = "MS-OT480",
  data = dflistna[["MS-OT480"]], isna = "impute")

c8_imputed_matrix_list[["MS-Astral"]] <- msmatrixpro(pro = prolist[["MS-Astral"]], cate = "MS-Astral",
  data = dflistna[["MS-Astral"]], isna = "impute")

c8_imputed_matrix_list[["Olink-Explore"]] <- olinkmatrixpro(pro = prolist[["Olink-Explore"]], cate = "Olink-Explore",
  data = dflistna[["Olink-Explore"]], isna = "impute")

c8_imputed_matrix_list[["Olink-HT"]] <- olinkmatrixpro(pro = prolist[["Olink-HT"]], cate = "Olink-HT",
  data = dflistna[["Olink-HT"]], isna = "impute")

c8_imputed_matrix_list[["Somascan"]] <- somamatrixpro(pro = prolist[["Somascan"]], cate = "Somascan",
  data = dflistna[["Somascan"]], isna = "impute")

# 02. Compare expression structure with t-SNE ---------------------------------

fig_id <- "C8_1"
catename <- "Within-platform tSNE of imputed protein matrices"

p_c8_1_tsne_list <- list()

for (i in names(c8_imputed_matrix_list)) {
  p_c8_1_tsne_list[[i]] <- tsnebind(c8_imputed_matrix_list[[i]], per = 7, til = i)
}

donor_display_map <- metadata %>% dplyr::distinct(SubjectID) %>%
  dplyr::mutate(subject_numeric = suppressWarnings(as.numeric(sub("^Subject", "", as.character(SubjectID))))) %>%
  dplyr::arrange(subject_numeric) %>%
  dplyr::mutate(DonorLabel = sprintf("Donor %02d", dplyr::row_number()))
tsne_donor_cols <- c("#D4C96E", "#807961", "#8DACB6", "#7D839C", "#CA9791",
  "#55794D", "#EFF3F9", "#924741", "#E0DC96", "#DA9F58")
for (i in names(p_c8_1_tsne_list)) {
  p_c8_1_tsne_list[[i]] <- p_c8_1_tsne_list[[i]] +
    scale_fill_manual(values = stats::setNames(tsne_donor_cols, donor_display_map$SubjectID),
    breaks = donor_display_map$SubjectID, labels = sub("^Donor ", "Individual ", donor_display_map$DonorLabel),
    name = "Individual", guide = guide_legend(override.aes = list(
      shape = 21, colour = "#525A73", size = 3, alpha = 1))) +
    scale_shape_manual(values = c(F = 24, M = 21), name = "Sex") +
    scale_size_manual(values = c(Day0 = 2.2, Day3 = 2.2, Day5 = 2.2), guide = "none") +
    labs(x = "t-SNE 1", y = "t-SNE 2")
}
p_c8_1_tsne_combined <- p_c8_1_tsne_list[[1]] +
  p_c8_1_tsne_list[[2]] + p_c8_1_tsne_list[[3]] + p_c8_1_tsne_list[[4]] + (p_c8_1_tsne_list[[5]] + themelegendright) +
  plot_layout(widths = c(1, 1, 1, 1, 1))

ggsave(paste0(wdplot, "F6a_within_platform_tsne_imputed_protein_matrices.pdf"), pdf_text_grob(p_c8_1_tsne_combined), width = 32,
  height = 6.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

p_c8_1_legend <-  p_c8_1_tsne_list[[1]] + p_c8_1_tsne_list[[2]] + p_c8_1_tsne_list[[3]] +themelegendtop+
  p_c8_1_tsne_list[[4]] + (p_c8_1_tsne_list[[5]] ) + plot_layout(widths = c(1, 1, 1, 1, 1))
ggsave(paste0(wdplot, "F6a_legend.pdf"), pdf_text_grob(p_c8_1_legend), width = 30, height = 6.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 03. Fit individual and sex association models -------------------------------

build_dep_protein_evidence <- function(dep_df, platform, fdr_cutoff = 0.05) {
  required_cols <- c("Protein", "somer", "Time_pval", "Sex_pval", "indi_pval", "Time_FDR", "Sex_FDR", "indi_FDR")
  missing_cols <- setdiff(required_cols, colnames(dep_df))
  feature_df <- dep_df %>%
    dplyr::transmute( platform = platform, Protein = as.character(Protein), Somer = dplyr::na_if(as.character(somer), ""),
    feature_id = ifelse(platform == "Somascan", as.character(somer), as.character(Protein)),
    Time_pval, Sex_pval, indi_pval, Time_FDR, Sex_FDR, indi_FDR
  ) %>% dplyr::filter(!is.na(Protein), nzchar(Protein), !is.na(feature_id), nzchar(feature_id))

  feature_long <- dplyr::bind_rows( dplyr::transmute(feature_df, platform, Protein, Somer, feature_id,
    effect = "Time", p_value = Time_pval, FDR = Time_FDR), dplyr::transmute(feature_df, platform, Protein, Somer, feature_id,
    effect = "Sex", p_value = Sex_pval, FDR = Sex_FDR), dplyr::transmute(feature_df, platform, Protein, Somer, feature_id,
    effect = "Individual", p_value = indi_pval, FDR = indi_FDR)
  )

  feature_long %>%
    dplyr::group_by(platform, effect, Protein) %>% dplyr::summarise( n_features_tested = dplyr::n_distinct(feature_id),
    n_fdr_significant = sum(is.finite(FDR) & FDR < fdr_cutoff),
    n_nominal_significant = sum(is.finite(p_value) & p_value < 0.05),
    min_p_value = if (all(!is.finite(p_value))) {
    NA_real_
  } else {
    min(p_value, na.rm = TRUE)
  },
    min_FDR = if (all(!is.finite(FDR))) {
    NA_real_
  } else {
    min(FDR, na.rm = TRUE)
  }, feature_ids = paste(sort(unique(feature_id)), collapse = ";"),
    significant_feature_ids = paste(sort(unique(feature_id[is.finite(FDR) & FDR < fdr_cutoff])), collapse = ";"),
    evidence_class = dplyr::case_when( n_fdr_significant == 0 ~ "no FDR-significant feature",
    n_features_tested == 1 ~ "single-feature evidence",
    n_fdr_significant == n_features_tested ~ "all tested features FDR-significant", TRUE ~ "partial multi-feature evidence"
  ), .groups = "drop"
  )
}
fig_id <- "C8_2"
catename <- "Mixed-effects model based DEP identification"

coefficient_or_na <- function(name) {
  if (name %in% names(coefficients)) {
    unname(coefficients[[name]])
  } else {
    NA_real_
  }
}
for (i in names(c8_imputed_matrix_list)) {
  c8_dep_matrix <- c8_imputed_matrix_list[[i]]
  c8_dep_long_df <- melt(as.matrix(c8_dep_matrix))
  colnames(c8_dep_long_df) <- c("Var1", "Var2", "value")
  c8_dep_long_df <- left_join(c8_dep_long_df, metadata, by = c("Var2" = "OriginalID"))
  colnames(c8_dep_long_df) <- c("Protein", "OriginalID", "Expression", "colnames(olink1)", "SubjectID", "Time", "Age",
    "Sex", "Individual", "age", "somaid")
  c8_dep_long_df$Time <- as.factor(c8_dep_long_df$Time)
  c8_dep_long_df$Sex <- as.factor(c8_dep_long_df$Sex)
  c8_dep_proteins <- unique(c8_dep_long_df$Protein)
  c8_dep_result_list <- list()
  c8_logfc_result_list <- list()

  for (j in c8_dep_proteins) {
    c8_dep_protein_df <- c8_dep_long_df %>% dplyr::filter(Protein == j)

    c8_dep_lmm <- lmer(Expression ~ Time + Sex + (1 | SubjectID),data = c8_dep_protein_df,REML = FALSE)

    c8_dep_random_anova <- ranova(c8_dep_lmm)
    c8_dep_fixed_anova <- anova(c8_dep_lmm, type = 3)
    coefficients <- lme4::fixef(c8_dep_lmm)
    c8_logfc_result_list[[j]] <- data.frame(
      log2FC_M_vs_F = coefficient_or_na("SexM"),
      log2FC_visit2_vs_visit1 = coefficient_or_na("Timevisit2"),
      log2FC_visit3_vs_visit1 = coefficient_or_na("Timevisit3"),
      Time_pval = c8_dep_fixed_anova["Time", "Pr(>F)"],
      Sex_pval = c8_dep_fixed_anova["Sex", "Pr(>F)"],
      indi_pval = c8_dep_random_anova["(1 | SubjectID)", "Pr(>Chisq)"], Protein = j)

    c8_dep_result_list[[j]] <- data.frame( Time_pval = c8_dep_fixed_anova["Time", "Pr(>F)"],
      Sex_pval = c8_dep_fixed_anova["Sex", "Pr(>F)"], indi_pval = c8_dep_random_anova["(1 | SubjectID)", "Pr(>Chisq)"],
      Protein = j
    )
  }
  c8_dep_result_df <- do.call(rbind, c8_dep_result_list)
  saveRDS(do.call(rbind, c8_logfc_result_list),
    paste0(wdresult, "C8_3_0_", i, "_dep_log2fc_results.rds"))

  saveRDS( c8_dep_result_df, file = paste0(wdresult, fig_id, "_", i, "_mixed_effects_dep_results.rds")
  )
}

# 2.1 DEP counts with FDR correction ------------------------------------------

fig_id <- "C8_2_1"
catename <- "DEP counts with FDR correction"

c8_dep_list_fdr <- list()
c8_dep_list_raw <- list()
c8_dep_evidence_list <- list()

for (i in names(c8_imputed_matrix_list)) {
  c8_dep_result_df <- readRDS(paste0(wdresult, "C8_2_", i, "_mixed_effects_dep_results.rds"))
  c8_dep_result_df <- c8_dep_result_df %>% dplyr::mutate(Time_FDR = p.adjust(Time_pval, method = "BH"),
    Sex_FDR  = p.adjust(Sex_pval,  method = "BH"),indi_FDR = p.adjust(indi_pval, method = "BH"))

  if (i == "Somascan") {
    colnames(c8_dep_result_df) <- c("Time_pval", "Sex_pval", "indi_pval", "somer","Time_FDR", "Sex_FDR", "indi_FDR")
    somer_protein2 <- somer_protein[, c("SeqId", "EntrezGeneSymbol")]
    colnames(somer_protein2) <- c("SeqId", "Protein")
    c8_dep_result_df <- left_join(c8_dep_result_df, somer_protein2,by = c("somer" = "SeqId"))
  } else {
    c8_dep_result_df$somer <- ""
  }
  c8_dep_evidence_list[[i]] <- build_dep_protein_evidence(c8_dep_result_df, platform = i)
}
c8_dep_protein_evidence <- dplyr::bind_rows(c8_dep_evidence_list)
for (i in names(c8_imputed_matrix_list)) {
  for (effect_name in c("Time", "Sex", "Individual")) {
    evidence_sub <- dplyr::filter(c8_dep_protein_evidence, platform == i, effect == effect_name)
    c8_dep_list_fdr[[paste0(i, "_", effect_name, "_DEP")]] <-
      unique(evidence_sub$Protein[evidence_sub$n_fdr_significant >= 1])
    c8_dep_list_raw[[paste0(i, "_", effect_name, "_DEP")]] <-
      unique(evidence_sub$Protein[evidence_sub$n_nominal_significant >= 1])
  }
}
saveRDS(c8_dep_protein_evidence, paste0(wdresult, "C8_2_1_dep_protein_evidence.rds"))

c8_dep_summary_df <- data.frame(Category = names(c8_dep_list_raw),GeneCount = lengths(c8_dep_list_raw),
  GeneCountFDR  = lengths(c8_dep_list_fdr),
  stringsAsFactors = FALSE)

c8_dep_summary_df$category1 <- sub("_.*", "", c8_dep_summary_df$Category)
c8_dep_summary_df$category2 <- sub("^[^_]*_", "", c8_dep_summary_df$Category)

c8_dep_summary_df$category1 <- factor(c8_dep_summary_df$category1,levels = unique(c8_dep_summary_df$category1))
c8_dep_summary_df$category2 <- factor(c8_dep_summary_df$category2,levels = unique(c8_dep_summary_df$category2))
# 2.2 DEP dot plot ------------------------------------------------------------

fig_id <- "C8_2_2"

catename <- "DEP dot plot"

c8_dep_bind_df <- NULL
for (i in names(c8_imputed_matrix_list)) {
  c8_dep_result_df <- readRDS(paste0(wdresult, "C8_2_", i, "_mixed_effects_dep_results.rds"))
  c8_dep_result_df$cate <- i
  c8_dep_result_df <- c8_dep_result_df %>% dplyr::mutate(Time_pval=p.adjust(Time_pval,"BH"), Sex_pval=p.adjust(Sex_pval,
    "BH"), indi_pval=p.adjust(indi_pval,"BH"))

  if (i == "Somascan") {
    colnames(c8_dep_result_df) <- c("Time_pval","Sex_pval","indi_pval","somer","cate")
    somer_protein2 <- somer_protein[,c("SeqId","EntrezGeneSymbol")]

    colnames(somer_protein2) <- c("SeqId","Protein")
    c8_dep_result_df <- left_join(c8_dep_result_df, somer_protein2, by=c("somer"="SeqId"))
  } else {
    c8_dep_result_df$somer <- ""
  }
  c8_dep_bind_df <- rbind(c8_dep_bind_df, c8_dep_result_df)
}
c8_dep_bind_df <- melt(c8_dep_bind_df)
c8_dep_bind_df <- dplyr::filter(c8_dep_bind_df, value < 0.05)
c8_dep_bind_df$value <- -log10(c8_dep_bind_df$value)

c8_dep_bind_df <- c8_dep_bind_df %>% dplyr::group_by(cate, variable) %>% dplyr::arrange(desc(value), .by_group = TRUE) %>%
  dplyr::mutate(top5_Protein = if_else(row_number() <= 5, as.character(Protein), "")) %>%
  dplyr::ungroup()

c8_dep_bind_df$cate <- factor(c8_dep_bind_df$cate, level=c("MS-OT480","MS-Astral","Olink-Explore","Olink-HT","Somascan"))
c8_dep_bind_df$variable <- factor(c8_dep_bind_df$variable, level=c("indi_pval","Time_pval","Sex_pval"))
c8_dep_display_dot_df <- c8_dep_bind_df %>% dplyr::filter(variable %in% db06_display_pvals)
# 2.3 DEP VENN ----------------------------------------------------------------

fig_id <- "C8_2_3"

catename <- "DEP overlap across platforms"

library(VennDiagram)

library(gridExtra)

c8_dep_bind_df$variable <- dplyr::case_when( c8_dep_bind_df$variable == "Time_pval" ~ "Time",
  c8_dep_bind_df$variable == "Sex_pval" ~ "Sex", TRUE ~ "Individual"
)
c8_dep_display_df <- c8_dep_bind_df %>% dplyr::filter(variable %in% db06_display_effects)
# 2.5 shared DEP lists --------------------------------------------------------

fig_id <- "C8_2_5"

catename <- "Shared DEP lists"

c8_dep_individual_df <- dplyr::filter(c8_dep_bind_df, variable == "Individual")
c8_dep_individual_list <- list()
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_individual_list[[j]] <- na.omit(unique(c8_dep_individual_df$Protein[c8_dep_individual_df$cate == j]))
}
c8_dep_individual_shared <- Reduce(intersect, c8_dep_individual_list[1:5])

c8_dep_time_df <- dplyr::filter(c8_dep_bind_df, variable == "Time")
c8_dep_time_list <- list()
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_time_list[[j]] <- na.omit(unique(c8_dep_time_df$Protein[c8_dep_time_df$cate == j]))
}
c8_dep_time_shared <- Reduce(intersect, c8_dep_time_list[1:5])

c8_dep_sex_df <- dplyr::filter(c8_dep_bind_df, variable == "Sex")
c8_dep_sex_list <- list()
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_sex_list[[j]] <- na.omit(unique(c8_dep_sex_df$Protein[c8_dep_sex_df$cate == j]))
}
c8_dep_sex_overlap_msot480_olinkexplore <- Reduce(intersect, list(c8_dep_sex_list[["MS-OT480"]], c8_dep_sex_list[["Olink-Explore"]]))
c8_dep_sex_overlap_somascan_olinkexplore <- Reduce(intersect, list(c8_dep_sex_list[["Somascan"]], c8_dep_sex_list[["Olink-Explore"]]))
c8_dep_sex_overlap_msot480_somascan <- Reduce(intersect, list(c8_dep_sex_list[["MS-OT480"]], c8_dep_sex_list[["Somascan"]]))

# 06. Compare individual and sex-associated protein overlaps ------------------

fig_id <- "C8_6_1"

catename <- "Overlap DEPs heatmap"

c8_dep_individual_heatmap_list <- list()
c8_dep_individual_df <- dplyr::filter(c8_dep_bind_df, variable == "Individual")
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_individual_heatmap_sub_df <- dplyr::filter(c8_dep_individual_df, cate == j)
  c8_dep_individual_heatmap_list[[paste0(j, "_Individual DEPs")]] <- na.omit(unique(c8_dep_individual_heatmap_sub_df$Protein))
}
c8_dep_time_heatmap_list <- list()
c8_dep_time_df <- dplyr::filter(c8_dep_bind_df, variable == "Time")
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_time_heatmap_sub_df <- dplyr::filter(c8_dep_time_df, cate == j)
  c8_dep_time_heatmap_list[[paste0(j, "_Time DEPs")]] <- na.omit(unique(c8_dep_time_heatmap_sub_df$Protein))
}
c8_dep_sex_heatmap_list <- list()
c8_dep_sex_df <- dplyr::filter(c8_dep_bind_df, variable == "Sex")
for (j in names(c8_imputed_matrix_list)) {
  c8_dep_sex_heatmap_sub_df <- dplyr::filter(c8_dep_sex_df, cate == j)
  c8_dep_sex_heatmap_list[[paste0(j, "_Sex DEPs")]] <- na.omit(unique(c8_dep_sex_heatmap_sub_df$Protein))
}
c8_dep_heatmap_list <- list()
for (j in c("Individual", "Time", "Sex")) {
  c8_dep_sub_df <- dplyr::filter(c8_dep_bind_df, variable == j)
  for (i in names(c8_imputed_matrix_list)) {
    c8_dep_heatmap_list[[paste0(i, "_", j, "_DEPs")]] <- na.omit(unique(c8_dep_sub_df$Protein[c8_dep_sub_df$cate == i]))
  }
}
# S9a uses raw P <= 0.05; shared FDR-selected lists remain unchanged.
c8_dep_display_heatmap_list <- list()
for (effect_name in db06_display_effects) {
  for (platform_name in names(c8_imputed_matrix_list)) {
    c8_dep_display_heatmap_list[[paste0(platform_name, "_", effect_name, "_DEPs")]] <- c8_dep_protein_evidence %>%
      dplyr::filter(platform == platform_name, effect == effect_name, is.finite(min_p_value), min_p_value <= 0.05) %>%
      dplyr::pull(Protein) %>% unique() %>% na.omit()
  }
}
p_c8_6_1 <- overlapheatmap(c8_dep_display_heatmap_list, "") + themelegendright + themetext30
c8_overlap_names <- names(c8_dep_display_heatmap_list)
c8_overlap_labels <- paste0(c8_overlap_names, " (", lengths(c8_dep_display_heatmap_list), ")")
c8_overlap_tested_list <- list()
for (group_name in c8_overlap_names) {
  platform_name <- sub("_(Individual|Sex)_DEPs$", "", group_name)
  effect_name <- sub("^.*_(Individual|Sex)_DEPs$", "\\1", group_name)
  c8_overlap_tested_list[[group_name]] <- c8_dep_protein_evidence %>%
    dplyr::filter(platform == platform_name, effect == effect_name, is.finite(min_p_value)) %>%
    dplyr::pull(Protein) %>% unique()
}
c8_overlap_plot_df <- p_c8_6_1$data
c8_overlap_plot_df$p_overlap <- NA_real_
for (i in seq_len(length(c8_overlap_names) - 1L)) {
  for (j in seq.int(i + 1L, length(c8_overlap_names))) {
    background <- intersect(c8_overlap_tested_list[[i]], c8_overlap_tested_list[[j]])
    set_a <- intersect(c8_dep_display_heatmap_list[[i]], background)
    set_b <- intersect(c8_dep_display_heatmap_list[[j]], background)
    if (length(background) && length(set_a) && length(set_b)) {
      p_overlap <- stats::phyper(length(intersect(set_a, set_b)) - 1L, length(set_a),
        length(background) - length(set_a), length(set_b), lower.tail = FALSE)
      pair_cells <- (c8_overlap_plot_df$List1 == c8_overlap_labels[i] &
        c8_overlap_plot_df$List2 == c8_overlap_labels[j]) |
        (c8_overlap_plot_df$List1 == c8_overlap_labels[j] & c8_overlap_plot_df$List2 == c8_overlap_labels[i])
      c8_overlap_plot_df$p_overlap[pair_cells] <- p_overlap
    }
  }
}
c8_overlap_plot_df$overlap_label <- as.character(c8_overlap_plot_df$Overlap)
c8_overlap_plot_df$overlap_star <- ifelse(
  is.finite(c8_overlap_plot_df$p_overlap) & c8_overlap_plot_df$p_overlap < 0.05, "*", "")
p_c8_6_1$data <- c8_overlap_plot_df
p_c8_6_1 <- p_c8_6_1 + geom_text(data = c8_overlap_plot_df, aes(label = overlap_star), nudge_y = 0.17,
  size = 2.5, colour = "black", vjust = 0.5, hjust = 0.5, show.legend = FALSE) +
  labs(caption = "* P < 0.05 (one-sided hypergeometric test; unadjusted)")

ggsave(paste0(wdplot, "S9a_individual_sex_overlap_deps_heatmap.pdf"), pdf_text_grob(p_c8_6_1), width = 20, height = 10, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
model_bundle <- list(c8_imputed_matrix_list = c8_imputed_matrix_list, metadata = metadata,
  somer_protein = somer_protein, c8_dep_individual_list = c8_dep_individual_list,
  c8_dep_sex_list = c8_dep_sex_list, c8_dep_protein_evidence = c8_dep_protein_evidence)
saveRDS(model_bundle, paste0(wdresult, "model_inputs.rds"))
# 08.1. Associated-protein counts, effects and overlaps -----------------------

platform_keys <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
platform_labels <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "SomaScan")
names(platform_labels) <- platform_keys
donor_colours <- c("#D4C96E", "#807961", "#8DACB6", "#7D839C", "#CA9791",
  "#55794D", "#EFF3F9", "#924741", "#E0DC96", "#DA9F58")
raw_p_cutoff <- 0.05

cache <- model_bundle
model_matrices <- cache$c8_imputed_matrix_list
metadata <- cache$metadata
subject_keys <- unique(as.character(metadata$SubjectID))
subject_numeric <- suppressWarnings(as.numeric(sub("^Subject", "", subject_keys)))
donor_display_map <- data.frame( subject_key = subject_keys, subject_numeric = subject_numeric, stringsAsFactors = FALSE)
donor_display_map <- donor_display_map[order(donor_display_map$subject_numeric), , drop = FALSE]
donor_display_map$DonorLabel <- sprintf("Donor %02d", seq_len(nrow(donor_display_map)))
donor_palette <- stats::setNames(donor_colours, donor_display_map$DonorLabel)

evidence_file <- file.path(in08, "C8_2_1_dep_protein_evidence.rds")
effect_files <- file.path(in08, paste0("C8_3_0_", platform_keys, "_dep_log2fc_results.rds"))
model_files <- file.path(in08, paste0("C8_2_", platform_keys, "_mixed_effects_dep_results.rds"))
evidence <- readRDS(evidence_file)
individual_model_results <- stats::setNames(lapply(model_files, readRDS), platform_keys)
required_columns <- c("platform", "effect", "Protein", "n_features_tested",
  "n_nominal_significant", "n_fdr_significant", "min_p_value", "min_FDR", "feature_ids")
horizontal_pair_theme <- theme( text = element_text(family = "sans"),
  axis.title.x = element_text(size = 8, face = "bold", colour = "black"),
  axis.title.y = element_text(size = 8, face = "bold", colour = "black"),
  axis.text.x = element_text(size = 8, angle = 0, vjust = 0, hjust = 0.5, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"), plot.title = element_text(size = 9, face = "plain"),
  plot.subtitle = element_text(size = 8), panel.background = element_rect(fill = NA, colour = NA),
  legend.background = element_rect(fill = NA, colour = NA), legend.title = element_text(size = 8),
  legend.text = element_text(size = 8), legend.key.size = grid::unit(0.5, "cm"),
  panel.grid.major = element_line(colour = "#EDEDED", linewidth = 0.5),
  panel.grid.minor = element_line(colour = "#EDEDED", linewidth = 0.35),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35), strip.background = element_blank(),
  axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.25),
  axis.ticks.length.x = grid::unit(0.08, "cm"), axis.ticks.length.y = grid::unit(0.08, "cm")
)

save_08_1_comparisons <- function(plot_list, pair_index, has_points, section, effect) {
  selected_keys <- c("MS-OT480__MS-Astral", "MS-Astral__Olink-HT",
    "MS-Astral__Somascan", "Olink-Explore__Olink-HT", "Olink-HT__Somascan")
  all_keys <- vapply(pair_index, paste, character(1), collapse = "__")
  selected <- match(selected_keys, all_keys)
  remaining <- setdiff(seq_along(pair_index), selected)
  stem <- if (effect == "Sex") {
    "donor_mean_sex_log2fc"
  } else {
    "donor_vs_other_9_log2fc"
  }
  selected_prefix <- if (effect == "Sex") {
    "F6c"
  } else {
    "F6b"
  }
  remaining_prefix <- if (effect == "Sex") {
    "S8e"
  } else {
    "S8d"
  }
  axis_separator <- if (effect == "Sex") {
    " individual-mean"
  } else {
    " log₂FC"
  }
  plot_list <- lapply(plot_list, function(p) {
    p + horizontal_pair_theme + labs(
      x = if (!is.null(p$labels$x)) {
      sub(axis_separator, paste0("\n", trimws(axis_separator)), p$labels$x, fixed = TRUE)
    } else {
      NULL
    },
      y = if (!is.null(p$labels$y)) {
      sub(axis_separator, paste0("\n", trimws(axis_separator)), p$labels$y, fixed = TRUE)
    } else {
      NULL
    }) + theme(legend.position = "none", plot.caption = element_blank(), plot.margin = margin(1, 1, 1, 1, unit = "mm"))
  })
  p_selected <- cowplot::plot_grid(plotlist = plot_list[selected], ncol = 5, nrow = 1, align = "hv", axis = "tblr")
  p_remaining <- cowplot::plot_grid(plotlist = plot_list[remaining], ncol = 5, nrow = 1, align = "hv", axis = "tblr")
  ggsave(file.path(wdplot, paste0(selected_prefix, "_selected_pairwise_", stem, "_1row_5col.pdf")),
    pdf_text_grob(p_selected), width = 32, height = 6.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  ggsave(file.path(wdplot, paste0(remaining_prefix, "_remaining_pairwise_", stem, "_1row_5col.pdf")),
    pdf_text_grob(p_remaining), width = 32, height = 6.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  example_index <- selected[has_points[selected]][1]
  p_example <- plot_list[[example_index]]
  if (effect == "Sex") {
    p_example <- p_example + scale_shape_manual( values = c("Both" = 16, "Platform 1" = 17, "Platform 2" = 15),
      labels = c("Both technologies", "Horizontal technology", "Vertical technology"), drop = FALSE) +
      labs(colour = "Direction on horizontal technology", shape = "P < 0.05 in") +
      guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 2.2, alpha = 1)),
      shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.2, alpha = 1)))
  } else {
    p_example <- p_example + labs(fill = "Individual", shape = "Individual P < 0.05 in") +
      guides(fill = guide_legend(order = 1, ncol = 5, byrow = TRUE, override.aes = list(shape = 21, size = 2.2, alpha = 1)),
      shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.2, alpha = 1)))
  }
  p_example <- p_example + theme(legend.position = "top", legend.box = "vertical", legend.key.size = grid::unit(0.4, "cm"),
    legend.title = element_text(size = 8, face = "plain"), legend.text = element_text(size = 8),
    legend.background = element_blank(), legend.key = element_blank())
  ggsave(file.path(wdplot, paste0(selected_prefix, "_pairwise_", stem, "_top_legend_example.pdf")),
    pdf_text_grob(p_example), width = 16, height = 12, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  invisible(list(selected = p_selected, remaining = p_remaining, legend_example = p_example))
}
# 02. Compare nominal-P and FDR counts ----------------------------------------

required_columns <- c("platform", "effect", "Protein", "n_features_tested", "n_nominal_significant", "n_fdr_significant",
  "min_p_value", "min_FDR")
platform_levels <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
platform_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral",
  "Olink-Explore" = "Olink-Explore", "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan"
)
association_levels <- c("Individual", "Sex")
association_labels <- c( "Individual" = "Individual-associated proteins", "Sex" = "Sex-associated proteins"
)
count_colours <- c( "Individual | nominal P < 0.05" = "#E2D26D", "Individual | BH-FDR < 0.05" = "#D39D2A",
  "Sex | nominal P < 0.05" = "#A6A77A", "Sex | BH-FDR < 0.05" = "#4F6B56"
)
count_labels <- c( "Individual | nominal P < 0.05" = "Individual, P < 0.05",
  "Individual | BH-FDR < 0.05" = "Individual, FDR < 0.05", "Sex | nominal P < 0.05" = "Sex, P < 0.05",
  "Sex | BH-FDR < 0.05" = "Sex, FDR < 0.05"
)

count_df <- evidence %>% dplyr::filter(effect %in% association_levels) %>% dplyr::transmute(
  platform = factor(platform, levels = platform_levels), association = factor(effect, levels = association_levels),
  Protein = as.character(Protein), nominal_associated = n_nominal_significant >= 1L,
  fdr_associated = n_fdr_significant >= 1L
) %>% dplyr::filter(!is.na(Protein), nzchar(Protein)) %>% dplyr::group_by(platform, association) %>%
  dplyr::summarise( tested_proteins = dplyr::n(), nominal_associated_proteins = sum(nominal_associated),
  fdr_associated_proteins = sum(fdr_associated),
  nominal_associated_percentage = 100 * nominal_associated_proteins / tested_proteins,
  fdr_associated_percentage = 100 * fdr_associated_proteins / tested_proteins, .groups = "drop"
)

sex_associated_proteins <- evidence %>%
  dplyr::filter(effect == "Sex") %>% dplyr::transmute( platform = as.character(platform), Protein = as.character(Protein),
  nominal_P_lt_0_05 = n_nominal_significant >= 1L, BH_FDR_lt_0_05 = n_fdr_significant >= 1L,
  n_features_tested = as.integer(n_features_tested), n_nominal_significant = as.integer(n_nominal_significant),
  n_fdr_significant = as.integer(n_fdr_significant), minimum_raw_P = min_p_value, minimum_BH_FDR = min_FDR
) %>% dplyr::filter( platform %in% platform_levels, !is.na(Protein), nzchar(Protein), nominal_P_lt_0_05 | BH_FDR_lt_0_05
) %>%
  dplyr::mutate( significance_label = dplyr::case_when(
  nominal_P_lt_0_05 & BH_FDR_lt_0_05 ~ "Nominal P < 0.05 and BH-FDR < 0.05", BH_FDR_lt_0_05 ~ "BH-FDR < 0.05",
  nominal_P_lt_0_05 ~ "Nominal P < 0.05"
)
) %>% dplyr::select( platform, Protein, significance_label, nominal_P_lt_0_05, BH_FDR_lt_0_05,
  n_features_tested, n_nominal_significant, n_fdr_significant, minimum_raw_P, minimum_BH_FDR
) %>% dplyr::arrange( factor(platform, levels = platform_levels), dplyr::desc(BH_FDR_lt_0_05),
  dplyr::desc(nominal_P_lt_0_05), Protein
)

plot_df <- bind_rows( count_df %>% dplyr::transmute(platform, association, tested_proteins, threshold = "nominal P < 0.05",
  associated_proteins = nominal_associated_proteins, associated_percentage = nominal_associated_percentage),
  count_df %>% dplyr::transmute(platform, association, tested_proteins, threshold = "BH-FDR < 0.05",
  associated_proteins = fdr_associated_proteins, associated_percentage = fdr_associated_percentage)
) %>%
  dplyr::mutate( measure = factor( paste(association, threshold, sep = " | "), levels = names(count_colours)
)
)

pdf_file <- file.path(wdplot, "S8a_nominal_p_fdr_individual_sex_associated_protein_counts_ggbreak.pdf")

p_p_fdr_counts <- ggplot( plot_df, aes(x = platform, y = associated_proteins, fill = measure)
) + geom_col( position = position_dodge2(width = 0.84, preserve = "single"),
  width = 0.78, colour = "black", linewidth = 0.3
) +
  geom_text( aes(label = associated_proteins), position = position_dodge2(width = 0.84, preserve = "single"),
  vjust = 0.5, hjust = -0.08, angle = 90, size = 2.1
) +
  scale_fill_manual( values = count_colours, labels = count_labels, name = "Association and threshold",
  guide = guide_legend(nrow = 2, byrow = TRUE)
) +
  scale_x_discrete(labels = platform_labels) + scale_y_break(breaks = c(3000, 7000), space = 0.35, scales = 0.9) +
  scale_y_continuous( limits = c(0, 8500), breaks = seq(0, 8000, by = 2000),
  labels = scales::label_number(big.mark = ","), expand = c(0, 0)
) + labs( x = "Platform", y = "Number of proteins",
  caption = "Mixed-effects model results; protein counted when at least one feature meets the named threshold."
) +
  theme_classic(base_size = 10) + theme( legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1),
  plot.caption = element_text(hjust = 0, size = 8)
)

ggsave(pdf_file, pdf_text_grob(p_p_fdr_counts), width = 12, height = 10, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
message("08.1 nominal-P/FDR count plot complete: ", basename(pdf_file))

counts_results <- list(counts = count_df, plot_data = plot_df, sex_associated_proteins = sex_associated_proteins, plot = p_p_fdr_counts)

# 03. Calculate donor-mean sex effects and pairwise comparisons ---------------

required_evidence <- c("platform", "effect", "Protein", "min_p_value",
  "n_features_tested", "n_nominal_significant", "feature_ids")
sex_evidence <- evidence[evidence$effect == "Sex" &
  evidence$platform %in% platform_keys, required_evidence, drop = FALSE]
sex_evidence$Protein <- trimws(as.character(sex_evidence$Protein))
feature_rows <- vector("list", length(platform_keys))
for (k in seq_along(platform_keys)) {
  platform <- platform_keys[k]
  effect <- readRDS(effect_files[k])
  required_effect <- c("Protein", "Sex_pval")
  feature_id <- as.character(effect$Protein)
  feature_rows[[k]] <- data.frame( platform = platform, feature_id = feature_id, raw_p = as.numeric(effect$Sex_pval),
    log2fc = NA_real_, stringsAsFactors = FALSE
  )
}
feature_df <- do.call(rbind, feature_rows)
rownames(feature_df) <- NULL

map_rows <- list()
for (i in seq_len(nrow(sex_evidence))) {
  ids <- strsplit(as.character(sex_evidence$feature_ids[i]), ";", fixed = TRUE)[[1]]
  ids <- trimws(ids)
  ids <- unique(ids[nzchar(ids)])
  map_rows[[i]] <- data.frame( platform = sex_evidence$platform[i], Protein = sex_evidence$Protein[i],
    feature_id = ids, stringsAsFactors = FALSE
  )
}
feature_map <- do.call(rbind, map_rows)
feature_df$Protein <- feature_map$Protein[match( paste(feature_df$platform, feature_df$feature_id),
  paste(feature_map$platform, feature_map$feature_id))]
unmapped <- is.na(feature_df$Protein)
feature_df <- feature_df[!unmapped, , drop = FALSE]

for (i in seq_len(nrow(sex_evidence))) {
  rows <- feature_df[feature_df$platform == sex_evidence$platform[i] &
    feature_df$Protein == sex_evidence$Protein[i], , drop = FALSE]
  finite_p <- rows$raw_p[is.finite(rows$raw_p)]
  observed_min <- if (length(finite_p)) {
    min(finite_p)
  } else {
    NA_real_
  }
  observed_n <- sum(is.finite(rows$raw_p) & rows$raw_p < raw_p_cutoff)
}
feature_df$nominal_feature <- is.finite(feature_df$raw_p) & feature_df$raw_p < raw_p_cutoff

for (platform in platform_keys) {
  model_matrix <- model_matrices[[platform]]
  sample_index <- match(as.character(colnames(model_matrix)), as.character(metadata$OriginalID))
  subject <- as.character(metadata$SubjectID[sample_index])
  sex <- as.character(metadata$Sex[sample_index])
  donor_ids <- unique(subject)
  donor_sex <- vapply(donor_ids, function(id) {
    values <- unique(sex[subject == id])
    values
  }, character(1))
  selected <- which(feature_df$platform == platform)
  feature_index <- match(feature_df$feature_id[selected], rownames(model_matrix))
  donor_means <- vapply(donor_ids, function(id) {
    rowMeans(model_matrix[feature_index, subject == id, drop = FALSE])
  }, numeric(length(selected)))
  donor_means <- matrix(donor_means, nrow = length(selected))
  feature_df$log2fc[selected] <-
    rowMeans(donor_means[, donor_sex == "M", drop = FALSE]) -
    rowMeans(donor_means[, donor_sex == "F", drop = FALSE])
}
platform_count_rows <- list()
for (k in seq_along(platform_keys)) {
  platform <- platform_keys[k]
  rows <- sex_evidence[sex_evidence$platform == platform, , drop = FALSE]
  platform_count_rows[[k]] <- data.frame( platform = platform_labels[[platform]], tested_proteins = nrow(rows),
    sex_associated_proteins_p_lt_0_05 = sum( is.finite(rows$min_p_value) & rows$min_p_value < raw_p_cutoff),
    stringsAsFactors = FALSE
  )
}
platform_counts <- do.call(rbind, platform_count_rows)
platform_counts <- rbind(platform_counts, data.frame( platform = "Five-platform unique union",
  tested_proteins = length(unique(sex_evidence$Protein)),
  sex_associated_proteins_p_lt_0_05 = length(unique(sex_evidence$Protein[
  is.finite(sex_evidence$min_p_value) & sex_evidence$min_p_value < raw_p_cutoff])), stringsAsFactors = FALSE
))

pair_index <- utils::combn(platform_keys, 2, simplify = FALSE)
pair_points <- list()
pair_summary_rows <- list()
for (pair in pair_index) {
  key <- paste(pair, collapse = "__")
  x_evidence <- sex_evidence[sex_evidence$platform == pair[1], , drop = FALSE]
  y_evidence <- sex_evidence[sex_evidence$platform == pair[2], , drop = FALSE]
  x_dep <- x_evidence$Protein[is.finite(x_evidence$min_p_value) &
    x_evidence$min_p_value < raw_p_cutoff]
  y_dep <- y_evidence$Protein[is.finite(y_evidence$min_p_value) &
    y_evidence$min_p_value < raw_p_cutoff]
  dep_union <- union(x_dep, y_dep)
  both_tested <- intersect(x_evidence$Protein, y_evidence$Protein)
  eligible <- intersect(dep_union, both_tested)

  x <- feature_df[feature_df$platform == pair[1] & feature_df$Protein %in% eligible,
    c("Protein", "feature_id", "raw_p", "log2fc", "nominal_feature"), drop = FALSE]
  y <- feature_df[feature_df$platform == pair[2] & feature_df$Protein %in% eligible,
    c("Protein", "feature_id", "raw_p", "log2fc", "nominal_feature"), drop = FALSE]
  names(x) <- c("Protein", "feature_x", "raw_p_x", "log2fc_x", "nominal_feature_x")
  names(y) <- c("Protein", "feature_y", "raw_p_y", "log2fc_y", "nominal_feature_y")
  points <- merge(x, y, by = "Protein", all = FALSE, sort = FALSE)

  excluded_non_dep_assay_points <- sum( !points$nominal_feature_x & !points$nominal_feature_y)
  points <- points[points$nominal_feature_x | points$nominal_feature_y, , drop = FALSE]
  points$platform_x <- rep(unname(platform_labels[[pair[1]]]), nrow(points))
  points$platform_y <- rep(unname(platform_labels[[pair[2]]]), nrow(points))
  points$nominal_protein_x <- points$Protein %in% x_dep
  points$nominal_protein_y <- points$Protein %in% y_dep
  points$finite_pair <- is.finite(points$log2fc_x) & is.finite(points$log2fc_y)
  points$direction_consistent <- points$finite_pair &
    sign(points$log2fc_x) == sign(points$log2fc_y) &
    points$log2fc_x != 0 & points$log2fc_y != 0
  if (pair[1] == "Somascan") {
    points$display_id <- paste0(points$Protein, " (", points$feature_x, ")")
  } else if (pair[2] == "Somascan") {
    points$display_id <- paste0(points$Protein, " (", points$feature_y, ")")
  } else {
    points$display_id <- points$Protein
  }
  points <- points[order(points$Protein, points$feature_x, points$feature_y), , drop = FALSE]
  pair_points[[key]] <- points

  finite <- points[points$finite_pair, , drop = FALSE]
  pearson_r <- spearman_rho <- NA_real_
  if (nrow(finite) >= 3L && stats::sd(finite$log2fc_x) > 0 &&
    stats::sd(finite$log2fc_y) > 0) {
    pearson_r <- stats::cor(finite$log2fc_x, finite$log2fc_y, method = "pearson")
    spearman_rho <- stats::cor(finite$log2fc_x, finite$log2fc_y, method = "spearman")
  }
  pair_summary_rows[[key]] <- data.frame( platform_x = unname(platform_labels[[pair[1]]]),
    platform_y = unname(platform_labels[[pair[2]]]), dep_union_all_proteins = length(dep_union),
    dep_union_tested_in_both = length(eligible), dep_in_both = length(intersect(x_dep, y_dep)),
    dep_x_only_within_tested_pair = length(setdiff(intersect(x_dep, both_tested), y_dep)),
    dep_y_only_within_tested_pair = length(setdiff(intersect(y_dep, both_tested), x_dep)),
    excluded_non_dep_assay_points = excluded_non_dep_assay_points, plotted_points_with_finite_log2fc = nrow(finite),
    plotted_distinct_proteins = length(unique(finite$Protein)),
    direction_consistent_points = sum(finite$direction_consistent),
    direction_discordant_points = sum(!finite$direction_consistent), pearson_r = pearson_r, spearman_rho = spearman_rho,
    stringsAsFactors = FALSE
  )
}
pair_summary <- do.call(rbind, pair_summary_rows)
all_pair_points <- do.call(rbind, pair_points)
names(all_pair_points)[names(all_pair_points) == "raw_p_x"] <- "p_x"
names(all_pair_points)[names(all_pair_points) == "raw_p_y"] <- "p_y"
all_pair_points$effect_estimator <- "mean_Donor_log2_expression_M_minus_F"
pair_summary$effect_estimator <- "mean_Donor_log2_expression_M_minus_F"

pair_plot_list <- vector("list", length(pair_index))
pair_plot_limits <- vector("list", length(pair_index))
for (j in seq_along(pair_index)) {
  pair <- pair_index[[j]]
  key <- paste(pair, collapse = "__")
  points <- pair_points[[key]]
  points <- points[points$finite_pair, , drop = FALSE]
  summary_row <- pair_summary[pair_summary$platform_x == platform_labels[[pair[1]]] &
    pair_summary$platform_y == platform_labels[[pair[2]]], , drop = FALSE]
  if (nrow(points)) {
    edge <- max(abs(c(points$log2fc_x, points$log2fc_y)), na.rm = TRUE)
    edge <- max(edge * 1.12, 0.1)
    pair_plot_limits[[j]] <- c(-edge * 1.12, edge * 1.12)
    points$x_direction <- ifelse(points$log2fc_x > 0, "Higher in males",
      ifelse(points$log2fc_x < 0, "Higher in females", "Zero difference"))
    points$significance <- ifelse(points$nominal_feature_x & points$nominal_feature_y,
      "Both", ifelse(points$nominal_feature_x, "Platform 1", "Platform 2"))
    points$significance <- factor(points$significance, levels = c("Both", "Platform 1", "Platform 2"))
    r_label <- paste0("Pearson r = ",
      if (is.finite(summary_row$pearson_r)) {
      sprintf("%.2f", summary_row$pearson_r)
    } else {
      "NA"
    }, "\nSpearman ρ = ",
      if (is.finite(summary_row$spearman_rho)) {
      sprintf("%.2f", summary_row$spearman_rho)
    } else {
      "NA"
    }, "\nn = ", nrow(points), " points")
    p_pair <- ggplot(points, aes(x = log2fc_x, y = log2fc_y)) +
      geom_vline(xintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.35) +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.35) +
      geom_point(aes(colour = x_direction, shape = significance), size = 2.1, alpha = 0.8) +
      annotate("label", x = -edge * 0.98, y = edge * 1.04, label = r_label, hjust = 0, vjust = 1, size = 3,
      label.size = 0.2, fill = "white") +
      scale_colour_manual(values = c("Higher in males" = "#3B76AF", "Higher in females" = "#C45A75",
      "Zero difference" = "#777777"), drop = FALSE) + scale_shape_manual(values = c("Both" = 16, "Platform 1" = 17,
      "Platform 2" = 15), labels = c(paste(platform_labels[[pair[1]]], "and", platform_labels[[pair[2]]]),
      platform_labels[[pair[1]]], platform_labels[[pair[2]]]), drop = FALSE) + guides(colour = guide_legend(order = 1),
      shape = guide_legend(order = 2)) + coord_equal(xlim = c(-edge * 1.12, edge * 1.12),
      ylim = c(-edge * 1.12, edge * 1.12), expand = FALSE) +
      labs(x = paste0(platform_labels[[pair[1]]], " individual-mean log₂FC (male vs female)"),
      y = paste0(platform_labels[[pair[2]]], " individual-mean log₂FC (male vs female)"),
      title = paste(platform_labels[[pair[1]]], "vs", platform_labels[[pair[2]]]),

      colour = paste(platform_labels[[pair[1]]], "direction"), shape = "P < 0.05 in") + theme_classic(base_size = 10) +
      theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 7),
      legend.position = "bottom", legend.box = "vertical")
  } else {
    p_pair <- ggplot() + theme_void(base_size = 10) +
      labs(title = paste(platform_labels[[pair[1]]], "vs", platform_labels[[pair[2]]]),
      subtitle = "No union DEP with finite Log₂FC in both technologies")
  }
  pair_plot_list[[j]] <- p_pair
}
pair_plot_grid <- vector("list", length(pair_plot_list))
for (j in seq_along(pair_plot_list)) {
  p_compact <- pair_plot_list[[j]] + theme( legend.position = "none", plot.caption = element_blank(),
    plot.title = element_text(size = 10), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
    plot.margin = margin(2, 3, 2, 3)
  )
  if (!is.null(pair_plot_limits[[j]])) {
    p_compact <- suppressMessages(p_compact + coord_cartesian(
      xlim = pair_plot_limits[[j]], ylim = pair_plot_limits[[j]], expand = FALSE))
  }
  pair_plot_grid[[j]] <- p_compact
}
sex_results <- list(feature_results = feature_df, platform_counts = platform_counts,
  pair_points = pair_points, pair_summary = pair_summary, points = all_pair_points,
  pair_index = pair_index, plots = pair_plot_grid, has_points = !vapply(pair_plot_limits, is.null, logical(1)))

# 04. Save selected and remaining sex-effect comparisons ----------------------

sex_plots <- save_08_1_comparisons(sex_results$plots, sex_results$pair_index,
  sex_results$has_points, section = "04", effect = "Sex")
# 05. Calculate donor-versus-other-nine effects and pairwise comparisons ------

required_evidence <- c("platform", "effect", "Protein", "n_features_tested",
  "n_nominal_significant", "min_p_value", "feature_ids")
individual_evidence <- evidence[evidence$effect == "Individual" &
  evidence$platform %in% platform_keys, required_evidence, drop = FALSE]
individual_evidence$Protein <- trimws(as.character(individual_evidence$Protein))
feature_rows <- vector("list", length(platform_keys))
for (k in seq_along(platform_keys)) {
  model_result <- individual_model_results[[k]]
  feature_id <- as.character(model_result$Protein)
  feature_rows[[k]] <- data.frame( platform = platform_keys[k], feature_id = feature_id,
    raw_p = as.numeric(model_result$indi_pval), stringsAsFactors = FALSE)
}
feature_df <- do.call(rbind, feature_rows)
rownames(feature_df) <- NULL

map_rows <- vector("list", nrow(individual_evidence))
for (i in seq_len(nrow(individual_evidence))) {
  ids <- strsplit(as.character(individual_evidence$feature_ids[i]), ";", fixed = TRUE)[[1]]
  ids <- unique(trimws(ids[nzchar(ids)]))
  map_rows[[i]] <- data.frame( platform = individual_evidence$platform[i],
    Protein = individual_evidence$Protein[i], feature_id = ids, stringsAsFactors = FALSE)
}
feature_map <- do.call(rbind, map_rows)
feature_df$Protein <- feature_map$Protein[match( paste(feature_df$platform, feature_df$feature_id),
  paste(feature_map$platform, feature_map$feature_id))]
unmapped <- is.na(feature_df$Protein)
feature_df <- feature_df[!unmapped, , drop = FALSE]
for (i in seq_len(nrow(individual_evidence))) {
  rows <- feature_df[feature_df$platform == individual_evidence$platform[i] &
    feature_df$Protein == individual_evidence$Protein[i], , drop = FALSE]
  finite_p <- rows$raw_p[is.finite(rows$raw_p)]
  observed_min <- if (length(finite_p)) {
    min(finite_p)
  } else {
    NA_real_
  }
  observed_n <- sum(is.finite(rows$raw_p) & rows$raw_p < raw_p_cutoff)
}
feature_df$nominal_feature <- is.finite(feature_df$raw_p) &
  feature_df$raw_p < raw_p_cutoff

platform_counts <- do.call(rbind, lapply(platform_keys, function(platform) {
  rows <- individual_evidence[individual_evidence$platform == platform, , drop = FALSE]
  data.frame(platform = unname(platform_labels[[platform]]), tested_proteins = nrow(rows),
    individual_associated_proteins_p_lt_0_05 = sum(rows$n_nominal_significant >= 1L), stringsAsFactors = FALSE)
}))
platform_counts <- rbind(platform_counts, data.frame( platform = "Five-platform unique union",
  tested_proteins = length(unique(individual_evidence$Protein)), individual_associated_proteins_p_lt_0_05 = length(unique(
  individual_evidence$Protein[individual_evidence$n_nominal_significant >= 1L])), stringsAsFactors = FALSE))

pair_index <- utils::combn(platform_keys, 2, simplify = FALSE)
pair_membership <- list()
needed_proteins <- stats::setNames(rep(list(character()), length(platform_keys)), platform_keys)
for (pair in pair_index) {
  key <- paste(pair, collapse = "__")
  x <- individual_evidence[individual_evidence$platform == pair[1], , drop = FALSE]
  y <- individual_evidence[individual_evidence$platform == pair[2], , drop = FALSE]
  x_dep <- x$Protein[x$n_nominal_significant >= 1L]
  y_dep <- y$Protein[y$n_nominal_significant >= 1L]
  dep_union <- union(x_dep, y_dep)
  eligible <- intersect(dep_union, intersect(x$Protein, y$Protein))
  pair_membership[[key]] <- list(dep_union = dep_union, eligible = eligible, x_dep = x_dep, y_dep = y_dep)
  for (platform in pair) {
    needed_proteins[[platform]] <-
      union(needed_proteins[[platform]], eligible)
  }
}
# 05.3. Calculate donor-versus-other-nine log2 fold changes -------------------

donor_fc_rows <- list()
for (platform in platform_keys) {
  model_matrix <- model_matrices[[platform]]
  sample_index <- match(as.character(colnames(model_matrix)), as.character(metadata$OriginalID))
  subject <- as.character(metadata$SubjectID[sample_index])
  selected <- feature_df[feature_df$platform == platform &
    feature_df$Protein %in% needed_proteins[[platform]], , drop = FALSE]
  for (i in seq_len(nrow(selected))) {
    expression <- as.numeric(model_matrix[selected$feature_id[i], ])
    donor_means <- vapply(donor_display_map$subject_key, function(id)
      mean(expression[subject == id]), numeric(1))
    donor_log2fc <- vapply(seq_along(donor_means), function(j)
      donor_means[j] - mean(donor_means[-j]), numeric(1))
    donor_fc_rows[[length(donor_fc_rows) + 1L]] <- data.frame( platform = platform, Protein = selected$Protein[i],
      feature_id = selected$feature_id[i], DonorLabel = donor_display_map$DonorLabel,
      donor_mean_log2_expression = as.numeric(donor_means), log2fc_vs_other_9 = as.numeric(donor_log2fc),
      raw_p = selected$raw_p[i], nominal_feature = selected$nominal_feature[i], stringsAsFactors = FALSE)
  }
}
donor_fc_df <- do.call(rbind, donor_fc_rows)
donor_fc_df$effect_estimator <- "mean_Donor_log2_minus_mean_other_9_Donors"
# 05.4. Summarise donor-stratified comparisons --------------------------------

pair_points <- list()
summary_rows <- list()
pooled_summary_rows <- list()
for (pair in pair_index) {
  key <- paste(pair, collapse = "__")
  member <- pair_membership[[key]]
  x <- donor_fc_df[donor_fc_df$platform == pair[1] &
    donor_fc_df$Protein %in% member$eligible, c("Protein", "feature_id", "DonorLabel", "donor_mean_log2_expression",
    "log2fc_vs_other_9", "raw_p", "nominal_feature"), drop = FALSE]
  y <- donor_fc_df[donor_fc_df$platform == pair[2] &
    donor_fc_df$Protein %in% member$eligible, c("Protein", "feature_id", "DonorLabel", "donor_mean_log2_expression",
    "log2fc_vs_other_9", "raw_p", "nominal_feature"), drop = FALSE]
  names(x) <- c("Protein", "feature_x", "DonorLabel", "donor_mean_log2_x", "log2fc_x", "p_x", "nominal_feature_x")
  names(y) <- c("Protein", "feature_y", "DonorLabel", "donor_mean_log2_y", "log2fc_y", "p_y", "nominal_feature_y")
  points <- merge(x, y, by = c("Protein", "DonorLabel"), all = FALSE, sort = FALSE)

  excluded_by_donor <- table(factor(points$DonorLabel[
    !points$nominal_feature_x & !points$nominal_feature_y], levels = donor_display_map$DonorLabel))
  points <- points[points$nominal_feature_x | points$nominal_feature_y, , drop = FALSE]
  points$platform_x <- rep(unname(platform_labels[[pair[1]]]), nrow(points))
  points$platform_y <- rep(unname(platform_labels[[pair[2]]]), nrow(points))
  points$nominal_protein_x <- points$Protein %in% member$x_dep
  points$nominal_protein_y <- points$Protein %in% member$y_dep
  points$finite_pair <- is.finite(points$log2fc_x) &
    is.finite(points$log2fc_y)
  points$significance <- ifelse(points$nominal_feature_x & points$nominal_feature_y,
    "Both", ifelse(points$nominal_feature_x, "Platform 1", "Platform 2"))
  points$significance <- factor(points$significance, levels = c("Both", "Platform 1", "Platform 2"))
  points <- points[order(points$DonorLabel, points$Protein, points$feature_x, points$feature_y), , drop = FALSE]
  pair_points[[key]] <- points
  pooled_points <- points[points$finite_pair, , drop = FALSE]
  pooled_pearson_r <- pooled_spearman_rho <- NA_real_
  if (nrow(pooled_points) >= 3L &&
    stats::sd(pooled_points$log2fc_x) > 0 &&
    stats::sd(pooled_points$log2fc_y) > 0) {
    pooled_pearson_r <- stats::cor(pooled_points$log2fc_x, pooled_points$log2fc_y, method = "pearson")
    pooled_spearman_rho <- stats::cor(pooled_points$log2fc_x, pooled_points$log2fc_y, method = "spearman")
  }
  pooled_summary_rows[[key]] <- data.frame( platform_x = unname(platform_labels[[pair[1]]]),
    platform_y = unname(platform_labels[[pair[2]]]), dep_union_all_proteins = length(member$dep_union),
    dep_union_tested_in_both = length(member$eligible), excluded_non_dep_assay_points = sum(excluded_by_donor),
    plotted_points_with_finite_log2fc = nrow(pooled_points),
    plotted_distinct_proteins = length(unique(pooled_points$Protein)),
    plotted_donors = length(unique(pooled_points$DonorLabel)),
    pearson_r = pooled_pearson_r, spearman_rho = pooled_spearman_rho, stringsAsFactors = FALSE)
  for (donor in donor_display_map$DonorLabel) {
    donor_points <- points[points$DonorLabel == donor & points$finite_pair, , drop = FALSE]
    pearson_r <- spearman_rho <- NA_real_
    if (nrow(donor_points) >= 3L &&
      stats::sd(donor_points$log2fc_x) > 0 &&
      stats::sd(donor_points$log2fc_y) > 0) {
      pearson_r <- stats::cor(donor_points$log2fc_x, donor_points$log2fc_y, method = "pearson")
      spearman_rho <- stats::cor(donor_points$log2fc_x, donor_points$log2fc_y, method = "spearman")
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame( platform_x = unname(platform_labels[[pair[1]]]),
      platform_y = unname(platform_labels[[pair[2]]]), DonorLabel = donor,
      dep_union_all_proteins = length(member$dep_union), dep_union_tested_in_both = length(member$eligible),
      dep_in_both = length(intersect(member$x_dep, member$y_dep)),
      excluded_non_dep_assay_points = as.integer(excluded_by_donor[donor]),
      plotted_points_with_finite_log2fc = nrow(donor_points),
      plotted_distinct_proteins = length(unique(donor_points$Protein)), pearson_r = pearson_r, spearman_rho = spearman_rho,
      stringsAsFactors = FALSE)
  }
}
pair_summary <- do.call(rbind, summary_rows)
pair_pooled_summary <- do.call(rbind, pooled_summary_rows)
all_pair_points <- do.call(rbind, pair_points)
all_pair_points$effect_estimator <- "mean_Donor_log2_minus_mean_other_9_Donors"
pair_summary$effect_estimator <- "mean_Donor_log2_minus_mean_other_9_Donors"
pair_pooled_summary$effect_estimator <- "mean_Donor_log2_minus_mean_other_9_Donors"

pair_overview_plots <- vector("list", length(pair_index))
pair_has_points <- logical(length(pair_index))
for (j in seq_along(pair_index)) {
  pair <- pair_index[[j]]
  key <- paste(pair, collapse = "__")
  points <- pair_points[[key]]
  points <- points[points$finite_pair, , drop = FALSE]
  points$DonorLabel <- factor(points$DonorLabel, levels = donor_display_map$DonorLabel)
  if (nrow(points)) {
    pair_has_points[j] <- TRUE
    pooled_row <- pair_pooled_summary[
      pair_pooled_summary$platform_x == platform_labels[[pair[1]]] &
      pair_pooled_summary$platform_y == platform_labels[[pair[2]]], , drop = FALSE]
    pooled_label <- paste0(
      "Pearson r = ", if (is.finite(pooled_row$pearson_r)) {
      sprintf("%.2f", pooled_row$pearson_r)
    } else {
      "NA"
    },
      "\nSpearman ρ = ", if (is.finite(pooled_row$spearman_rho)) {
      sprintf("%.2f", pooled_row$spearman_rho)
    } else {
      "NA"
    }, "\nn = ", pooled_row$plotted_points_with_finite_log2fc, " points")
    pair_overview_plots[[j]] <- ggplot(points, aes(x = log2fc_x, y = log2fc_y)) +
      geom_vline(xintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.25) +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.25) +
      geom_point(aes(fill = DonorLabel, shape = significance), colour = "grey35", stroke = 0.2, size = 1.2, alpha = 0.8) +
      annotate("label", x = -Inf, y = Inf, label = pooled_label, hjust = -0.05, vjust = 1.1, size = 2.25,
      label.size = 0.15, fill = "white") + scale_fill_manual(values = donor_palette, drop = FALSE,
      labels = function(x) sub("^Donor ", "Individual ", x), name = "Individual") +
      scale_shape_manual(values = c("Both" = 21, "Platform 1" = 24, "Platform 2" = 22),
      labels = c("Both", "Horizontal platform", "Vertical platform"), drop = FALSE, name = "Individual P < 0.05 in") +
      labs(x = paste0(platform_labels[[pair[1]]], " log₂FC vs other 9"),
      y = paste0(platform_labels[[pair[2]]], " log₂FC vs other 9"), title = paste(platform_labels[[pair[1]]], "vs",
      platform_labels[[pair[2]]])) +
      theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9),
      axis.title = element_text(size = 7), axis.text = element_text(size = 6), legend.position = "bottom")
  } else {
    pair_overview_plots[[j]] <- ggplot() + theme_void(base_size = 8) +
      labs(title = paste(platform_labels[[pair[1]]], "vs", platform_labels[[pair[2]]]), subtitle = "No paired points")
  }
}
individual_results <- list(feature_results = feature_df, donor_fc = donor_fc_df, platform_counts = platform_counts, pair_points = pair_points,
  pair_summary = pair_summary, pooled_summary = pair_pooled_summary, points = all_pair_points,
  pair_index = pair_index, plots = pair_overview_plots, has_points = pair_has_points)

# 06. Save selected and remaining individual-effect comparisons ---------------

individual_plots <- save_08_1_comparisons(individual_results$plots,
  individual_results$pair_index, individual_results$has_points, section = "06", effect = "Individual")
# 07. Compare BH-FDR associated-protein sets ----------------------------------

somer_protein <- cache$somer_protein
dep_platform_keys <- platform_keys
col_platform_internal <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
  "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF", "Overlap" = "#D1C1E1FF"
)
platform_display_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral",
  "Olink-Explore" = "Olink-Explore", "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan", "Overlap" = "Overlap"
)

soma_map_df <- somer_protein %>% dplyr::transmute( somamer_id = as.character(.data$SeqId),
  protein = as.character(.data$EntrezGeneSymbol)
) %>%
  dplyr::filter( !is.na(.data$somamer_id), .data$somamer_id != "", !is.na(.data$protein), .data$protein != ""
) %>% dplyr::distinct(.data$somamer_id, .keep_all = TRUE)

dep_long_records <- list()
for (platform_key in dep_platform_keys) {
  platform_dep_df <- individual_model_results[[platform_key]]
  required_model_columns <- c("Time_pval", "Sex_pval", "indi_pval", "Protein")
  platform_dep_df <- platform_dep_df %>% dplyr::mutate( Time_pval = stats::p.adjust(.data$Time_pval, method = "BH"),
    Sex_pval = stats::p.adjust(.data$Sex_pval, method = "BH"), indi_pval = stats::p.adjust(.data$indi_pval, method = "BH"),
    platform_key = platform_key
  )
  if (platform_key == "Somascan") {
    platform_dep_df <- platform_dep_df %>% dplyr::transmute( platform_key = .data$platform_key,
      Time_pval = .data$Time_pval, Sex_pval = .data$Sex_pval, indi_pval = .data$indi_pval,
      somamer_id = as.character(.data$Protein)
    ) %>% dplyr::left_join(soma_map_df, by = "somamer_id")
  } else {
    platform_dep_df <- platform_dep_df %>%
      dplyr::transmute( platform_key = .data$platform_key, Time_pval = .data$Time_pval, Sex_pval = .data$Sex_pval,
      indi_pval = .data$indi_pval, protein = as.character(.data$Protein)
    )
  }
  dep_long_records[[platform_key]] <- platform_dep_df %>% tidyr::pivot_longer(
    cols = c("Time_pval", "Sex_pval", "indi_pval"), names_to = "effect_key", values_to = "fdr_value"
  ) %>% dplyr::filter(is.finite(.data$fdr_value), .data$fdr_value < 0.05) %>% dplyr::mutate( effect = dplyr::recode(
    .data$effect_key, "indi_pval" = "Individual", "Sex_pval" = "Sex", "Time_pval" = "Time"
  )
  )
}
dep_display_df <- dplyr::bind_rows(dep_long_records) %>%
  dplyr::filter(.data$effect %in% c("Individual", "Sex"), !is.na(.data$protein))

p_association_venn_list <- list()
for (effect_name in c("Individual", "Sex")) {
  effect_df <- dep_display_df %>% dplyr::filter(.data$effect == effect_name)
  effect_sets <- lapply(dep_platform_keys, function(platform_key) {
    unique(effect_df$protein[effect_df$platform_key == platform_key])
  })
  names(effect_sets) <- unname(platform_display_labels[dep_platform_keys])
  effect_colours <- unname(col_platform_internal[dep_platform_keys])
  p_association_venn_list[[effect_name]] <- VennDiagram::venn.diagram( x = effect_sets, filename = NULL, col = "white",
    fill = effect_colours, alpha = 0.6, cex = 0.45, cat.cex = 0.45, cat.col = "black", lwd = 0.5, cat.lwd = 0.5,
    margin = 0.1, cat.dist = rep(0.2, length(effect_sets)), main = paste0(effect_name, " (FDR < 0.05)"), main.cex = 0.5
  )
}
p_association_venn <- do.call( gridExtra::arrangeGrob, c(p_association_venn_list, ncol = 2)
)

venn_pdf_file <- file.path(wdplot, "S9b_individual_sex_dep_venn_platform_colours.pdf")

ggplot2::ggsave(venn_pdf_file, pdf_text_grob(p_association_venn), width = 12, height = 5.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

venn_sets <- lapply(c("Individual", "Sex"), function(effect_name) {
  effect_df <- dep_display_df[dep_display_df$effect == effect_name, , drop = FALSE]
  stats::setNames(lapply(dep_platform_keys, function(key)
    unique(effect_df$protein[effect_df$platform_key == key])), unname(platform_display_labels[dep_platform_keys]))
})
names(venn_sets) <- c("Individual", "Sex")

venn_results <- list(sets = venn_sets, feature_evidence = dep_display_df, plot = p_association_venn)

# 08. Compare nominal-P associated-protein sets -------------------------------

somer_protein <- cache$somer_protein
dep_platform_keys <- platform_keys
col_platform_internal <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
  "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF", "Overlap" = "#D1C1E1FF"
)
platform_display_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral",
  "Olink-Explore" = "Olink-Explore", "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan", "Overlap" = "Overlap"
)

soma_map_df <- somer_protein %>%
  dplyr::transmute( somamer_id = as.character(.data$SeqId), protein = as.character(.data$EntrezGeneSymbol)
) %>% dplyr::filter( !is.na(.data$somamer_id), .data$somamer_id != "", !is.na(.data$protein), .data$protein != ""
) %>%
  dplyr::distinct(.data$somamer_id, .keep_all = TRUE)

dep_long_records <- list()
for (platform_key in dep_platform_keys) {
  platform_dep_df <- individual_model_results[[platform_key]]
  required_model_columns <- c("Sex_pval", "indi_pval", "Protein")
  platform_dep_df <- platform_dep_df %>% dplyr::mutate( platform_key = platform_key
  )
  if (platform_key == "Somascan") {
    platform_dep_df <- platform_dep_df %>% dplyr::transmute( platform_key = .data$platform_key, Sex_pval = .data$Sex_pval,
      indi_pval = .data$indi_pval, somamer_id = as.character(.data$Protein)
    ) %>%
      dplyr::left_join(soma_map_df, by = "somamer_id")
  } else {
    platform_dep_df <- platform_dep_df %>%
      dplyr::transmute( platform_key = .data$platform_key, Sex_pval = .data$Sex_pval, indi_pval = .data$indi_pval,
      protein = as.character(.data$Protein)
    )
  }
  dep_long_records[[platform_key]] <- platform_dep_df %>% tidyr::pivot_longer( cols = c("Sex_pval", "indi_pval"),
    names_to = "effect_key", values_to = "p_value"
  ) %>% dplyr::filter(is.finite(.data$p_value), .data$p_value < 0.05) %>% dplyr::mutate( effect = dplyr::recode(
    .data$effect_key, "indi_pval" = "Individual", "Sex_pval" = "Sex"
  )
  )
}
dep_display_df <- dplyr::bind_rows(dep_long_records) %>%
  dplyr::filter(.data$effect %in% c("Individual", "Sex"), !is.na(.data$protein))

p_association_venn_list <- list()
for (effect_name in c("Individual", "Sex")) {
  effect_df <- dep_display_df %>% dplyr::filter(.data$effect == effect_name)
  effect_sets <- lapply(dep_platform_keys, function(platform_key) {
    unique(effect_df$protein[effect_df$platform_key == platform_key])
  })
  names(effect_sets) <- unname(platform_display_labels[dep_platform_keys])
  effect_colours <- unname(col_platform_internal[dep_platform_keys])
  p_association_venn_list[[effect_name]] <- VennDiagram::venn.diagram( x = effect_sets, filename = NULL, col = "white",
    fill = effect_colours, alpha = 0.6, cex = 0.45, cat.cex = 0.45, cat.col = "black", lwd = 0.5, cat.lwd = 0.5,
    margin = 0.1, cat.dist = rep(0.2, length(effect_sets)), main = paste0(effect_name, " (P < 0.05)"), main.cex = 0.5
  )
}
p_association_venn <- do.call( gridExtra::arrangeGrob, c(p_association_venn_list, ncol = 2)
)

venn_pdf_file <- file.path(wdplot, "S9b_nominal_p_individual_sex_dep_venn_platform_colours.pdf")

ggplot2::ggsave(venn_pdf_file, pdf_text_grob(p_association_venn), width = 12, height = 5.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

venn_sets <- lapply(c("Individual", "Sex"), function(effect_name) {
  effect_df <- dep_display_df[dep_display_df$effect == effect_name, , drop = FALSE]
  stats::setNames(lapply(dep_platform_keys, function(key)
    unique(effect_df$protein[effect_df$platform_key == key])), unname(platform_display_labels[dep_platform_keys]))
})
names(venn_sets) <- c("Individual", "Sex")

nominal_venn_results <- list(sets = venn_sets, feature_evidence = dep_display_df, plot = p_association_venn)

saveRDS(list(counts = counts_results, sex = sex_results, individual = individual_results,
  venn_fdr = venn_results, venn_nominal = nominal_venn_results),
  paste0(wdresult, "association_comparisons.rds"))
# 08.2. Olink associations across detection sets ------------------------------
# 08.2 ------------------------------------------------------------------------

fit_olink_models <- function(expression_matrix, platform) {
  metadata_model <- metadata[match(colnames(expression_matrix), metadata$OriginalID),
    c("OriginalID", "SubjectID", "Time", "Sex"), drop = FALSE]
  metadata_model$SubjectID <- factor(metadata_model$SubjectID)
  metadata_model$Time <- factor(metadata_model$Time, levels = c("Day0", "Day3", "Day5"))
  metadata_model$Sex <- factor(metadata_model$Sex)

  protein_ids <- rownames(expression_matrix)
  result_rows <- vector("list", length(protein_ids))

  for (i in seq_along(protein_ids)) {
    model_df <- data.frame( Expression = as.numeric(expression_matrix[i, ]), SubjectID = metadata_model$SubjectID,
      Time = metadata_model$Time, Sex = metadata_model$Sex
    )
    model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]

    fit_result <- tryCatch( {
      model <- suppressMessages(suppressWarnings(lmerTest::lmer( Expression ~ Time + Sex + (1 | SubjectID), data = model_df,
        REML = FALSE
      )))
      fixed_anova <- stats::anova(model, type = 3)
      random_anova <- lmerTest::ranova(model)
      data.frame( Time_pval = fixed_anova["Time", "Pr(>F)"], Sex_pval = fixed_anova["Sex", "Pr(>F)"],
        Individual_pval = random_anova["(1 | SubjectID)", "Pr(>Chisq)"], n_observations = nrow(model_df), model_ok = TRUE,
        singular_fit = lme4::isSingular(model, tol = 1e-4), stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame( Time_pval = NA_real_, Sex_pval = NA_real_, Individual_pval = NA_real_, n_observations = nrow(model_df),
        model_ok = FALSE, singular_fit = NA, stringsAsFactors = FALSE
      )
    })

    fit_result$Protein <- protein_ids[[i]]
    fit_result$Platform <- platform
    result_rows[[i]] <- fit_result

    if (i %% 250 == 0 || i == length(protein_ids)) {
      message(platform, ": fitted ", i, "/", length(protein_ids), " protein models")
    }
  }
  dplyr::bind_rows(result_rows) %>%
    dplyr::select(Platform, Protein, Time_pval, Sex_pval, Individual_pval, n_observations, model_ok, singular_fit)
}

load_registered_assays(wddata)
metadata <- read_processed_metadata(wddata)[1:30, ]
detection_sets <- readRDS(file.path(repo_dir, "result", "Fig2",
  "01-01_detection_classification_sets.rds"
))

rownames(olink3072)[rownames(olink3072) == "NTproBNP"] <- "NT-proBNP"
rownames(olink3072)[rownames(olink3072) == "WARS"] <- "WARS1"
lod3072$Assay[lod3072$Assay == "NTproBNP"] <- "NT-proBNP"
lod3072$Assay[lod3072$Assay == "WARS"] <- "WARS1"

lod3072$Assay[lod3072$Assay == "CERT1"] <- "CERT"

metadata$Time <- dplyr::recode( metadata$Time, "visit1" = "Day0", "visit2" = "Day3", "visit3" = "Day5"
)
metadata$Time <- factor(metadata$Time, levels = c("Day0", "Day3", "Day5"))
metadata$SubjectID <- factor(metadata$SubjectID, levels = unique(metadata$SubjectID))

predefined_sets <- list( "Olink-Explore" = intersect(unique(rownames(olink3072)), unique(lod3072$Assay)),
  "Olink-HT" = intersect(unique(rownames(olinkht)), unique(lodht$Assay))
)
protein_sets <- list( "Predefined protein targets" = predefined_sets,
  "Detected proteins" = detection_sets[["detected"]][c("Olink-Explore", "Olink-HT")],
  "Individual detected proteins" = detection_sets[["individual_detected"]][c("Olink-Explore", "Olink-HT")],
  "Broadly detected proteins" = detection_sets[["broad_detected"]][c("Olink-Explore", "Olink-HT")]
)

dflistna <- list("Olink-Explore" = olink3072, "Olink-HT" = olinkht)
expression_matrices <- list( "Olink-Explore" = olinkmatrixpro(
  pro = predefined_sets[["Olink-Explore"]], cate = "Olink-Explore", data = olink3072, isna = "impute", keepfilter = FALSE
), "Olink-HT" = olinkmatrixpro( pro = predefined_sets[["Olink-HT"]], cate = "Olink-HT",
  data = olinkht, isna = "impute", keepfilter = FALSE
)
)

model_results <- dplyr::bind_rows( fit_olink_models(expression_matrices[["Olink-Explore"]], "Olink-Explore"),
  fit_olink_models(expression_matrices[["Olink-HT"]], "Olink-HT")
)

# 02. Associated-protein counts and platform-specific plots -------------------

protein_set_levels <- c( "Predefined protein targets", "Detected proteins", "Individual detected proteins",
  "Broadly detected proteins"
)
association_levels <- c("Individual-associated proteins", "Sex-associated proteins")

summary_rows <- list()
evidence_rows <- list()
row_index <- 1L
for (platform in c("Olink-Explore", "Olink-HT")) {
  platform_results <- dplyr::filter(model_results, Platform == platform)
  for (protein_set in protein_set_levels) {
    members <- unique(protein_sets[[protein_set]][[platform]])
    scope_results <- platform_results %>% dplyr::filter(Protein %in% members)

    for (association in association_levels) {
      p_column <- if (association == "Individual-associated proteins") {
        "Individual_pval"
      } else {
        "Sex_pval"
      }
      evidence <- scope_results %>% dplyr::transmute( Platform = platform, Protein_set = protein_set,
        Association = association, Protein = Protein, P_value = .data[[p_column]],
        FDR = stats::p.adjust(.data[[p_column]], method = "BH"), nominal_significant = !is.na(P_value) & P_value < 0.05,
        FDR_significant = !is.na(FDR) & FDR < 0.05, model_ok = model_ok, singular_fit = singular_fit
      )
      evidence_rows[[row_index]] <- evidence
      summary_rows[[row_index]] <- data.frame( Platform = platform, Protein_set = protein_set, Association = association,
        protein_set_size = length(members), proteins_tested = sum(!is.na(evidence$P_value)),
        nominal_associated_proteins = sum(evidence$nominal_significant),
        fdr_associated_proteins = sum(evidence$FDR_significant),

        associated_proteins = sum(evidence$FDR_significant), stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}
association_evidence <- dplyr::bind_rows(evidence_rows)
association_summary <- dplyr::bind_rows(summary_rows) %>% dplyr::mutate(
  Protein_set = factor(Protein_set, levels = protein_set_levels),
  Association = factor(Association, levels = association_levels),
  x_category = interaction(Protein_set, Association, sep = " — ", lex.order = TRUE)
) %>% dplyr::arrange(Platform, Protein_set, Association)

saveRDS( association_evidence, file.path(wdresult, "08.2-02_olink_association_by_protein_set_evidence.rds")
)

platform_colors <- c("Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF")
threshold_levels <- c("nominal P < 0.05", "BH-FDR < 0.05")

association_plot_data <- dplyr::bind_rows( association_summary %>% dplyr::transmute( Platform, Protein_set, Association,
  threshold = "nominal P < 0.05", associated_proteins = nominal_associated_proteins),
  association_summary %>% dplyr::transmute( Platform, Protein_set, Association,
  threshold = "BH-FDR < 0.05", associated_proteins = fdr_associated_proteins)
) %>%
  dplyr::mutate(threshold = factor(threshold, levels = threshold_levels))
shared_y_max <- max(association_plot_data$associated_proteins, na.rm = TRUE)

plot_one_platform <- function(platform) {
  plot_df <- association_plot_data %>% dplyr::filter(Platform == platform) %>% dplyr::mutate( x_category = factor(
    paste(as.character(Protein_set), as.character(Association), sep = "\n"), levels = paste(
    rep(protein_set_levels, each = length(association_levels)), rep(association_levels, times = length(protein_set_levels)),
    sep = "\n"
  )
  )
  )
  label_offset <- max(1, shared_y_max * 0.015)
  bar_position <- position_dodge(width = 0.80)

  ggplot(plot_df, aes(x = x_category, y = associated_proteins, group = threshold)) +
    geom_col( aes(alpha = threshold), position = bar_position, fill = unname(platform_colors[[platform]]),
    color = "black", linewidth = 0.35, width = 0.72
  ) +
    geom_text( aes(y = associated_proteins + label_offset, label = associated_proteins),
    position = bar_position, size = 2.6, vjust = 0
  ) +
    facet_wrap(~Platform) + scale_x_discrete(labels = function(x)
    gsub("Broadly detected proteins", "Broad detected proteins", x, fixed = TRUE)) + scale_alpha_manual(
    values = c("nominal P < 0.05" = 0.45, "BH-FDR < 0.05" = 1),
    breaks = threshold_levels, labels = c("P < 0.05", "FDR < 0.05"), drop = FALSE, name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(fill = unname(platform_colors[[platform]])))
  ) +
    scale_y_continuous( limits = c(0, shared_y_max * 1.10 + label_offset), expand = expansion(mult = c(0, 0))
  ) + labs(x = NULL, y = "Number of associated proteins") + theme2_nobackplus + theme( legend.position = "bottom",
    legend.title = element_text(size = 9), legend.text = element_text(size = 9), strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10), plot.margin = margin(t = 8, r = 12, b = 10, l = 24, unit = "pt")
  )
}
p_olink_explore <- plot_one_platform("Olink-Explore")
p_olink_ht <- plot_one_platform("Olink-HT")

ggsave( file.path(wdplot, "S8b_Olink-Explore_associated_counts_by_protein_set.pdf"),
  pdf_text_grob(p_olink_explore), width = 15, height = 12, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)
ggsave( file.path(wdplot, "S8c_Olink-HT_associated_counts_by_protein_set.pdf"),
  pdf_text_grob(p_olink_ht), width = 15, height = 12, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
)

message("08.2 completed: two platform-specific nominal-P / BH-FDR associated-protein count panels were generated.")

saveRDS(list(models = model_results, protein_sets = protein_sets, summary = association_summary),
  paste0(wdresult, "olink_association_by_detection_set.rds"))
# 08.10. Selected individual and sex-associated protein profiles --------------

platform_keys <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
platform_labels <- setNames(c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "SomaScan"), platform_keys)
y_labels <- setNames(c("log2 intensity", "log2 intensity", "NPX", "NPX", "log2 RFU"), platform_keys)
col_donor <- c("#D4C96E", "#807961", "#8DACB6", "#7D839C", "#CA9791", "#55794D", "#EFF3F9", "#924741", "#E0DC96", "#DA9F58")

col_sex <- c("Female" = "#CC79A7", "Male" = "#56B4E9")
page_ncol <- 4L
page_width_cm <- 32

trajectory_theme <- theme_classic(base_size = 8) + theme(axis.line = element_line(colour = "black", linewidth = 0.3),
  axis.ticks = element_line(linewidth = 0.3), axis.title = element_text(size = 8), axis.text = element_text(size = 8),
  plot.title = element_text(size = 9, face = "plain"),
  legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.key.size = grid::unit(0.4, "cm"),
  plot.margin = margin(5, 7, 5, 5))
cache <- model_bundle
expression_matrices <- cache$c8_imputed_matrix_list
metadata <- cache$metadata
somer_protein <- cache$somer_protein
individual_protein_sets <- cache$c8_dep_individual_list
sex_protein_sets <- cache$c8_dep_sex_list
dep_evidence <- cache$c8_dep_protein_evidence
required_evidence_columns <- c("platform", "effect", "Protein", "feature_ids",
  "significant_feature_ids", "n_fdr_significant")
soma_dep_evidence <- dep_evidence[dep_evidence$platform == "Somascan" &
  dep_evidence$effect %in% c("Individual", "Sex"), , drop = FALSE]
soma_dep_evidence$display_feature_ids <- soma_dep_evidence$significant_feature_ids
soma_dep_evidence$n_display_significant <- soma_dep_evidence$n_fdr_significant
soma_dep_evidence$display_criterion <- "BH-FDR < 0.05"
soma_map <- unique(somer_protein[, c("SeqId", "EntrezGeneSymbol")])
soma_map <- soma_map[!is.na(soma_map$SeqId) & !is.na(soma_map$EntrezGeneSymbol) &
  nzchar(soma_map$SeqId) & nzchar(soma_map$EntrezGeneSymbol), , drop = FALSE]
# 03. Align Donor and visit labels --------------------------------------------

subject_keys <- unique(as.character(metadata$SubjectID))
subject_numeric <- suppressWarnings(as.numeric(sub("^Subject", "", subject_keys)))
donor_display_map <- data.frame(subject_key = subject_keys, subject_numeric = subject_numeric)
donor_display_map <- donor_display_map[order(donor_display_map$subject_numeric), ]
donor_display_map$DonorLabel <- sprintf("Donor %02d", seq_len(10))
col_donor <- setNames(col_donor, donor_display_map$DonorLabel)

plot_metadata <- data.frame( sample_key = as.character(metadata$OriginalID),
  donor = donor_display_map$DonorLabel[match(as.character(metadata$SubjectID), donor_display_map$subject_key)],
  day = unname(c(Day0 = 0, Day3 = 3, Day5 = 5)[as.character(metadata$Time)]),
  sex = unname(c(F = "Female", M = "Male")[as.character(metadata$Sex)]))
plot_metadata$donor <- factor(plot_metadata$donor, levels = names(col_donor))
plot_metadata$sex <- factor(plot_metadata$sex, levels = names(col_sex))
plot_metadata <- plot_metadata[order(plot_metadata$donor, plot_metadata$day), ]
for (platform in platform_keys) {
  mat <- expression_matrices[[platform]]
}
# 04. Select proteins for the current displays --------------------------------

individual_candidates <- sort(unique(Reduce(intersect, individual_protein_sets[platform_keys])))
sex_nominal_sets <- setNames(vector("list", length(platform_keys)), platform_keys)
for (platform in platform_keys) {
  rows <- dep_evidence[dep_evidence$platform == platform & dep_evidence$effect == "Sex", ]
  sex_nominal_sets[[platform]] <- sort(unique(as.character(rows$Protein[
    is.finite(rows$min_p_value) & rows$min_p_value < 0.05])))
}
sex_candidates <- sort(Reduce(intersect, sex_nominal_sets))

individual_examples <- individual_candidates
sex_examples <- sex_candidates

# 04.1. Select nominal-P sex-associated SOMAmers ------------------------------

soma_sex_model_file <- file.path(in08, "C8_2_Somascan_mixed_effects_dep_results.rds")
soma_sex_model <- readRDS(soma_sex_model_file)
set_sex_nominal_display_ids <- function(proteins, evidence = soma_dep_evidence, model = soma_sex_model) {
  for (protein in proteins) {
    idx <- which(evidence$Protein == protein & evidence$effect == "Sex")
    if (!length(idx)) {
      next
    }
    ids <- strsplit(as.character(evidence$feature_ids[idx]), ";", fixed = TRUE)[[1]]
    model_idx <- match(ids, as.character(model$Protein))
    p_values <- as.numeric(model$Sex_pval[model_idx])
    finite_p <- p_values[is.finite(p_values)]
    observed_min <- if (length(finite_p)) {
      min(finite_p)
    } else {
      NA_real_
    }
    nominal_ids <- sort(unique(ids[is.finite(p_values) & p_values < 0.05]))
    evidence$display_feature_ids[idx] <- paste(nominal_ids, collapse = ";")
    evidence$n_display_significant[idx] <- length(nominal_ids)
    evidence$display_criterion[idx] <- "Raw P < 0.05"
  }
  evidence
}
soma_dep_evidence <- set_sex_nominal_display_ids(sex_examples)
# 05. Prepare effect-specific protein displays --------------------------------

get_significant_somamers <- function(protein, effect = c("Individual", "Sex"),
  matrices = expression_matrices, mapping = soma_map,
  evidence = soma_dep_evidence) {
  effect <- match.arg(effect)
  row <- evidence[evidence$Protein == protein & evidence$effect == effect, , drop = FALSE]
  if (!nrow(row)) {
    return(character())
  }
  ids <- strsplit(as.character(row$display_feature_ids), ";", fixed = TRUE)[[1]]
  ids <- sort(unique(trimws(ids[nzchar(trimws(ids))])))
  tested_ids <- strsplit(as.character(row$feature_ids), ";", fixed = TRUE)[[1]]
  mapped_ids <- unique(as.character(mapping$SeqId[mapping$EntrezGeneSymbol == protein]))
  ids
}
# 05.1. Standardise expression across the 30 matched samples ------------------

prepare_protein_donor_zscore <- function(protein, matrices = expression_matrices,
  sample_metadata = plot_metadata, mapping = soma_map,
  effect = c("Individual", "Sex"), evidence = soma_dep_evidence) {
  effect <- match.arg(effect)
  somamer_ids <- get_significant_somamers(protein, effect, matrices, mapping, evidence)
  platforms <- c(platform_keys[1:4], rep("Somascan", max(1L, length(somamer_ids))))
  features <- c(rep(protein, 4), if (length(somamer_ids)) {
    somamer_ids
  } else {
    NA_character_
  })
  technology_labels <- c(unname(platform_labels[platform_keys[1:4]]),
    if (length(somamer_ids)) {
    paste0("SomaScan\n", somamer_ids)
  } else {
    "SomaScan"
  })
  donor_order <- levels(sample_metadata$donor)
  donor_rows <- match(donor_order, as.character(sample_metadata$donor))
  donor_df_list <- list()
  visit_df_list <- list()
  unavailable <- character()

  for (i in seq_along(platforms)) {
    platform <- platforms[i]
    feature <- features[i]
    technology <- technology_labels[i]
    if (is.na(feature) || !feature %in% rownames(matrices[[platform]])) {
      unavailable[technology] <- if (platform == "Somascan" && is.na(feature)) {
        "No significant\nSOMAmer"
      } else {
        "Not available"
      }
      next
    }
    values <- as.numeric(matrices[[platform]][feature, sample_metadata$sample_key])
    if (!any(is.finite(values))) {
      unavailable[technology] <- "Not available"
      next
    }
    sample_mean <- mean(values)
    sample_sd <- stats::sd(values)
    if (!is.finite(sample_sd) || sample_sd == 0) {
      unavailable[technology] <- "No variation\nz-score undefined"
      next
    }
    sample_z <- (values - sample_mean) / sample_sd
    donor_means <- vapply(donor_order, function(donor) {
      mean(values[as.character(sample_metadata$donor) == donor])
    }, numeric(1))
    donor_z_means <- vapply(donor_order, function(donor) {
      mean(sample_z[as.character(sample_metadata$donor) == donor])
    }, numeric(1))
    donor_df_list[[technology]] <- data.frame( technology = technology, donor = donor_order,
      sex = sample_metadata$sex[donor_rows], donor_mean = unname(donor_means), z_score = unname(donor_z_means))

    visit_df_list[[technology]] <- data.frame(
      technology = technology, donor = sample_metadata$donor, day = sample_metadata$day, z_score = sample_z)
  }
  donor_df <- do.call(rbind, donor_df_list)
  rownames(donor_df) <- NULL
  donor_df$technology <- factor(donor_df$technology, levels = technology_labels)
  donor_df$donor <- factor(donor_df$donor, levels = donor_order)
  donor_df$sex <- factor(donor_df$sex, levels = names(col_sex))
  visit_df <- do.call(rbind, visit_df_list)
  rownames(visit_df) <- NULL
  visit_df$technology <- factor(visit_df$technology, levels = technology_labels)
  visit_df$donor <- factor(visit_df$donor, levels = donor_order)
  visit_df <- visit_df[order(visit_df$technology, visit_df$donor, visit_df$day), ]
  list(data = donor_df, visits = visit_df, technology_labels = technology_labels, unavailable = unavailable)
}
# 05.2. Compare sex effects across platforms ----------------------------------

plot_protein_sex <- function(protein, matrices = expression_matrices, sample_metadata = plot_metadata, mapping = soma_map,
  evidence = soma_dep_evidence, drop_unavailable = FALSE) {
  sex_display <- prepare_protein_donor_zscore(protein, matrices, sample_metadata, mapping,
    effect = "Sex", evidence = evidence)
  if (drop_unavailable) {
    sex_display$data$technology <- droplevels(sex_display$data$technology)
    sex_display$technology_labels <- levels(sex_display$data$technology)
    sex_display$unavailable <- character()
  }
  p_sex <- ggplot(sex_display$data, aes(x = technology, y = z_score, fill = sex)) +
    geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(aes(group = interaction(technology, sex)), width = 0.6,
    position = position_dodge(width = 0.72), outlier.shape = NA, colour = "grey35", linewidth = 0.3, alpha = 0.45) +
    geom_point(aes(group = sex), shape = 21, size = 2, stroke = 0.25, colour = "grey30",
    position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.72, seed = 81)) +
    scale_fill_manual(values = col_sex, limits = names(col_sex), drop = FALSE, name = "Sex") +
    scale_x_discrete(limits = sex_display$technology_labels, drop = FALSE) +
    labs(title = protein, subtitle = "Five-platform shared Sex-associated protein (raw P < 0.05)",
    x = "Technology / SOMAmer", y = "Mean of three sample z-scores") + trajectory_theme + theme(legend.position = "bottom")
  if (length(sex_display$unavailable)) {
    missing_df <- data.frame(technology = names(sex_display$unavailable), label = unname(sex_display$unavailable))
    p_sex <- p_sex + geom_text(data = missing_df, aes(x = technology, y = 0, label = label),
      inherit.aes = FALSE, size = 2.4, colour = "grey40")
  }
  p_sex
}
# 05.3. Compare visits within each Donor --------------------------------------

plot_protein_individual <- function(protein, matrices = expression_matrices,
  sample_metadata = plot_metadata, mapping = soma_map,
  day_labels = c("Day 0", "Day 3", "Day 5")) {
  individual_display <- prepare_protein_donor_zscore(protein, matrices, sample_metadata, mapping, effect = "Individual")
  p_individual <- ggplot(individual_display$visits,
    aes(x = factor(day, levels = c(0, 3, 5)), y = z_score, colour = donor, group = donor)) +
    geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed", linewidth = 0.3) +
    geom_line(linewidth = 0.45, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = col_donor, limits = names(col_donor), drop = FALSE,
    labels = function(x) sub("^Donor ", "Individual ", x), name = "Individual") +
    scale_x_discrete(limits = c("0", "3", "5"), labels = day_labels, drop = FALSE, expand = expansion(add = 0.35)) +
    facet_grid(cols = vars(technology), drop = FALSE, switch = "x") +
    labs(title = protein, x = "Sampling time within technology / SOMAmer", y = "Standardized expression (sample z-score)") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    trajectory_theme + theme(legend.position = "bottom", strip.placement = "outside",
    strip.background = element_blank(), strip.text.x = element_text(size = 8), panel.spacing.x = grid::unit(0.5, "cm"))
  if (length(individual_display$unavailable)) {
    missing_df <- data.frame(technology = names(individual_display$unavailable),
      label = unname(individual_display$unavailable))
    p_individual <- p_individual + geom_text(data = missing_df,
      aes(x = "3", y = 0, label = label), inherit.aes = FALSE, size = 2.4, colour = "grey40")
  }
  p_individual
}
# 06. Prepare the selected protein plot lists ---------------------------------

p_individual_list <- setNames(vector("list", length(individual_examples)), individual_examples)
for (protein in individual_examples) {
  p_individual_list[[protein]] <- plot_protein_individual(protein)
}
p_sex_list <- setNames(vector("list", length(sex_examples)), sex_examples)
for (protein in sex_examples) {
  p_sex_list[[protein]] <- plot_protein_sex(protein)
}
# 10. Combine the selected individual-associated examples ---------------------

individual_protein_dir <- wdplot
dir.create(individual_protein_dir, recursive = TRUE, showWarnings = FALSE)

selected_individual_proteins <- c("LPA", "RNASE3")
individual_example_plots <- lapply(selected_individual_proteins, function(protein) {
  plot_protein_individual(protein) + labs(title = protein, subtitle = NULL)
})
p_individual_examples <- wrap_plots(individual_example_plots, ncol = 2, nrow = 1, widths = c(5, 6), guides = "collect") &
  theme(legend.position = "bottom")
individual_examples_pdf <- file.path(individual_protein_dir, "F6d_LPA_RNASE3.pdf")

ggsave(individual_examples_pdf, pdf_text_grob(p_individual_examples), width = 28, height = 10, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 11. Save the selected sex-associated example --------------------------------

sex_single_examples <- unique(c(sex_examples, "INSL3", "KLK3", "SERPINF1", "SHBG", "KLK14"))
sex_nominal_union <- unique(unlist(sex_nominal_sets, use.names = FALSE))
sex_single_evidence <- set_sex_nominal_display_ids(sex_single_examples)
p_sex_single_list <- setNames(vector("list", length(sex_single_examples)), sex_single_examples)
for (protein in sex_single_examples) {
  p_sex_single_list[[protein]] <- plot_protein_sex(protein, evidence = sex_single_evidence, drop_unavailable = TRUE)
}
sex_protein_dir <- wdplot
dir.create(sex_protein_dir, recursive = TRUE, showWarnings = FALSE)

for (protein in "SHBG") {
  protein_pdf <- file.path(sex_protein_dir, paste0("F6e_", protein, ".pdf"))

  p_protein <- p_sex_single_list[[protein]] + labs(title = protein, subtitle = NULL)
  ggsave(protein_pdf, pdf_text_grob(p_protein), width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
}
# 12. Combine the remaining sex-associated proteins ---------------------------

sex_combined_proteins <- setdiff(names(p_sex_single_list), "SHBG")
p_sex_combined_list <- setNames(vector("list", 6), sex_combined_proteins)
for (protein in sex_combined_proteins) {
  p_sex_combined_list[[protein]] <- p_sex_single_list[[protein]] + labs(title = protein, subtitle = NULL)
}
p_sex_combined <- wrap_plots(p_sex_combined_list, ncol = 3, nrow = 2, byrow = TRUE, guides = "collect") &
  theme(legend.position = "bottom")

sex_combined_pdf <- file.path(wdplot, "S11_sex_associated_proteins_without_shbg_2row_3col.pdf")

dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)
ggsave(sex_combined_pdf, pdf_text_grob(p_sex_combined), width = 36, height = 16, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 13.1. Combine the 25 selected individual-associated proteins ----------------

individual_25_proteins <- c( "ADH4", "AHSP", "ANG", "BDNF", "CCL18", "CCL5", "CD5L", "COL6A3", "CRH",
  "CXCL3", "F12", "F7", "GHRL", "GSR", "HLA-A", "IGFBP2", "LEP", "LPA",
  "PF4", "PPBP", "RNASE3", "SHBG", "SPINK5", "TGFBI", "TIMP3"
)
p_individual_25_list <- setNames(vector("list", length(individual_25_proteins)), individual_25_proteins)
for (protein in individual_25_proteins) {
  p_individual_25_list[[protein]] <- suppressMessages( plot_protein_individual(protein) +
    scale_x_discrete(limits = c("0", "3", "5"), labels = c("0", "3", "5"), drop = FALSE, expand = expansion(add = 0.35)) +
    labs(title = protein, subtitle = NULL, x = "Day", y = "Z-score")
  )
}
individual_25_ncol <- 3L
individual_25_nrow <- ceiling(length(individual_25_proteins) / individual_25_ncol)
p_individual_25_combined <- wrap_plots(p_individual_25_list, ncol = individual_25_ncol, nrow = individual_25_nrow,
  byrow = TRUE, guides = "collect") & theme(legend.position = "bottom")

# 13.1.1. Save the selected individual-associated assembly --------------------

individual_25_pdf <- file.path(wdplot, "S10_individual_associated_proteins_selected_25_3col.pdf")

dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

ggsave(individual_25_pdf, pdf_text_grob(p_individual_25_combined), width = 35, height = 44,
  units = "cm", limitsize = FALSE, device = grDevices::cairo_pdf, family = "Helvetica")
# 14. Export the five-platform associated-protein union -----------------------

build_associated_union_table <- function(effect = c("Sex", "Individual"),
  selection = c("P", "FDR")) {
  effect <- match.arg(effect)
  selection <- match.arg(selection)
  required_columns <- c("platform", "effect", "Protein", "min_p_value", "min_FDR",
    "feature_ids", "n_features_tested", "n_nominal_significant", "n_fdr_significant")
  evidence <- dep_evidence[dep_evidence$effect == effect &
    dep_evidence$platform %in% platform_keys, , drop = FALSE]
  selected_values <- if (selection == "P") {
    evidence$min_p_value
  } else {
    evidence$min_FDR
  }
  proteins <- sort(unique(as.character(evidence$Protein[
    is.finite(selected_values) & selected_values < 0.05])))
  union_df <- data.frame(Protein = proteins, Selection = paste0(effect, ": ", selection, " < 0.05"),
    stringsAsFactors = FALSE, check.names = FALSE)

  p_column <- if (effect == "Sex") {
    "Sex_pval"
  } else {
    "indi_pval"
  }
  soma_model_df <- data.frame(SOMAmer = as.character(soma_sex_model$Protein),
    P = as.numeric(soma_sex_model[[p_column]]), stringsAsFactors = FALSE)
  soma_model_df$FDR <- stats::p.adjust(soma_model_df$P, method = "BH")

  for (statistic in c("P", "FDR")) {
    for (platform in platform_keys) {
      rows <- evidence[evidence$platform == platform, , drop = FALSE]
      column <- paste0(platform_labels[[platform]], if (platform == "Somascan") {
        "_min_"
      } else {
        "_"
      }, statistic)
      source_column <- if (statistic == "P") {
        "min_p_value"
      } else {
        "min_FDR"
      }
      union_df[[column]] <- as.numeric(rows[[source_column]][match(proteins, rows$Protein)])
    }
  }
  for (statistic in c("P", "FDR")) {
    union_df[[paste0("Associated_platforms_", statistic, "_lt_0.05")]] <- vapply(proteins, function(protein) {
      rows <- evidence[evidence$Protein == protein, , drop = FALSE]
      values <- if (statistic == "P") {
        rows$min_p_value
      } else {
        rows$min_FDR
      }
      paste(unname(platform_labels[as.character(rows$platform[
        is.finite(values) & values < 0.05])]), collapse = ";")
    }, character(1))
  }
  union_df$SomaScan_SOMAmer_ids <- NA_character_
  union_df$SomaScan_min_P_SOMAmer <- NA_character_
  union_df$SomaScan_min_FDR_SOMAmer <- NA_character_
  soma_rows <- list()

  for (i in seq_along(proteins)) {
    protein <- proteins[i]
    row <- evidence[evidence$platform == "Somascan" & evidence$Protein == protein, , drop = FALSE]
    if (!nrow(row)) {
      next
    }
    ids <- strsplit(as.character(row$feature_ids), ";", fixed = TRUE)[[1]]
    ids <- sort(unique(trimws(ids[nzchar(trimws(ids))])))
    model_idx <- match(ids, soma_model_df$SOMAmer)
    mapped_ids <- as.character(soma_map$SeqId[soma_map$EntrezGeneSymbol == protein])
    feature_df <- soma_model_df[model_idx, , drop = FALSE]
    finite_p <- feature_df$P[is.finite(feature_df$P)]
    finite_fdr <- feature_df$FDR[is.finite(feature_df$FDR)]
    observed_p <- if (length(finite_p)) {
      min(finite_p)
    } else {
      NA_real_
    }
    observed_fdr <- if (length(finite_fdr)) {
      min(finite_fdr)
    } else {
      NA_real_
    }
    union_df$SomaScan_SOMAmer_ids[i] <- paste(ids, collapse = ";")
    union_df$SomaScan_min_P_SOMAmer[i] <- paste(ids[is.finite(feature_df$P) &
      feature_df$P == observed_p], collapse = ";")
    union_df$SomaScan_min_FDR_SOMAmer[i] <- paste(ids[is.finite(feature_df$FDR) &
      feature_df$FDR == observed_fdr], collapse = ";")
    soma_rows[[protein]] <- data.frame(Protein = protein, Effect = effect,
      SOMAmer = ids, P = feature_df$P, FDR = feature_df$FDR, P_lt_0.05 = is.finite(feature_df$P) & feature_df$P < 0.05,
      FDR_lt_0.05 = is.finite(feature_df$FDR) & feature_df$FDR < 0.05, stringsAsFactors = FALSE, check.names = FALSE)
  }
  soma_detail_df <- if (length(soma_rows)) {
    do.call(rbind, soma_rows)
  } else {
    data.frame(Protein = character(), Effect = character(), SOMAmer = character(),
      P = numeric(), FDR = numeric(), P_lt_0.05 = logical(), FDR_lt_0.05 = logical(), check.names = FALSE)
  }
  rownames(soma_detail_df) <- NULL
  list(proteins = union_df, somamers = soma_detail_df)
}
# 15. Export nominal-P sex-associated proteins --------------------------------

sex_union_tables <- build_associated_union_table("Sex", selection = "P")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
sex_union_file <- file.path(wdresult, paste0(plot_prefix, "15_sex_associated_protein_union_P_FDR.csv"))

write.csv(sex_union_tables$proteins, sex_union_file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

# 16. Export BH-FDR individual-associated proteins ----------------------------

individual_union_tables <- build_associated_union_table("Individual", selection = "FDR")
individual_union_file <- file.path(wdresult, paste0(plot_prefix, "16_individual_associated_protein_union_P_FDR.csv"))

write.csv(individual_union_tables$proteins, individual_union_file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

saveRDS(list(protein_evidence = model_bundle$c8_dep_protein_evidence,
  sex = sex_union_tables,
  individual = individual_union_tables,
  individual_profiles = p_individual_list,
  sex_profiles = p_sex_list),
  paste0(wdresult, "DEP_results.rds"))
