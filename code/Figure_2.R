# Figure 2: complete analysis and manuscript figures
# Run from the repository root.
source(file.path("code", "setup.R"), encoding = "UTF-8")
wdresult <- paste0(file.path(repo_dir, "result", "Fig2"), "/")
wdplot <- paste0(file.path(repo_dir, "plot", "Fig2"), "/")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

plot_prefix <- "01-"
platform_levels <- c("MS-OT480", "MS-Astral","Olink-Explore", "Olink-HT","SomaScan")

# 01. Load the registered detection inputs prepared from the parent project --

f2_input_file <- file.path(wddata, "processed", "Figure_2_detection_inputs.rds")
f2_inputs <- readRDS(f2_input_file)
list2env(f2_inputs, envir = environment())
metadata <- read_processed_metadata(wddata)
metadata <- metadata[1:30, , drop = FALSE]
metadata$somaid <- paste0("S", metadata$OriginalID)
metadata$Time <- dplyr::recode(metadata$Time, "visit1" = "Day0", "visit2" = "Day3", "visit3" = "Day5")
metadata$Time <- factor(metadata$Time, levels = c("Day0", "Day3", "Day5"))
metadata$SubjectID <- factor(metadata$SubjectID, levels = unique(metadata$SubjectID))
saveRDS(detection_classification_sets, paste0(wdresult, plot_prefix, "01_detection_classification_sets.rds"))

# 01-01. Plot protein detection counts ----------------------------------------

protein_category_order <- c( "Predefined protein targets", "Detected proteins", "Individual detected proteins",
  "Broadly detected proteins", "Undetected predefined protein"
)

df_plot <- tibble::tibble( platform = rep(platform_levels, each = length(protein_category_order)),
  protein_set = rep(protein_category_order, times = length(platform_levels)), n_protein = c(rbind(
  f2a_predefined_counts[platform_levels], f2a_detected_counts[platform_levels],
  f2a_individual_detected_counts[platform_levels], f2a_broad_detected_counts[platform_levels],
  f2a_undetected_counts[platform_levels]
))
) %>%
  dplyr::mutate(platform = factor(platform, levels = platform_levels),
  protein_set = factor(protein_set, levels = protein_category_order))

p_protein_counts <- ggplot(df_plot, aes(x = platform, y = n_protein, fill = protein_set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_text(aes(label = n_protein), position = position_dodge(width = 0.7), vjust = 0,hjust = 0,angle = 45,
  size = 2.5) + scale_fill_manual(values = c("Predefined protein targets" = "#439EB7",
  "Detected proteins" = "#475272", "Individual detected proteins" = "#AEC7E8FF",
  "Broadly detected proteins" = "#D28770", "Undetected predefined protein" = "#D3D3D3"),name = "") +
  theme2_nobacknoerror + labs(x = NULL, y = "Number of proteins") + theme(axis.text.y = element_text(angle = 90))

ggsave(paste0(wdplot, "F2a_protein_count_barplot.pdf"),pdf_text_grob(p_protein_counts), width = 15, height = 5.5, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

affinity_venn_file <- paste0(wdplot,"S1b_affinity_detection_classification_venn.pdf")
grDevices::cairo_pdf(affinity_venn_file, width = 18 / 2.54, height = 7.5 / 2.54, family = "Helvetica")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(
  0.86,0.14), "null")),clip = "off"))
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,3), layout.pos.row = 1, clip = "off"))

for (i in seq_len(nrow(f2a_detection_split_summary))) {
  summary_row <- f2a_detection_split_summary[i, ,drop = FALSE]
  grid::pushViewport(grid::viewport(layout.pos.col = i,clip = "off"))
  venn_grob <- VennDiagram::draw.pairwise.venn(area1 = as.numeric(summary_row$predefined_n),
    area2 = as.numeric(summary_row$detected_n), cross.area = as.numeric(summary_row$predefined_detected_overlap_n),
    category = c("", ""), fill = c("#439EB7","#475272"), alpha = c(0.8, 0.8),cat.cex = 0,cex = c(1.15,0, 1.15),
    margin = 0.18, ind = FALSE)
  undetected_n <- as.numeric(summary_row$undetected_n)
  numeric_label_nodes <- which(vapply( seq_along(venn_grob),
    function(i) {
    inherits(venn_grob[[i]],"text") &&
      grepl("^[0-9]+$", venn_grob[[i]][["label"]])
  }, logical(1)
  ))

  if (undetected_n > 0 && length(numeric_label_nodes) >= 2) {
    undetected_node <- numeric_label_nodes[which.max(vapply( venn_grob[numeric_label_nodes],
      function(node) as.numeric(node[["x"]]), numeric(1)
    ))]
    venn_grob[[undetected_node]][["label"]] <- paste0(undetected_n, " Undetected")
  }
  grid::grid.draw(pdf_text_grob(venn_grob))
  grid::grid.text(summary_row$platform, x = grid::unit(0.5, "npc"), y = grid::unit(0.98, "npc"), just = "top",
    gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica", col = "#000000"))
  if (as.numeric(summary_row$undetected_n) == 0) {
    grid::grid.text("0 Undetected", x = grid::unit(0.5,"npc"), y = grid::unit(0.08, "npc"), just = "bottom",
      gp = grid::gpar(fontsize = 7, fontfamily = "Helvetica", col = "#000000"))
  }
  grid::popViewport()
}
grid::popViewport()
grid::pushViewport(grid::viewport(layout.pos.row = 2, clip = "off"))
grid::grid.rect(x = grid::unit(0.38, "npc"), y = grid::unit(0.55, "npc"), width = grid::unit(0.025, "npc"),
  height = grid::unit(0.16, "npc"), gp = grid::gpar(fill = "#439EB7", col = NA))
grid::grid.text("Predefined protein targets", x = grid::unit(0.40, "npc"),y = grid::unit(0.55, "npc"),
  just = "left", gp = grid::gpar(fontsize = 7, fontfamily = "Helvetica", col = "#000000"))
grid::grid.rect(x = grid::unit(0.62,"npc"),y = grid::unit(0.55,"npc"),width = grid::unit(0.025,"npc"),
  height = grid::unit(0.16,"npc"), gp = grid::gpar(fill = "#475272",col = NA))
grid::grid.text("Detected proteins", x = grid::unit(0.64, "npc"), y = grid::unit(0.55,"npc"),just = "left",
  gp = grid::gpar(fontsize = 7, fontfamily = "Helvetica", col = "#000000"))
grid::popViewport()
grid::popViewport()
dev.off()

# 01-02. Plot detection frequency across visits -------------------------------

summarise_soma_visit_detection <- function(soma_df, protein_set = NULL) {
  if (!is.null(protein_set)) {
    soma_df <- soma_df %>% dplyr::filter(EntrezGeneSymbol %in% protein_set)
  }
  somamer_df <- soma_df %>% dplyr::filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    dplyr::group_by(SubjectID,SeqId, EntrezGeneSymbol) %>%
    dplyr::summarise(n_valid_visits = sum(valid), n_detected_visits = sum(detected),
    somer_detection_rate = dplyr::if_else(n_valid_visits > 0, n_detected_visits / n_valid_visits,NA_real_), .groups = "drop")

  protein_visit_df <- soma_df %>% dplyr::filter(!is.na(EntrezGeneSymbol),EntrezGeneSymbol != "") %>%
    dplyr::group_by(SubjectID,Time, EntrezGeneSymbol) %>%
    dplyr::summarise(n_somers = dplyr::n_distinct(SeqId), n_valid_somers = sum(valid),
    n_detected_somers = sum(detected), protein_visit_valid = n_valid_somers > 0,
    protein_visit_detected = n_detected_somers > 0, .groups = "drop")

  somamer_stats <- somamer_df %>%
    dplyr::group_by(SubjectID, EntrezGeneSymbol) %>% dplyr::summarise( n_somers = dplyr::n(),
    min_somer_detection_rate = if (any(!is.na(somer_detection_rate))) {
    min(somer_detection_rate,na.rm = TRUE)
  } else {
    NA_real_
  },

    max_somer_detection_rate = if (any(!is.na(somer_detection_rate))) {
    max(somer_detection_rate,na.rm = TRUE)
  } else {
    NA_real_
  },

    sd_somer_detection_rate = if (sum(!is.na(somer_detection_rate)) > 1) {
    stats::sd(somer_detection_rate, na.rm = TRUE)
  } else {
    NA_real_
  },

    n_detected_somers = sum(n_detected_visits > 0), protein_detected_any = n_detected_somers > 0, .groups = "drop"
  )

  protein_df <- protein_visit_df %>% dplyr::group_by(SubjectID, EntrezGeneSymbol) %>%
    dplyr::summarise(n_valid_visits = sum(protein_visit_valid), n_detected_visits = sum(protein_visit_detected),
    n_somers = max(n_somers), protein_detected_any = any(protein_visit_detected),.groups = "drop") %>%
    dplyr::mutate(protein_detection_score = dplyr::if_else(n_valid_visits > 0, n_detected_visits / n_valid_visits, NA_real_)) %>%

    dplyr::left_join(somamer_stats %>% dplyr::select(-n_somers, -protein_detected_any),
    by = c("SubjectID", "EntrezGeneSymbol"))

  list(somamer = somamer_df,protein = protein_df)
}
summarise_donor_detection_frequency <- function(score_df, score_column) {
  score_df %>% dplyr::transmute(SubjectID,number = .data[[score_column]]) %>% dplyr::filter(!is.na(number)) %>%
    dplyr::count(SubjectID, number, name = "count") %>% dplyr::group_by(SubjectID) %>%
    dplyr::mutate(percentage = count / sum(count)) %>%
    dplyr::ungroup()
}

det_freq_ot480 <- indi_prot(ms_ot480_uniquepep2[covered_protein_sets[["MS-OT480"]],], cate = "ms",metadata)
det_freq_astral <- indi_prot(ms_astral_uniquepep2[covered_protein_sets[["MS-Astral"]],], cate = "ms", metadata)
det_freq_explore <- indi_prot(lod3072 %>% dplyr::filter(Assay %in% covered_protein_sets[["Olink-Explore"]]),
  cate = "olink",metadata)
det_freq_ht <- indi_prot(lodht %>% dplyr::filter(Assay %in% covered_protein_sets[["Olink-HT"]]),cate = "olink", metadata)

soma_visit_full <- summarise_soma_visit_detection(soma_detection_master, covered_protein_sets[["SomaScan"]])

det_freq_soma <- summarise_donor_detection_frequency(soma_visit_full$protein, "n_detected_visits") %>%
  dplyr::group_by(SubjectID) %>%
  tidyr::complete(number = c(0, 1,2, 3), fill = list(count = 0, percentage = 0)) %>% dplyr::ungroup()

plot_individual_detection <- function(freq_df, tit) {
  p <- plot_indipro(freq_df, tit = tit)
  colour_scale <- p$scales$get_scales("colour")
  colour_scale$name <- "Individual"
  colour_scale$labels <- sprintf("Individual %02d", seq_along(colour_scale$breaks))
  p
}
detection_percentage_scale <- ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1,by = 0.2),
  labels = function(x) x * 100,expand = c(0, 0))
p_detection_frequency_ms_ot480 <- plot_individual_detection(det_freq_ot480, tit = "MS-OT480") + detection_percentage_scale +
  labs(x = "", y = "Percentage of detected proteins (%)")
p_detection_frequency_ms_astral <- plot_individual_detection(det_freq_astral,tit = "MS-Astral") + detection_percentage_scale +
  labs(x = "", y = "") + themeynone
p_detection_frequency_olink_explore <- plot_individual_detection(det_freq_explore, tit = "Olink-Explore") +
  detection_percentage_scale + labs(x = "Detection frequency across visits",y = "") + themeynone
p_detection_frequency_olink_ht <- plot_individual_detection(det_freq_ht, tit = "Olink-HT") + detection_percentage_scale +
  labs(x = "", y = "") + themeynone
p_detection_frequency_soma <- plot_individual_detection(det_freq_soma, tit = "SomaScan") + detection_percentage_scale +
  labs(x = "", y = "") + themeynone

grDevices::cairo_pdf(paste0(wdplot, "F2d_individual_detection_freq.pdf"), width = 11 / 2.54, height = 5.5 / 2.54, family = "Helvetica")
grid::grid.draw(pdf_text_grob(p_detection_frequency_ms_ot480 + p_detection_frequency_ms_astral + p_detection_frequency_olink_explore +
  p_detection_frequency_olink_ht + p_detection_frequency_soma + plot_layout(widths = c(1, 1, 1, 1, 1))))
dev.off()

grDevices::cairo_pdf(paste0(wdplot, "F2d_individual_detection_freq_legend.pdf"), width = 15 / 2.54, height = 5.5 / 2.54, family = "Helvetica")
grid::grid.draw(pdf_text_grob(p_detection_frequency_ms_ot480 + themelegendtop))
dev.off()

det_freq_ot480$cate <- "MS-OT480"
det_freq_astral$cate <- "MS-Astral"
det_freq_explore$cate <- "Olink-Explore"
det_freq_ht$cate <- "Olink-HT"
det_freq_soma$cate <- "SomaScan"
det_freq_bind <- dplyr::bind_rows(det_freq_ot480, det_freq_astral, det_freq_explore,det_freq_ht, det_freq_soma)

det_pct_levels <- c(3)
det_pct_ot480 <- det_freq_ot480 %>% dplyr::filter(number %in% det_pct_levels) %>% dplyr::group_by(SubjectID) %>%
  dplyr::summarise(total_percentage = sum(percentage,na.rm = TRUE))
det_pct_astral <- det_freq_astral %>% dplyr::filter(number %in% det_pct_levels) %>%
  dplyr::group_by(SubjectID) %>% dplyr::summarise(total_percentage = sum(percentage,na.rm = TRUE))
det_pct_explore <- det_freq_explore %>% dplyr::filter(number %in% det_pct_levels) %>%
  dplyr::group_by(SubjectID) %>% dplyr::summarise(total_percentage = sum(percentage, na.rm = TRUE))
det_pct_ht <- det_freq_ht %>% dplyr::filter(number %in% det_pct_levels) %>% dplyr::group_by(SubjectID) %>%
  dplyr::summarise(total_percentage = sum(percentage, na.rm = TRUE))
det_pct_soma <- det_freq_soma %>% dplyr::filter(number %in% det_pct_levels) %>% dplyr::group_by(SubjectID) %>%
  dplyr::summarise(total_percentage = sum(percentage,na.rm = TRUE))
median(det_pct_ot480$total_percentage)
median(det_pct_astral$total_percentage)
median(det_pct_explore$total_percentage)
median(det_pct_ht$total_percentage)
median(det_pct_soma$total_percentage)

# 01-03. Summarise donor support ----------------------------------------------

build_ms_donor_detection <- function(matrix_data, platform_name, proteins,metadata_df) {
  sample_ids <- as.character(metadata_df$OriginalID)

  dplyr::bind_rows(lapply(proteins, function(protein_name) {
    values <- matrix_data[protein_name, sample_ids, drop = TRUE]
    tibble::tibble( platform = platform_name, protein = protein_name,
      OriginalID = sample_ids,detected = !is.na(values) & values != 0
    )
  })) %>%
    dplyr::left_join(metadata_df %>% dplyr::select(OriginalID, SubjectID), by = "OriginalID")
}
build_olink_donor_detection <- function(lod_data,platform_name,proteins,metadata_df) {
  lod_data %>% dplyr::filter(Assay %in% proteins) %>%
    dplyr::transmute(platform = platform_name, protein = as.character(Assay),
    OriginalID = as.character(SampleID), detected = !is.na(NPX) & !is.na(LOD) & NPX >= LOD) %>%
    dplyr::inner_join(metadata_df %>% dplyr::select(OriginalID,SubjectID),by = "OriginalID") %>%
    dplyr::group_by(platform,protein, OriginalID,SubjectID) %>% dplyr::summarise(detected = any(detected),.groups = "drop")
}
soma_donor_detection_df <- soma_detection_master %>% dplyr::filter(EntrezGeneSymbol %in% detected_protein_sets[["SomaScan"]]) %>%
  dplyr::group_by(EntrezGeneSymbol, Sample,SubjectID) %>% dplyr::summarise(detected = any(detected), .groups = "drop") %>%
  dplyr::transmute(platform = "SomaScan", protein = EntrezGeneSymbol, SubjectID, detected)

protein_visit_detection_df <- dplyr::bind_rows(
  build_ms_donor_detection(ms_ot480_uniquepep2,"MS-OT480", detected_protein_sets[["MS-OT480"]], metadata),
  build_ms_donor_detection(ms_astral_uniquepep2,"MS-Astral", detected_protein_sets[["MS-Astral"]],metadata),
  build_olink_donor_detection(lod3072, "Olink-Explore", detected_protein_sets[["Olink-Explore"]],metadata),
  build_olink_donor_detection(lodht, "Olink-HT",detected_protein_sets[["Olink-HT"]], metadata), soma_donor_detection_df
)

protein_donor_detection_df <- protein_visit_detection_df %>% dplyr::group_by(platform, protein,SubjectID) %>%
  dplyr::summarise(n_detected_visits = sum(detected),individual_level_detection = n_detected_visits >= 2, .groups = "drop")

individual_detected_targets_df <- dplyr::bind_rows(lapply(platform_levels, function(platform_name) {
  tibble::tibble(platform = platform_name, protein = individual_detected_protein_sets[[platform_name]])
}))

donor_support_distribution_03 <- protein_donor_detection_df %>% dplyr::inner_join(individual_detected_targets_df,
  by = c("platform","protein")) %>%
  dplyr::group_by(platform, protein) %>%
    dplyr::summarise(n_donors_with_individual_detection = sum(individual_level_detection),.groups = "drop") %>%
  dplyr::count(platform, n_donors_with_individual_detection,name = "n_detected_proteins") %>%
  dplyr::group_by(platform) %>% tidyr::complete(n_donors_with_individual_detection = 1:10, fill = list(n_detected_proteins = 0L)) %>%
  dplyr::mutate(percentage = n_detected_proteins / sum(n_detected_proteins)) %>% dplyr::ungroup() %>%
  dplyr::mutate(platform = factor(platform,levels = platform_levels),
  n_donors_with_individual_detection = factor(n_donors_with_individual_detection,levels = 1:10)) %>%
  dplyr::arrange(platform, n_donors_with_individual_detection)

plot_donor_support <- function(platform_name, input_df,y_column, y_label) {
  plot_df <- input_df %>% dplyr::filter(platform == platform_name)
  ggplot2::ggplot(plot_df,ggplot2::aes(x = n_donors_with_individual_detection,y = .data[[y_column]])) +
    ggplot2::geom_col(width = 0.68,fill = "#AEC7E8FF",colour = "white", linewidth = 0.2) +
    ggplot2::labs(x = "Number of individuals with individual-level detection",y = y_label, title = platform_name) +
    theme1_noback + themey90
}
donor_percentage_plot_list <- lapply(platform_levels, plot_donor_support, input_df = donor_support_distribution_03,
  y_column = "percentage", y_label = "Percentage of individual detected proteins (%)"
)
donor_percentage_plot_list <- lapply(donor_percentage_plot_list, function(p) {
  p + ggplot2::scale_y_continuous(breaks = seq(0,1, by = 0.2),limits = c(0, 1), labels = function(x) x * 100,
    expand = ggplot2::expansion(mult = c(0,0.03)))
})

for (i in seq_along(donor_percentage_plot_list)) {
  if (i != 3) {
    donor_percentage_plot_list[[i]] <- donor_percentage_plot_list[[i]] + labs(x = NULL)
  }
}
for (i in 2:length(donor_percentage_plot_list)) {
  donor_percentage_plot_list[[i]] <- donor_percentage_plot_list[[i]] + themeynone
}
donor_percentage_plot_list[[3]] <- donor_percentage_plot_list[[3]] +
  labs(x = "Number of individuals with individual-level detection")

grDevices::cairo_pdf(paste0(wdplot, "F2e_individual_detected_protein_donor_support_percentage.pdf"), width = 11 / 2.54, height = 5.5 / 2.54, family = "Helvetica")
grid::grid.draw(pdf_text_grob(donor_percentage_plot_list[[1]] + donor_percentage_plot_list[[2]] + donor_percentage_plot_list[[3]] +
  donor_percentage_plot_list[[4]] + donor_percentage_plot_list[[5]] +
  plot_layout(widths = c(1, 1, 1, 1, 1))))
dev.off()
# 01-04. Plot donor support outside broadly detected sets ---------------------

donor_subset_platform_levels <- c("MS-OT480", "MS-Astral","Olink-Explore", "Olink-HT")

individual_not_broad_targets_df <- dplyr::bind_rows(lapply( donor_subset_platform_levels,
  function(platform_name) {
  tibble::tibble( platform = platform_name,
    protein = setdiff(individual_detected_protein_sets[[platform_name]], broad_detected_protein_sets[[
    platform_name]])
  )
}))

donor_support_distribution_03b <- protein_visit_detection_df %>% dplyr::inner_join(individual_not_broad_targets_df,
  by = c("platform", "protein")) %>%
  dplyr::group_by(platform, protein,SubjectID) %>% dplyr::summarise(n_detected_visits = sum(detected), .groups = "drop") %>%
  dplyr::group_by(platform, protein) %>%
  dplyr::summarise(n_individuals_with_individual_detection = sum(n_detected_visits >= 2),.groups = "drop") %>%
  dplyr::count(platform,n_individuals_with_individual_detection,name = "n_individual_not_broad_proteins") %>%
  dplyr::group_by(platform) %>% tidyr::complete(n_individuals_with_individual_detection = 1:10,
    fill = list(n_individual_not_broad_proteins = 0L)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(platform = factor(platform, levels = donor_subset_platform_levels),
  n_individuals_with_individual_detection = factor(n_individuals_with_individual_detection, levels = 1:10)) %>%
  dplyr::arrange(platform, n_individuals_with_individual_detection)

donor_subset_count_y_max <- max(donor_support_distribution_03b$n_individual_not_broad_proteins)
donor_subset_colours <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E", "Olink-Explore" = "#264653",
  "Olink-HT" = "#AEC7E8FF"
)
donor_subset_plot_list <- vector("list", length(donor_subset_platform_levels))

for (platform_index in seq_along(donor_subset_platform_levels)) {
  platform_name <- donor_subset_platform_levels[[platform_index]]
  plot_df <- donor_support_distribution_03b %>% dplyr::filter(platform == platform_name)
  donor_subset_plot_list[[platform_index]] <- ggplot2::ggplot(plot_df,
    ggplot2::aes(x = n_individuals_with_individual_detection, y = n_individual_not_broad_proteins)) +
    ggplot2::geom_col(width = 0.68, fill = donor_subset_colours[[platform_name]], colour = "#475272", linewidth = 0.3) +

    ggplot2::geom_text(ggplot2::aes(label = n_individual_not_broad_proteins), vjust = -0.35, size = 3) +
    ggplot2::scale_y_continuous(limits = c(0, donor_subset_count_y_max * 1.12),
    expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "Number of individuals with individual-level detection",
    y = "Number of individual detected - broadly detected proteins", title = platform_name) + theme2_nobackplus +
    themetext30
}
for (i in seq_along(donor_subset_plot_list)) {
  if (i != 3) {
    donor_subset_plot_list[[i]] <- donor_subset_plot_list[[i]] + labs(x = NULL)
  }
}
for (i in 2:length(donor_subset_plot_list)) {
  donor_subset_plot_list[[i]] <- donor_subset_plot_list[[i]] + themeynone
}
grDevices::cairo_pdf(paste0(wdplot, "S1c_individual_not_broad_detection_support_count.pdf"), width = 24 / 2.54, height = 6 / 2.54, family = "Helvetica")
grid::grid.draw(pdf_text_grob(donor_subset_plot_list[[1]] + donor_subset_plot_list[[2]] + donor_subset_plot_list[[3]] + donor_subset_plot_list[[4]] +
  plot_layout(widths = c(1, 1, 1,1))))
dev.off()
# 01-05. Plot the detected-protein UpSet --------------------------------------

library(VennDetail)
ven_input <- covered_protein_sets
ven <- venndetail(ven_input)
wide_ven <- result(ven, wide = TRUE)

library(ComplexUpset)
wide_ven <- wide_ven[, -7]
wide_ven <- na.omit(wide_ven)
rownames(wide_ven) <- NULL
wide_ven <- column_to_rownames(wide_ven, var = "Detail")

p_detected_upset <- ComplexUpset::upset(wide_ven, intersect = colnames(wide_ven),
  matrix = (intersection_matrix(geom = geom_point(shape = "square", size = 3),
  segment = geom_segment(linetype = "dotted", size = 0.2),outline_color = list(active = "#475272", inactive = "grey70")) +
  scale_color_manual(values = c("TRUE" = "#475272", "FALSE" = "#E8ECFBFF"), labels = c("TRUE" = "","FALSE" = ""),
  breaks = c("TRUE", "FALSE"),) + scale_y_discrete(position = "left") +
  annotate(geom = "text",label = "Look here →", x = "Comedy-Drama", y = "Drama", size = 5,hjust = 1)),
  base_annotations = list("Intersection size" = intersection_size(text = list(angle = 45,hjust = 0,vjust = 0.2,size = 3),
  mapping = aes(fill = "bars_color")) + scale_fill_manual(values = c("bars_color" = "#475272"), guide = "none")),
  set_sizes = (upset_set_size(geom = geom_bar(aes(x = group), width = 0.5,color = "#475272", fill = "#475272")) +
  theme(axis.ticks.x = element_line()) + geom_text(aes(label = ..count..),hjust = 0.2, size = 3, stat = "count") +
  annotate(geom = "text", label = "@",x = "Drama",y = 850,color = "white",size = 3) + expand_limits(y = 1100) +
  theme(axis.text.x = element_text(angle = 30))))

p_detected_upset_layout <- p_detected_upset + plot_layout(heights = c(1, 1)) + plot_layout(widths = c(3, 10))

ggsave(paste0(wdplot, "F2b_platform_overlap_upset.pdf"),pdf_text_grob(p_detected_upset_layout), width = 20, height = 9, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 01-06. Plot detection frequency in the five-platform overlap ----------------

det_freq_overlap_ot480 <- indi_prot(ms_ot480_uniquepep2[detected_overlap_5, ], cate = "ms", metadata)
det_freq_overlap_astral <- indi_prot(ms_astral_uniquepep2[detected_overlap_5, ], cate = "ms", metadata)
det_freq_overlap_explore <- indi_prot(lod3072 %>% dplyr::filter(Assay %in% detected_overlap_5), cate = "olink", metadata)
det_freq_overlap_ht <- indi_prot(lodht %>% dplyr::filter(Assay %in% detected_overlap_5), cate = "olink", metadata)

soma_visit_detected_overlap <- summarise_soma_visit_detection(soma_detection_master, detected_overlap_5)

det_freq_overlap_soma <- summarise_donor_detection_frequency(soma_visit_detected_overlap$protein, "n_detected_visits") %>%

  dplyr::group_by(SubjectID) %>% tidyr::complete(number = c(0,1, 2,3), fill = list(count = 0,percentage = 0)) %>%
  dplyr::ungroup()

det_freq_overlap_ot480$cate <- "MS-OT480"
det_freq_overlap_astral$cate <- "MS-Astral"
det_freq_overlap_explore$cate <- "Olink-Explore"
det_freq_overlap_ht$cate <- "Olink-HT"
det_freq_overlap_soma$cate <- "SomaScan"
det_freq_overlap_bind <- dplyr::bind_rows(det_freq_overlap_ot480, det_freq_overlap_astral,
  det_freq_overlap_explore,det_freq_overlap_ht, det_freq_overlap_soma)

# 01-07. Plot donor support in the five-platform overlap ----------------------

individual_detected_overlap_targets_df <- dplyr::bind_rows(lapply(platform_levels, function(platform_name) {
  tibble::tibble(platform = platform_name, protein = individual_detected_overlap_5)
}))

individual_detected_overlap_support_distribution_06 <- protein_visit_detection_df %>%
  dplyr::inner_join(individual_detected_overlap_targets_df,by = c("platform","protein")) %>%
  dplyr::group_by(platform,protein,SubjectID) %>%
  dplyr::summarise(n_detected_visits = sum(detected), individual_level_detection = n_detected_visits >= 2,
  .groups = "drop") %>%
  dplyr::group_by(platform, protein) %>%
    dplyr::summarise(n_donors_with_individual_detection = sum(individual_level_detection), .groups = "drop") %>%
  dplyr::count(platform, n_donors_with_individual_detection, name = "n_individual_detected_overlap_proteins") %>%
  dplyr::group_by(platform) %>% tidyr::complete(n_donors_with_individual_detection = 1:10,
    fill = list(n_individual_detected_overlap_proteins = 0L)) %>%
  dplyr::mutate(percentage = n_individual_detected_overlap_proteins / sum(n_individual_detected_overlap_proteins)) %>%

  dplyr::ungroup() %>% dplyr::mutate(platform = factor(platform, levels = platform_levels),
  n_donors_with_individual_detection = factor(n_donors_with_individual_detection,levels = 1:10)) %>%
  dplyr::arrange(platform, n_donors_with_individual_detection)

overlap_donor_percentage_plot_list <- lapply(platform_levels,plot_donor_support,
  input_df = individual_detected_overlap_support_distribution_06,
  y_column = "percentage", y_label = "Percentage of individual detected-overlap proteins"
)
overlap_donor_percentage_plot_list <- lapply(overlap_donor_percentage_plot_list,function(p) {
  p + ggplot2::scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0, 1), labels = function(x) x * 100,
    expand = ggplot2::expansion(mult = c(0, 0.03)))
})

for (i in seq_along(overlap_donor_percentage_plot_list)) {
  if (i != 3) {
    overlap_donor_percentage_plot_list[[i]] <- overlap_donor_percentage_plot_list[[i]] + labs(x = NULL)
  }
}
for (i in 2:length(overlap_donor_percentage_plot_list)) {
  overlap_donor_percentage_plot_list[[i]] <- overlap_donor_percentage_plot_list[[i]] + themeynone
}
overlap_donor_percentage_plot_list[[3]] <- overlap_donor_percentage_plot_list[[3]] +
  labs(x = "Number of individuals with individual-level detection")
grDevices::cairo_pdf(paste0(wdplot,"S1d_individual_detected_overlap_support_percentage.pdf"), width = 11 / 2.54, height = 5.5 / 2.54, family = "Helvetica")
grid::grid.draw(pdf_text_grob(overlap_donor_percentage_plot_list[[1]] + overlap_donor_percentage_plot_list[[2]] +
  overlap_donor_percentage_plot_list[[3]] + overlap_donor_percentage_plot_list[[4]] +
  overlap_donor_percentage_plot_list[[5]] + plot_layout(widths = c(1, 1, 1, 1,1))))
dev.off()
# 01-08. Plot detected-protein overlap ----------------------------------------

overlap_protein_sets <- detected_protein_sets

n_overlap_platforms <- length(overlap_protein_sets)
overlap_matrix <- matrix(0, nrow = n_overlap_platforms, ncol = n_overlap_platforms)
jaccard_matrix <- matrix(0,nrow = n_overlap_platforms,ncol = n_overlap_platforms)
overlap_protein_counts <- sapply(overlap_protein_sets,length)
overlap_platform_labels <- paste0(names(overlap_protein_sets)," (", overlap_protein_counts, ")")
rownames(overlap_matrix) <- overlap_platform_labels
colnames(overlap_matrix) <- overlap_platform_labels
rownames(jaccard_matrix) <- overlap_platform_labels
colnames(jaccard_matrix) <- overlap_platform_labels

for (row_index in 1:n_overlap_platforms) {
  for (column_index in 1:n_overlap_platforms) {
    overlap_count <- length(intersect(overlap_protein_sets[[row_index]], overlap_protein_sets[[column_index]]))
    overlap_matrix[row_index,column_index] <- overlap_count
    jaccard_matrix[row_index,
      column_index] <- jaccard_index(overlap_protein_sets[[row_index]], overlap_protein_sets[[column_index]])
  }
}
overlap_counts_df <- reshape2::melt(overlap_matrix, varnames = c("List1", "List2"), value.name = "Overlap")
overlap_jaccard_df <- reshape2::melt(jaccard_matrix, varnames = c("List1", "List2"), value.name = "Jaccard")
df <- merge(overlap_counts_df, overlap_jaccard_df,by = c("List1", "List2"))

p_detect_overlap_heatmap <- ggplot(df, aes(x = List1, y = List2)) + geom_tile(aes(fill = Jaccard), color = "white",size = 0.7) +
  geom_text(aes(label = Overlap),size = 2.5,color = "white",vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(name = "Jaccard index", low = "#F5F3B1", mid = "#C55548", high = "#1F1243", midpoint = 0.5,
  na.value = "#E6E6E6", guide = "colorbar") + labs(title = "", x = "",y = "") + themelegendtop + themetext30

ggsave(paste0(wdplot,"F2c_detected_overlap_heatmap.pdf"), pdf_text_grob(p_detect_overlap_heatmap), width = 11, height = 9, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 01-09. Summarise and plot protein missingness -------------------------------

summarise_soma_missing <- function(soma_df, protein_set = NULL) {
  if (!is.null(protein_set)) {
    soma_df <- soma_df %>% dplyr::filter(EntrezGeneSymbol %in% protein_set)
  }
  somamer_df <- soma_df %>% dplyr::filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    dplyr::group_by(SeqId, EntrezGeneSymbol) %>%
    dplyr::summarise(n_valid = sum(valid), n_detected = sum(detected),n_miss = sum(valid & !detected),
    detection_rate = dplyr::if_else(n_valid > 0,n_detected / n_valid, NA_real_),.groups = "drop")

  somamer_stats <- somamer_df %>% dplyr::group_by(EntrezGeneSymbol) %>%
    dplyr::summarise(n_somers = dplyr::n(), n_detected_somers = sum(n_detected > 0),
    n_miss_min = if (any(!is.na(n_miss))) {
    min(n_miss, na.rm = TRUE)
  } else {
    NA_real_
  },
    n_miss_max = if (any(!is.na(n_miss))) {
    max(n_miss, na.rm = TRUE)
  } else {
    NA_real_
  },
    n_miss_sd = if (sum(!is.na(n_miss)) > 1) {
    stats::sd(n_miss, na.rm = TRUE)
  } else {
    NA_real_
  }, .groups = "drop")

  protein_sample_df <- soma_df %>% dplyr::filter(!is.na(EntrezGeneSymbol),EntrezGeneSymbol != "") %>%
    dplyr::group_by(Sample, EntrezGeneSymbol) %>%
    dplyr::summarise(n_somers = dplyr::n_distinct(SeqId), n_valid_somers = sum(valid),
    n_detected_somers = sum(detected), protein_sample_valid = n_valid_somers > 0,
    protein_sample_detected = n_detected_somers > 0, .groups = "drop")

  protein_df <- protein_sample_df %>% dplyr::group_by(EntrezGeneSymbol) %>%
    dplyr::summarise(n_valid = sum(protein_sample_valid), n_detected = sum(protein_sample_detected),
    n_miss = sum(protein_sample_valid & !protein_sample_detected), n_somers = max(n_somers), .groups = "drop") %>%
    dplyr::mutate(detection_rate = dplyr::if_else(n_valid > 0,n_detected / n_valid, NA_real_),
    protein_detected_any = n_detected > 0) %>%
    dplyr::left_join(somamer_stats %>% dplyr::select(-n_somers), by = "EntrezGeneSymbol") %>%
    dplyr::mutate(prop_detected_somers = n_detected_somers / n_somers)

  list(somamer = somamer_df,protein = protein_df)
}
missingness_colours <- c( "MS-OT480"      = "#475272", "MS-Astral"     = "#68A39D", "Olink-Explore" = "#D9B381",
  "Olink-HT"      = "#ADC482", "SomaScan"      = "#B5252E"
)

build_ms_missingness <- function(matrix_data,platform_name) {
  data.frame( platform = platform_name, protein = rownames(matrix_data),
    n_valid = apply(matrix_data, 1, function(x) sum(!is.na(x))),
    n_detected = apply(matrix_data, 1, function(x) sum(!is.na(x) & x != 0)),
    n_miss = apply(matrix_data, 1,function(x) sum(!is.na(x) & x == 0)), n_somers = 1L
  ) %>%
    dplyr::mutate(detection_rate = dplyr::if_else(n_valid > 0, n_detected / n_valid, NA_real_))
}
build_olink_missingness <- function(lod_df, platform_name, protein_set = NULL) {
  if (!is.null(protein_set)) {
    lod_df <- lod_df %>% dplyr::filter(Assay %in% protein_set)
  }
  lod_df %>% dplyr::mutate(valid = !is.na(NPX) & !is.na(LOD),detected = valid & NPX >= LOD) %>% dplyr::group_by(Assay) %>%
    dplyr::summarise(n_valid = sum(valid), n_detected = sum(detected), n_miss = sum(valid & !detected),
    n_somers = 1L, .groups = "drop") %>%
    dplyr::rename(protein = Assay) %>%
    dplyr::mutate(detection_rate = dplyr::if_else(n_valid > 0, n_detected / n_valid, NA_real_), platform = platform_name)
}

miss_per_protein_ot480 <- build_ms_missingness(ms_ot480_uniquepep2[covered_protein_sets[[
  "MS-OT480"]],],"MS-OT480")
miss_per_protein_astral <- build_ms_missingness(ms_astral_uniquepep2[covered_protein_sets[[
  "MS-Astral"]],],"MS-Astral")
miss_per_protein_explore <- build_olink_missingness(lod3072,"Olink-Explore", covered_protein_sets[[
  "Olink-Explore"]])
miss_per_protein_ht <- build_olink_missingness(lodht, "Olink-HT", covered_protein_sets[["Olink-HT"]])

soma_missing_full <- summarise_soma_missing(soma_detection_master, covered_protein_sets[["SomaScan"]])

miss_per_protein_soma <- soma_missing_full$protein %>% dplyr::rename(protein = EntrezGeneSymbol) %>%
  dplyr::mutate(platform = "SomaScan")

miss_per_sample_all <- bind_rows(miss_per_protein_ot480,miss_per_protein_astral, miss_per_protein_explore,
  miss_per_protein_ht,miss_per_protein_soma) %>% dplyr::mutate(platform = factor(platform,levels = platform_levels))

missingness_y_max <- 30

plot_missingness <- function(df_sub, fill_col, y_max,ylab = "", title = "") {
  ggplot(df_sub,aes(x = platform, y = n_miss,fill = platform)) + geom_violin(trim = TRUE, alpha = 0.85, width = 0.9,color = NA) +
    scale_fill_manual(values = fill_col) +
    scale_y_continuous(limits = c(0,y_max), breaks = seq(0, y_max, by = 5)) + labs(x = NULL, y = ylab, title = title) +
    theme2_nobacknoerror + themelegendnone
}
p_missingness_ms_ot480 <- plot_missingness(dplyr::filter(miss_per_sample_all, platform == "MS-OT480"),
  missingness_colours,missingness_y_max, ylab = "Missing samples per protein (n)", title = "MS-OT480")
p_missingness_ms_astral <- plot_missingness(dplyr::filter(miss_per_sample_all, platform == "MS-Astral"),
  missingness_colours,missingness_y_max, title = "MS-Astral") + themeynone
p_missingness_olink_explore <- plot_missingness(dplyr::filter(miss_per_sample_all,platform == "Olink-Explore"),
  missingness_colours, missingness_y_max,title = "Olink-Explore") + themeynone
p_missingness_olink_ht <- plot_missingness(dplyr::filter(miss_per_sample_all, platform == "Olink-HT"),
  missingness_colours,missingness_y_max, title = "Olink-HT") + themeynone
p_missingness_soma <- plot_missingness(dplyr::filter(miss_per_sample_all, platform == "SomaScan"),
  missingness_colours,missingness_y_max,title = "SomaScan") + themeynone

p_missingness <- p_missingness_ms_ot480 + p_missingness_ms_astral + p_missingness_olink_explore +
  p_missingness_olink_ht + p_missingness_soma +
  plot_layout(nrow = 1,widths = c(1.2,1,1, 1, 1))

ggsave(paste0(wdplot, "F2f_missing_value_per_sample_violin.pdf"),pdf_text_grob(p_missingness), width = 12,height = 6, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 01-10. Plot missingness in the individual-detected overlap ------------------

missingness_overlap_proteins <- individual_detected_overlap_5

miss_ov_ot480 <- build_ms_missingness(ms_ot480_uniquepep2[missingness_overlap_proteins, ], "MS-OT480")
miss_ov_astral <- build_ms_missingness(ms_astral_uniquepep2[missingness_overlap_proteins, ],"MS-Astral")
miss_ov_explore <- build_olink_missingness(lod3072, "Olink-Explore",missingness_overlap_proteins)
miss_ov_ht <- build_olink_missingness(lodht, "Olink-HT",missingness_overlap_proteins)

soma_missing_overlap <- summarise_soma_missing(soma_detection_master, missingness_overlap_proteins)

miss_ov_soma <- soma_missing_overlap$protein %>% dplyr::rename(protein = EntrezGeneSymbol) %>%
  dplyr::mutate(platform = "SomaScan")

miss_ov_all <- bind_rows(miss_ov_ot480, miss_ov_astral,miss_ov_explore, miss_ov_ht,miss_ov_soma) %>%
  dplyr::mutate(platform = factor(platform,levels = platform_levels))

p_overlap_missingness_ms_ot480 <- plot_missingness(dplyr::filter(miss_ov_all,platform == "MS-OT480"),
  missingness_colours, missingness_y_max, ylab = "Missing samples per protein (n)", title = "MS-OT480")
p_overlap_missingness_ms_astral <- plot_missingness(dplyr::filter(miss_ov_all, platform == "MS-Astral"),
  missingness_colours,missingness_y_max,title = "MS-Astral") + themeynone
p_overlap_missingness_olink_explore <- plot_missingness(dplyr::filter(miss_ov_all,platform == "Olink-Explore"),
  missingness_colours, missingness_y_max, title = "Olink-Explore") + themeynone
p_overlap_missingness_olink_ht <- plot_missingness(dplyr::filter(miss_ov_all, platform == "Olink-HT"),
  missingness_colours, missingness_y_max, title = "Olink-HT") + themeynone
p_overlap_missingness_soma <- plot_missingness(dplyr::filter(miss_ov_all, platform == "SomaScan"),
  missingness_colours, missingness_y_max, title = "SomaScan") + themeynone

p_overlap_missingness <- p_overlap_missingness_ms_ot480 + p_overlap_missingness_ms_astral +
  p_overlap_missingness_olink_explore + p_overlap_missingness_olink_ht + p_overlap_missingness_soma +
  plot_layout(nrow = 1, widths = c(1.2, 1, 1,1,1))

ggsave(paste0(wdplot,"S1e_individual_detected_overlap_missing_violin.pdf"),pdf_text_grob(p_overlap_missingness),
  width = 12, height = 6, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
# 16. Plot platform-coloured detection-overlap Venn diagrams ------------------

detection_sets_file <- file.path(repo_dir, "result","Fig2","01-01_detection_classification_sets.rds")
detection_sets <- readRDS(detection_sets_file)
detection_platform_keys <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT","Somascan")
col_platform_internal <- c( "MS-OT480" = "#D9B381","MS-Astral" = "#B5252E",
  "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF","Overlap" = "#D1C1E1FF"
)
platform_display_labels <- c( "MS-OT480" = "MS-OT480", "MS-Astral" = "MS-Astral",
  "Olink-Explore" = "Olink-Explore", "Olink-HT" = "Olink-HT", "Somascan" = "SomaScan", "Overlap" = "Overlap"
)

detection_scope_keys <- c("detected", "individual_detected","broad_detected")
detection_scope_labels <- c( "detected" = "Detected proteins", "individual_detected" = "Individual detected proteins",
  "broad_detected" = "Broadly detected proteins"
)

p_detection_platform_colour_venn_list <- list()

for (scope_key in detection_scope_keys) {
  platform_sets <- detection_sets[[scope_key]][detection_platform_keys]
  names(platform_sets) <- unname(platform_display_labels[detection_platform_keys])
  p_detection_platform_colour_venn_list[[scope_key]] <- VennDiagram::venn.diagram(x = platform_sets,
    filename = NULL,col = "white", fill = unname(col_platform_internal[detection_platform_keys]),alpha = 0.6,
    cex = 0.8,cat.cex = 0.75, cat.col = "black",lwd = 0.5,cat.lwd = 0.5, margin = 0.1, cat.dist = rep(0.1, length(platform_sets)),
    main = unname(detection_scope_labels[scope_key]),main.cex = 0.8)
}
p_detection_platform_colour_venn <- do.call( gridExtra::arrangeGrob, c(p_detection_platform_colour_venn_list,ncol = 3)
)
detection_venn_pdf_file <- file.path(wdplot, "S1a_five_platform_detection_overlap_venn_platform_colours.pdf")
ggplot2::ggsave(detection_venn_pdf_file, pdf_text_grob(p_detection_platform_colour_venn),width = 27, height = 9,units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

saveRDS(list(protein_sets = detection_classification_sets,
  undetected = undetected_protein_sets,
  protein_counts = df_plot,
  detection_frequency = list(MS_OT480 = det_freq_ot480, MS_Astral = det_freq_astral, Olink_Explore = det_freq_explore, Olink_HT = det_freq_ht, SomaScan = det_freq_soma),
  individual_support = donor_support_distribution_03,
  overlap_individual_support = individual_detected_overlap_support_distribution_06,
  missing_per_sample = miss_per_sample_all,
  missing_overlap = miss_ov_all),
  paste0(wdresult, "detection_results.rds"))
