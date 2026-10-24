# Figure 5: complete analysis and manuscript figures
# Run from the repository root.
source(file.path("code", "setup.R"), encoding = "UTF-8")
wdresult <- paste0(file.path(repo_dir, "result", "Fig5"), "/")
wdplot <- paste0(file.path(repo_dir, "plot", "Fig5"), "/")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

scope_key <- "detected"

plot_prefix <- "05-"

load_registered_assays(wddata)
ot480 <- (diaori1)
astral <- (astori1)


rownames(olink3072)[rownames(olink3072) == "NTproBNP"] <- "NT-proBNP"
rownames(olink3072)[rownames(olink3072) == "WARS"] <- "WARS1"

detection_sets_file <- file.path(repo_dir, "result", "Fig2", "01-01_detection_classification_sets.rds")
detection_sets <- readRDS(detection_sets_file)
prolist_by_platform <- detection_sets[[scope_key]]
overlap5 <- Reduce(intersect, unname(prolist_by_platform))

prolist <- stats::setNames(rep(list(overlap5), length(prolist_by_platform)), names(prolist_by_platform))
scope_label <- switch( scope_key, detected = "Detected proteins", individual_detected = "Individual detected proteins",
  broad_detected = "Broadly detected proteins"
)
overlap_scope_label <- switch( scope_key, detected = "Five-platform shared Detected proteins",
  individual_detected = "Five-platform shared Individual detected proteins",
  broad_detected = "Five-platform shared Broadly detected proteins"
)
dflistna <- read_processed_matrices(wddata)
names(dflistna)[names(dflistna) == "Olink-3072"] <- "Olink-Explore"

metadata <- read_processed_metadata(wddata)
metadata <- metadata[1:30,]
metadata$somaid <- paste0("S",metadata$OriginalID)
metadata$Time <- dplyr::recode(metadata$Time, "visit1" = "Day0", "visit2" = "Day3", "visit3" = "Day5")
metadata$Time <- factor(metadata$Time, levels = c("Day0", "Day3", "Day5"))
metadata$SubjectID <- factor(metadata$SubjectID ,levels= unique(metadata$SubjectID ))

# 09. Compare positively correlated protein sets ------------------------------

catename <- paste0("Overlap of positively correlated proteins across platforms: ", scope_label)

poscorr_matrix <- list()

poscorr_matrix[["MS-OT480"]] <- msmatrixpro(pro = prolist[["MS-OT480"]], cate = "MS-OT480",
  data = dflistna[["MS-OT480"]], isna = "NA")

poscorr_matrix[["MS-Astral"]] <- msmatrixpro(pro = prolist[["MS-Astral"]], cate = "MS-Astral",
  data = dflistna[["MS-Astral"]], isna = "NA")

poscorr_matrix[["Olink-Explore"]] <- olinkmatrixpro(pro = prolist[["Olink-Explore"]], cate = "Olink-Explore",
  data = dflistna[["Olink-Explore"]], isna = FALSE, keepfilter = FALSE)

poscorr_matrix[["Olink-HT"]] <- olinkmatrixpro(pro = prolist[["Olink-HT"]], cate = "Olink-HT",
  data = dflistna[["Olink-HT"]], isna = FALSE, keepfilter = FALSE)

poscorr_matrix[["Somascan"]] <- somamatrixpro(pro = prolist[["Somascan"]], cate = "Somascan",
  data = dflistna[["Somascan"]], isna = FALSE, keepfilter = FALSE)

poscorr_corr <- procorr(poscorr_matrix)

poscorr_pro_list <- list()
poscorr_tested_pro_list <- list()
for(i in unique(poscorr_corr$cate)) {
  poscorr_pro <- poscorr_corr %>% dplyr::filter(cate == i, Correlation > 0.4, p_value < 0.05) %>%
    dplyr::distinct(Protein) %>% dplyr::pull(Protein)
  poscorr_pro <- unique(poscorr_pro)
  poscorr_pro_list[[i]] <- poscorr_pro

  poscorr_tested_pro <- poscorr_corr %>% dplyr::filter(cate == i, !is.na(Correlation)) %>% dplyr::distinct(Protein) %>%
    dplyr::pull(Protein) %>% unique()
  poscorr_tested_pro_list[[i]] <- poscorr_tested_pro
}
poscorr_inputlist <- list()
for(i in names(poscorr_pro_list)) {
  poscorr_inputlist[[paste0(i)]] <- poscorr_pro_list[[i]]
}
poscorr_n <- length(poscorr_inputlist)
poscorr_overlap_matrix <- matrix(0, nrow = poscorr_n, ncol = poscorr_n)
poscorr_jaccard_matrix <- matrix(0, nrow = poscorr_n, ncol = poscorr_n)
poscorr_p_matrix <- matrix(NA_real_, nrow = poscorr_n, ncol = poscorr_n)
poscorr_lens <- sapply(poscorr_inputlist, length)
poscorr_new_names <- paste0(names(poscorr_inputlist), " (", poscorr_lens, ")")
rownames(poscorr_overlap_matrix) <- poscorr_new_names
colnames(poscorr_overlap_matrix) <- poscorr_new_names
rownames(poscorr_jaccard_matrix) <- poscorr_new_names
colnames(poscorr_jaccard_matrix) <- poscorr_new_names
rownames(poscorr_p_matrix) <- poscorr_new_names
colnames(poscorr_p_matrix) <- poscorr_new_names

for(i in 1:poscorr_n) {
  for(j in 1:poscorr_n) {
    poscorr_overlap_count <- length(intersect(poscorr_inputlist[[i]], poscorr_inputlist[[j]]))
    poscorr_overlap_matrix[i, j] <- poscorr_overlap_count
    poscorr_jaccard_matrix[i, j] <- jaccard_index(poscorr_inputlist[[i]], poscorr_inputlist[[j]])
    if (i != j) {
      bg_ij <- intersect(poscorr_tested_pro_list[[names(poscorr_inputlist)[i]]],
        poscorr_tested_pro_list[[names(poscorr_inputlist)[j]]])
      a_ij <- intersect(poscorr_inputlist[[i]], bg_ij)
      b_ij <- intersect(poscorr_inputlist[[j]], bg_ij)
      x_ij <- length(intersect(a_ij, b_ij))
      if (length(bg_ij) > 0 && length(a_ij) > 0 && length(b_ij) > 0) {
        poscorr_p_matrix[i, j] <- phyper(q = x_ij - 1, m = length(a_ij), n = length(bg_ij) - length(a_ij),
          k = length(b_ij), lower.tail = FALSE)
      }
    }
  }
}
poscorr_p_adj_matrix <- matrix(NA_real_, nrow = poscorr_n, ncol = poscorr_n, dimnames = dimnames(poscorr_p_matrix))
poscorr_unique_pair_idx <- which(upper.tri(poscorr_p_matrix) & !is.na(poscorr_p_matrix))
poscorr_p_adj_matrix[poscorr_unique_pair_idx] <- p.adjust( poscorr_p_matrix[poscorr_unique_pair_idx], method = "BH"
)
poscorr_p_adj_matrix[lower.tri(poscorr_p_adj_matrix)] <-
  t(poscorr_p_adj_matrix)[lower.tri(poscorr_p_adj_matrix)]

poscorr_df_overlap <- reshape2::melt(poscorr_overlap_matrix, varnames = c("List1", "List2"), value.name = "Overlap")
poscorr_df_jaccard <- reshape2::melt(poscorr_jaccard_matrix, varnames = c("List1", "List2"), value.name = "Jaccard")
poscorr_df_p <- reshape2::melt(poscorr_p_matrix, varnames = c("List1", "List2"), value.name = "p_overlap")
poscorr_df_p_adj <- reshape2::melt(poscorr_p_adj_matrix, varnames = c("List1", "List2"), value.name = "p_adj_overlap")
poscorr_df <- merge(poscorr_df_overlap, poscorr_df_jaccard, by = c("List1", "List2"))
poscorr_df <- poscorr_df %>% dplyr::left_join(poscorr_df_p, by = c("List1", "List2")) %>%
  dplyr::left_join(poscorr_df_p_adj, by = c("List1", "List2")) %>% dplyr::mutate(

  overlap_sig = dplyr::case_when( is.na(p_adj_overlap) ~ "", p_adj_overlap < 0.05 ~ "*", TRUE ~ ""
), overlap_label = ifelse(overlap_sig == "", as.character(Overlap), paste0(Overlap, "\n", overlap_sig))
)

poscorr_jaccard_summary <- as.data.frame(which(upper.tri(poscorr_overlap_matrix), arr.ind = TRUE)) %>%
  dplyr::transmute( correlation_comparison_1 = names(poscorr_inputlist)[row],
  correlation_comparison_2 = names(poscorr_inputlist)[col], n_positive_corr_proteins_1 = as.integer(poscorr_lens[row]),
  n_positive_corr_proteins_2 = as.integer(poscorr_lens[col]),
  n_overlapping_positive_corr_proteins = as.integer(poscorr_overlap_matrix[cbind(row, col)]),
  jaccard_index = as.numeric(poscorr_jaccard_matrix[cbind(row, col)]),
  raw_p_overlap = as.numeric(poscorr_p_matrix[cbind(row, col)]),
  fdr_bh_overlap = as.numeric(poscorr_p_adj_matrix[cbind(row, col)]),
  fdr_bh_lt_0_05 = ifelse(fdr_bh_overlap < 0.05, "*", "")
)

p_poscorr_heatmap <- ggplot(poscorr_df, aes(x = List1, y = List2)) +
  geom_tile(aes(fill = Jaccard), color = "white", size = 0.7) +
  geom_text(aes(label = overlap_label), size = 3.5, color = "white", vjust = 0.5, hjust = 0.5, lineheight = 0.8) +
  scale_size_continuous(range = c(7, 19)) + scale_color_gradient(low = "white", high = "steelblue") +
  scale_fill_gradient2(name = "Jaccard index", low = "#F5F3B1", mid = "#C55548", high = "#1F1243", midpoint = 0.5,
    na.value = "#E6E6E6",
  guide = guide_colorbar(title.position = "left",title.hjust = 0.5,barwidth = unit(4, "cm"),barheight = unit(0.45, "cm"),
  ticks = TRUE  ) ) + labs(title = catename, x = "", y = "") + themetext30+themelegendtop

ggsave(paste0(wdplot, "S6_positive_corr_overlap_heatmap.pdf"), pdf_text_grob(p_poscorr_heatmap), width = 28, height = 13, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

# 05.11. Shared-protein correlation matrices and undetected Olink comparisons ---

# 01. Align input matrices and shared proteins --------------------------------

input_files <- c( detection = file.path(repo_dir, "result", "Fig2",
  "01-01_detection_classification_sets.rds"), ms_ot480 = file.path(wddata, "processed", "MS-OT480_qualified_log2.rds"),
  ms_astral = file.path(wddata, "processed", "MS-Astral_qualified_log2.rds"),
  metadata = file.path(wddata, "processed", "metadata.rds")
)
detection_sets <- readRDS(input_files[["detection"]])
platforms <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
shared_proteins <- Reduce(intersect, detection_sets[["detected"]][platforms])
metadata <- readRDS(input_files[["metadata"]])
metadata <- metadata[seq_len(30L), , drop = FALSE]
sample_ids <- as.character(metadata$OriginalID)
time_labels <- dplyr::recode(as.character(metadata$Time), visit1 = "Day0", visit2 = "Day3", visit3 = "Day5")
sample_info <- data.frame(sample = sample_ids, subject = as.character(metadata$SubjectID),
  time = time_labels, stringsAsFactors = FALSE)
load_registered_assays(wddata)
# 02. Prepare the platform expression matrices --------------------------------

prepare_ms <- function(platform, input_file, protein_scope = shared_proteins) {
  raw_table <- readRDS(input_file)
  required <- c("Genes", sample_ids)
  raw_table <- raw_table[raw_table$Genes %in% protein_scope, , drop = FALSE]
  absent <- setdiff(protein_scope, raw_table$Genes)
  duplicate_genes <- unique(raw_table$Genes[duplicated(raw_table$Genes)])
  raw <- as.matrix(raw_table[match(protein_scope, raw_table$Genes), sample_ids, drop = FALSE])
  rownames(raw) <- protein_scope
  colnames(raw) <- sample_ids
  raw <- apply(raw, 2, function(x) suppressWarnings(as.numeric(x)))
  raw <- matrix(raw, nrow = length(protein_scope), ncol = length(sample_ids), dimnames = list(protein_scope, sample_ids))
  observed_n <- rowSums(is.finite(raw))
  log2_raw <- raw
  list(matrix = log2_raw, log2_raw = log2_raw, audit = data.frame(platform = platform, protein = rownames(raw),
    observed_n = observed_n, remaining_na = rowSums(!is.finite(log2_raw))))
}
prepare_olink <- function(platform, lod_table, protein_scope = shared_proteins,
  replace_below_lod = TRUE, use_raw_npx = FALSE) {
  required <- c("Assay", "SampleID", "NPX", "LOD")
  long <- lod_table %>% dplyr::filter(Assay %in% protein_scope, SampleID %in% sample_ids) %>%
    dplyr::select(dplyr::all_of(required))
  long <- long %>% dplyr::group_by(Assay) %>%
    dplyr::mutate(minimum_npx = if (any(is.finite(NPX))) {
    min(NPX[is.finite(NPX)])
  } else {
    NA_real_
  },

    below_lod = is.finite(NPX) & is.finite(LOD) & NPX < LOD, missing_npx = !is.finite(NPX), missing_lod = !is.finite(LOD),
    value = if (use_raw_npx) {
    NPX
  } else {
    dplyr::if_else((replace_below_lod & below_lod) | missing_npx, minimum_npx, NPX)
  }) %>% dplyr::ungroup()
  expected_n <- length(unique(long$Assay)) * length(sample_ids)
  wide <- long %>% dplyr::select(dplyr::all_of(c("Assay", "SampleID", "value"))) %>%
    tidyr::pivot_wider(names_from = "SampleID", values_from = "value")
  mat <- as.matrix(wide[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- wide$Assay
  audit <- long %>% dplyr::group_by(Assay) %>% dplyr::summarise(below_lod_n = sum(below_lod),
    missing_npx_n = sum(missing_npx), missing_lod_n = sum(missing_lod), .groups = "drop") %>%
    dplyr::transmute(platform = platform, protein = Assay, below_lod_n, missing_npx_n, missing_lod_n)
  list(matrix = mat, audit = audit)
}
ms_ot <- prepare_ms("MS-OT480", input_files[["ms_ot480"]])
ms_as <- prepare_ms("MS-Astral", input_files[["ms_astral"]])
ol_exp <- prepare_olink("Olink-Explore", lod3072)
ol_ht <- prepare_olink("Olink-HT", lodht)

map <- somer_protein %>%
  dplyr::transmute(protein = trimws(as.character(EntrezGeneSymbol)), somer = as.character(SeqId)) %>%
  dplyr::filter(!is.na(protein), nzchar(protein), !is.na(somer), nzchar(somer), protein %in% shared_proteins) %>%
  dplyr::distinct()
map <- map %>% dplyr::group_by(somer) %>% dplyr::filter(dplyr::n_distinct(protein) == 1L) %>% dplyr::ungroup()
soma_raw <- as.matrix(soma)
map <- map %>% dplyr::filter(somer %in% rownames(soma_raw))
soma_matrix <- soma_raw[map$somer, sample_ids, drop = FALSE]
storage.mode(soma_matrix) <- "numeric"
soma_matrix[!is.finite(soma_matrix) | soma_matrix <= 0] <- NA_real_
soma_matrix <- log2(soma_matrix)
matrices <- list(`MS-OT480` = ms_ot$matrix, `MS-Astral` = ms_as$matrix,
  `Olink-Explore` = ol_exp$matrix, `Olink-HT` = ol_ht$matrix, Somascan = soma_matrix)
# 03. Calculate paired Spearman correlations ----------------------------------

safe_rho <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
    length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}
weighted_spearman <- function(x, y, protein) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]

  y <- y[keep]

  protein <- protein[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  weight <- 1 / as.numeric(table(protein)[protein])
  weight <- weight / sum(weight)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  mx <- sum(weight * rx)

  my <- sum(weight * ry)
  vx <- sum(weight * (rx - mx)^2)

  vy <- sum(weight * (ry - my)^2)
  if (vx <= 0 || vy <= 0) {
    return(NA_real_)
  }
  sum(weight * (rx - mx) * (ry - my)) / sqrt(vx * vy)
}
pair_map <- function(a, b) {
  if (b == "Somascan" && !soma_protein_level) {
    map %>% dplyr::filter(protein %in% rownames(matrices[[a]])) %>% dplyr::transmute(row_a = protein, row_b = somer,
      protein = protein)
  } else {
    common <- intersect(rownames(matrices[[a]]), rownames(matrices[[b]]))
    data.frame(row_a = common, row_b = common, protein = common)
  }
}
soma_protein_level <- FALSE
sample_rho <- function(a, b, rows, sample_a, sample_b) {
  x <- matrices[[a]][rows$row_a, sample_a]
  y <- matrices[[b]][rows$row_b, sample_b]
  if (b == "Somascan") {
    weighted_spearman(x, y, rows$protein)
  } else {
    safe_rho(x, y)
  }
}
sample_n <- function(a, b, rows, sample_a, sample_b) {
  sum(is.finite(matrices[[a]][rows$row_a, sample_a]) &
    is.finite(matrices[[b]][rows$row_b, sample_b]))
}
# 04. Prepare correlation panels and output paths -----------------------------

corr_level_cols <- c("Per-protein cross-platform correlation" = "#4E79A7FF",
  "Same-sample cross-platform correlation" = "#F28E2B")
corr_matrix_display <- list( upper_triangle = "Per-protein cross-platform correlation",
  lower_triangle = "Same-sample cross-platform correlation", histogram_bins = 16L, colours = corr_level_cols,
  median_line = "dashed_vertical_line", protein_counts = "title_only_five_platform_shared_membership")
display_platforms <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "SomaScan")
pair_label <- function(a, b) paste(a, "vs", b)
header_plot <- function(label, vertical = FALSE) {
  ggplot() + annotate("text", x = .5, y = .5, label = label,
    angle = if (vertical) {
    90
  } else {
    0
  }, fontface = "bold", size = 4) + xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(panel.background = element_rect(fill = "#D9D9D9"), plot.margin = margin(1, 1, 1, 1, "pt"))
}
empty_plot <- function(diagonal = FALSE) {
  ggplot() + theme_void() + theme(
    panel.background = element_rect(fill = if (diagonal) {
    "#F2F2F2"
  } else {
    "white"
  },
    colour = if (diagonal) {
    "#D9D9D9"
  } else {
    NA
  }), plot.margin = margin(1, 1, 1, 1, "pt"))
}
cell_plot <- function(inner) {
  background <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "#D9D9D9"), plot.margin = margin(1, 1, 1, 1, "pt"))
  background + patchwork::inset_element(inner, left = .12, bottom = .12, right = .84, top = .84)
}
correlation_histogram_cell <- function(df, fill, show_median_line = TRUE) {
  if (!nrow(df)) {
    return(empty_plot(TRUE))
  }
  median_rho <- median(df$value)
  inner <- ggplot(df, aes(x = value, y = after_stat(density))) + geom_histogram(bins = corr_matrix_display$histogram_bins,
    fill = fill, colour = "white", linewidth = .15, boundary = -1) +
    (if (show_median_line) {
    geom_vline(xintercept = median_rho, linetype = "dashed", colour = "grey20", linewidth = .4, show.legend = FALSE)
  }) + annotate("text", x = .95, y = Inf, label = sprintf("%.2f", median_rho),
    hjust = 1, vjust = 1.2, size = 2.7, colour = "grey20") +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) + scale_y_continuous(expand = expansion(mult = c(0, .18))) +
    labs(x = "Spearman correlation coefficient", y = "Density") + theme_classic(base_size = 9) +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7), plot.margin = margin(0, 0, 0, 0, "pt"))
  cell_plot(inner)
}
upper_cell <- function(pair, show_median_line = TRUE) {
  df <- upper_df %>% dplyr::filter(cate == pair)
  if (!nrow(df)) {
    return(empty_plot(TRUE))
  }
  correlation_histogram_cell(df, fill = unname(corr_level_cols[["Per-protein cross-platform correlation"]]),
    show_median_line = show_median_line)
}
lower_cell <- function(pair) {
  df <- lower_df %>% dplyr::filter(cate == pair)
  if (!nrow(df)) {
    return(empty_plot(TRUE))
  }
  correlation_histogram_cell(df, fill = unname(corr_level_cols[["Same-sample cross-platform correlation"]]))
}
legend_df <- data.frame(type = factor(names(corr_level_cols), levels = names(corr_level_cols)), x = 1:2, y = 1)
p_corr_legend <- ggplot(legend_df, aes(x = x, y = y, fill = type)) + geom_point(shape = 22, size = 4) +
  scale_fill_manual(values = corr_level_cols, name = "Correlation level") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(shape = 22, size = 4))) +
  theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 10))

outputs <- c( zero_matrix = file.path(wdplot, "F5a_ms_na_zero_cross_platform_corr_matrix.pdf"),
  zero_legend = file.path(wdplot, "F5a_ms_na_zero_cross_platform_corr_matrix_legend.pdf"),
  zero_evidence = file.path(wdresult, "05.11-05_ms_na_zero_correlation_evidence.rds"),
  zero_audit = file.path(wdresult, "05.11-05_ms_na_zero_audit.rds"),

  undetected_violin = file.path(wdplot, "F5b_shared_olink_undetected_cross_platform_correlation_violin.pdf"),
  undetected_evidence = file.path(wdresult, "05.11-07_shared_olink_undetected_cross_platform_correlation.rds"),

  individual_matrix = file.path(wdplot, "S7a_individual_detected_cross_platform_corr_matrix.pdf"),
  individual_legend = file.path(wdplot, "S7a_individual_detected_cross_platform_corr_matrix_legend.pdf"),
  individual_evidence = file.path(wdresult, "05.11-09_individual_detected_correlation_evidence.rds"),
  individual_audit = file.path(wdresult, "05.11-09_individual_detected_audit.rds"),

  broad_matrix = file.path(wdplot, "S7b_broad_detected_cross_platform_corr_matrix.pdf"),
  broad_legend = file.path(wdplot, "S7b_broad_detected_cross_platform_corr_matrix_legend.pdf"),
  broad_evidence = file.path(wdresult, "05.11-10_broad_detected_correlation_evidence.rds"),
  broad_audit = file.path(wdresult, "05.11-10_broad_detected_audit.rds")
)

dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)
# 05. Compare shared detected proteins with MS missing values set to zero -----

matrices[["MS-OT480"]] <- ms_ot$log2_raw
matrices[["MS-Astral"]] <- ms_as$log2_raw
for (platform in c("MS-OT480", "MS-Astral")) {
  matrices[[platform]][!is.finite(matrices[[platform]])] <- 0
}
soma_protein_ids <- unique(map$protein)
soma_protein_matrix <- vapply(soma_protein_ids, function(protein_id) {
  assay_rows <- which(map$protein == protein_id)
  apply(soma_matrix[assay_rows, , drop = FALSE], 2, function(values) {
    values <- values[is.finite(values)]
    if (length(values)) {
      stats::median(values)
    } else {
      NA_real_
    }
  })
}, numeric(length(sample_ids)))
soma_protein_matrix <- t(soma_protein_matrix)
dimnames(soma_protein_matrix) <- list(soma_protein_ids, sample_ids)
matrices[["Somascan"]] <- soma_protein_matrix
soma_protein_level <- TRUE

same_sample_list <- list()

within_individual_list <- list()

protein_list <- list()
k_sample <- 1L

k_individual <- 1L

k_protein <- 1L
for (i in seq_len(length(platforms) - 1L)) {
  for (j in seq.int(i + 1L, length(platforms))) {
    a <- platforms[i]

    b <- platforms[j]
    pair <- paste(a, "vs", if (b == "Somascan") {
      "SomaScan"
    } else {
      b
    })
    rows <- pair_map(a, b)
    for (sample_id in sample_ids) {
      same_sample_list[[k_sample]] <- data.frame( cate = pair, sample = sample_id,
        Correlation = sample_rho(a, b, rows, sample_id, sample_id),
        n_assay_pairs = sample_n(a, b, rows, sample_id, sample_id))
      k_sample <- k_sample + 1L
    }
    for (subject_id in unique(sample_info$subject)) {
      subject_samples <- sample_info[sample_info$subject == subject_id, ]
      for (time_a in c("Day0", "Day3", "Day5")) {
        for (time_b in setdiff(c("Day0", "Day3", "Day5"), time_a)) {
          sample_a <- subject_samples$sample[subject_samples$time == time_a]
          sample_b <- subject_samples$sample[subject_samples$time == time_b]
          within_individual_list[[k_individual]] <- data.frame(
            cate = pair, subject = subject_id, rep_a = time_a, rep_b = time_b,
            corr = sample_rho(a, b, rows, sample_a, sample_b), n_assay_pairs = sample_n(a, b, rows, sample_a, sample_b))
          k_individual <- k_individual + 1L
        }
      }
    }
    for (r in seq_len(nrow(rows))) {
      x <- matrices[[a]][rows$row_a[r], sample_ids]
      y <- matrices[[b]][rows$row_b[r], sample_ids]
      protein_list[[k_protein]] <- data.frame( cate = pair, Protein = rows$protein[r],
        Somer = if (b == "Somascan" && !soma_protein_level) {
        rows$row_b[r]
      } else {
        NA_character_
      }, Correlation = safe_rho(x, y), n_paired_samples = sum(is.finite(x) & is.finite(y)))
      k_protein <- k_protein + 1L
    }
  }
}
per_sample_corr_df <- dplyr::bind_rows(same_sample_list)
per_individual_corr_df <- dplyr::bind_rows(within_individual_list)
per_protein_corr_df <- dplyr::bind_rows(protein_list)

lower_df <- per_sample_corr_df %>% dplyr::transmute(cate, value = Correlation) %>% dplyr::filter(is.finite(value))
upper_df <- per_protein_corr_df %>%
  dplyr::transmute(cate, value = Correlation) %>% dplyr::filter(is.finite(value))
plot_list <- c(list(empty_plot()), lapply(display_platforms, header_plot))
for (i in seq_along(display_platforms)) {
  plot_list[[length(plot_list) + 1L]] <- header_plot(display_platforms[i], TRUE)
  for (j in seq_along(display_platforms)) {
    if (i == j) {
      cell <- empty_plot(TRUE)
    } else if (i > j) {
      cell <- lower_cell(pair_label(display_platforms[j], display_platforms[i]))
    } else {
      cell <- upper_cell(pair_label(display_platforms[i], display_platforms[j]), show_median_line = TRUE)
    }
    plot_list[[length(plot_list) + 1L]] <- cell
  }
}
p_corr_matrix <- wrap_plots(plotlist = plot_list, ncol = 6, byrow = TRUE) +
  plot_layout(widths = c(.28, rep(1, 5)), heights = c(.28, rep(1, 5))) +
  plot_annotation(title = sprintf("Five-platform shared Detected proteins (N = %d); MS NA = 0",
  length(unique(shared_proteins))), subtitle = "SomaScan: median log2 RFU across SOMAmer per protein and sample")

saveRDS(list(ms = dplyr::bind_rows( ms_ot$audit %>% dplyr::mutate(zero_filled_n = rowSums(!is.finite(ms_ot$log2_raw))),
  ms_as$audit %>% dplyr::mutate(zero_filled_n = rowSums(!is.finite(ms_as$log2_raw)))),
  olink = dplyr::bind_rows(ol_exp$audit, ol_ht$audit)), outputs[["zero_audit"]])
saveRDS(list(same_sample = per_sample_corr_df, within_individual = per_individual_corr_df,
  per_protein = per_protein_corr_df), outputs[["zero_evidence"]])
ggsave(outputs[["zero_matrix"]], pdf_text_grob(p_corr_matrix), width = 24, height = 17, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
ggsave(outputs[["zero_legend"]], pdf_text_grob(p_corr_legend), width = 24, height = 18, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")

# 07. Compare undetected Olink targets with other-platform detected proteins ---

focal_olink_platforms <- c("Olink-Explore", "Olink-HT")
olink_tables <- list("Olink-Explore" = lod3072, "Olink-HT" = lodht)
undetected_by_olink <- lapply(focal_olink_platforms, function(platform) {
  predefined <- unique(as.character(olink_tables[[platform]]$Assay))
  setdiff(predefined[!is.na(predefined) & nzchar(predefined)], as.character(detection_sets[["detected"]][[platform]]))
})
names(undetected_by_olink) <- focal_olink_platforms
undetected_union <- unique(unlist(undetected_by_olink, use.names = FALSE))
focal_undetected_matrices <- list()
for (platform in focal_olink_platforms) {
  target_ids <- undetected_by_olink[[platform]]
  if (length(target_ids)) {
    focal_undetected_matrices[[platform]] <- prepare_olink( platform, olink_tables[[platform]], protein_scope = target_ids,
      use_raw_npx = TRUE)$matrix
  }
}
external_undetected_matrices <- list()
for (platform in setdiff(platforms, "Somascan")) {
  target_ids <- intersect(undetected_union,
    if (platform %in% focal_olink_platforms) {
    unique(as.character(olink_tables[[platform]]$Assay))
  } else {
    undetected_union
  })
  if (platform %in% focal_olink_platforms) {
    target_ids <- intersect(target_ids, unique(as.character(olink_tables[[platform]]$Assay)))
    if (length(target_ids)) {
      external_undetected_matrices[[platform]] <-
        prepare_olink(platform, olink_tables[[platform]], protein_scope = target_ids, use_raw_npx = TRUE)$matrix
    }
  } else {
    raw_table <- readRDS(input_files[[if (platform == "MS-OT480") {
      "ms_ot480"
    } else {
      "ms_astral"
    }]])
    external_targets <- intersect(target_ids, as.character(detection_sets[["detected"]][[platform]]))
    raw_table <- raw_table[raw_table$Genes %in% external_targets, , drop = FALSE]
    duplicate_targets <- unique(raw_table$Genes[duplicated(raw_table$Genes)])
    if (nrow(raw_table)) {
      ext <- as.matrix(raw_table[, sample_ids, drop = FALSE])
      storage.mode(ext) <- "numeric"
      rownames(ext) <- raw_table$Genes
      ext[!is.finite(ext)] <- 0
      external_undetected_matrices[[platform]] <- ext
    }
  }
}
map_all_somers <- somer_protein %>% dplyr::transmute(protein = trimws(as.character(EntrezGeneSymbol)),
  somer = as.character(SeqId)) %>% dplyr::filter(!is.na(protein), nzchar(protein), !is.na(somer), nzchar(somer)) %>%
  dplyr::distinct()
map_all_somers <- map_all_somers %>% dplyr::group_by(somer) %>%
  dplyr::filter(dplyr::n_distinct(protein) == 1L) %>% dplyr::ungroup()
soma_undetected_rows <- map_all_somers %>%
  dplyr::filter(somer %in% rownames(soma_raw), protein %in% undetected_union,
  protein %in% as.character(detection_sets[["detected"]][["Somascan"]]))
soma_full_numeric <- soma_raw
storage.mode(soma_full_numeric) <- "numeric"
soma_full_numeric[!is.finite(soma_full_numeric) | soma_full_numeric <= 0] <- NA_real_
soma_full_numeric <- log2(soma_full_numeric)

undetected_pair_records <- list()

undetected_pair_index <- 1L
undetected_coverage_records <- list()

undetected_coverage_index <- 1L
for (focal in focal_olink_platforms) {
  if (is.null(focal_undetected_matrices[[focal]])) {
    next
  }
  focal_matrix <- focal_undetected_matrices[[focal]]
  for (external in setdiff(platforms, focal)) {
    external_detected <- as.character(detection_sets[["detected"]][[external]])
    available_external <- if (external == "Somascan") {
      unique(soma_undetected_rows$protein)
    } else {
      rownames(external_undetected_matrices[[external]])
    }
    overlap <- intersect(undetected_by_olink[[focal]], intersect(external_detected, available_external))
    undetected_coverage_records[[undetected_coverage_index]] <- data.frame(
      focal_olink = focal, external_platform = external, n_shared_undetected_focal_proteins = length(overlap),
      stringsAsFactors = FALSE)
    undetected_coverage_index <- undetected_coverage_index + 1L
    for (protein_id in overlap) {
      if (external == "Somascan") {
        somer_ids <- soma_undetected_rows$somer[soma_undetected_rows$protein == protein_id]
        for (somer_id in somer_ids) {
          undetected_pair_records[[undetected_pair_index]] <- data.frame(
            focal_olink = focal, external_platform = external, protein = protein_id, somamer_id = somer_id,
            rho = safe_rho(focal_matrix[protein_id, sample_ids], soma_full_numeric[somer_id, sample_ids]),
            stringsAsFactors = FALSE)
          undetected_pair_index <- undetected_pair_index + 1L
        }
      } else {
        ext_matrix <- external_undetected_matrices[[external]]
        if (!is.null(ext_matrix) && protein_id %in% rownames(ext_matrix)) {
          undetected_pair_records[[undetected_pair_index]] <- data.frame(
            focal_olink = focal, external_platform = external, protein = protein_id, somamer_id = NA_character_,
            rho = safe_rho(focal_matrix[protein_id, sample_ids], ext_matrix[protein_id, sample_ids]),
            stringsAsFactors = FALSE)
          undetected_pair_index <- undetected_pair_index + 1L
        }
      }
    }
  }
}
undetected_pair_df <- dplyr::bind_rows(undetected_pair_records)
undetected_coverage_df <- dplyr::bind_rows(undetected_coverage_records)
undetected_protein_rho_df <- undetected_pair_df %>% dplyr::group_by(focal_olink, external_platform, protein) %>%
  dplyr::summarise(rho = if (any(is.finite(rho))) {
  median(rho[is.finite(rho)])
} else {
  NA_real_
}, n_somamers = sum(!is.na(somamer_id)), n_valid_somamer_rho = sum(!is.na(somamer_id) & is.finite(rho)),
  .groups = "drop") %>% dplyr::filter(is.finite(rho)) %>% dplyr::mutate(
  focal_olink = factor(focal_olink, levels = focal_olink_platforms),
  external_platform = factor(external_platform, levels = platforms))
undetected_median_df <- undetected_protein_rho_df %>% dplyr::group_by(focal_olink, external_platform) %>%
  dplyr::summarise(median_rho = median(rho), n_valid_proteins = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(undetected_coverage_df,
  by = c("focal_olink" = "focal_olink", "external_platform" = "external_platform")) %>%
  dplyr::mutate(label = sprintf("%.2f", median_rho), y = 1.08)
platform_colours <- c("MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
  "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF")
undetected_fill_cols <- c(platform_colours, "Somascan" = "#69B190FF")
undetected_plot <- ggplot2::ggplot(undetected_protein_rho_df,
  ggplot2::aes(x = external_platform, y = rho, fill = external_platform)) +
  ggplot2::geom_violin(trim = FALSE, colour = "transparent", alpha = .85, width = .82) +
  ggplot2::geom_boxplot(width = .18, fill = "white", colour = "white", outlier.shape = NA, alpha = .85) +
  ggplot2::geom_text(data = undetected_median_df, ggplot2::aes(x = external_platform, y = y, label = label),
  inherit.aes = FALSE, size = 3, colour = "black") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::scale_fill_manual(values = undetected_fill_cols, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = c("Somascan" = "SomaScan")) +
  ggplot2::coord_cartesian(ylim = c(-1, 1.15), clip = "off") +
  ggplot2::facet_wrap(~ focal_olink, nrow = 1, scales = "free_x") +
  ggplot2::labs(title = sprintf(paste0("Olink Undetected proteins shared with other platforms\n",
  "(N = %d unique proteins across comparisons)"), dplyr::n_distinct(undetected_pair_df$protein)),

  x = NULL, y = "Spearman correlation coefficient") + ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
  plot.margin = ggplot2::margin(t = 14, r = 8, b = 8, l = 8, unit = "pt"))
ggplot2::ggsave(outputs[["undetected_violin"]], pdf_text_grob(undetected_plot), width = 22, height = 7, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
saveRDS(list(method = list( focal_olink = "raw_NPX_including_below_LOD_and_NA; no_imputation",
  external_olink = "raw_NPX_including_below_LOD_and_NA; no_imputation",
  external_ms = "log2_positive_intensity_then_NA_to_zero",
  somascan = "per_SOMAmer_rho_then_median_per_protein_and_platform"), somamer_level = undetected_pair_df,
  protein_level = undetected_protein_rho_df, overlap_coverage = undetected_coverage_df,
  platform_medians = undetected_median_df), outputs[["undetected_evidence"]])

run_zero_scope_matrix <- function(scope_key, scope_label, plot_title, output_prefix) {
  per_platform_sets <- detection_sets[[scope_key]]
  scope_proteins <- Reduce(intersect, per_platform_sets[platforms])
  ms_ot_scope <- prepare_ms("MS-OT480", input_files[["ms_ot480"]], protein_scope = scope_proteins)
  ms_as_scope <- prepare_ms("MS-Astral", input_files[["ms_astral"]], protein_scope = scope_proteins)
  ms_ot_zero <- ms_ot_scope$log2_raw
  ms_as_zero <- ms_as_scope$log2_raw
  ms_ot_zero[!is.finite(ms_ot_zero)] <- 0
  ms_as_zero[!is.finite(ms_as_zero)] <- 0
  ol_exp_scope <- prepare_olink("Olink-Explore", lod3072, protein_scope = scope_proteins)
  ol_ht_scope <- prepare_olink("Olink-HT", lodht, protein_scope = scope_proteins)

  scope_map <- somer_protein %>% dplyr::transmute(protein = trimws(as.character(EntrezGeneSymbol)),
    somer = as.character(SeqId)) %>%
    dplyr::filter(!is.na(protein), nzchar(protein), protein %in% scope_proteins,
    !is.na(somer), nzchar(somer), somer %in% rownames(soma_raw)) %>% dplyr::distinct()
  scope_map <- scope_map %>% dplyr::group_by(somer) %>%
    dplyr::filter(dplyr::n_distinct(protein) == 1L) %>% dplyr::ungroup()
  soma_scope_numeric <- soma_raw[scope_map$somer, sample_ids, drop = FALSE]
  storage.mode(soma_scope_numeric) <- "numeric"
  soma_scope_numeric[!is.finite(soma_scope_numeric) | soma_scope_numeric <= 0] <- NA_real_
  soma_scope_numeric <- log2(soma_scope_numeric)
  soma_scope_proteins <- unique(scope_map$protein)
  soma_scope_matrix <- matrix(NA_real_, nrow = length(soma_scope_proteins),
    ncol = length(sample_ids), dimnames = list(soma_scope_proteins, sample_ids))
  for (protein_id in soma_scope_proteins) {
    rows <- which(scope_map$protein == protein_id)
    soma_scope_matrix[protein_id, ] <- apply(
      soma_scope_numeric[rows, , drop = FALSE], 2, function(values) {
      values <- values[is.finite(values)]
      if (length(values)) {
        stats::median(values)
      } else {
        NA_real_
      }
    })
  }
  matrices <<- list(`MS-OT480` = ms_ot_zero, `MS-Astral` = ms_as_zero,
    `Olink-Explore` = ol_exp_scope$matrix, `Olink-HT` = ol_ht_scope$matrix, Somascan = soma_scope_matrix)
  soma_protein_level <<- TRUE

  same_sample_records <- list()

  individual_records <- list()

  protein_records <- list()
  idx_sample <- 1L

  idx_individual <- 1L

  idx_protein <- 1L
  for (i in seq_len(length(platforms) - 1L)) {
    for (j in seq.int(i + 1L, length(platforms))) {
      a <- platforms[i]

      b <- platforms[j]
      pair <- paste(a, "vs", if (b == "Somascan") {
        "SomaScan"
      } else {
        b
      })
      rows <- pair_map(a, b)
      for (sample_id in sample_ids) {
        same_sample_records[[idx_sample]] <- data.frame( cate = pair, sample = sample_id,
          Correlation = sample_rho(a, b, rows, sample_id, sample_id),
          n_assay_pairs = sample_n(a, b, rows, sample_id, sample_id))
        idx_sample <- idx_sample + 1L
      }
      for (subject_id in unique(sample_info$subject)) {
        subject_samples <- sample_info[sample_info$subject == subject_id, , drop = FALSE]
        for (time_a in c("Day0", "Day3", "Day5")) {
          for (time_b in setdiff(c("Day0", "Day3", "Day5"), time_a)) {
            sample_a <- subject_samples$sample[subject_samples$time == time_a]
            sample_b <- subject_samples$sample[subject_samples$time == time_b]
            individual_records[[idx_individual]] <- data.frame(
              cate = pair, subject = subject_id, rep_a = time_a, rep_b = time_b,
              corr = sample_rho(a, b, rows, sample_a, sample_b), n_assay_pairs = sample_n(a, b, rows, sample_a, sample_b))
            idx_individual <- idx_individual + 1L
          }
        }
      }
      for (r in seq_len(nrow(rows))) {
        x <- matrices[[a]][rows$row_a[r], sample_ids]
        y <- matrices[[b]][rows$row_b[r], sample_ids]
        protein_records[[idx_protein]] <- data.frame( cate = pair, Protein = rows$protein[r], Somer = NA_character_,
          Correlation = safe_rho(x, y), n_paired_samples = sum(is.finite(x) & is.finite(y)))
        idx_protein <- idx_protein + 1L
      }
    }
  }
  same_sample_df <- dplyr::bind_rows(same_sample_records)
  within_individual_df <- dplyr::bind_rows(individual_records)
  per_protein_df <- dplyr::bind_rows(protein_records)
  lower_df <<- same_sample_df %>% dplyr::transmute(cate, value = Correlation) %>% dplyr::filter(is.finite(value))
  upper_df <<- per_protein_df %>%
    dplyr::transmute(cate, value = Correlation) %>% dplyr::filter(is.finite(value))
  plot_list_scope <- c(list(empty_plot()), lapply(display_platforms, header_plot))
  for (i in seq_along(display_platforms)) {
    plot_list_scope[[length(plot_list_scope) + 1L]] <- header_plot(display_platforms[i], TRUE)
    for (j in seq_along(display_platforms)) {
      if (i == j) {
        cell <- empty_plot(TRUE)
      } else if (i > j) {
        cell <- lower_cell(pair_label(display_platforms[j], display_platforms[i]))
      } else {
        cell <- upper_cell(pair_label(display_platforms[i], display_platforms[j]), show_median_line = TRUE)
      }
      plot_list_scope[[length(plot_list_scope) + 1L]] <- cell
    }
  }
  p_scope <- wrap_plots(plotlist = plot_list_scope, ncol = 6, byrow = TRUE) +
    plot_layout(widths = c(.28, rep(1, 5)), heights = c(.28, rep(1, 5))) +
    plot_annotation(title = sprintf("%s (N = %d)", plot_title, length(unique(scope_proteins))))

  output_matrix <- outputs[[paste0(output_prefix, "_matrix")]]
  output_legend <- outputs[[paste0(output_prefix, "_legend")]]
  output_evidence <- outputs[[paste0(output_prefix, "_evidence")]]
  output_audit <- outputs[[paste0(output_prefix, "_audit")]]

  ggplot2::ggsave(output_matrix, pdf_text_grob(p_scope), width = 24, height = 17, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  ggplot2::ggsave(output_legend, pdf_text_grob(p_corr_legend), width = 24, height = 18, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  saveRDS(list(same_sample = same_sample_df, within_individual = within_individual_df, per_protein = per_protein_df),
    output_evidence)
  saveRDS(list(ms = list( `MS-OT480` = data.frame(protein = rownames(ms_ot_zero),
    observed_n = rowSums(is.finite(ms_ot_scope$log2_raw)), zero_filled_n = rowSums(!is.finite(ms_ot_scope$log2_raw))),
    `MS-Astral` = data.frame(protein = rownames(ms_as_zero), observed_n = rowSums(is.finite(ms_as_scope$log2_raw)),
    zero_filled_n = rowSums(!is.finite(ms_as_scope$log2_raw)))),
    olink = dplyr::bind_rows(ol_exp_scope$audit, ol_ht_scope$audit),
    soma_mapping = dplyr::count(scope_map, protein, name = "n_somamers")), output_audit)

  invisible(list(matrix = p_scope, evidence = per_protein_df))
}
run_zero_scope_matrix("individual_detected", "Individual detected", "Overlap individual detected proteins", "individual")
run_zero_scope_matrix("broad_detected", "Broad detected", "Overlap broad detected proteins", "broad")

saveRDS(list(correlations = poscorr_corr,
  positive_proteins = poscorr_pro_list,
  overlap = poscorr_df),
  paste0(wdresult, "positive_correlation_overlap.rds"))
