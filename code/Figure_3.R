# Figure 3: complete analysis and manuscript figures
# Run from the repository root.
source(file.path("code", "setup.R"), encoding = "UTF-8")
wdresult <- paste0(file.path(repo_dir, "result", "Fig3"), "/")
wdplot <- paste0(file.path(repo_dir, "plot", "Fig3"), "/")
dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
dir.create(wdplot, recursive = TRUE, showWarnings = FALSE)

protein_scopes <- c("detected", "individual_detected", "broad_detected")
for (scope_key in protein_scopes) {
  is_07_1_wrapper_run <- scope_key != "detected"

  wdresult <- paste0(file.path(repo_dir, "result", "Fig3", scope_key), "/")
  dir.create(wdresult, recursive = TRUE, showWarnings = FALSE)
  plot_prefix <- switch(scope_key, detected = "07-", individual_detected = "07.1-individual-",
    broad_detected = "07.1-broadly-")

  detection_sets_file <- file.path(repo_dir, "result", "Fig2", "01-01_detection_classification_sets.rds")
  detection_sets <- readRDS(detection_sets_file)
  prolist <- detection_sets[[scope_key]]
  scope_label <- switch( scope_key, detected = "Detected proteins", individual_detected = "Individual detected proteins",
    broad_detected = "Broadly detected proteins"
  )
  scope_short_label <- switch( scope_key, detected = "detected", individual_detected = "individual detected",
    broad_detected = "broadly detected"
  )
  hpa_proteinatlas_file <- file.path(wddata, "reference", "HPA", "proteinatlas.tsv")
  hpa_mapping_file <- file.path(wddata, "reference", "HPA", "hpa-ann-v2.csv")

  HPA_cite <- read.delim(hpa_proteinatlas_file, sep = "\t", header = TRUE, fill = TRUE)
  HPA_cite <- HPA_cite[,c(1,3,5,19,16,9,59,58,70,71,57,20,21,22,23)]
  HPA_cite$HPAgene <- HPA_cite$Gene

  hpaannbind <- read_csv(hpa_mapping_file)
  hpaannbindpro <- hpaannbind[,c(1,16)]
  hpaannbindpro <- na.omit(hpaannbindpro)
  hpaannbindpro <- hpaannbindpro[!duplicated(hpaannbindpro$protein1) & !is.na(hpaannbindpro$protein1),]

  HPA_cite <- left_join(HPA_cite, hpaannbindpro, by = "HPAgene")
  HPA_cite2 <- HPA_cite %>% dplyr::mutate(protein1 = if_else(is.na(protein1), HPAgene, protein1))
  HPA_cite2 <- HPA_cite2[,colnames(hpaannbind)]

  backprotein <- Reduce(union, prolist[1:length(names(prolist))])
  HPA_cite2 <- dplyr::filter(HPA_cite2, protein1 %in% backprotein)

  ann_list_detect <- prolist

  tissue_gene_detect <- tissuelistfunc(inputdir = HPA_cite2, backprotein = backprotein, outputdir = paste0(wdresult,
    "07-01_tissue_gene.RDS"))
  sc_gene_detect <- sclistfunc(inputdir = HPA_cite2, backprotein = backprotein, outputdir = paste0(wdresult,
    "07-01_sc_gene.RDS"))

  tissue_fe_detect <- calc_fold_enrichment(query_list = ann_list_detect, ref_list = tissue_gene_detect, background = backprotein)

  sc_fe_detect <- calc_fold_enrichment(query_list = ann_list_detect, ref_list = sc_gene_detect, background = backprotein)

  inputlist <- ann_list_detect
  cate <- gsub(" ", "_", scope_short_label)

  secretome_reference <- secretomefunc(backprotein)
  secregenelo <- secretome_reference$location
  secregenefun <- secretome_reference$function_sets
  secretome_detect_loc_fe <- calc_fold_enrichment(inputlist <- ann_list_detect, secregenelo, backprotein)
  secretome_detect_fun_fe <- calc_fold_enrichment(inputlist, secregenefun, backprotein)
  secretome_detect_res <- list(loc = secretome_detect_loc_fe, fun = secretome_detect_fun_fe)
  # 03. Compare tissue and cell-type enrichment ---------------------------------

  add_enrichment_odds_ratio <- function(fe_table) {
    required_cols <- c("query_size", "ref_size", "background_size", "overlap", "p_value")
    out <- fe_table
    a <- as.numeric(out$overlap)
    b <- as.numeric(out$query_size) - a
    c <- as.numeric(out$ref_size) - a
    d <- as.numeric(out$background_size) - as.numeric(out$query_size) - as.numeric(out$ref_size) + a
    zero_cell <- a == 0 | b == 0 | c == 0 | d == 0
    out$odds_ratio <- ifelse( zero_cell, ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)), (a * d) / (b * c)
    )
    out$odds_ratio_method <- ifelse(zero_cell, "Haldane–Anscombe 0.5 corrected", "Uncorrected")
    out
  }
  enrichment_raw_p_stars <- function(p_value) {
    dplyr::case_when( is.na(p_value) ~ "", p_value < 0.05 ~ "*", TRUE ~ ""
    )
  }
  plot_enrichment_lollipop <- function(fe_table, out_pdf, top_n, scope_label, reference_label) {
    plot_df <- add_enrichment_odds_ratio(fe_table) %>%
      dplyr::filter(!is.na(overlap), !is.na(ref_size), !is.na(odds_ratio), is.finite(odds_ratio))
    significant_ref_sets <- plot_df %>% dplyr::group_by(ref_set) %>% dplyr::summarise(
      any_significant = any(!is.na(p_value) & p_value < 0.05), .groups = "drop"
    ) %>% dplyr::filter(any_significant) %>% dplyr::pull(ref_set)
    plot_df <- dplyr::filter(plot_df, ref_set %in% significant_ref_sets)
    platform_order <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
    platform_colors <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
      "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF"
    )
    present_platforms <- unique(as.character(plot_df$query_set))
    platform_levels <- c( platform_order[platform_order %in% present_platforms], setdiff(present_platforms, platform_order)
    )
    platform_colors <- c(platform_colors, stats::setNames(rep("grey50", length(setdiff(present_platforms,
      names(platform_colors)))), setdiff(present_platforms, names(platform_colors))))

    ref_rank <- plot_df %>%
      dplyr::group_by(ref_set) %>% dplyr::summarise(total_overlap = sum(overlap, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total_overlap), ref_set) %>% dplyr::slice_head(n = top_n)
    ref_levels <- rev(ref_rank$ref_set)
    plot_df <- plot_df %>% dplyr::filter(ref_set %in% ref_rank$ref_set) %>% dplyr::mutate(
      ref_set = factor(ref_set, levels = ref_levels), query_set = factor(query_set, levels = platform_levels),
      y_base = as.numeric(ref_set), y_offset = (as.numeric(query_set) - (length(platform_levels) + 1) / 2) * 0.16,
      y_plot = y_base + y_offset, significance = enrichment_raw_p_stars(p_value),
      star_x = overlap + max(0.5, max(overlap, na.rm = TRUE) * 0.03)
    )
    significant_df <- dplyr::filter(plot_df, significance != "")

    p <- ggplot(plot_df) +
      geom_segment( aes(x = 0, xend = overlap, y = y_plot, yend = y_plot, colour = query_set),
      linewidth = 0.35, alpha = 1, show.legend = FALSE
    ) +

      geom_point(aes(x = overlap, y = y_plot, colour = query_set, size = odds_ratio), alpha = 1) + geom_text(
      data = significant_df, aes(x = star_x, y = y_plot, label = significance),
      inherit.aes = FALSE, size = 2.1, fontface = "bold", colour = "black", show.legend = FALSE
    ) + scale_colour_manual( values = platform_colors[platform_levels], breaks = platform_levels, name = "Platform"
    ) +
      scale_size_continuous( name = "Odds Ratio", range = c(0.9, 3.2), trans = "log10"
    ) +
      scale_y_continuous(breaks = seq_along(ref_levels), labels = ref_levels, expand = expansion(add = 0.45)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs( x = "Proteins overlapping reference set",
      y = reference_label, title = paste0(scope_label, " enrichment")
    ) +
      theme2_nobacknoerror + themetext30 + themelegendright + theme( axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9), plot.title = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7, hjust = 0), legend.key.size = grid::unit(4, "mm")
    )

    invisible(list(plot = p, plot_data = plot_df, reference_order = ref_rank))
  }
  tissue_lollipop_result <- plot_enrichment_lollipop( tissue_fe_detect,
    out_pdf = paste0(wdplot, "F3b-07-03a_tissue_lollipop_enrichment.pdf"), top_n = 15, scope_label = scope_label,
    reference_label = "Tissue-specific reference set"
  )
  sc_lollipop_result <- plot_enrichment_lollipop( sc_fe_detect,
    out_pdf = paste0(wdplot, "F3c-07-03a_sc_lollipop_enrichment.pdf"), top_n = 18, scope_label = scope_label,
    reference_label = "Cell-type-specific reference set"
  )

  lollipop_combined <- patchwork::wrap_plots( tissue_lollipop_result$plot, sc_lollipop_result$plot, nrow = 1, ncol = 2,
    widths = c(1, 1)
  )
  if (!is_07_1_wrapper_run) {
    ggsave( paste0(wdplot, "F3bc_tissue_sc_lollipop_enrichment.pdf"), pdf_text_grob(lollipop_combined), width = 32, height = max( 12,
      0.65 * max(length(tissue_lollipop_result$reference_order), length(sc_lollipop_result$reference_order)) + 3
    ), units = "cm", device = grDevices::cairo_pdf, family = "Helvetica"
    )
  }
  # 03.1 Summarise significant tissue and cell-type enrichment ------------------

  plot_enrichment_lollipop_significant_only <- function(fe_table, out_pdf, top_n,
    scope_label, reference_label) {
    valid_df <- add_enrichment_odds_ratio(fe_table) %>%
      dplyr::filter( !is.na(overlap), !is.na(ref_size), !is.na(odds_ratio), is.finite(odds_ratio), ref_size >= 5
    )
    significant_ref_sets <- valid_df %>%
      dplyr::group_by(ref_set) %>% dplyr::summarise( any_significant = any(!is.na(p_value) & p_value < 0.05),
      .groups = "drop"
    ) %>%
      dplyr::filter(any_significant) %>% dplyr::pull(ref_set)
    platform_order <- c("MS-OT480", "MS-Astral", "Olink-Explore", "Olink-HT", "Somascan")
    platform_colors <- c( "MS-OT480" = "#D9B381", "MS-Astral" = "#B5252E",
      "Olink-Explore" = "#264653", "Olink-HT" = "#AEC7E8FF", "Somascan" = "#69B190FF"
    )
    present_platforms <- unique(as.character(valid_df$query_set))
    platform_levels <- c( platform_order[platform_order %in% present_platforms], setdiff(present_platforms, platform_order)
    )
    extra_platforms <- setdiff(present_platforms, names(platform_colors))
    if (length(extra_platforms) > 0) {
      platform_colors <- c( platform_colors, stats::setNames(rep("grey50", length(extra_platforms)), extra_platforms)
      )
    }
    ref_rank <- valid_df %>% dplyr::filter(ref_set %in% significant_ref_sets) %>% dplyr::group_by(ref_set) %>%
      dplyr::summarise(total_overlap = sum(overlap, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total_overlap), ref_set) %>% dplyr::slice_head(n = top_n)
    ref_levels <- rev(ref_rank$ref_set)

    plot_df <- valid_df %>% dplyr::filter(ref_set %in% ref_rank$ref_set, !is.na(p_value), p_value < 0.05) %>%
      dplyr::mutate( ref_set = factor(ref_set, levels = ref_levels), query_set = factor(query_set, levels = platform_levels),
      y_base = as.numeric(ref_set), y_offset = (as.numeric(query_set) - (length(platform_levels) + 1) / 2) * 0.16,
      y_plot = y_base + y_offset
    )
    p <- ggplot(plot_df) +
      geom_segment( aes(x = 0, xend = overlap, y = y_plot, yend = y_plot, colour = query_set),
      linewidth = 0.35, alpha = 1, show.legend = FALSE
    ) +
      geom_point( aes(x = overlap, y = y_plot, colour = query_set, size = odds_ratio), alpha = 1
    ) + scale_colour_manual( values = platform_colors[platform_levels], breaks = platform_levels, name = "Platform"
    ) +
      scale_size_continuous( name = "Odds Ratio", range = c(0.9, 3.2), trans = "log10"
    ) + scale_y_continuous( breaks = seq_along(ref_levels), labels = ref_levels, expand = expansion(add = 0.45)
    ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs( x = "Proteins overlapping reference set",
      y = reference_label, title = paste0(scope_label, " enrichment")
    ) +
      theme2_nobacknoerror + themetext30 + themelegendright + theme( axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9), plot.title = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7, hjust = 0), legend.key.size = grid::unit(4, "mm")
    )

    invisible(list(plot = p, plot_data = plot_df, reference_order = ref_rank))
  }
  tissue_lollipop_significant_result <- plot_enrichment_lollipop_significant_only( tissue_fe_detect,
    out_pdf = paste0(wdplot, plot_prefix, "03b_tissue_lollipop_significant.pdf"), top_n = 15, scope_label = scope_label,
    reference_label = "Tissue-specific reference set"
  )
  sc_lollipop_significant_result <- plot_enrichment_lollipop_significant_only( sc_fe_detect,
    out_pdf = paste0(wdplot, plot_prefix, "03b_sc_lollipop_significant.pdf"), top_n = 18, scope_label = scope_label,
    reference_label = "Cell-type-specific reference set"
  )

  significant_lollipop_combined <- patchwork::wrap_plots( (tissue_lollipop_significant_result$plot+themelegendnone),
    sc_lollipop_significant_result$plot, nrow = 1, ncol = 2, widths = c(1, 1)
  )

  if (is_07_1_wrapper_run) {
    figure_prefix <- if (scope_key == "individual_detected") {
      "S2a"
    } else {
      "S2b"
    }
    ggsave(paste0(wdplot, figure_prefix, "_tissue_sc_lollipop_significant.pdf"), pdf_text_grob(significant_lollipop_combined),
      width = 26, height = max(8.5, 0.5 * max(length(tissue_lollipop_significant_result$reference_order),
      length(sc_lollipop_significant_result$reference_order)) + 2.5), units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
  }
  if (scope_key == "detected") {
    library(webshot)
    library(networkD3)

    p1 <- secretome_coverage(secregenelo, inputlist)
    p2 <- secretome_coverage(secregenefun, inputlist)
    df1 <- p1[["heatmap"]][["data"]]
    df2 <- p2[["heatmap"]][["data"]]

    df1 <- df1 %>% dplyr::mutate( X1_number = as.integer(str_extract(X1, "(?<=\\()[0-9]+(?=\\))")),
      X2_number = as.integer(str_extract(X2, "(?<=\\()[0-9]+(?=\\))")), X1 = str_replace(X1, " \\(.*\\)$", ""),
      X2 = str_replace(X2, " \\(.*\\)$", ""), number = X2_number*value
    )%>% dplyr::select(source = X2, target = X1, value = number)
    df1$target <- gsub(" ", "_", df1$target)
    df1$source <- gsub(" ", "_", df1$source)

    dfbind <- rbind(df1)

    df <- dfbind %>% dplyr::filter(!is.na(value)) %>% dplyr::mutate(source = as.character(source),
      target = as.character(target))
    nodes_vec <- unique(c(df$source, df$target))
    node2 <- data.frame(name = nodes_vec, group = nodes_vec)
    df <- df %>% dplyr::mutate(IDsource = match(source, node2$name) - 1, IDtarget = match(target, node2$name) - 1,
      group = target)

    col <- c("MS-OT480"="#D9B381", "MS-Astral"="#B5252E" ,"Olink-Explore" = "#264653", "Olink-HT"= "#AEC7E8","Somascan"="#69B190",
      "Immunoglobulin genes" ="#C27877","Intracellular and membrane"="#EDC4BB","Secreted - unknown location"="grey",
      "Secreted in brain"="#D4C0AC","Secreted in female reproductive system"="#9E7A64",
        "Secreted in male reproductive system"="#7D576F" ,
      "Secreted in other tissues"="#C3A8D1" ,
      "Secreted to blood"="#BC5D50","Secreted to digestive system"="#C97231","Secreted to extracellular matrix"="#E1C367"  )
    names(col) <- gsub(" ", "_", names(col))

    color_scale <- paste0( "d3.scaleOrdinal()", ".domain([\"", paste(names(col), collapse = "\", \""), "\"])",
      ".range([\"", paste(col, collapse = "\", \""), "\"])"
    )

    df$source <- factor(df$source,levels = unique(df$source))

    p <- sankeyNetwork( Links = df, Nodes = node2, Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "name",
      NodeGroup = "group", LinkGroup = "group", colourScale = color_scale, fontSize = 12, fontFamily = "Helvetica",
      nodeWidth = 20, nodePadding = 19
    )
    p <- htmlwidgets::onRender(p, '
      function(el) {
        d3.select(el).selectAll("text")
          .style("font-family", "Helvetica")
          .style("fill", "#000000")
          .style("fill-opacity", 1)
          .style("opacity", 1)
          .style("text-shadow", "none")
          .text(function() { return this.textContent.replace(/[\\u2010-\\u2014\\u2212]/g, "-"); });
      }
    ')

    sankey_html_file <- file.path(wdresult, paste0(plot_prefix, "04_secretome_sankey.html"))
    saveNetwork(p, file = sankey_html_file, selfcontained = FALSE)

    if (!is_07_1_wrapper_run) {
      # PhantomJS needs an absolute, percent-encoded file URL for non-ASCII paths.
      sankey_html_path <- normalizePath(sankey_html_file, winslash = "/", mustWork = TRUE)
      sankey_html_url <- paste0("file://", gsub("%2F", "/",
        utils::URLencode(enc2utf8(sankey_html_path), reserved = TRUE), fixed = TRUE))
      sankey_pdf_file <- file.path(normalizePath(wdplot, winslash = "/", mustWork = TRUE), "F3a_secretome_sankey.pdf")
      webshot::webshot(sankey_html_url, file = sankey_pdf_file, vwidth = 350, vheight = 600)
    }

    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(DOSE)

    library(RColorBrewer)

    backgroundgene <- Reduce(union, prolist[1:5])
    go_enrich_function <- function( gene , backgroundgene) {
      OrgDb <- "org.Hs.eg.db"
      go_up <- enrichGO(gene=gene, OrgDb="org.Hs.eg.db", universe = backgroundgene, keyType = "SYMBOL", ont="CC",
        pAdjustMethod = "BH", qvalueCutoff=0.05, pvalueCutoff=0.05,
      )
      r_go_up <- go_up@result
      r_go_up <- dplyr::filter(r_go_up, pvalue < 0.05)
      return(r_go_up)
    }
    custom_colors <- c(
      "Overlap protein" =  "#D1C1E1FF", "MS-OT480"="#D9B381", "MS-Astral"="#B5252E" ,"Olink-Explore" = "#264653",
      "Olink-HT"= "#AEC7E8","Somascan"="#69B190"

    )

    result1 <- NULL

    p <- list()
    i <-  "MS-OT480"
    for (i in names(prolist)) {
      proteininput <- prolist[[i]]

      gore <- go_enrich_function(proteininput,backgroundgene)
      gore$cate <- i
      gore <- gore %>% head(n = 5)
      result1 <- rbind(result1,gore)

      p[[i]] <- ggplot(data = gore, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)), fill = cate)) +
        geom_bar(stat = "identity", width = 0.7, position = position_dodge(preserve = "single")) +
        labs(x = "−log10(P value)",fill="") +
        scale_y_discrete(position = "right") + theme2_nobackplus + scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
        scale_fill_manual(values = custom_colors,name="")+ geom_text(data = gore, aes(x = 0.1,
        y = reorder(Description, -log10(pvalue)), label = Description), color = "black", size = 2.5, vjust = 0.5, hjust = 0) +
        theme( axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y.right = element_blank(),
        axis.ticks.y.left = element_line(), axis.line.y.right = element_blank(), axis.line.y.left = element_line(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      ) + themelegendleft
    }
    saveRDS(result1, paste0(wdresult, "go_cc_results.rds"))
    combined_plot <- do.call(gridExtra::grid.arrange, c(p, ncol = 1))

    if (!is_07_1_wrapper_run) {
      ggsave(paste0(wdplot, "F3a_go_cc_all.pdf"),pdf_text_grob(combined_plot),width = 8,  height = 12, units = "cm", device = grDevices::cairo_pdf, family = "Helvetica")
    }

  }
  saveRDS(list(protein_sets = prolist, tissue_reference = tissue_gene_detect,
    cell_reference = sc_gene_detect, tissue_enrichment = tissue_fe_detect,
    cell_enrichment = sc_fe_detect, secretome_enrichment = secretome_detect_res),
    paste0(wdresult, "enrichment_results.rds"))
}
