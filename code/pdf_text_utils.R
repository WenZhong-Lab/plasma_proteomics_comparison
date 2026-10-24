# Prepare export copies only; retain analysis objects and non-text graphics.
configure_pdf_text <- function() {
  if (!capabilities("cairo")) {
    stop("Cairo is required to export PDFs with embedded editable fonts.")
  }
  if (isNamespaceLoaded("showtext")) {
    showtext::showtext_auto(FALSE)
  }
  ggplot2::theme_set(ggplot2::theme_get() +
    ggplot2::theme(text = ggplot2::element_text(family = "Helvetica", colour = "#000000")))
}

pdf_text_grob <- function(plot) {
  if (!inherits(plot, c("ggplot", "patchwork", "gList")) && !grid::is.grob(plot)) {
    stop("PDF text preparation requires a ggplot, patchwork or grid grob.")
  }
  # Build on the destination device at its original export dimensions.
  grid::gTree(source_plot = plot, cl = "pdf_text_export")
}

makeContent.pdf_text_export <- function(x) {
  plot <- x$source_plot
  # Let ggbreak assemble its original broken-axis layout before editing text.
  if (inherits(plot, "ggbreak")) {
    draw_break <- utils::getS3method("grid.draw", "ggbreak")
    plot <- draw_break(plot, recording = FALSE)
  }
  if (inherits(plot, "patchwork")) {
    plot <- patchwork::patchworkGrob(plot)
  } else if (inherits(plot, "ggplot")) {
    plot <- ggplot2::ggplotGrob(plot)
  } else if (inherits(plot, "gList")) {
    plot <- grid::gTree(children = plot)
  }
  if (!grid::is.grob(plot)) {
    stop("PDF text preparation requires a ggplot, patchwork or grid grob.")
  }
  # Materialise deferred axis/annotation grobs so their text is included.
  plot <- grid::grid.force(plot)
  normalise_text <- function(grob) {
    if (inherits(grob, "text")) {
      if (is.null(grob$gp)) {
        grob$gp <- grid::gpar()
      }
      grob$gp$fontfamily <- "Helvetica"
      grob$gp$col <- "#000000"
      grob$gp$alpha <- 1
      if (is.character(grob$label)) {
        grob$label <- gsub("[\u2010-\u2014\u2212]", "-", enc2utf8(grob$label), perl = TRUE)
      }
    }
    if (length(grob$grobs)) {
      grob$grobs <- lapply(grob$grobs, normalise_text)
    }
    if (length(grob$children)) {
      grob$children <- do.call(grid::gList, lapply(grob$children, normalise_text))
    }
    grob
  }
  grid::setChildren(x, grid::gList(normalise_text(plot)))
}
