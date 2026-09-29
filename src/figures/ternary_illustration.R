# Figure 1, 4: compact colocalization schematic.
# Run from the repository root:
#   Rscript ternary_illustration_revised.R
#
# Scientific content is unchanged:
#   original T = P3 (left vertex; Independent)
#   original R = P4 (top vertex; COLOC)
#   original L = 1 - P3 - P4 (right vertex)
#   classification boundaries remain P3 = 0.7 and P4 = 0.7.
#
# The original polygons are projected explicitly to an equilateral triangle.
# This uses ggplot2 rather than ggtern so the publication labels and small-font
# layout are deterministic and are not affected by ternary-axis clipping.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 first: install.packages('ggplot2')", call. = FALSE)
}
if (utils::packageVersion("ggplot2") < "3.4.0") {
  stop("This script requires ggplot2 >= 3.4.0 (linewidth support).", call. = FALSE)
}

# ----------------------- Settings to adjust ------------------------------
FINAL_WIDTH_MM  <- 28     # Width of this SMALL PANEL in the final figure,
                          # including its white margins; NOT the full Fig. 1.
EXPORT_WIDTH_IN <- 3      # Keep the original 3 x 3 inch export aspect ratio.
FONT_FAMILY     <- "Arial"

# These sizes are the intended sizes AFTER insertion at FINAL_WIDTH_MM.
REGION_PT <- 4.5          # ~3/4 of previous size; COLOC / Independent / Undetermined
AXIS_PT   <- 4.2          # ~3/4 of previous size; P3 / P4 / 1 - P3 - P4
TICK_PT   <- 4.0          # ~3/4 of previous size; threshold labels

THRESHOLD <- 0.7          # Scientific parameter: unchanged from the original.
OUTPUT_DIR <- "results/figures"
OUTPUT_STEM <- "ternary_illustration_fig1"  # Does not overwrite the old PDF.
SAVE_ACTUAL_SIZE_PDF <- FALSE               # Do not write the actual-size proof PDF.
SAVE_PNG <- FALSE                          # Do not write PNG preview.
# -----------------------------------------------------------------------

stopifnot(
  is.finite(FINAL_WIDTH_MM), FINAL_WIDTH_MM > 0,
  is.finite(EXPORT_WIDTH_IN), EXPORT_WIDTH_IN > 0,
  is.finite(THRESHOLD), THRESHOLD > 0.5, THRESHOLD < 1,
  all(is.finite(c(REGION_PT, AXIS_PT, TICK_PT))),
  all(c(REGION_PT, AXIS_PT, TICK_PT) > 0)
)
if (FINAL_WIDTH_MM < 26) {
  warning("This layout was designed for approximately 28 mm. Check label spacing at smaller widths.")
}

# Check fonts when systemfonts is available; do not install packages or fonts.
resolve_font <- function(requested) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    message("systemfonts is not installed; exact font availability cannot be checked. ",
            "The graphics device may substitute ", requested, ".")
    return(requested)
  }
  available <- unique(systemfonts::system_fonts()$family)
  candidates <- c(requested, "Helvetica", "Liberation Sans", "Nimbus Sans")
  hit <- match(tolower(candidates), tolower(available))
  found <- which(!is.na(hit))
  chosen <- if (length(found)) available[hit[found[1]]] else "sans"
  if (!identical(tolower(chosen), tolower(requested))) {
    warning(sprintf("Font '%s' is unavailable; using '%s'.", requested, chosen))
  }
  chosen
}

USE_CAIRO <- isTRUE(capabilities("cairo"))
if (USE_CAIRO) {
  FONT_FAMILY <- resolve_font(FONT_FAMILY)
} else {
  # The standard PDF device supports Helvetica without extra font setup.
  FONT_FAMILY <- "Helvetica"
  warning("Cairo is unavailable: PDFs will use standard Helvetica; PNG export is skipped.")
}

# The original vertices, parameterised without changing their values at 0.7.
tau <- THRESHOLD
q <- 1 - tau
plot_data <- data.frame(
  T = c(1, tau, tau, 0, q, 0, tau, tau, 0, 0, q),
  R = c(0, q, 0, 1, tau, tau, q, 0, 0, tau, tau),
  L = c(0, 0, q, 0, 0, q, 0, q, 1, q, 0),
  Label = factor(c("T", "T", "T", "R", "R", "R", "L", "L", "L", "L", "L"),
                 levels = c("L", "R", "T"))
)
stopifnot(all(abs(rowSums(plot_data[c("T", "R", "L")]) - 1) < 1e-10))

# Fix the original ggplot2 discrete colours explicitly, including alpha = 0.75.
REGION_COLOURS <- c(L = "#F8766D", R = "#00BA38", T = "#619CFF")

# Coordinates are on a fixed square canvas; margins are included in its width.
# T=P3, R=P4, L=residual. Only translation and a uniform scale are applied.
ternary_xy <- function(T, R, L) {
  data.frame(x = 0.13 + 0.76 * (L + 0.5 * R),
             y = 0.21 + 0.76 * sqrt(3) / 2 * R)
}
polygon_data <- cbind(plot_data, with(plot_data, ternary_xy(T, R, L)))
outline <- ternary_xy(c(1, 0, 0, 1), c(0, 1, 0, 0), c(0, 0, 1, 0))

segment_between <- function(a, b) {
  data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y)
}
# Draw each internal boundary once (rather than stacking polygon outlines).
boundaries <- rbind(
  segment_between(ternary_xy(q, tau, 0), ternary_xy(0, tau, q)),
  segment_between(ternary_xy(tau, q, 0), ternary_xy(tau, 0, q))
)

# Increasing-probability directions, matching the original axis orientation.
axis_point <- function(axis, value, offset = 0) {
  p <- switch(axis,
    P3 = ternary_xy(value, 1 - value, 0),
    P4 = ternary_xy(0, value, 1 - value),
    residual = ternary_xy(1 - value, 0, value),
    stop("Unknown axis: ", axis)
  )
  normal <- switch(axis,
    P3 = c(-sqrt(3) / 2, 0.5),
    P4 = c(sqrt(3) / 2, 0.5),
    residual = c(0, -1)
  )
  p$x <- p$x + offset * normal[1]
  p$y <- p$y + offset * normal[2]
  p
}
axis_arrows <- rbind(
  segment_between(axis_point("P3", 0.12, 0.060), axis_point("P3", 0.62, 0.060)),
  segment_between(axis_point("P4", 0.12, 0.060), axis_point("P4", 0.62, 0.060)),
  segment_between(axis_point("residual", 0.55, 0.035),
                  axis_point("residual", 0.95, 0.035))
)
axis_labels <- rbind(
  cbind(axis_point("P3", 0.42, 0.098), label = "P3", angle = 60),
  cbind(axis_point("P4", 0.42, 0.098), label = "P4", angle = -60),
  cbind(axis_point("residual", 0.74, 0.115), label = "1 - P3 - P4", angle = 0)
)

# A schematic needs the decision thresholds, not fifteen small tick labels.
# All probabilities now use the 0--1 scale, so 0.7 is the original 70%.
ticks <- do.call(rbind, lapply(c("P3", "P4"), function(axis) {
  segment_between(axis_point(axis, tau), axis_point(axis, tau, 0.018))
}))
tick_labels <- rbind(
  cbind(axis_point("P3", tau, 0.034), label = format(tau, trim = TRUE), hjust = 1),
  cbind(axis_point("P4", tau, 0.034), label = format(tau, trim = TRUE), hjust = 0)
)

# Keep the black-outline / white-fill labels. Place them so they partially cover
# the coloured regions, without leader lines, matching the user's preferred style.
region_labels <- data.frame(
  x = c(0.575, 0.235, 0.475),
  y = c(0.765, 0.200, 0.395),
  label = c("COLOC", "Independent", "Undetermined")
)

plot_ternary <- function(export_width_in = EXPORT_WIDTH_IN) {
  enlargement <- export_width_in * 25.4 / FINAL_WIDTH_MM
  # ggplot2 geom_text size is in mm by default; theme text size is in points.
  # Keep the conversion explicit for ggplot2 3.4/3.5 and later releases.
  pt_per_mm <- 72.27 / 25.4
  text_size <- function(final_pt) final_pt * enlargement / pt_per_mm
  line_size <- function(final_pt) final_pt * enlargement / pt_per_mm

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = polygon_data, ggplot2::aes(x, y, group = Label, fill = Label),
      colour = NA, alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(values = REGION_COLOURS, guide = "none") +
    ggplot2::geom_path(
      data = outline, ggplot2::aes(x, y),
      colour = "black", linewidth = line_size(0.45), linejoin = "mitre"
    ) +
    ggplot2::geom_segment(
      data = boundaries, ggplot2::aes(x, y, xend = xend, yend = yend),
      colour = "black", linewidth = line_size(0.35)
    ) +
    ggplot2::geom_segment(
      data = axis_arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
      colour = "black", linewidth = line_size(0.30),
      arrow = grid::arrow(length = grid::unit(0.55 * enlargement, "mm"),
                          type = "open", angle = 25)
    ) +
    ggplot2::geom_segment(
      data = ticks, ggplot2::aes(x, y, xend = xend, yend = yend),
      colour = "black", linewidth = line_size(0.30)
    ) +
    ggplot2::geom_text(
      data = tick_labels, ggplot2::aes(x, y, label = label, hjust = hjust),
      family = FONT_FAMILY, size = text_size(TICK_PT), colour = "black"
    ) +
    ggplot2::geom_text(
      data = axis_labels, ggplot2::aes(x, y, label = label, angle = angle),
      family = FONT_FAMILY, size = text_size(AXIS_PT), colour = "black"
    ) +
    ggplot2::geom_label(
      data = region_labels, ggplot2::aes(x, y, label = label),
      family = FONT_FAMILY, size = text_size(REGION_PT), colour = "black",
      fill = "white", label.size = line_size(0.30), label.padding = grid::unit(0.10, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1),
                        expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = FONT_FAMILY) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

if (!dir.exists(OUTPUT_DIR) && !dir.create(OUTPUT_DIR, recursive = TRUE)) {
  stop("Cannot create output directory: ", OUTPUT_DIR, call. = FALSE)
}

save_pdf <- function(filename, width_in) {
  if (USE_CAIRO) {
    ggplot2::ggsave(filename, plot = plot_ternary(width_in),
                   device = grDevices::cairo_pdf, family = FONT_FAMILY,
                   width = width_in, height = width_in, units = "in", bg = "white")
  } else {
    ggplot2::ggsave(filename, plot = plot_ternary(width_in),
                   device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE,
                   width = width_in, height = width_in, units = "in", bg = "white")
  }
  message("Saved: ", filename)
}

main_pdf <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".pdf"))
save_pdf(main_pdf, EXPORT_WIDTH_IN)
if (SAVE_ACTUAL_SIZE_PDF) {
  save_pdf(file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, "_actual_size.pdf")),
           FINAL_WIDTH_MM / 25.4)
}
if (SAVE_PNG && USE_CAIRO) {
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".png")),
    plot = plot_ternary(), device = "png", type = "cairo",
    width = EXPORT_WIDTH_IN, height = EXPORT_WIDTH_IN, units = "in",
    dpi = 600, bg = "transparent"
  )
}

message(sprintf(
  "Panel width: %.1f mm; font: %s; final text sizes: regions %.1f pt, axes %.1f pt, ticks %.1f pt.",
  FINAL_WIDTH_MM, FONT_FAMILY, REGION_PT, AXIS_PT, TICK_PT
))
message("Insert the complete square image without cropping its margins; otherwise the font-size conversion changes.")
