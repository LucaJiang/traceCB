# illustration of colocalization analysis
library(ggtern)
library(ggplot2)

plot_data <- data.frame(
  "T" = c(1, 0.7, 0.7, 0, 0.3, 0, 0.7, 0.7, 0, 0, 0.3),
  "R" = c(0, 0.3, 0, 1, 0.7, 0.7, 0.3, 0, 0, 0.7, 0.7),
  "L" = c(0, 0, 0.3, 0, 0, 0.3, 0, 0.3, 1, 0.3, 0),
  "Label" = as.factor(c(
    "T", "T", "T", "R", "R", "R", "L", "L", "L", "L", "L"
  ))
)
dfLabels <- plyr::ddply(plot_data, "Label", function(df) {
  if (df$Label[1] == "T") {
    label <- "Independent"
  } else if (df$Label[1] == "R") {
    label <- "COLOC"
  } else if (df$Label[1] == "L") {
    label <- "Undetermined"
  } else {
    label <- "Unknown"
  }
  means <- colMeans(df[setdiff(colnames(df), "Label")])
  result <- data.frame(
    T = means["T"],
    R = means["R"],
    L = means["L"],
    Label = label
  )
  return(result)
})
plot_ternary <- function(data) {
  result <- ggtern(
    data = data,
    aes(
      x = T,
      y = R,
      z = L,
      color = Label
    )
  ) +
    geom_polygon(
      mapping = aes(fill = Label),
      alpha = 0.75,
      size = 0.5,
      color = "black"
    ) +
    geom_label(
      data = dfLabels,
      mapping = aes(label = Label),
      size = 2.6,
      color = "black",
      fill = "white",
      alpha = 1
    ) +
    theme_showgrid_minor() +
    theme_showarrows() +
    guides(color = "none", fill = "none") +
    Tlab("", "P4 = P(colocalized signals)") +
    Llab("", "P3 = P(independent signals)") +
    Rlab("", "1 - P3 - P4") +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    )

  return(result)
}
ggsave(
  "/Users/lucajiang/learn/CityU/traceCB/data/img/ternary_illu.pdf",
  plot = plot_ternary(plot_data),
  width = 3,
  height = 3
)
