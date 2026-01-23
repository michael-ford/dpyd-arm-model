# =============================================================================
# R/visualizations/plot_network.R
# Figure 2: Network Diagram (Star Topology)
# =============================================================================
# Generates a network diagram showing the structure of the DPYD NMA.
# All variant treatments connect to Wild-type (WT) in a star topology,
# reflecting the arm-based model structure.
# =============================================================================

# Source constants for colors, labels, and theme (if not already loaded)
if (!exists("VARIANT_COLORS")) {
  source("R/visualizations/constants.R")
}

library(igraph)
library(ggraph)

#' Generate Network Diagram for DPYD NMA
#'
#' Creates a star-topology network diagram showing the relationships between
#' DPYD variants in the network meta-analysis. Node sizes are proportional
#' to the number of study arms for each treatment.
#'
#' @param nma_data Data frame containing NMA data with columns: s.id, t.id, r, n
#' @param treatment_names Character vector of treatment names corresponding to t.id values
#' @param output_dir Output directory path (default: "output/wt_unified")
#'
#' @return A ggplot2/ggraph object containing the network diagram
#'
#' @examples
#' nma_data <- readRDS("output/wt_unified/wt_unified_nma_data.rds")
#' treatment_names <- readRDS("output/wt_unified/wt_unified_treatment_names.rds")
#' p <- plot_network(nma_data, treatment_names)

plot_network <- function(nma_data, treatment_names, output_dir = "output/wt_unified") {

  # 1. Count arms (rows) per treatment
  arm_counts <- table(treatment_names[nma_data$t.id])

  # 2. Create edge list (all variants connect to WT - star topology)
  # In arm-based model, all treatments are compared through WT
  edges <- data.frame(
    from = rep("WT", 4),
    to = c("HapB3", "2846hetho", "2Ahetho", "13hetho"),
    stringsAsFactors = FALSE
  )

  # 3. Create graph object
  g <- graph_from_data_frame(edges, directed = FALSE)

  # 4. Add node attributes
  V(g)$size <- as.numeric(arm_counts[V(g)$name])
  V(g)$color <- VARIANT_COLORS[V(g)$name]
  V(g)$label <- TREATMENT_LABELS[V(g)$name]

  # 5. Create plot with ggraph
  p <- ggraph(g, layout = "star") +
    geom_edge_link(
      width = 1.5,
      color = "gray40",
      alpha = 0.8
    ) +
    geom_node_point(
      aes(size = size, color = name),
      show.legend = FALSE
    ) +
    geom_node_text(
      aes(label = label),
      repel = TRUE,
      size = 4,
      fontface = "bold"
    ) +
    scale_color_manual(values = VARIANT_COLORS) +
    scale_size_continuous(range = c(8, 20)) +
    theme_void() +
    labs(
      title = "Network Structure: DPYD Variants",
      subtitle = "Node size proportional to number of study arms"
    )

  # 6. Create figures directory if it doesn't exist
  fig_dir <- file.path(output_dir, "figures")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
  }

  # 7. Save figure
  ggsave(
    filename = file.path(fig_dir, "figure2_network.png"),
    plot = p,
    width = FIG_DIMS$network["width"],
    height = FIG_DIMS$network["height"],
    dpi = FIG_DPI
  )

  cat("Saved: figures/figure2_network.png\n")

  return(p)
}
