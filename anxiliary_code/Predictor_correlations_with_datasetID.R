#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                 Predictor_correlations_with_datasetID.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#
rm(list=ls())
#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(skimr)
library(sf)
#fix_windows_histograms()

source(here::here("R/00_Configuration.R"))

predictors <- readRDS(here("Data/output/1_data/2_predictors.rds"))


#----------------------------------------------------------#
# Test statistic difference in predictor for datasetID -----
#----------------------------------------------------------#

dd_meta <- predictors %>%
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  skim_tee()

library(rstatix)

test_results <- predictors %>%
  filter(samplingPeriodID== 1) %>%
  select(-atlas_relCirc) %>%
  pivot_longer(cols = where(is.numeric), names_to = "variable", values_to = "vals") %>%
  select(datasetID, variable, vals) %>%
  group_by(variable) %>%
  dunn_test(vals ~ datasetID, p.adjust.method = "BH")

# Show results
library(flextable)
test_results %>% kableExtra::kable()
test_results %>% flextable() %>%
  colformat_double(digits = 3) %>%
  autofit() %>%
  save_as_docx(path = here("Documentation/Influence_datasetID_dunn_test.docx"))

library(ggcorrplot)
library(corrr)

cor_matrix<- predictors %>%
  filter(samplingPeriodID== 1) %>% select(where(is.numeric)) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs", diagonal = 0)

# Replace perfect correlations with NA
cor_matrix[cor_matrix == 1 | cor_matrix == -1] <- 0

# Plot the cleaned correlation network
network_plot(cor_matrix, min_cor = 0.3)


# custom Plotly network:
library(corrr)
library(dplyr)
library(igraph)
library(plotly)

# Compute Spearman correlation matrix
cor_matrix <- predictors %>%
  filter(samplingPeriodID == 1) %>%
  select(where(is.numeric)) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter significant correlations (|r| > 0.3)
edges <- cor_matrix %>%
  stretch() %>%
  filter(abs(r) >= 0.3, x != y)  # Remove self-correlations

# Create igraph object
graph <- graph_from_data_frame(edges, directed = FALSE)

# Get node positions using force-directed layout
layout <- layout_with_fr(graph)
layout_df <- as.data.frame(layout)
colnames(layout_df) <- c("x", "y")

# Get edge list for plotting
edges_df <- get.data.frame(graph)

# Assign edge colors based on node connections
edges_df$color <- ifelse(edges_df$from == "AOO" | edges_df$to == "AOO", "red",
                         ifelse(edges_df$from == "rel_occ_Ncells" | edges_df$to == "rel_occ_Ncells", "blue", "#FFEB3B33"))

# Create Plotly network graph
p <- plot_ly(type = "scatter", mode = "markers+text")

# Add nodes
p <- p %>%
  add_trace(
    x = layout_df$x, y = layout_df$y,
    text = names(V(graph)), textposition = "top center",
    marker = list(size = 10, color = "lightgrey"),
    hoverinfo = "text"
  )

# Add edges individually (fixing the color issue)
for (i in 1:nrow(edges_df)) {
  p <- p %>%
    add_segments(
      x = layout_df[match(edges_df$from[i], names(V(graph))), "x"],
      y = layout_df[match(edges_df$from[i], names(V(graph))), "y"],
      xend = layout_df[match(edges_df$to[i], names(V(graph))), "x"],
      yend = layout_df[match(edges_df$to[i], names(V(graph))), "y"],
      line = list(color = edges_df$color[i], width = abs(edges_df$r[i]) * 5),
      hoverinfo = "none"
    )
}

# Final layout
p <- p %>%
  layout(
    title = "Correlation Network Plot (Plotly)",
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
  )

p

library(corrr)
library(dplyr)
library(igraph)
library(plotly)

# Compute Spearman correlation matrix
cor_matrix <- predictors %>%
  filter(samplingPeriodID == 1) %>%
  select(where(is.numeric)) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter significant correlations (|r| > 0.3)
edges <- cor_matrix %>%
  stretch() %>%
  filter(abs(r) >= 0.3, x != y)  # Remove self-correlations

# Extract red network (AOO connections)
red_edges <- edges %>%
  filter(x == "AOO" | y == "AOO")

# Extract blue network (rel_occ_Ncells connections)
blue_edges <- edges %>%
  filter(x == "rel_occ_Ncells" | y == "rel_occ_Ncells")

# Find all nodes that appear in the blue network
blue_nodes <- unique(c(blue_edges$x, blue_edges$y))

# Remove nodes from red network if they appear in the blue network
filtered_red_edges <- red_edges %>%
  filter(!(x %in% blue_nodes | y %in% blue_nodes))

# Create igraph object for this subset
graph_subset <- graph_from_data_frame(filtered_red_edges, directed = FALSE)

# Use force-directed layout with increased spacing
layout <- layout_with_fr(graph_subset, repulserad = 3)
layout_df <- as.data.frame(layout)
colnames(layout_df) <- c("x", "y")

# Get edge list for plotting
edges_df <- get.data.frame(graph_subset)

# Create Plotly network graph
p <- plot_ly(type = "scatter", mode = "markers+text")

# Add nodes
p <- p %>%
  add_trace(
    x = layout_df$x, y = layout_df$y,
    text = names(V(graph_subset)), textposition = "top center",
    marker = list(size = 16, color = "black", opacity = 1),
    hoverinfo = "text"
  )

# Add only red edges
for (i in 1:nrow(edges_df)) {
  corr_value <- round(edges_df$r[i], 2)  # Round correlation value
  hover_text <- paste0("Correlation: ", corr_value)

  p <- p %>%
    add_segments(
      x = layout_df[match(edges_df$from[i], names(V(graph_subset))), "x"],
      y = layout_df[match(edges_df$from[i], names(V(graph_subset))), "y"],
      xend = layout_df[match(edges_df$to[i], names(V(graph_subset))), "x"],
      yend = layout_df[match(edges_df$to[i], names(V(graph_subset))), "y"],
      line = list(color = "red", width = abs(edges_df$r[i]) * 8, opacity = 0.6),
      hoverinfo = "text",
      text = hover_text
    )
}

# Final layout adjustments
p <- p %>%
  layout(
    title = "Filtered Red Network (No Blue Nodes)",
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
  )

p


### For all variables separately:
library(corrr)
library(dplyr)
library(igraph)
library(plotly)

# Compute Spearman correlation matrix
cor_matrix <- predictors %>%
  filter(samplingPeriodID == 1) %>%
  select(where(is.numeric)) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter significant correlations (|r| > 0.3)
edges <- cor_matrix %>%
  stretch() %>%
  filter(abs(r) >= 0.3, x != y)  # Remove self-correlations

# Get a list of all unique variables in the network
unique_variables <- unique(c(edges$x, edges$y))

# Create a list to store plots
plot_list <- list()

# Loop through each variable
for (var in unique_variables) {

  # Filter edges that involve the current variable
  edges_subset <- edges %>%
    filter(x == var | y == var)

  # Skip if no edges exist for this variable
  if (nrow(edges_subset) == 0) next

  # Create igraph object for this subset
  graph_subset <- graph_from_data_frame(edges_subset, directed = FALSE)

  # Use a better layout for node spreading
  layout <- layout_with_kk(graph_subset)  # Kamada-Kawai layout for better separation
  layout_df <- as.data.frame(layout)
  colnames(layout_df) <- c("x", "y")

  # Get edge list for plotting
  edges_df <- get.data.frame(graph_subset)

  # Assign edge colors
  edges_df$color <- ifelse(edges_df$from == var | edges_df$to == var, "red", "gray")

  # Create Plotly network graph
  p <- plot_ly(type = "scatter", mode = "markers+text")

  # Add nodes
  p <- p %>%
    add_trace(
      x = layout_df$x, y = layout_df$y,
      text = names(V(graph_subset)), textposition = "top center",
      marker = list(size = 12, color = "blue"),
      hoverinfo = "text"
    )

  # Add edges individually (to ensure correct colors and display correlation values)
  for (i in 1:nrow(edges_df)) {
    corr_value <- round(edges_df$r[i], 2)  # Round correlation value to 2 decimal places
    hover_text <- paste0("Correlation: ", corr_value)

    p <- p %>%
      add_segments(
        x = layout_df[match(edges_df$from[i], names(V(graph_subset))), "x"],
        y = layout_df[match(edges_df$from[i], names(V(graph_subset))), "y"],
        xend = layout_df[match(edges_df$to[i], names(V(graph_subset))), "x"],
        yend = layout_df[match(edges_df$to[i], names(V(graph_subset))), "y"],
        line = list(color = edges_df$color[i], width = abs(edges_df$r[i]) * 5, opacity = 0.3),  # Adjust opacity
        hoverinfo = "text",
        text = hover_text  # Display correlation value on hover
      )
  }

  # Final layout
  p <- p %>%
    layout(
      title = paste("Network for:", var),
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
    )

  # Store in the plot list
  plot_list[[var]] <- p
}

library(plotly)

# Convert plot list to a subplot layout (arrange in a grid)
grid_plot <- subplot(plot_list, nrows = ceiling(sqrt(length(plot_list))),
                     titleX = TRUE, titleY = TRUE,
                     margin = 0.05)

# Show the grid plot
grid_plot


plot_list[[4]]

















#############################
############################

# Compute Spearman correlation matrix
cor_matrix <- predictors %>%
  filter(samplingPeriodID == 1) %>%
  select(where(is.numeric)) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter significant correlations (|r| > 0.3)
edges <- cor_matrix %>%
  stretch() %>%
  filter( x != y)  # Remove self-correlations

# Extract red network (AOO connections)
red_edges <- edges %>%
  filter(x == "AOO" | y == "AOO")

# Extract blue network (rel_occ_Ncells connections)
blue_edges <- edges %>%
  filter(x == "rel_occ_Ncells" | y == "rel_occ_Ncells")


# Create a dataframe with all unique variables from the 'y' column
common_df <- full_join(
  red_edges %>% select(y, r) %>% rename(r_AOO = r) %>% filter(y != "AOO"),
  blue_edges %>% select(y, r) %>% rename(r_rel_occ_Ncells = r) %>% filter(y != "rel_occ_Ncells"),
  by = "y"
) %>%
  arrange(desc(abs(r_AOO)), desc(abs(r_rel_occ_Ncells))) %>%   # Sort by absolute correlation strength
  distinct()
# %>%
#  mutate(
#    strong_cor_AOO = case_when(abs(r_AOO) >= 0.5 ~ 1,
#                               .default = 0),
#
#    strong_cor_prev = case_when(abs(r_rel_occ_Ncells) >= 0.5 ~ 1,
#                                .default = 0)
#  )

# Print the dataframe
print(common_df, n = 45)

common_df %>%
  write.csv(here("Documentation/META_relative_absolute_AOO_correlation_comparison.csv"))





common_df_long <- common_df %>%
  pivot_longer(cols = c(r_AOO, r_rel_occ_Ncells),
               names_to = "variable",
               values_to = "r") %>%
  mutate(strong_cor = case_when(abs(r) >= 0.5 ~ 1,
                                .default = 0))
p_plotly <- ggplotly(
  ggplot(common_df_long, aes(y = y, x = r, color = factor(variable), shape = factor(strong_cor))) +
    geom_point(size = 3, aes(fill = factor(variable))) +  # Fill for filled circles
    scale_shape_manual(values = c("0" = 1, "1" = 21)) +  # Open = 1, Closed = 21
    scale_fill_manual(values = c("red", "blue")) +  # Adjust colors as needed
    scale_color_manual(values = c("red", "blue"))+
    theme_minimal()+
    ggplot2::geom_vline(xintercept = 0, col = "grey"))


library(htmlwidgets)

# Save interactive Plotly plot
saveWidget(p_plotly, "plot.html", selfcontained = TRUE)


library(ggplot2)

pp <- ggplot(common_df_long) +
  aes(x = r, y = y, fill = variable, colour = variable) +
  geom_point(size = 2L,
             shape = "circle small") +
  scale_fill_manual(values = c(r_AOO = "#FF004D", r_rel_occ_Ncells = "#0C7DFF"
  )) +
  scale_color_manual(values = c(r_AOO = "#FF004D", r_rel_occ_Ncells = "#0C7DFF")) +
  theme_minimal()

pp
saveWidget(pp, here("figures/1_data/corr_plot_rel_abs_AOO.html"), selfcontained = TRUE)

common_df_long %>% filter(strong_cor == 1) %>% print(n = 45)

kableExtra::kable(common_df_long)
kableExtra::kable(common_df)
