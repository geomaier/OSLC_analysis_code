# Plot the aggregate BRAF-status counts for submitted Figure 6.2.
# Set OSLC_BRAF_STATUS_FILE to BRAF_STATUS_TISSUE.csv and optionally OSLC_BRAF_FIGURE_OUTPUT.

library(ggplot2)
library(dplyr)
library(readr)
library(grid)

# Use Arial on Windows
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

# Input and separate output folder
infile <- Sys.getenv("OSLC_BRAF_STATUS_FILE")
if (!nzchar(infile) || !file.exists(infile)) {
  stop("Set OSLC_BRAF_STATUS_FILE to the existing BRAF_STATUS_TISSUE.csv file.",
       call. = FALSE)
}
output_parent <- Sys.getenv("OSLC_BRAF_FIGURE_OUTPUT", dirname(infile))
outdir <- file.path(output_parent, "Figure_6_2_Antrum",
                    paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (file.exists(outdir) || !dir.create(outdir, recursive = TRUE)) {
  stop("Could not create a new output directory: ", outdir, call. = FALSE)
}

# Paper palette
paper_cols <- c(
  "SSL_HP"       = "#FF0000",
  "H_ADJ_SSL_HP" = "#E9C61D",
  "H_COL"        = "#1313EC",
  "CA"           = "#825CA6DE",
  "ANTRUM"       = "#F155D0"
)

# Font settings
font_family <- "Arial"
font_pt <- 7
x_font_pt <- 6
pt_to_gg <- function(x) x / 2.845276

# Read data
df <- read_csv(infile, show_col_types = FALSE)
colnames(df) <- c("Sample", "V600E", "WT")

df <- df %>%
  mutate(
    Sample = toupper(trimws(Sample)),
    V600E  = as.numeric(V600E),
    WT     = as.numeric(WT)
  )

# Extract SSL and HP separately
ssl <- df %>% filter(Sample == "SSL")
hp  <- df %>% filter(Sample == "HP")

ssl_pos   <- ssl$V600E
hp_pos    <- hp$V600E
ssl_total <- ssl$V600E + ssl$WT
hp_total  <- hp$V600E + hp$WT

# Build plotting table
# Combine SSL + HP into SSL_HP
# Label both ANTRUM and legacy STOMACH entries as ANTRUM
plot_df <- df %>%
  mutate(
    Group = case_when(
      Sample %in% c("SSL", "HP")        ~ "SSL_HP",
      Sample == "H_ADJ_SSL_HP"          ~ "H_ADJ_SSL_HP",
      Sample == "H_COL"                 ~ "H_COL",
      Sample == "CA"                    ~ "CA",
      Sample %in% c("ANTRUM", "STOMACH") ~ "ANTRUM",
      TRUE                              ~ Sample
    )
  ) %>%
  group_by(Group) %>%
  summarise(
    V600E = sum(V600E, na.rm = TRUE),
    WT    = sum(WT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total   = V600E + WT,
    Percent = ifelse(Total > 0, round(100 * V600E / Total), 0),
    TopLabel = ifelse(
      V600E > 0,
      paste0(V600E, "/", Total, "\n", Percent, "%"),
      paste0(V600E, "/", Total)
    ),
    LabelY = ifelse(V600E > 0, V600E + 0.6, 0.35)
  )

# Group order for the submitted figure
plot_df$Group <- factor(
  plot_df$Group,
  levels = c("SSL_HP", "H_ADJ_SSL_HP", "H_COL", "CA", "ANTRUM")
)

plot_df <- plot_df %>% arrange(Group)

# Inside labels for SSL and HP within the SSL_HP bar
# Show counts only, no percentages
inside_labels <- data.frame(
  x = 1,
  y = c(
    ssl_pos / 2,
    ssl_pos + hp_pos / 2
  ),
  label = c(
    paste0("SSL\n", ssl_pos, "/", ssl_total),
    paste0("HP\n", hp_pos, "/", hp_total)
  )
)

# White divider between SSL and HP sections
divider_df <- data.frame(
  x = 1 - 0.26,
  xend = 1 + 0.26,
  y = ssl_pos,
  yend = ssl_pos
)

# Plot
p <- ggplot(plot_df, aes(x = Group, y = V600E, fill = Group)) +

  geom_col(
    width = 0.52,
    colour = "black",
    linewidth = 0.35
  ) +

  geom_segment(
    data = divider_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 0.5
  ) +

  geom_text(
    aes(y = LabelY, label = TopLabel),
    family = font_family,
    size = pt_to_gg(font_pt),
    colour = "black",
    vjust = 0,
    lineheight = 0.9
  ) +

  geom_text(
    data = inside_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = font_family,
    size = pt_to_gg(font_pt),
    colour = "white",
    lineheight = 0.9
  ) +

  scale_fill_manual(values = paper_cols) +

  scale_x_discrete(
    labels = c(
      "SSL_HP"       = "SSL_HP",
      "H_ADJ_SSL_HP" = "H_ADJ_SSL_HP",
      "H_COL"        = "H_COL",
      "CA"           = "CA",
      "ANTRUM"       = "ANTRUM"
    ),
    expand = expansion(add = c(0.65, 0.35))
  ) +

  scale_y_continuous(
    breaks = c(0, 10, 20, 30),
    expand = expansion(mult = c(0, 0))
  ) +

  coord_cartesian(
    ylim = c(0, 34),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = expression(italic(BRAF) ~ " V600E+ samples (n)")
  ) +

  theme_classic(base_family = font_family, base_size = font_pt) +
  theme(
    legend.position = "none",

    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),

    axis.text.x = element_text(
      family = font_family,
      size = x_font_pt,
      colour = "black",
      angle = 20,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 2)
    ),

    axis.text.y = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),

    axis.title.y = element_text(
      family = font_family,
      size = font_pt,
      colour = "black",
      margin = margin(r = 4)
    ),

    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.ticks.length = unit(1.2, "mm"),

    plot.margin = margin(4, 5, 3, 3, "mm")
  )

print(p)

# Save outputs
ggsave(
  filename = file.path(outdir, "Panel_C_BRAF_STATUS_TISSUE_ANTRUM.pdf"),
  plot = p,
  width = 85,
  height = 58,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Panel_C_BRAF_STATUS_TISSUE_ANTRUM.jpg"),
  plot = p,
  width = 85,
  height = 58,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# Show the saved files and open their folder on Windows.
message("Plots saved in: ", normalizePath(outdir, winslash = "/"))
if (.Platform$OS.type == "windows") {
  shell.exec(normalizePath(outdir))
}
