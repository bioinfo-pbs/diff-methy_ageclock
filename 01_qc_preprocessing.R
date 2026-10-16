# 01_qc_preprocessing.R

idat_dir <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Data Bu Tiwi/209648820170/"
metadata_file <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Hasil/metadata.csv"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("sesame", quietly = TRUE))
  BiocManager::install("sesame", ask = FALSE, update = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

library(sesame)
library(ggplot2)
library(ggrepel)

# Palet Pandu Biosains
pandu_red    <- "#D24645"
pandu_orange <- "#FE9232"
pandu_navy   <- "#00002E"
pandu_gray   <- "#4F4F67"
pandu_black  <- "#1F1A17"
pandu_white  <- "#FFFFFF"

# Persiapan SeSAMe
dir.create(
  path.expand("~/.cache/R/ExperimentHub"),
  recursive = TRUE,
  showWarnings = FALSE
)

sesameDataCache()

idats <- searchIDATprefixes(idat_dir)

# Generate probe-level beta value 
beta_probe <- openSesame(
  idats,
  prep = "QCDPB",
  collapseToPfx = FALSE
)

colnames(beta_probe) <- basename(idats)

out_dir <- file.path(dirname(metadata_file), "1_qc_preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  beta_probe,
  file.path(out_dir, "beta_probe.rds")
)

# Generate CpG-level beta value
beta <- betasCollapseToPfx(beta_probe)
colnames(beta) <- basename(idats)

saveRDS(
  beta,
  file.path(out_dir, "beta_cpg.rds")
)

# Metadata
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (anyDuplicated(metadata$sample))
  stop("Ada ID sample yang duplikat pada metadata.")

i <- match(colnames(beta), metadata$sample)

if (anyNA(i))
  stop("Ada sampel yang tidak ditemukan di metadata.")

metadata <- metadata[i, , drop = FALSE]

sample_label <- metadata$sample

if ("sample_name(opsional)" %in% names(metadata)) {
  x <- trimws(metadata[["sample_name(opsional)"]])
  use_name <- !is.na(x) & nzchar(x)
  sample_label[use_name] <- x[use_name]
}

sample_colors <- setNames(
  grDevices::colorRampPalette(
    c(
      pandu_navy,
      pandu_gray,
      pandu_red,
      pandu_orange
    )
  )(ncol(beta)),
  sample_label
)

# Sample QC
pvals <- openSesame(
  idats,
  prep = "QCD",
  func = pOOBAH,
  return.pval = TRUE
)

call_rate <- colMeans(
  !is.na(pvals) & pvals <= 0.05
) * 100

qc <- data.frame(
  sample = sample_label,
  call_rate = call_rate,
  valid_cpg = colSums(!is.na(beta)),
  total_cpg = nrow(beta)
)

qc$percent <- qc$valid_cpg / qc$total_cpg * 100
qc$row <- rev(seq_len(nrow(qc)))

qc$call_rate_label <- sprintf("%.2f", qc$call_rate)

qc$valid_label <- sprintf(
  "%s / %s",
  format(qc$valid_cpg, big.mark = ",", trim = TRUE),
  format(qc$total_cpg, big.mark = ",", trim = TRUE)
)

qc$percent_label <- sprintf("%.2f", qc$percent)

p <- ggplot(qc) +
  geom_hline(
    yintercept = qc$row - 0.5,
    colour = grDevices::adjustcolor(pandu_gray, alpha.f = 0.18),
    linewidth = 0.35
  ) +
  geom_text(
    aes(0.03, row, label = sample),
    hjust = 0,
    size = 4,
    colour = pandu_black
  ) +
  geom_text(
    aes(0.45, row, label = call_rate_label),
    size = 4,
    colour = pandu_gray
  ) +
  geom_text(
    aes(0.70, row, label = valid_label),
    size = 4,
    colour = pandu_gray
  ) +
  geom_text(
    aes(0.92, row, label = percent_label),
    hjust = 0.5,
    size = 4,
    fontface = "bold",
    colour = pandu_gray
  ) +
  annotate(
    "text",
    x = 0.03,
    y = max(qc$row) + 1,
    label = "Sample",
    hjust = 0,
    fontface = "bold",
    size = 4,
    colour = pandu_navy
  ) +
  annotate(
    "text",
    x = 0.45,
    y = max(qc$row) + 1,
    label = "Call Rate (%)",
    fontface = "bold",
    size = 4,
    colour = pandu_navy
  ) +
  annotate(
    "text",
    x = 0.70,
    y = max(qc$row) + 1,
    label = "Valid CpGs",
    fontface = "bold",
    size = 4,
    colour = pandu_navy
  ) +
  annotate(
    "text",
    x = 0.92,
    y = max(qc$row) + 1,
    label = "Valid CpGs (%)",
    hjust = 0.5,
    fontface = "bold",
    size = 4,
    colour = pandu_navy
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0.5, max(qc$row) + 1.4)
  ) +
  labs(title = "Sample QC") +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13.2,
      colour = pandu_navy,
      margin = margin(b = 8)
    ),
    plot.margin = margin(8, 10, 8, 10)
  )

ggsave(
  file.path(out_dir, "01_sample_qc.png"),
  p,
  width = 10,
  height = 3.6,
  dpi = 300,
  bg = pandu_white
)


# Beta-value distribution
density_data <- do.call(
  rbind,
  lapply(seq_len(ncol(beta)), function(i) {
    d <- density(
      beta[, i],
      na.rm = TRUE,
      from = 0,
      to = 1
    )
    
    data.frame(
      beta = d$x,
      density = d$y,
      sample = sample_label[i]
    )
  })
)

density_data$sample <- factor(
  density_data$sample,
  levels = sample_label
)

p <- ggplot(
  density_data,
  aes(beta, density, colour = sample)
) +
  geom_line(
    linewidth = 0.9,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = sample_colors
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Beta value",
    y = "Density",
    title = "Beta-value Distribution",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(
      face = "bold",
      colour = pandu_navy
    )
  ) +
  guides(
    colour = guide_legend(
      ncol = 4,
      byrow = TRUE
    )
  )

ggsave(
  file.path(out_dir, "02_beta_distribution.png"),
  p,
  width = 10,
  height = 6.2,
  dpi = 300,
  bg = pandu_white
)


# PCA. 10.000 CpG paling bervariasi
beta_pca <- beta[
  rowSums(is.na(beta)) == 0,
  ,
  drop = FALSE
]

variance <- apply(beta_pca, 1, var)

top <- order(
  variance,
  decreasing = TRUE
)[seq_len(min(10000, length(variance)))]

pca <- prcomp(
  t(beta_pca[top, , drop = FALSE]),
  scale. = FALSE
)

explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- data.frame(
  sample = sample_label,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

p <- ggplot(
  pca_data,
  aes(PC1, PC2, label = sample)
) +
  geom_point(
    size = 4,
    colour = pandu_red
  ) +
  geom_text_repel(
    size = 4,
    colour = pandu_navy,
    box.padding = 0.6,
    point.padding = 0.45,
    segment.colour = pandu_gray,
    segment.linewidth = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", explained[1]),
    y = sprintf("PC2 (%.1f%%)", explained[2]),
    title = "Principal Component Analysis"
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.18)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.18)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      colour = pandu_navy
    ),
    plot.margin = margin(15, 30, 15, 20)
  )

ggsave(
  file.path(out_dir, "03_pca.png"),
  p,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = pandu_white
)

# Software Versions
software_dir <- file.path(dirname(metadata_file), "software_versions")
dir.create(software_dir, recursive = TRUE, showWarnings = FALSE)

software_versions <- data.frame(
  Software = c("R", "Bioconductor", "SeSAMe", "sesameData", "ggplot2", "ggrepel"),
  Version = c(
    as.character(getRversion()),
    as.character(BiocManager::version()),
    as.character(packageVersion("sesame")),
    as.character(packageVersion("sesameData")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("ggrepel"))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  software_versions,
  file.path(software_dir, "software_versions.csv"),
  row.names = FALSE
)

# Ringkasan akhir Console
cat("Jumlah sampel:", ncol(beta), "\n")
cat("Jumlah CpG:", nrow(beta), "\n")
cat("Selesai.\n")