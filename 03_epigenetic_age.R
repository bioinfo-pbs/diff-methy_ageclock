# 03_epigenetic_age.R

metadata_file <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Hasil/metadata.csv"

if (!requireNamespace("methylclock", quietly = TRUE))
  BiocManager::install("methylclock", ask = FALSE, update = FALSE)

library(methylclock)
library(ggplot2)
library(ggrepel)

# Palet Pandu Biosains
pandu_red    <- "#D24645"
pandu_orange <- "#FE9232"
pandu_navy   <- "#00002E"
pandu_gray   <- "#4F4F67"
pandu_white  <- "#FFFFFF"

# Data hasil preprocessing Script 1
beta <- readRDS(file.path(dirname(metadata_file), "1_qc_preprocessing", "beta_cpg.rds"))

# Metadata
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[
  match(colnames(beta), metadata$sample),
  ,
  drop = FALSE
]

sample_label <- metadata$sample

if ("sample_name(opsional)" %in% names(metadata)) {
  x <- trimws(metadata[["sample_name(opsional)"]])
  use_name <- !is.na(x) & nzchar(x)
  sample_label[use_name] <- x[use_name]
}

# CpG Horvath Pan-Tissue Methylation Clock
horvath_coef <- methylclockData::get_coefHorvath()

horvath_cpg <- setdiff(
  as.character(horvath_coef$CpGmarker),
  "(Intercept)"
)

beta_clock <- beta[
  intersect(horvath_cpg, rownames(beta)),
  ,
  drop = FALSE
]

no_signal <- rownames(beta_clock)[
  rowSums(!is.na(beta_clock)) == 0
]

beta_clock <- beta_clock[
  !rownames(beta_clock) %in% no_signal,
  ,
  drop = FALSE
]

missing_cpg <- setdiff(
  horvath_cpg,
  rownames(beta_clock)
)

available_cpg <- nrow(beta_clock)
valid_per_sample <- colSums(!is.na(beta_clock))

if (length(missing_cpg) > 0) {
  
  reference_beta <- horvath_coef$medianByCpG[
    match(missing_cpg, horvath_coef$CpGmarker)
  ]
  
  reference_matrix <- matrix(
    reference_beta,
    nrow = length(missing_cpg),
    ncol = ncol(beta_clock),
    dimnames = list(missing_cpg, colnames(beta_clock))
  )
  
  beta_clock <- rbind(beta_clock, reference_matrix)
}

beta_clock <- beta_clock[
  horvath_cpg,
  ,
  drop = FALSE
]

clock_result <- DNAmAge(
  beta_clock,
  clocks = "Horvath",
  fastImp = FALSE,
  normalize = FALSE,
  cell.count = FALSE
)

dnam_age <- clock_result$Horvath

# CpG Horvath Pan-Tissue Age Difference
results <- data.frame(
  sample = sample_label,
  chronological_age = metadata$age,
  dnam_age = dnam_age
)

results$age_difference <-
  results$dnam_age - results$chronological_age

clock_table <- data.frame(
  sample = sample_label,
  actual_age_years = results$chronological_age,
  dnam_age_years = results$dnam_age,
  age_difference = results$age_difference
)

# Methylation Clock + Age Difference csv file
out_dir <- file.path(dirname(metadata_file), "3_epigenetic_age")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  clock_table,
  file.path(out_dir, "07_MethylationClock_AgeDifference.csv"),
  row.names = FALSE
)

# Methylation Clock Plot
age_limits <- range(
  c(results$chronological_age, results$dnam_age),
  na.rm = TRUE
)

p <- ggplot(
  results,
  aes(
    chronological_age,
    dnam_age,
    label = sample
  )
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
  scale_x_continuous(
    breaks = seq(
      floor(age_limits[1]),
      ceiling(age_limits[2]),
      by = 1
    ),
    labels = as.character
  ) +
  coord_fixed(
    ratio = 1,
    xlim = age_limits,
    ylim = age_limits
  ) +
  labs(
    x = "Actual Age (years)",
    y = "DNAm Age (years)",
    title = "DNA Methylation Clock",
    subtitle = "Horvath Pan-Tissue Clock"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      colour = pandu_navy
    ),
    plot.subtitle = element_text(
      colour = pandu_gray
    )
  )

ggsave(
  file.path(out_dir, "08_methylation_clock.png"),
  p,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = pandu_white
)


# Age Difference Plot
results$direction <- ifelse(
  results$age_difference > 0,
  "Older than Actual Age",
  ifelse(
    results$age_difference < 0,
    "Younger than Actual Age",
    "Same as Actual Age"
  )
)

results$sample <- factor(
  results$sample,
  levels = rev(sample_label)
)

p <- ggplot(
  results,
  aes(sample, age_difference, fill = direction)
) +
  geom_hline(
    yintercept = 0,
    colour = pandu_gray,
    linewidth = 0.5
  ) +
  geom_col(
    width = 0.65
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Younger than Actual Age" = pandu_orange,
      "Older than Actual Age" = pandu_red,
      "Same as Actual Age" = pandu_gray
    ),
    breaks = c(
      "Younger than Actual Age",
      "Older than Actual Age",
      "Same as Actual Age"
    )
  ) +
  labs(
    x = NULL,
    y = "DNAm Age - Actual Age (years)",
    title = "Epigenetic Age Difference",
    subtitle = "Horvath Pan-Tissue Clock",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      colour = pandu_navy
    ),
    plot.subtitle = element_text(
      colour = pandu_gray
    )
  )

ggsave(
  file.path(out_dir, "09_age_difference.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = pandu_white
)

# Software Versions
software_dir <- file.path(dirname(metadata_file), "software_versions")
version_file <- file.path(software_dir, "software_versions.csv")

software_versions <- read.csv(
  version_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!identical(names(software_versions), c("Software", "Version")))
  stop("Format software_versions.csv tidak sesuai.")

new_versions <- data.frame(
  Software = c("methylclock", "methylclockData"),
  Version = c(
    as.character(packageVersion("methylclock")),
    as.character(packageVersion("methylclockData"))
  ),
  stringsAsFactors = FALSE
)

software_versions <- software_versions[
  !software_versions$Software %in% new_versions$Software,
  ,
  drop = FALSE
]
software_versions <- rbind(software_versions, new_versions)

write.csv(software_versions, version_file, row.names = FALSE)

# Ringkasan akhir Console
cat("\nCpG Horvath tersedia:", available_cpg, "/", length(horvath_cpg), "\n")
cat("CpG menggunakan referensi Horvath:", length(missing_cpg), "\n")
cat("Kelengkapan CpG Horvath per sampel (sebelum imputasi):\n")
for (i in seq_along(sample_label)) {
  cat(sample_label[i], ":", valid_per_sample[i], "/", length(horvath_cpg), "\n")
}
print(results)
cat("Jumlah sampel:", nrow(results), "\n")
cat("Selesai.\n")