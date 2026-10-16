# 02_differential_methylation.R

beta_probe_file <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Hasil/1_qc_preprocessing/beta_probe.rds"
beta_cpg_file <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Hasil/1_qc_preprocessing/beta_cpg.rds"
metadata_file <- "/home/farah/Fikri/6_Array_Methy_ageclock_Stunting/Hasil/metadata.csv"

reference_group  <- "Control"
comparison_group <- "Case"

fdr_cutoff <- 0.05

# Set Bioconductor agar jawaban otomatis yes tanpa tunggu konfirmasi
options(
  ANNOTATION_HUB_ASK = FALSE,
  EXPERIMENT_HUB_ASK = FALSE
)

if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", ask = FALSE, update = FALSE)

if (!requireNamespace("DMRcate", quietly = TRUE))
  BiocManager::install("DMRcate", ask = FALSE, update = FALSE)

if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38", ask = FALSE, update = FALSE)

library(limma)
library(DMRcate)
library(ggplot2)

# Palet Pandu Biosains
pandu_red    <- "#D24645"
pandu_orange <- "#FE9232"
pandu_navy   <- "#00002E"
pandu_gray   <- "#4F4F67"
pandu_white  <- "#FFFFFF"

# Persiapan IlluminaHumanMethylationEPICv2anno.20a1.hg38
dir.create(
  path.expand("~/.cache/R/AnnotationHub"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Data hasil preprocessing Script 1
beta_probe <- readRDS(beta_probe_file)
beta <- readRDS(beta_cpg_file)

# Metadata
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"sex" %in% names(metadata))
  stop("Kolom sex wajib tersedia di metadata.")

sex <- toupper(trimws(as.character(metadata$sex)))

if (anyNA(sex) || any(!sex %in% c("M", "F")))
  stop("Kolom sex wajib diisi M atau F untuk seluruh sampel.")

metadata <- metadata[
  match(colnames(beta), metadata$sample),
  ,
  drop = FALSE
]

use_sample <- metadata$group %in% c(reference_group, comparison_group)

metadata <- metadata[use_sample, , drop = FALSE]
beta <- beta[, use_sample, drop = FALSE]
beta_probe <- beta_probe[, use_sample, drop = FALSE]

if (
  sum(metadata$group == reference_group) < 2 ||
  sum(metadata$group == comparison_group) < 2
)
  stop("Masing-masing group harus memiliki minimal 2 sampel.")

if (!is.numeric(metadata$age) || anyNA(metadata$age))
  stop("Kolom age wajib berisi angka usia untuk seluruh sampel.")

metadata$group <- factor(
  metadata$group,
  levels = c(reference_group, comparison_group)
)

# Menentukan apakah Age bisa dimasukkan sebagai kovariat dalam model statistik
model_messages <- character()
model_terms <- "group"
candidate <- model.matrix(~ group + age, data = metadata)

if (qr(candidate)$rank == ncol(candidate) &&
    nrow(candidate) > ncol(candidate)) {
  model_terms <- c(model_terms, "age")
} else {
  model_messages <- c(model_messages, "Age tidak dimasukkan: tidak dapat diestimasi bersama group.")
}

# Menentukan apakah Jenis Kelamin bisa dimasukkan sebagai kovariat dalam model statistik
metadata$sex <- factor(toupper(trimws(as.character(metadata$sex))))

if (nlevels(metadata$sex) > 1) {
  candidate <- model.matrix(
    reformulate(c(model_terms, "sex")),
    data = metadata
  )
  
  if (qr(candidate)$rank == ncol(candidate) &&
      nrow(candidate) > ncol(candidate)) {
    model_terms <- c(model_terms, "sex")
  } else {
    model_messages <- c(model_messages, "Sex tidak dimasukkan: model tidak dapat mengestimasi kovariat ini.")
  }
} else {
  model_messages <- c(model_messages, "Sex tidak dimasukkan: seluruh sampel memiliki jenis kelamin yang sama.")
}

# Menentukan apakah Batch Ekstraksi bisa dimasukkan sebagai kovariat dalam model statistik
if ("batch(opsional)" %in% names(metadata)) {
  batch <- trimws(as.character(metadata[["batch(opsional)"]]))
  
  if (all(!is.na(batch) & nzchar(batch))) {
    metadata$batch <- factor(batch)
    
    if (nlevels(metadata$batch) > 1) {
      candidate <- model.matrix(
        reformulate(c(model_terms, "batch")),
        data = metadata
      )
      
      if (qr(candidate)$rank == ncol(candidate) &&
          nrow(candidate) > ncol(candidate)) {
        model_terms <- c(model_terms, "batch")
      } else {
        model_messages <- c(model_messages, "Extraction batch tidak dimasukkan: model tidak dapat mengestimasi kovariat ini.")
      }
    }
  }
}

design <- model.matrix(
  reformulate(model_terms),
  data = metadata
)

# DMP
m_value <- log2(beta / (1 - beta))
m_value[!is.finite(m_value)] <- NA_real_

reference_sample <- metadata$group == reference_group
comparison_sample <- metadata$group == comparison_group

keep <-
  rowSums(!is.na(m_value[, reference_sample, drop = FALSE])) >= 2 &
  rowSums(!is.na(m_value[, comparison_sample, drop = FALSE])) >= 2

valid_patterns <- apply(
  1L * !is.na(m_value[keep, , drop = FALSE]),
  1,
  paste0,
  collapse = ""
)

unique_patterns <- unique(valid_patterns)

pattern_ok <- vapply(unique_patterns, function(pattern) {
  valid <- strsplit(pattern, "", fixed = TRUE)[[1]] == "1"
  x <- design[valid, , drop = FALSE]
  nrow(x) > ncol(x) && qr(x)$rank == ncol(x)
}, logical(1))

keep[keep] <- pattern_ok[match(valid_patterns, unique_patterns)]

m_value <- m_value[keep, , drop = FALSE]
beta_dmp <- beta[keep, , drop = FALSE]

fit <- eBayes(lmFit(m_value, design))

dmp_raw <- topTable(
  fit,
  coef = 2,
  number = Inf,
  sort.by = "P"
)

mean_reference <- rowMeans(
  beta_dmp[, reference_sample, drop = FALSE],
  na.rm = TRUE
)

mean_comparison <- rowMeans(
  beta_dmp[, comparison_sample, drop = FALSE],
  na.rm = TRUE
)

dmp <- data.frame(
  CpG = rownames(dmp_raw),
  mean_beta_reference = mean_reference[rownames(dmp_raw)],
  mean_beta_comparison = mean_comparison[rownames(dmp_raw)],
  delta_beta =
    mean_comparison[rownames(dmp_raw)] -
    mean_reference[rownames(dmp_raw)],
  logFC_M = dmp_raw$logFC,
  p_value = dmp_raw$P.Value,
  FDR = dmp_raw$adj.P.Val
)

# Anotasi genom dan gen EPIC v2 (hg38) untuk DMP
anno_env <- new.env(parent = baseenv())
data(
  "Locations",
  package = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  envir = anno_env
)
data(
  "Other",
  package = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  envir = anno_env
)

anno_locations <- as.data.frame(anno_env$Locations)
anno_other <- as.data.frame(anno_env$Other)

if (!all(c("UCSC_RefGene_Name", "UCSC_RefGene_Group") %in% names(anno_other)))
  stop("Kolom anotasi gen/region RefGene tidak tersedia dalam package EPIC v2.")

location_idx <- match(
  dmp$CpG,
  sub("_.*$", "", rownames(anno_locations))
)
gene_idx <- match(
  dmp$CpG,
  sub("_.*$", "", rownames(anno_other))
)

collapse_gene_annotations <- function(gene_text, region_text, cpg) {
  if (is.na(gene_text) || !nzchar(trimws(gene_text)))
    return(c(associated.genes = NA_character_, gene.regions = NA_character_))
  
  genes <- trimws(strsplit(gene_text, ";", fixed = TRUE)[[1]])
  genes_unique <- unique(genes[nzchar(genes)])
  
  if (length(genes_unique) == 0)
    return(c(associated.genes = NA_character_, gene.regions = NA_character_))
  
  associated <- paste(genes_unique, collapse = ";")
  
  if (is.na(region_text) || !nzchar(trimws(region_text)))
    return(c(associated.genes = associated, gene.regions = NA_character_))
  
  regions <- trimws(strsplit(region_text, ";", fixed = TRUE)[[1]])
  
  if (length(genes) != length(regions))
    stop("Jumlah entri Gene dan Gene Region tidak cocok untuk CpG: ", cpg)
  
  region_per_gene <- vapply(genes_unique, function(g) {
    x <- unique(regions[genes == g])
    x <- x[nzchar(x)]
    if (length(x) == 0)
      return(NA_character_)
    region_label <- paste(x, collapse = "|")
    if (length(genes_unique) == 1)
      region_label
    else
      paste0(g, ":", region_label)
  }, character(1))
  
  region_per_gene <- region_per_gene[!is.na(region_per_gene)]
  region_result <- if (length(region_per_gene) > 0) {
    paste(region_per_gene, collapse = ";")
  } else {
    NA_character_
  }
  
  c(associated.genes = associated, gene.regions = region_result)
}

gene_annotation <- vapply(seq_along(gene_idx), function(i) {
  collapse_gene_annotations(
    as.character(anno_other$UCSC_RefGene_Name[gene_idx[i]]),
    as.character(anno_other$UCSC_RefGene_Group[gene_idx[i]]),
    dmp$CpG[i]
  )
}, c(associated.genes = NA_character_, gene.regions = NA_character_))

dmp <- data.frame(
  CpG = dmp$CpG,
  seqnames = anno_locations$chr[location_idx],
  position = anno_locations$pos[location_idx],
  dmp[, -1, drop = FALSE],
  associated.genes = gene_annotation["associated.genes", ],
  gene.regions = gene_annotation["gene.regions", ]
)

# DMP csv file
out_dir <- file.path(dirname(metadata_file), "2_diff_methylation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  dmp,
  file.path(out_dir, "04_dmp_results.csv"),
  row.names = FALSE
)

# DMP volcano plot
significant <- !is.na(dmp$FDR) & dmp$FDR < fdr_cutoff

higher_label <- paste("Higher in", comparison_group)
lower_label <- paste("Lower in", comparison_group)

dmp$status <- "Not significant"
dmp$status[significant & dmp$logFC_M > 0] <- higher_label
dmp$status[significant & dmp$logFC_M < 0] <- lower_label

dmp$neg_log10_fdr <- -log10(
  pmax(dmp$FDR, .Machine$double.xmin)
)

volcano_colors <- c(
  "Not significant" =
    grDevices::adjustcolor(pandu_gray, alpha.f = 0.35)
)

volcano_colors[higher_label] <- pandu_red
volcano_colors[lower_label] <- pandu_orange

p <- ggplot(
  dmp,
  aes(logFC_M, neg_log10_fdr, colour = status)
) +
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    colour = pandu_gray
  ) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed",
    linewidth = 0.5,
    colour = pandu_gray
  ) +
  scale_colour_manual(
    values = volcano_colors,
    breaks = c(lower_label, "Not significant", higher_label)
  ) +
  labs(
    x = "logFC M",
    y = "-log10(FDR)",
    title = "Differentially Methylated Positions",
    subtitle = paste(comparison_group, "vs", reference_group),
    colour = NULL
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
  file.path(out_dir, "05_dmp_volcano.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = pandu_white
)

# DMR
dmr_annotation <- cpg.annotate(
  datatype = "array",
  object = beta_probe,
  what = "Beta",
  arraytype = "EPICv2",
  epicv2Remap = TRUE,
  epicv2Filter = "mean",
  analysis.type = "differential",
  design = design,
  coef = 2,
  fdr = fdr_cutoff
)

no_sig_dmr <- sum(dmr_annotation@ranges$is.sig, na.rm = TRUE) == 0

if (no_sig_dmr) {
  
  dmr <- data.frame(
    DMR = character(),
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    no.cpgs = integer(),
    min_smoothed_fdr = numeric(),
    Stouffer = numeric(),
    HMFDR = numeric(),
    Fisher = numeric(),
    maxdiff = numeric(),
    meandiff = numeric(),
    overlapping.genes = character()
  )
  
} else {
  
  dmr_fit <- dmrcate(
    dmr_annotation,
    lambda = 1000,
    C = 2
  )
  
  dmr <- as.data.frame(
    extractRanges(
      dmr_fit,
      genome = "hg38"
    )
  )
  
  dmr$DMR <- paste0("DMR", seq_len(nrow(dmr)))
  
  dmr <- dmr[
    ,
    c("DMR", setdiff(colnames(dmr), "DMR")),
    drop = FALSE
  ]
}

write.csv(
  dmr,
  file.path(out_dir, "06_dmr_results.csv"),
  row.names = FALSE
)

# Info arah comparison untuk dipakai di HTML report
write.csv(
  data.frame(
    reference_group = reference_group,
    comparison_group = comparison_group
  ),
  file.path(out_dir, "comparison_info.csv"),
  row.names = FALSE
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
  Software = c("limma", "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "DMRcate"),
  Version = c(
    as.character(packageVersion("limma")),
    as.character(packageVersion("IlluminaHumanMethylationEPICv2anno.20a1.hg38")),
    as.character(packageVersion("DMRcate"))
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
cat("\n")
if (length(model_messages) > 0)
  cat(paste(model_messages, collapse = "\n"), "\n")
cat("Model:", paste(model_terms, collapse = " + "), "\n")
cat("Comparison:", comparison_group, "vs", reference_group, "\n")
cat("DMP FDR <", fdr_cutoff, ":", sum(significant), "\n")
if (no_sig_dmr)
  cat("Tidak ditemukan DMR pada FDR ", fdr_cutoff, ".\n", sep = "")
cat("DMR:", nrow(dmr), "\n")
cat("Selesai.\n")