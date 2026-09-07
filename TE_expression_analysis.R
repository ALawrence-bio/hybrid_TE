# =============================================================================
# TE Expression Analysis Pipeline — Ovaries & Testes
# DESeq2 differential expression + chromosome categories + inheritance patterns
#
# Authors: Alex Lemopoulos, Robin Ségurel
# Repository: https://github.com/<your-username>/<your-repo>
#
# USAGE:
#   1. Fill in the CONFIGURATION section below with your own file paths and
#      sample metadata. All paths marked "# <-- EDIT" must be updated.
#   2. Run the entire script.
#   3. Per-tissue results (plots + data) are stored in `results$ova` and
#      `results$tes` and can be accessed after the pipeline completes.
#   4. Panel figures are assembled at the bottom of the script.
#
# INPUT FILES REQUIRED (per tissue):
#   - Collapsed TE count matrix (.tsv), from TElocal
#   - Chromosome category mapping file (.tsv)
#   - Loci-specific count matrix with chromosome category (.tsv)
#   - (Optional) Gene count matrix (.tsv) for gene-level dds objects
#
# OUTPUT:
#   - Differential expression tables (.tab)
#   - Normalised count tables (.csv)
#   - Plots (individual .png files + panel figures)
#   - (Optional) Saved DESeq2 dds objects (.rds)
#
# DEPENDENCIES:
#   Install required packages with:
#   install.packages(c("ggplot2","dplyr","scales","tidyr","stringr",
#                      "readr","patchwork","cowplot","gplots","pheatmap",
#                      "RColorBrewer","magrittr"))
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
# =============================================================================

# -- Libraries ----------------------------------------------------------------
library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(stringr)
library(readr)
library(patchwork)   # for panel figure assembly at the end

# =============================================================================
# CONFIGURATION
# Edit these two blocks — everything else runs automatically
# =============================================================================
# Set your working directory — all output files will be saved here
setwd("/path/to/your/working/directory")  # <-- EDIT

tissue_configs <- list(

  ova = list(
    label        = "Ovaries",
    short        = "ova",

    # Input files                                        # <-- EDIT paths below
    counts_collapsed = "path/to/ovaries_TE_counts_collapsed.tsv",
    chrom_map        = "path/to/ovaries_chromosome_mapping.tsv",
    counts_loci      = "path/to/ovaries_loci_counts_with_chrom_category.tsv",

    # Output files                                       # <-- EDIT paths below
    out_diff_expr    = "ova_diff_express.tab",
    out_diff_sig     = "ova_diff_express_sig05.tab",

    # Sample metadata — must match column order in count matrix exactly
    # Replace with your own sample names and species assignments
    sample_names   = c("sample1", "sample2", "sample3"),  # <-- EDIT
    sample_species = c(rep("pDom", 1), rep("pIta", 1), rep("pHis", 1))  # <-- EDIT
  ),

  tes = list(
    label        = "Testes",
    short        = "testes",

    # Input files                                        # <-- EDIT paths below
    counts_collapsed = "path/to/testes_TE_counts_collapsed.tsv",
    chrom_map        = "path/to/testes_chromosome_mapping.tsv",
    counts_loci      = "path/to/testes_loci_counts_with_chrom_category.tsv",

    # Output files                                       # <-- EDIT paths below
    out_diff_expr    = "testes_diff_express.tab",
    out_diff_sig     = "testes_diff_express_sig05.tab",

    # Sample metadata — must match column order in count matrix exactly
    # Replace with your own sample names and species assignments
    sample_names   = c("sample1", "sample2", "sample3"),  # <-- EDIT
    sample_species = c(rep("pDom", 1), rep("pIta", 1), rep("pHis", 1))  # <-- EDIT
  )

)

# -- Shared thresholds (apply to both tissues) --------------------------------
log2FC_threshold <- 0.32
padj_threshold   <- 0.05

# -- Shared colour palettes ---------------------------------------------------
cat_levels <- c("A", "W", "Z", "A+W", "A+Z", "W+Z", "A+W+Z")

# Scatter plot: over- and underdominant share the same yellow
inheritance_colors_scatter <- c(
  "Conserved"    = "grey",
  "Additive"     = "purple",
  "Overdominant" = "#FFD700",
  "Underdominant"= "orange",
  "pDom-dominant"= "blue",
  "pHis-dominant"= "red"
)

# Bar chart: underdominant distinguished by orange
inheritance_colors_bar <- c(
  "Conserved"    = "grey",
  "Additive"     = "purple",
  "Overdominant" = "#FFD700",
  "Underdominant"= "orange",
  "pDom-dominant"= "blue",
  "pHis-dominant"= "red"
)

# Factor levels for scatterplot colour legend
inheritance_factor_levels <- c(
  "Conserved","pDom-dominant","Additive","pHis-dominant","Overdominant","Underdominant"
)

# Factor levels for stacked bar chart (bottom-to-top order)
inheritance_factor_levels_bar <- c(
  "Conserved","pDom-dominant","Additive","pHis-dominant","Overdominant","Underdominant"
)

# -- Shared species colour scale ----------------------------------------------
# Locks in consistent colours across both tissues regardless of sample order.
# Adjust hex codes here if you want different colours.
species_colours <- c("pDom" = "blue", "pHis" = "red", "pIta" = "#FFD700")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Normalise chromosome category strings to canonical A/W/Z order
normalize_category <- function(x) {
  parts <- unlist(strsplit(gsub("\\s+", "", x), "\\+"))
  parts <- parts[parts %in% c("A", "W", "Z")]
  paste(parts[order(match(parts, c("A","W","Z")))], collapse = "+")
}

# Build a volcano plot from a DESeq2 results object.
# Significance is encoded as a labelled factor so the legend reads exactly the
# desired text rather than TRUE/FALSE.
make_volcano_plot <- function(res, title) {
  down_trt <- paste0("up in\n",str_sub(title,1,4))
  up_trt <- paste0("up in\n",str_sub(title,9,12))
  # species_names <- c("pDom" = "house",
  #                    "pIta" = "Italian",
  #                    "pHis" = "Spanish")
  #  sp1 <- str_sub(title, 1, 4)
  #  sp2 <- str_sub(title, 9, 12)
  # down_trt <- paste0("up in\n", species_names[sp1])
  # up_trt   <- paste0("up in\n", species_names[sp2])
  
  sig_label     <- "padj < 0.05 & abs(log2FoldChange) > 0.32"
  nonsig_label  <- "Not significant"
  as.data.frame(res) %>%
    mutate(
      TE_id      = rownames(res),
      significant = factor(
        ifelse(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold,
               sig_label, nonsig_label),
        levels = c(nonsig_label, sig_label)
      )
    ) %>%
    na.omit() %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
      values = c("grey", "red"),
      name   = NULL
    ) +
    annotate(geom = "text", x = -Inf, y = Inf, label= down_trt, hjust =-0.1, vjust = 1.1) +
    annotate(geom = "text", x = Inf, y = Inf, label= up_trt, hjust =1.1, vjust = 1.1) +
    ylim(0, 20) + xlim(-8, 8) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
}

# Extract TE_id, LFC, padj, pvalue from a DESeq2 results object
make_results_summary <- function(res) {
  data.frame(
    TE_id          = rownames(res),
    log2FoldChange = res$log2FoldChange,
    padj           = res$padj,
    pvalue         = res$pvalue
  ) %>% na.omit()
}

# IDs to reclassify as Unknown (manually curated — parsed class unreliable)
ids_to_unknown <- c(
  "rnd-5_family-6382",   # was rnd-5_family-6382#DNA/Dada
  "rnd-5_family-4100",   # was rnd-5_family-4100#PLE/Chlamys
  "rnd-5_family-2053",   # was rnd-5_family-2053#SINE/MIR
  "rnd-5_family-32486",  # was rnd-5_family-32486#SINE/ID
  "rnd-5_family-72847",  # was rnd-5_family-72847#SINE/ID
  "rnd-5_family-4428",   # was rnd-5_family-4428#RC/Helitron
  "rnd-5_family-5246",   # was rnd-5_family-5246#RC/Helitron
  "rnd-5_family-8403",   # was rnd-5_family-8403#RC/Helitron
  "rnd-5_family-1398",   # was rnd-5_family-1398#RC/Helitron
  "rnd-5_family-2412"    # was rnd-5_family-2412#RC/Helitron
)

# Derive broad TE class from TE_id string (e.g. "rnd-5_family-775/LINE/CR1")
add_te_class <- function(df) {
  df %>%
    mutate(
      subclass_raw   = sapply(strsplit(TE_id, "/"), `[`, 2),
      superfamily    = sub("-.*", "", sapply(strsplit(TE_id, "/"), `[`, 3)),
      subclass_broad = ifelse(superfamily == "Penelope", "PLE", subclass_raw),
      TE_class       = ifelse(
        subclass_broad %in% c("DNA","LINE","LTR","PLE","RC","SINE","Unknown"),
        subclass_broad, "Unknown"
      )
    )
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# Runs the full pipeline for one tissue and returns all plots + data objects
# =============================================================================

run_te_analysis <- function(cfg) {

  label <- cfg$label
  short <- cfg$short
  cat("\n\n", strrep("=", 70), "\n")
  cat("  Running analysis:", label, "\n")
  cat(strrep("=", 70), "\n\n")

  # ---------------------------------------------------------------------------
  # 1. LOAD COUNT MATRIX & BUILD SAMPLE TABLE
  # ---------------------------------------------------------------------------
  count_mat <- as.matrix(
    read.csv(cfg$counts_collapsed, sep = "\t", header = TRUE, row.names = "TE_id")
  )
  cat("Count matrix dimensions:", nrow(count_mat), "TEs x", ncol(count_mat), "samples\n")

  sample_table <- data.frame(
    species   = factor(cfg$sample_species),
    row.names = cfg$sample_names
  )

  # ---------------------------------------------------------------------------
  # 2. CHROMOSOME CATEGORY BAR CHART (collapsed TEs)
  # ---------------------------------------------------------------------------
  te_map <- read_tsv(cfg$chrom_map, show_col_types = FALSE)
  stopifnot("TE_id" %in% names(te_map), "category" %in% names(te_map))

  te_map <- te_map %>%
    mutate(across(c(TE_id, category), trimws),
           category = vapply(category, normalize_category, character(1)))

  annot_chr <- data.frame(TE_id = trimws(rownames(count_mat)),
                          stringsAsFactors = FALSE) %>%
    left_join(te_map[, c("TE_id", "category")], by = "TE_id")

  cat("TEs in count matrix:", nrow(annot_chr), "\n")
  cat("Mapped to chromosome category:", sum(!is.na(annot_chr$category)), "\n")
  cat("Unmapped:", sum(is.na(annot_chr$category)), "\n")
  print(table(annot_chr$category, useNA = "ifany"))

  counts_by_cat <- annot_chr %>%
    filter(!is.na(category)) %>%
    mutate(category = factor(category, levels = cat_levels)) %>%
    dplyr::count(category) %>%
    complete(category = factor(cat_levels, levels = cat_levels), fill = list(n = 0)) %>%
    arrange(category)

  plot_chrom_bar <- ggplot(counts_by_cat, aes(x = category, y = n)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Chromosome category",
         y = "Number of TEs",
         title = "TEs by chromosome category (A/W/Z)")

  # ---------------------------------------------------------------------------
  # 3. LOCI-SPECIFIC CHROMOSOME CATEGORY PLOTS
  # ---------------------------------------------------------------------------
  df_chrom <- read.delim(cfg$counts_loci, stringsAsFactors = FALSE)

  # Use only sample columns that are present in the loci count file
  loci_sample_cols <- intersect(cfg$sample_names, colnames(df_chrom))

  loci_summary <- df_chrom %>%
    group_by(category) %>%
    summarise(
      loci         = dplyr::n(),
      total_counts = sum(across(all_of(loci_sample_cols))),
      .groups      = "drop"
    )
  print(loci_summary)

  plot_loci_all <- ggplot(df_chrom, aes(x = category)) +
    geom_bar() + theme_bw() +
    labs(x = "Chromosome category", y = "Number of TE loci",
         title = "All TE loci by chromosome category")

  df_chrom_expr <- df_chrom %>%
    filter(rowSums(across(all_of(loci_sample_cols))) > 0)

  plot_loci_expr <- ggplot(df_chrom_expr, aes(x = category)) +
    geom_bar() + theme_bw() +
    labs(x = "Chromosome category", y = "Number of TE loci",
         title = paste0(label, " — expressed TE loci by chromosome category"))

  # ---------------------------------------------------------------------------
  # 4. DESeq2: SETUP, NORMALISATION, QC
  # ---------------------------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = sample_table,
    design    = ~ species
  )

  dds$species <- factor(dds$species, levels = c("pDom","pHis","pIta"))
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- estimateSizeFactors(dds)

  write.csv(data.frame(SizeFactor = sizeFactors(dds)),
            paste0(short, "_size_factors.csv"), row.names = TRUE)

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  write.csv(norm_counts, paste0(short, "_normalized_counts.csv"), row.names = TRUE)

  rld <- rlog(dds, blind = TRUE)

  # PCA
  plot_pca <- DESeq2::plotPCA(rld, intgroup = "species", ntop = 500) +
    scale_color_manual(values = species_colours, name = "Species") +
    ggtitle(paste0("PCA — pDom, pHis, pIta (", label, ", top 500 variable TEs)"))

  # Sample distance heatmap (base graphics — printed directly, not stored as ggplot)
  dist_mat <- as.matrix(dist(t(assay(rld))))
  rownames(dist_mat) <- colnames(dist_mat) <- with(
    colData(dds), paste(species, cfg$sample_names, sep = " : ")
  )
  hc    <- hclust(as.dist(dist_mat))
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(dist_mat,
            Rowv = as.dendrogram(hc), symm = TRUE,
            trace = "none", col = rev(hmcol), margin = c(10, 10),
            main = paste(label, "— sample distances"))

  # ---------------------------------------------------------------------------
  # 5. DESeq2: DIFFERENTIAL EXPRESSION (pDom vs pHis, parent-parent)
  # ---------------------------------------------------------------------------
  dds <- DESeq(dds)
  cat("\nAvailable contrasts:\n")
  print(DESeq2::resultsNames(dds))

  dds$species <- factor(dds$species, levels = c("pDom","pHis","pIta"))
  res_table <- results(dds, contrast = c("species", "pDom", "pHis"))
  res_table <- res_table[order(res_table$padj), ]
  res_table <- cbind(TE_id = rownames(res_table),
                     as.data.frame(res_table, row.names = NULL))

  write.table(res_table, file = cfg$out_diff_expr,
              sep = "\t", quote = FALSE, row.names = FALSE)

  resSig <- subset(res_table, padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold)
  write.table(resSig, file = cfg$out_diff_sig,
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\nDE summary (pDom vs pHis):\n")
  cat("  Upregulated in pDom:", sum(resSig$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("  Downregulated in pDom:", sum(resSig$log2FoldChange < 0, na.rm = TRUE), "\n")

  
  plotMA(dds, ylim = c(-5, 5),
         main = paste0("DEx TEs in pHis vs pDom (", label, ")"), alpha = padj_threshold)

  # Top 50 DE TEs heatmap
  top50 <- res_table$TE_id[
    !is.na(res_table$padj) &
      res_table$padj < padj_threshold &
      abs(res_table$log2FoldChange) > log2FC_threshold
  ][1:50]
  counts_top50 <- counts(dds, normalized = TRUE)[top50, ]
  heatmap.2(counts_top50,
            col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
            Rowv = FALSE, Colv = FALSE, scale = "none",
            dendrogram = "none", trace = "none",
            main = paste(label, "— top 50 DE TEs"))
  write.csv(counts_top50, paste0(short, "_counts_top50.csv"), row.names = TRUE)

  # ---------------------------------------------------------------------------
  # 6. PAIRWISE RESULTS & VOLCANO PLOTS
  # ---------------------------------------------------------------------------
  res_ip <- results(dds, contrast = c("species", "pDom", "pIta"))
                                    # positive = elevated in pIta relative to pDom
  res_ih <- results(dds, contrast = c("species", "pHis", "pIta"))
                                    #positive = elevated in pIta relative to pHis
  res_ph <- results(dds, contrast = c("species", "pDom", "pHis"))
                                    # positive = elevated in pDom relative to pHis
  cat("\nDE summary (pIta vs pHis):\n")
  cat("  Upregulated in pIta:", sum(res_ih$log2FoldChange > 0&res_ih$padj < 0.05, na.rm = TRUE), "\n")
  cat("  Downregulated in pIta:", sum(res_ih$log2FoldChange < 0 & res_ih$padj < 0.05, na.rm = TRUE), "\n")
  
  write.csv(make_results_summary(res_ip),
            paste0(short, "_summary_pIta_vs_pDom.csv"), row.names = FALSE)
  write.csv(make_results_summary(res_ih),
            paste0(short, "_summary_pIta_vs_pHis.csv"), row.names = FALSE)
  write.csv(make_results_summary(res_ph),
            paste0(short, "_summary_pDom_vs_pHis.csv"), row.names = FALSE)

  plot_volcano_ip <- make_volcano_plot(res_ip, paste0("pDom vs pIta (", label, ")"))
  plot_volcano_ih <- make_volcano_plot(res_ih, paste0("pHis vs pIta (", label, ")"))
  plot_volcano_ph <- make_volcano_plot(res_ph, paste0("pHis vs pDom (", label, ")"))

  # ---------------------------------------------------------------------------
  # 7. INHERITANCE CLASSIFICATION
  # ---------------------------------------------------------------------------
  df <- data.frame(
    TE_id   = rownames(dds),
    LFC_ih  = res_ih$log2FoldChange,
    LFC_ip  = res_ip$log2FoldChange,
    LFC_ph = res_ph$log2FoldChange,
    padj_ih = res_ih$padj,
    padj_ip = res_ip$padj,
    padj_ph = res_ph$padj
  ) %>% na.omit()

  #df$significant <- df$padj_ip < padj_threshold | df$padj_ih < padj_threshold
df$DE_ih <- df$padj_ih < padj_threshold &  abs(df$LFC_ih) > log2FC_threshold
df$DE_ip <- df$padj_ip < padj_threshold &  abs(df$LFC_ip) > log2FC_threshold
df$DE_ph <- df$padj_ph < padj_threshold &  abs(df$LFC_ph) > log2FC_threshold

  
# #  df <- df %>%
#     mutate(inheritance = case_when(
#       # Additive: hybrid is intermediate (both LFCs within threshold)
#       significant &
#         LFC_ip >= -log2FC_threshold & LFC_ip <= log2FC_threshold &
#         LFC_ih >= -log2FC_threshold & LFC_ih <= log2FC_threshold  ~ "Additive",
#       # Overdominant: hybrid significantly exceeds both parents
#       significant &
#         LFC_ip >  log2FC_threshold  & LFC_ih >  log2FC_threshold  ~ "Overdominant",
#       # Underdominant: hybrid is significantly below both parents
#       significant &
#         LFC_ip < -log2FC_threshold  & LFC_ih < -log2FC_threshold  ~ "Underdominant",
#       # pDom-dominant: hybrid is close to pDom (LFC_ip within threshold)
#       # but significantly different from pHis in either direction (abs LFC_ih)
#       significant &
#         abs(LFC_ip) <= log2FC_threshold &
#         abs(LFC_ih) >  log2FC_threshold                           ~ "pDom-dominant",
#       # pHis-dominant: hybrid is close to pHis (LFC_ih within threshold)
#       # but significantly different from pDom in either direction (abs LFC_ip)
#       significant &
#         abs(LFC_ih) <= log2FC_threshold &
#         abs(LFC_ip) >  log2FC_threshold                           ~ "pHis-dominant",
#       TRUE                                                         ~ "Conserved"
#     )) %>%
#     add_te_class() %>%
#     mutate(inheritance = factor(inheritance, levels = inheritance_factor_levels))

df <- df %>% 
      mutate(inheritance = case_when(
        # Additive: hybrid is intermediate (both LFCs within threshold)
          DE_ph == TRUE & DE_ih == FALSE & DE_ip == FALSE  ~ "Additive",
        # Overdominant: hybrid significantly exceeds both parents
        padj_ip <0.05 & padj_ih <0.05 &
          LFC_ip >  log2FC_threshold  & LFC_ih >  log2FC_threshold  ~ "Overdominant",
        # Underdominant: hybrid is significantly below both parents
        padj_ip <0.05 & padj_ih <0.05 &
          LFC_ip < -log2FC_threshold  & LFC_ih < -log2FC_threshold  ~ "Underdominant",
        # pDom-dominant: hybrid is close to pDom (LFC_ip within threshold)
        # but significantly different from pHis in either direction (abs LFC_ih)
        padj_ih <0.05 &
          abs(LFC_ip) <= log2FC_threshold &
          abs(LFC_ih) >  log2FC_threshold                           ~ "pDom-dominant",
        # pHis-dominant: hybrid is close to pHis (LFC_ih within threshold)
        # but significantly different from pDom in either direction (abs LFC_ip)
        padj_ip <0.05 &
          abs(LFC_ih) <= log2FC_threshold &
          abs(LFC_ip) >  log2FC_threshold                           ~ "pHis-dominant",
        TRUE                                                         ~ "Conserved"
      )) %>%
      add_te_class() %>%
      mutate(
        # Override TE_class for manually reclassified IDs.
        # TE_id format: "rnd-5_family-4428/RC/Helitron"
        # Extract base ID before "/" and match exactly (case-insensitive).
        # The %in% exact match prevents partial hits (e.g. 1398 vs 13985).
        TE_id_base = sub("/.*", "", TE_id),
        TE_class   = ifelse(tolower(TE_id_base) %in% tolower(ids_to_unknown),
                            "Unknown", TE_class),
        inheritance = factor(inheritance, levels = inheritance_factor_levels)
      )

  # Save inheritance subsets
  for (cat_name in c("Additive","Overdominant","Underdominant",
                     "pDom-dominant","pHis-dominant","Conserved")) {
    write.csv(
      df[df$inheritance == cat_name, ],
      paste0(short, "_filt_", gsub("-", "_", tolower(cat_name)), "_TEs.csv"),
      row.names = FALSE
    )
  }
  write.csv(df[df$inheritance %in% c("Overdominant","Underdominant"), ],
            paste0(short, "_filt_transgressive_TEs.csv"), row.names = FALSE)
  write.csv(df, paste0(short, "_filt_all_TEs_inheritance.csv"), row.names = FALSE)

  # Inheritance scatterplot
  plot_inheritance_scatter <- ggplot() +
    geom_point(data = df[df$inheritance=="Conserved",], aes(x = LFC_ip, y = LFC_ih, color = inheritance), alpha = 0.7) +
    geom_point(data = df[df$inheritance!="Conserved",], aes(x = LFC_ip, y = LFC_ih, color = inheritance), alpha = 0.7) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dotted") +
    geom_hline(yintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dotted") +
    scale_color_manual(values = inheritance_colors_scatter) +
    theme_minimal() +
    labs(x = "LFC: Italian vs house", y = "LFC: Italian vs Spanish",
         title = paste0("Inheritance patterns — ", label),
         color = "Inheritance Pattern")

  # Inheritance proportion bar chart
  inheritance_summary <- as.data.frame(table(df$inheritance)) %>%
    setNames(c("inheritance_category", "count")) %>%
    mutate(
      proportion           = count / sum(count),
      inheritance_category = factor(inheritance_category,
                                    levels = c("Additive","Conserved","Overdominant",
                                               "Underdominant","pDom-dominant","pHis-dominant"))
    )

  plot_inheritance_bar <- ggplot(inheritance_summary,
                                 aes(x = inheritance_category, y = proportion,
                                     fill = inheritance_category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = inheritance_colors_bar,
                      name   = "Inheritance Category",
                      labels = c("pDom-dominant" = "house-dominant",
                                 "pHis-dominant" = "Spanish-dominant")) +
    scale_x_discrete(labels = c("pDom-dominant" = "house-dominant",
                                "pHis-dominant" = "Spanish-dominant")) +    scale_y_continuous(labels = percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste0("Inheritance categories — ", label),
         x = "Inheritance Category", y = "Proportion of TEs")
    
  # ---------------------------------------------------------------------------
  # 8. INHERITANCE BY TE CLASS (stacked proportional bar chart)
  # ---------------------------------------------------------------------------
  # Factor levels are applied to df BEFORE summarising so they are preserved
  # through group_by/summarise. Rows are then arranged in factor level order
  # (bottom-to-top) since geom_bar stacks in row order, not factor level order.
  df$inheritance <- factor(df$inheritance, levels = inheritance_factor_levels_bar)

  # Actual counts per TE_class x inheritance combination
  te_counts <- df %>%
    group_by(TE_class, inheritance, .drop = FALSE) %>%
    summarise(count = dplyr::n(), .groups = "drop")

  # Full grid ensuring all combinations are present (including zeros)
  te_counts_grid <- expand.grid(
    TE_class    = sort(unique(df$TE_class)),
    inheritance = inheritance_factor_levels_bar,
    stringsAsFactors = FALSE
  ) %>%
    left_join(te_counts, by = c("TE_class", "inheritance")) %>%
    mutate(count = ifelse(is.na(count), 0L, count)) %>%
    group_by(TE_class) %>%
    mutate(total_count = sum(count)) %>%
    ungroup()

  # "All TEs" bar: collapse across TE classes
  all_te_summary <- te_counts_grid %>%
    group_by(inheritance) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(TE_class    = "All TEs",
           total_count = sum(count))

  # Combine, compute proportions, and arrange rows in factor level order
  te_inh_summary <- bind_rows(te_counts_grid, all_te_summary) %>%
    mutate(
      proportion  = ifelse(total_count > 0, count / total_count, 0),
      inheritance = factor(inheritance, levels = inheritance_factor_levels_bar)
    ) %>%
    arrange(TE_class, inheritance)   # row order = bottom-to-top stack order

  class_levels <- intersect(
    c("All TEs","DNA","LINE","LTR","RC","SINE","Unknown"),
    unique(te_inh_summary$TE_class)
  )
  te_inh_summary$TE_class <- factor(te_inh_summary$TE_class, levels = class_levels)

  x_labels <- te_inh_summary %>%
    group_by(TE_class) %>%
    summarise(total_count = dplyr::first(total_count), .groups = "drop") %>%
    { setNames(paste0(.$TE_class, "\n(n = ", .$total_count, ")"), .$TE_class) }

  plot_te_class_bar <- ggplot(te_inh_summary,
                              aes(x = TE_class, y = proportion, fill = inheritance)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = inheritance_colors_bar,
      breaks = inheritance_factor_levels_bar,   # legend order
      drop   = FALSE
    ) +
    scale_y_continuous(labels = percent_format()) +
    scale_x_discrete(labels = x_labels) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    labs(title = paste0("Inheritance by TE class — ", label),
         x = "TE class", y = "Proportion of TEs",
         fill = "Inheritance")

  # ---------------------------------------------------------------------------
  # Print all ggplots for this tissue
  # ---------------------------------------------------------------------------
  print(plot_chrom_bar)
  print(plot_loci_all)
  print(plot_loci_expr)
  print(plot_pca)
  print(plot_volcano_ip)
  print(plot_volcano_ih)
  print(plot_volcano_ph)
  print(plot_inheritance_scatter)
  print(plot_inheritance_bar)
  print(plot_te_class_bar)

  # ---------------------------------------------------------------------------
  # Return all plots and key data objects
  # ---------------------------------------------------------------------------
  list(
    plots = list(
      chrom_bar           = plot_chrom_bar,
      loci_all            = plot_loci_all,
      loci_expr           = plot_loci_expr,
      pca                 = plot_pca,
      volcano_ip          = plot_volcano_ip,
      volcano_ih          = plot_volcano_ih,
      volcano_ph          = plot_volcano_ph,
      inheritance_scatter = plot_inheritance_scatter,
      inheritance_bar     = plot_inheritance_bar,
      te_class_bar        = plot_te_class_bar
    ),
    data = list(
      dds            = dds,
      rld            = rld,
      df_inheritance = df,
      res_ip         = res_ip,
      res_ih         = res_ih,
      res_ph         = res_ph,
      counts_by_cat  = counts_by_cat,
      te_inh_summary = te_inh_summary,
      df_chrom       = df_chrom,
      df_chrom_expr  = df_chrom_expr
    )
  )
}

# =============================================================================
# RUN PIPELINE FOR BOTH TISSUES
# =============================================================================

results <- lapply(tissue_configs, run_te_analysis)

# All plots and data are now accessible as:
#   results$ova$plots$pca
#   results$tes$plots$inheritance_scatter
#   results$ova$data$dds   etc.

# =============================================================================
# PANEL FIGURE ASSEMBLY
# =============================================================================

# Helper: takes a stored PCA plot, removes its title and legend, applies
# the shared species colour scale, and locks the aspect ratio so both
# panels have identical shapes regardless of their data ranges.
prepare_pca <- function(pca_plot) {
  pca_plot +
    scale_color_manual(values = species_colours, name = "Species") +
    theme(legend.position    = "none",
          plot.title         = element_blank(),
          panel.background   = element_blank(),
          panel.border       = element_rect(colour = "black", fill = NA),
          aspect.ratio       = 1)   # square panel regardless of data ranges
}

# =============================================================================
# FIGURE 1: PCA — Testes (A) + Ovaries (B), shared legend to the right
# =============================================================================

# prepare_pca() strips the legend, so we re-add it to one plot and set
# legend.position = "none" on the other. patchwork then collects the
# single surviving legend into guide_area().
pca_tes_clean <- prepare_pca(results$tes$plots$pca) +
  scale_color_manual(values = species_colours, name = "Species") +
  theme(legend.position = "right")   # keep legend on testes plot for collection

pca_ova_clean <- prepare_pca(results$ova$plots$pca)  # legend already stripped

figure1_pca <- (pca_tes_clean | pca_ova_clean | guide_area()) +
  plot_layout(
    widths  = c(1, 1, 0.2),   # two equal plot columns + narrow legend column
    guides  = "collect"        # move the surviving legend into guide_area()
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag          = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.05, 1.0)
  )

print(figure1_pca)

ggsave("figure1_pca.pdf", figure1_pca, width = 10, height = 4)
ggsave("figure1_pca.png", figure1_pca, width = 10, height = 4, dpi = 300)

# =============================================================================
# FIGURE 2: Volcano plots — 3 rows x 2 columns (cowplot)
# Row 1: pDom vs pHis  (A = testes, B = ovaries)
# Row 2: pIta vs pHis  (C = testes, D = ovaries)
# Row 3: pIta vs pDom  (E = testes, F = ovaries)
# Shared legend above the grid, bold A-F tags inside each plot
# =============================================================================
library(cowplot)

# Helper: removes title and legend, keeps white background with border
prepare_volcano <- function(v_plot) {
  v_plot +
    theme(plot.title       = element_blank(),
          legend.position  = "none",
          panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA))
}

# Extract legend from a copy of one plot with legend enabled.
# get_plot_component() is used for compatibility with newer cowplot versions;
# falls back to get_legend() if needed.
volcano_legend_source <- prepare_volcano(results$tes$plots$volcano_ph) +
  theme(legend.position = "top",
        legend.text     = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"))

volcano_legend <- tryCatch(
  cowplot::get_plot_component(volcano_legend_source, "guide-box-top", return_all = TRUE),
  error = function(e) cowplot::get_legend(volcano_legend_source)
)

# Prepare all six plots
v_A <- prepare_volcano(results$tes$plots$volcano_ph)   # A: Testes,  pDom vs pHis
v_B <- prepare_volcano(results$ova$plots$volcano_ph)   # B: Ovaries, pDom vs pHis
v_C <- prepare_volcano(results$tes$plots$volcano_ih)   # C: Testes,  pIta vs pHis
v_D <- prepare_volcano(results$ova$plots$volcano_ih)   # D: Ovaries, pIta vs pHis
v_E <- prepare_volcano(results$tes$plots$volcano_ip)   # E: Testes,  pIta vs pDom
v_F <- prepare_volcano(results$ova$plots$volcano_ip)   # F: Ovaries, pIta vs pDom

# Assemble 3x2 grid
plot_grid_6 <- cowplot::plot_grid(
  v_A, v_B,
  v_C, v_D,
  v_E, v_F,
  ncol           = 2,
  labels         = c("A", "B", "C", "D", "E", "F"),
  label_size     = 14,
  label_fontface = "bold",
  label_x        = 0.05,
  label_y        = 1.05
)

# Place legend above the grid
figure2_volcano <- cowplot::plot_grid(
  volcano_legend,
  plot_grid_6,
  ncol        = 1,
  rel_heights = c(0.06, 1)   # increase 0.06 if legend is clipped
)

print(figure2_volcano)

ggsave("figure2_volcano.pdf", figure2_volcano, width = 8, height = 10)
ggsave("figure2_volcano.png", figure2_volcano, width = 8, height = 10, dpi = 300)

# =============================================================================
# FIGURE 3: Inheritance scatterplots (row 1) + bar charts (row 2)
# A = testes scatter, B = ovaries scatter
# C = testes bar,     D = ovaries bar
# Shared bar chart legend to the right, vertically centred between rows
# =============================================================================

# Helper: remove legend and title from scatter plots
prepare_scatter <- function(p) {
  p + theme(legend.position = "none",
            plot.title      = element_blank())
}

# Helper: remove legend and title from bar charts; expand y-axis upper limit
# so geom_text percentage labels above the tallest bar are not clipped
prepare_bar <- function(p) {
  p + theme(legend.position = "none",
            plot.title      = element_blank()) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL)
}

# Extract legend from bar chart (inheritance_colors_bar palette)
inh_legend <- cowplot::get_plot_component(
  results$tes$plots$inheritance_bar +
    theme(legend.position = "right"),
  "guide-box-right",
  return_all = TRUE
)

# Prepare the four plots
s_A <- prepare_scatter(results$tes$plots$inheritance_scatter)   # A: Testes scatter
s_B <- prepare_scatter(results$ova$plots$inheritance_scatter)   # B: Ovaries scatter
b_C <- prepare_bar(results$tes$plots$inheritance_bar)           # C: Testes bar
b_D <- prepare_bar(results$ova$plots$inheritance_bar)           # D: Ovaries bar

# Build the 2x2 plot grid with a blank spacer row between the two rows.
# rel_heights = c(1, gap, 1): increase the middle value to add more spacing.
row_gap <- 0.08

grid_4 <- cowplot::plot_grid(
  s_A, s_B,
  NULL, NULL,        # spacer row
  b_C, b_D,
  ncol        = 2,
  rel_heights = c(1, row_gap, 1)
)

# Combine grid and legend side by side, legend centred vertically
figure3_inh <- cowplot::plot_grid(
  grid_4,
  inh_legend,
  ncol        = 2,
  rel_widths  = c(1, 0.25)   # widened to prevent legend being clipped
)

# Overlay panel labels. With the spacer, each plot row occupies
# 1/(2 + row_gap) of the total height. Labels are nudged above the top
# edge of each panel by using y values just above those boundaries.
row_h    <- 1 / (2 + row_gap)   # fractional height of one plot row
# The figure is shrunk slightly within the canvas (ymax = 0.97) to leave a
# top margin where A/B labels can sit without being clipped.
top_margin <- 0.03   # fraction of figure height reserved at top for A/B labels
top_y      <- 1.0 - (top_margin * 0.3)      # A/B: sits in the top margin
bottom_y   <- 1.0 - top_margin - row_h - (row_gap * 0.1)  # C/D: top of spacer row

figure3_inh <- ggdraw() +
  draw_plot(ggplot() +                              # white background
              theme(panel.background = element_rect(fill = "white", colour = NA),
                    plot.background  = element_rect(fill = "white", colour = NA))) +
  draw_plot(figure3_inh,                            # shrink plot down to leave top margin
            x = 0, y = 0, width = 1, height = 1 - top_margin) +
  draw_label("A", x = 0.02,  y = top_y,    fontface = "bold", size = 14) +
  draw_label("B", x = 0.42, y = top_y,    fontface = "bold", size = 14) +
  draw_label("C", x = 0.02,  y = bottom_y, fontface = "bold", size = 14) +
  draw_label("D", x = 0.42, y = bottom_y, fontface = "bold", size = 14)

print(figure3_inh)

ggsave("figure3_inheritance.pdf", figure3_inh, width = 10, height = 8)
ggsave("figure3_inheritance.png", figure3_inh, width = 10, height = 8, dpi = 300)

# =============================================================================
# COMBINED TE CLASS BAR CHART — Testes and Ovaries side by side
# For each TE class the testes bar appears on the left and ovaries on the right,
# touching with a thin white dividing line. Each bar has its own n= sub-label.
# Individual tissue plots (plot_te_class_bar in section 8) are kept separate.
# =============================================================================

# Pull the per-tissue summary tables and drop any NA TE classes
te_combined <- bind_rows(
  results$tes$data$te_inh_summary %>% mutate(tissue = "Testes"),
  results$ova$data$te_inh_summary %>% mutate(tissue = "Ovaries")
) %>%
  filter(!is.na(TE_class)) %>%
  mutate(tissue = factor(tissue, levels = c("Testes", "Ovaries")))

# Canonical class order (same as individual plots)
class_levels_combined <- intersect(
  c("All TEs","DNA","LINE","LTR","RC","SINE","Unknown"),
  unique(as.character(te_combined$TE_class))
)
te_combined$TE_class <- factor(te_combined$TE_class, levels = class_levels_combined)

# Build the interaction factor with explicitly ordered levels so bars are
# interleaved: All TEs.Testes, All TEs.Ovaries, DNA.Testes, DNA.Ovaries, ...
correct_bar_levels <- as.vector(rbind(
  paste0(class_levels_combined, ".Testes"),
  paste0(class_levels_combined, ".Ovaries")
))

te_combined <- te_combined %>%
  mutate(bar_id = factor(
    paste0(TE_class, ".", tissue),
    levels = correct_bar_levels
  ))

# n= totals per tissue per class for x-axis sub-labels
totals_combined <- te_combined %>%
  group_by(TE_class, tissue) %>%
  summarise(total_count = dplyr::first(total_count), .groups = "drop")

# Build named label vector: one entry per bar_id level (format: "Class.Tissue").
# Testes bars get the class name + both n= values; Ovaries bars get a blank.
bar_labels <- levels(te_combined$bar_id) %>%
  { setNames(., .) } %>%
  sapply(function(bid) {
    parts <- strsplit(bid, ".", fixed = TRUE)[[1]]
    cls   <- parts[1]
    tis   <- parts[length(parts)]   # last element is tissue
    # class name may contain spaces but not dots, so cls is always parts[1]
    if (tis == "Testes") {
      n_tes <- totals_combined$total_count[
        totals_combined$TE_class == cls & totals_combined$tissue == "Testes"]
      n_ova <- totals_combined$total_count[
        totals_combined$TE_class == cls & totals_combined$tissue == "Ovaries"]
      paste0(cls, "
Testes n=", n_tes, "
Ovaries n=", n_ova)
    } else {
      ""   # blank label for Ovaries bar; class label shown on Testes bar
    }
  })

# Number of TE classes — used to place the white dividing lines
n_classes <- length(class_levels_combined)

# Use a numeric x axis so we can manually control bar positions:
# within each TE class the two bars are placed 0.45 apart (touching);
# between TE classes the gap is 1.5 units.
inner_gap  <- 0.45   # distance between testes and ovaries bar centres
outer_gap  <- 1.5    # distance between the centres of adjacent TE class pairs

# Compute bar centre x positions
pair_centres <- seq(0, by = outer_gap, length.out = n_classes)
tes_x <- pair_centres - inner_gap / 2
ova_x <- pair_centres + inner_gap / 2

# Add numeric x positions to the data
te_combined <- te_combined %>%
  mutate(x_pos = ifelse(tissue == "Testes",
                        tes_x[as.integer(TE_class)],
                        ova_x[as.integer(TE_class)]))

# Dividing line x positions: midpoint between testes and ovaries bars
divider_x <- pair_centres   # sits exactly between the two bars in each pair

# Testes/Ovaries tick label positions and text
tick_x      <- c(rbind(tes_x, ova_x))
tick_labels <- rep(c("Testes", "Ovaries"), times = n_classes)
tick_n      <- c(rbind(
  totals_combined$total_count[totals_combined$tissue == "Testes"],
  totals_combined$total_count[totals_combined$tissue == "Ovaries"]
))
tick_text <- paste0(tick_labels, "
n=", tick_n)

# TE class label positions: midpoint of each pair, placed below tick labels
class_label_x <- pair_centres

plot_te_class_combined <- ggplot(
    te_combined,
    aes(x = x_pos, y = proportion, fill = inheritance, group = interaction(x_pos, inheritance))
  ) +
  geom_col(position = "stack", width = inner_gap - 0.02, colour = NA) +
  # White dividing line between testes and ovaries within each TE class
  geom_vline(
    data        = data.frame(xintercept = divider_x),
    aes(xintercept = xintercept),
    colour      = "white",
    linewidth   = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = inheritance_colors_bar,
    breaks = inheritance_factor_levels_bar,
    name   = "Inheritance Category",
    drop   = FALSE
  ) +
  scale_y_continuous(labels = percent_format()) +
  # Numeric x axis: testes/ovaries tick labels at bar positions
  scale_x_continuous(
    breaks = tick_x,
    labels = tick_text,
    expand = expansion(add = 0.6)
  ) +
  # TE class labels spanning each pair, drawn below the tick labels
  annotate("text",
           x        = class_label_x,
           y        = -0.085,
           label    = class_levels_combined,
           size     = 3,
           fontface = "bold",
           hjust    = 0.5) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top",
    plot.margin     = margin(t = 5, r = 10, b = 40, l = 10)  # extra bottom for class labels
  ) +
  labs(
    title = "Inheritance by TE class — Testes and Ovaries",
    x     = "TE class",
    y     = "Proportion of TEs"
  )

print(plot_te_class_combined)

ggsave("figure4_te_class_bar.pdf", plot_te_class_combined, width = 12, height = 6)
ggsave("figure4_te_class_bar.png", plot_te_class_combined, width = 12, height = 6, dpi = 300)

# =============================================================================
# FIGURE 3+4 COMBINED: Inheritance scatterplots + bar charts (rows 1-2)
#                      + Combined TE class bar chart (row 3)
# A-D as in figure 3; E = combined TE class bar chart
# Single shared legend below the panel
# =============================================================================

# Strip legend and title from combined TE class bar chart for panel use
te_class_panel <- plot_te_class_combined +
  theme(legend.position = "none",
        plot.title      = element_blank())

# Extract legend from bar chart for shared bottom legend
bottom_legend <- cowplot::get_plot_component(
  results$tes$plots$inheritance_bar +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 9),
          legend.key.size  = unit(0.4, "cm")) +
    guides(fill = guide_legend(nrow = 1)),
  "guide-box-bottom",
  return_all = TRUE
)

# Stack the three rows: 2x2 grid, then combined TE class bar
# rel_heights: top two rows (equal) + bottom bar chart row
combined_grid <- cowplot::plot_grid(
  grid_4,
  te_class_panel,
  ncol        = 1,
  rel_heights = c(1, 0.6, 0.7)
)

# Add legend below
combined_with_legend <- cowplot::plot_grid(
  combined_grid,
  bottom_legend,
  ncol        = 1,
  rel_heights = c(1, 0.03)
)

# Add top margin for A/B labels and overlay all panel labels
top_margin_comb <- 0.03
row_h_comb      <- 1 / (2 + row_gap)   # height of one row within the top 2x2

# y positions are in the combined figure's coordinate space.
# The top 2x2 grid occupies rel_height 1/(1+0.7) = ~0.588 of combined_grid,
# and the legend strip takes ~0.06/(1.06) of the total figure height.
plot_frac   <- 1 / 1.06          # fraction of total height that is plots (vs legend)
top_grid_h  <- plot_frac * (1 / (1 + 0.7))   # fraction of total = top 2x2 portion
row_h_abs   <- top_grid_h * row_h_comb        # absolute height of one scatter/bar row

top_y_comb    <- plot_frac - (top_margin_comb * 0.3)
bottom_y_comb <- plot_frac - top_margin_comb - row_h_abs - (row_gap * top_grid_h * 0.1)
bar_y_comb    <- top_grid_h - 0.02            # E label: just above TE class bar

figure3_4_combined <- ggdraw() +
  draw_plot(ggplot() +
              theme(panel.background = element_rect(fill = "white", colour = NA),
                    plot.background  = element_rect(fill = "white", colour = NA))) +
  draw_plot(combined_with_legend,
            x = 0, y = 0, width = 1, height = 1 - top_margin_comb) +
  draw_label("A", x = 0.02,  y = 0.98,    fontface = "bold", size = 14) +
  draw_label("B", x = 0.52,  y = 0.98,    fontface = "bold", size = 14) +
  draw_label("C", x = 0.02,  y = 0.67, fontface = "bold", size = 14) +
  draw_label("D", x = 0.52,  y = 0.67, fontface = "bold", size = 14) +
  draw_label("E", x = 0.02,  y = 0.38,    fontface = "bold", size = 14)

print(figure3_4_combined)

ggsave("figure3_4_combined.pdf", figure3_4_combined, width = 10, height = 12)
ggsave("figure3_4_combined.png", figure3_4_combined, width = 10, height = 12, dpi = 300)



# =============================================================================
# FIGURE 5: Chromosome category panel — 2 rows x 3 columns
#
# Row 1 (all TEs — categories based on where copies are found):
#   A = TEs by all 7 chromosome categories (A/W/Z combinations) — shared
#   B = TEs by all 7 categories, expressed in testes
#   C = TEs by all 7 categories, expressed in ovaries
#
# Row 2 (loci-specific — single copy so only A, W, or Z):
#   D = All TE loci by chromosome category (A/W/Z only) — shared
#   E = Expressed TE loci by chromosome category — testes
#   F = Expressed TE loci by chromosome category — ovaries
# =============================================================================

# Helper: remove title and legend
prepare_chrom <- function(p) {
  p + theme(plot.title      = element_blank(),
            legend.position = "none")
}

# --- Row 1: collapsed TE counts by all 7 categories --------------------------

# Helper: filter counts_by_cat to only TEs present in a tissue's dds object
# (i.e. passed the low-count filter) then recount by chromosome category.
# annot_chr is rebuilt from the stored chrom map data via counts_by_cat.
make_chrom_bar_expressed <- function(tissue_key, tissue_label) {
  # TE IDs that passed filtering in this tissue's dds
  expressed_ids <- rownames(results[[tissue_key]]$data$dds)

  # Reuse the full chrom map annotation stored via counts_by_cat by
  # re-reading it from the chrom_map file — but filter to expressed IDs only
  te_map_tissue <- read_tsv(tissue_configs[[tissue_key]]$chrom_map,
                             show_col_types = FALSE) %>%
    mutate(across(c(TE_id, category), trimws),
           category = vapply(category, normalize_category, character(1)))

  counts_expr <- te_map_tissue %>%
    filter(TE_id %in% expressed_ids, !is.na(category)) %>%
    mutate(category = factor(category, levels = cat_levels)) %>%
    dplyr::count(category) %>%
    complete(category = factor(cat_levels, levels = cat_levels), fill = list(n = 0)) %>%
    arrange(category)

  ggplot(counts_expr, aes(x = category, y = n)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Chromosome category",
         y = paste0("Number of TEs (expressed in ", tissue_label, ")"),
         title = paste0("Expressed TEs by chromosome category (", tissue_label, ")"))
}

# Row 1 plots
ch_A <- prepare_chrom(results$tes$plots$chrom_bar)                        # A: all TEs, shared
ch_B <- prepare_chrom(make_chrom_bar_expressed("tes", "testes"))           # B: expressed, testes
ch_C <- prepare_chrom(make_chrom_bar_expressed("ova", "ovaries"))          # C: expressed, ovaries

# --- Row 2: loci-specific A/W/Z only -----------------------------------------
ch_D <- prepare_chrom(results$tes$plots$loci_all)                          # D: all loci, shared

ch_E <- prepare_chrom(results$tes$plots$loci_expr) +
  labs(y = "Number of TE loci (expressed)")                                # E: expressed loci, testes

ch_F <- prepare_chrom(results$ova$plots$loci_expr) +
  labs(y = "Number of TE loci (expressed)") +
  scale_y_continuous(labels = scales::scientific)                           # F: expressed loci, ovaries

# Assemble 2x3 grid
chrom_grid <- cowplot::plot_grid(
  ch_A, ch_B, ch_C,
  ch_D, ch_E, ch_F,
  ncol           = 3,
  labels         = c("A", "B", "C", "D", "E", "F"),
  label_size     = 14,
  label_fontface = "bold",
  label_x        = 0.05,
  label_y        = 1.01
)

# Wrap in ggdraw with a top margin so top row labels have room
top_margin_fig5 <- 0.04

figure5_chromosome <- ggdraw() +
  draw_plot(ggplot() +
              theme(panel.background = element_rect(fill = "white", colour = NA),
                    plot.background  = element_rect(fill = "white", colour = NA))) +
  draw_plot(chrom_grid, x = 0, y = 0, width = 1, height = 1 - top_margin_fig5)

print(figure5_chromosome)

ggsave("figure5_chromosome.pdf", figure5_chromosome, width = 14, height = 8)
ggsave("figure5_chromosome.png", figure5_chromosome, width = 14, height = 8, dpi = 300)

# =============================================================================
# SAVE DDS OBJECTS
# TE dds objects are extracted from results; gene dds objects are built fresh
# from the gene count files using the same sample metadata as the TE pipeline.
# All objects saved as .rds to the R_data_objects directory.
# =============================================================================

rdata_dir <- "path/to/R_data_objects"  # <-- EDIT: directory to save .rds files

# --- Save filtered TE dds objects --------------------------------------------
saveRDS(results$tes$data$dds,
        file.path(rdata_dir, "dds_TE_testes.rds"))
saveRDS(results$ova$data$dds,
        file.path(rdata_dir, "dds_TE_ovaries.rds"))

cat("Saved TE dds objects\n")

# --- Build and save gene dds objects -----------------------------------------
gene_counts_dir <- "path/to/gene_count_files"  # <-- EDIT: directory containing gene count files

# Testes gene dds
testes_gene_counts <- read.delim(
  file.path(gene_counts_dir, "testes_TElocal_counts_genes.tsv"),
  row.names = 1, stringsAsFactors = FALSE
)
testes_gene_sample_table <- data.frame(
  row.names = tissue_configs$tes$sample_names,
  species   = factor(tissue_configs$tes$sample_species)
)
# Ensure column order matches sample table
testes_gene_counts <- testes_gene_counts[, tissue_configs$tes$sample_names]

dds_gene_testes <- DESeqDataSetFromMatrix(
  countData = testes_gene_counts,
  colData   = testes_gene_sample_table,
  design    = ~ species
)
dds_gene_testes <- dds_gene_testes[rowSums(counts(dds_gene_testes)) > 1, ]
dds_gene_testes <- DESeq(dds_gene_testes)

saveRDS(dds_gene_testes, file.path(rdata_dir, "dds_gene_testes.rds"))
cat("Saved gene dds: testes\n")

# Ovaries gene dds
ovaries_gene_counts <- read.delim(
  file.path(gene_counts_dir, "ovaries_TElocal_counts_genes.tsv"),
  row.names = 1, stringsAsFactors = FALSE
)
ovaries_gene_sample_table <- data.frame(
  row.names = tissue_configs$ova$sample_names,
  species   = factor(tissue_configs$ova$sample_species)
)
ovaries_gene_counts <- ovaries_gene_counts[, tissue_configs$ova$sample_names]

dds_gene_ovaries <- DESeqDataSetFromMatrix(
  countData = ovaries_gene_counts,
  colData   = ovaries_gene_sample_table,
  design    = ~ species
)
dds_gene_ovaries <- dds_gene_ovaries[rowSums(counts(dds_gene_ovaries)) > 1, ]
dds_gene_ovaries <- DESeq(dds_gene_ovaries)

saveRDS(dds_gene_ovaries, file.path(rdata_dir, "dds_gene_ovaries.rds"))
cat("Saved gene dds: ovaries\n")

# =============================================================================
# PAIRWISE DE SUMMARY — numbers and percentages for results text
# Thresholds: padj < 0.05 AND abs(log2FoldChange) > 0.32 (matching volcano plots)
# Comparisons:
#   res_ph: pDom vs pHis  (positive LFC = higher in pDom)
#   res_ip: pIta vs pDom  (positive LFC = higher in pIta)
#   res_ih: pIta vs pHis  (positive LFC = higher in pIta)
# =============================================================================

summarise_de <- function(res, comparison_label, tissue_label) {
  df <- as.data.frame(res) %>%
    filter(!is.na(padj), !is.na(log2FoldChange))
  
  total    <- nrow(df)
  sig      <- df %>% filter(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold)
  n_up     <- sum(sig$log2FoldChange > 0)
  n_down   <- sum(sig$log2FoldChange < 0)
  n_sig    <- nrow(sig)
  pct_up   <- round(n_up   / total * 100, 1)
  pct_down <- round(n_down / total * 100, 1)
  pct_sig  <- round(n_sig  / total * 100, 1)
  pct_cons <- round((total - n_sig) / total * 100, 1)
  
  cat("\n---", tissue_label, "|", comparison_label, "---\n")
  cat("  Total TEs tested:        ", total, "\n")
  cat("  Significant (total):     ", n_sig,  " (", pct_sig,  "%)\n", sep = "")
  cat("  Upregulated (pos LFC):   ", n_up,   " (", pct_up,   "%)\n", sep = "")
  cat("  Downregulated (neg LFC): ", n_down, " (", pct_down, "%)\n", sep = "")
  cat("  Conserved:               ", total - n_sig, " (", pct_cons, "%)\n", sep = "")
}

# Testes
summarise_de(results$tes$data$res_ph, "pDom vs pHis",  "Testes")
summarise_de(results$tes$data$res_ip, "pIta vs pDom",  "Testes")
summarise_de(results$tes$data$res_ih, "pIta vs pHis",  "Testes")

# Ovaries
summarise_de(results$ova$data$res_ph, "pDom vs pHis",  "Ovaries")
summarise_de(results$ova$data$res_ip, "pIta vs pDom",  "Ovaries")
summarise_de(results$ova$data$res_ih, "pIta vs pHis",  "Ovaries")

# Capture all output and write to file
sink("pairwise_DE_summary.txt")
summarise_de(results$tes$data$res_ph, "pDom vs pHis",  "Testes")
summarise_de(results$tes$data$res_ip, "pIta vs pDom",  "Testes")
summarise_de(results$tes$data$res_ih, "pIta vs pHis",  "Testes")
summarise_de(results$ova$data$res_ph, "pDom vs pHis",  "Ovaries")
summarise_de(results$ova$data$res_ip, "pIta vs pDom",  "Ovaries")
summarise_de(results$ova$data$res_ih, "pIta vs pHis",  "Ovaries")
sink()

# =============================================================================
# CHROMOSOME CATEGORY SUMMARIES — saved to 4 separate .txt files
#
# File 1: TE counts by all 7 chromosome categories (A/W/Z combinations)
#         Source: counts_by_cat (collapsed TEs from chrom map, figure A)
# File 2: TE loci counts for A, W, Z only
#         Source: df_chrom (all loci, figure B)
# File 3: Expressed TE loci counts for A, W, Z — testes
#         Source: df_chrom_expr (testes)
# File 4: Expressed TE loci counts for A, W, Z — ovaries
#         Source: df_chrom_expr (ovaries)
# =============================================================================

# Helper: summarise a data frame's category column into n and percentage
summarise_chrom <- function(df, category_col, title) {
  cats  <- df[[category_col]]
  cats  <- cats[!is.na(cats) & cats != ""]
  total <- length(cats)
  tbl   <- table(cats)
  out   <- data.frame(
    Category   = names(tbl),
    N          = as.integer(tbl),
    Percentage = round(as.integer(tbl) / total * 100, 2),
    stringsAsFactors = FALSE
  )
  list(title = title, total = total, data = out)
}

write_summary <- function(summary_list, filepath) {
  sink(filepath)
  cat(summary_list$title, "\n")
  cat(strrep("-", nchar(summary_list$title)), "\n")
  cat("Total:", summary_list$total, "\n\n")
  print(summary_list$data, row.names = FALSE)
  sink()
  cat("Written:", filepath, "\n")
}

# --- File 1: all 7 categories from counts_by_cat (figure A) -----------------
# counts_by_cat already has n per category; rebuild as raw vector for helper
counts_by_cat_tes <- results$tes$data$counts_by_cat
chrom_all7 <- rep(counts_by_cat_tes$category, counts_by_cat_tes$n)
chrom_all7_df <- data.frame(category = as.character(chrom_all7))
s1 <- summarise_chrom(chrom_all7_df, "category",
                      "TEs by chromosome category — all 7 groups (A/W/Z combinations)")
# Ensure all 7 levels present and in order
s1$data <- s1$data[match(cat_levels, s1$data$Category), ]
s1$data <- s1$data[!is.na(s1$data$Category), ]
write_summary(s1, "chrom_summary_all7_categories.txt")

# --- File 2: A, W, Z only from df_chrom loci (figure B) --------------------
df_loci_awz <- results$tes$data$df_chrom %>%
  filter(category %in% c("A", "W", "Z"))
s2 <- summarise_chrom(df_loci_awz, "category",
                      "TE loci by chromosome category — A, W, Z only (all loci, figure B)")
write_summary(s2, "chrom_summary_loci_AWZ.txt")

# --- File 3: expressed TE loci A/W/Z — testes -------------------------------
df_expr_tes_awz <- results$tes$data$df_chrom_expr %>%
  filter(category %in% c("A", "W", "Z"))
s3 <- summarise_chrom(df_expr_tes_awz, "category",
                      "Expressed TE loci by chromosome category — A, W, Z (Testes)")
write_summary(s3, "chrom_summary_expressed_AWZ_testes.txt")

# --- File 4: expressed TE loci A/W/Z — ovaries ------------------------------
df_expr_ova_awz <- results$ova$data$df_chrom_expr %>%
  filter(category %in% c("A", "W", "Z"))
s4 <- summarise_chrom(df_expr_ova_awz, "category",
                      "Expressed TE loci by chromosome category — A, W, Z (Ovaries)")
write_summary(s4, "chrom_summary_expressed_AWZ_ovaries.txt")
