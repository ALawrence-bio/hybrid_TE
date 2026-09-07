suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggtext))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))

###################################################################################
########################## Repeat landscape plot generation  ######################
################################ for all passer species ###########################
###################################################################################

# Set variables
out_directory="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/repeat_landscapes"
# for pIta
species_name_pIta="pIta"
in_gff_pIta="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/pIta/earlGrey_pIta_EarlGrey/earlGrey_pIta_summaryFiles/earlGrey_pIta.filteredRepeats.gff"
out_directory_pIta="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/repeat_landscapes" #/repeat_landscape_recreation_test

# for pDom
species_name_pDom="pDom"
in_gff_pDom="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/pDom/earlGrey_pDom_EarlGrey/earlGrey_pDom_summaryFiles/earlGrey_pDom.filteredRepeats.gff"
out_directory_pDom="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/repeat_landscapes" #/repeat_landscape_recreation_test

# for pHis
species_name_pHis="pHis"
in_gff_pHis="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/pHis/earlGrey_pHis_EarlGrey/earlGrey_pHis_summaryFiles/earlGrey_pHis.filteredRepeats.gff"
out_directory_pHis="C:/Users/alexi/Documents/Lund_University_HybridExpress_project/Transposable_Elements/earlGrey/repeat_landscapes" #/repeat_landscape_recreation_test

# for all passer species (for TE subclass across species plots)
subclass_name_DNA="DNA transposon"


# Create plot titles for passer species
plot_title_pIta <- paste0("Repeat landscape of *", gsub("_", " ", species_name_pIta), "*")
title_plot_pIta <- ggplot() + labs(title = plot_title_pIta) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())
plot_title_pDom <- paste0("Repeat landscape of *", gsub("_", " ", species_name_pDom), "*")
title_plot_pDom <- ggplot() + labs(title = plot_title_pDom) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())
plot_title_pHis <- paste0("Repeat landscape of *", gsub("_", " ", species_name_pHis), "*")
title_plot_pHis <- ggplot() + labs(title = plot_title_pHis) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())

# Create plot titles for TE subclasses
plot_title_DNA <- paste0("Repeat landscape of *", gsub("_", " ", subclass_name_DNA), "*")
title_plot_DNA <- ggplot() + labs(title = plot_title_DNA) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())

# Read in data
divergence_eg_gff_pIta <- read_gff(in_gff_pIta)
divergence_eg_gff_pDom <- read_gff(in_gff_pDom)
divergence_eg_gff_pHis <- read_gff(in_gff_pHis)

# Breakdown classification of repeats
divergence_eg_tes_gff_pIta <- divergence_eg_gff_pIta %>%
  dplyr::mutate(subclass = sub("/.*", "", type),
                superfamily = sub("-.*", "", sub(".*/", "", type)))
divergence_eg_tes_gff_pDom <- divergence_eg_gff_pDom %>%
  dplyr::mutate(subclass = sub("/.*", "", type),
                superfamily = sub("-.*", "", sub(".*/", "", type)))
divergence_eg_tes_gff_pHis <- divergence_eg_gff_pHis %>%
  dplyr::mutate(subclass = sub("/.*", "", type),
                superfamily = sub("-.*", "", sub(".*/", "", type)))

# Fix Penelopes
divergence_eg_tes_gff_pIta <- divergence_eg_tes_gff_pIta %>%
  dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass)) %>%
  dplyr::mutate(subclass = ifelse(subclass %in% c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Unknown"), subclass, "Other")) %>%
  dplyr::mutate(named_subclass = case_when(subclass == "DNA" ~ "DNA Transposon",
                                           subclass == "LTR" ~ "LTR Retrotransposon",
                                           subclass == "PLE" ~ "Penelope",
                                           subclass == "RC" ~ "Rolling Circle",
                                           .default = subclass))
divergence_eg_tes_gff_pDom <- divergence_eg_tes_gff_pDom %>%
  dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass)) %>%
  dplyr::mutate(subclass = ifelse(subclass %in% c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Unknown"), subclass, "Other")) %>%
  dplyr::mutate(named_subclass = case_when(subclass == "DNA" ~ "DNA Transposon",
                                           subclass == "LTR" ~ "LTR Retrotransposon",
                                           subclass == "PLE" ~ "Penelope",
                                           subclass == "RC" ~ "Rolling Circle",
                                           .default = subclass))
divergence_eg_tes_gff_pHis <- divergence_eg_tes_gff_pHis %>%
  dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass)) %>%
  dplyr::mutate(subclass = ifelse(subclass %in% c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Unknown"), subclass, "Other")) %>%
  dplyr::mutate(named_subclass = case_when(subclass == "DNA" ~ "DNA Transposon",
                                           subclass == "LTR" ~ "LTR Retrotransposon",
                                           subclass == "PLE" ~ "Penelope",
                                           subclass == "RC" ~ "Rolling Circle",
                                           
                                           .default = subclass))
# Reclassify specific TEs into the Unknown category
# These families are being moved to Unknown regardless of their original annotation.
# In the GFF, these are identified via the ID attribute in column 9, which read_gff
# parses into a column called 'ID'. The IDs are in uppercase e.g. RND-5_FAMILY-6382.
ids_to_unknown <- c(
  "RND-5_FAMILY-6382",   # was rnd-5_family-6382#DNA/Dada
  "RND-5_FAMILY-4100",   # was rnd-5_family-4100#PLE/Chlamys
  "RND-5_FAMILY-2053",   # was rnd-5_family-2053#SINE/MIR
  "RND-5_FAMILY-32486",  # was rnd-5_family-32486#SINE/ID
  "RND-5_FAMILY-72847",  # was rnd-5_family-72847#SINE/ID
  "RND-5_FAMILY-4428",   # was rnd-5_family-4428#RC/Helitron
  "RND-5_FAMILY-5246",   # was rnd-5_family-5246#RC/Helitron
  "RND-5_FAMILY-8403",   # was rnd-5_family-8403#RC/Helitron
  "RND-5_FAMILY-1398",   # was rnd-5_family-1398#RC/Helitron
  "RND-5_FAMILY-2412"    # was rnd-5_family-2412#RC/Helitron
)

reclassify_to_unknown <- function(df) {
  df %>%
    dplyr::mutate(
      subclass       = dplyr::if_else(ID %in% ids_to_unknown, "Unknown", subclass),
      named_subclass = dplyr::if_else(ID %in% ids_to_unknown, "Unknown", named_subclass),
      superfamily    = dplyr::if_else(ID %in% ids_to_unknown, "Unknown", superfamily)
    )
}

divergence_eg_tes_gff_pIta <- reclassify_to_unknown(divergence_eg_tes_gff_pIta)
divergence_eg_tes_gff_pDom <- reclassify_to_unknown(divergence_eg_tes_gff_pDom)
divergence_eg_tes_gff_pHis <- reclassify_to_unknown(divergence_eg_tes_gff_pHis)

#Reclassify manually curated TEs from "Other" into correct categories based on annotation
# valid subclasses to allow
valid_subclasses <- c("DNA","LINE","LTR","SINE","PLE","RC")

# Reclassify curated entries: if subclass == "Other" and ID ends with _PASHIS/_PASHISC/_PASDOM,
# set subclass from the prefix of superfamily (before _, -, or /).
# Also fix the special "unspecified" superfamily for GGERVL etc. using the family/ID.
# for pIta
divergence_eg_tes_gff_pIta <- divergence_eg_tes_gff_pIta %>%
  dplyr::mutate(
    # identify curated families by suffix
    is_curated = stringr::str_detect(ID, "(_PASHIS|_PASHISC|_PASDOM)$"),
    
    # clean superfamily text
    superfamily_trim = stringr::str_trim(superfamily),
    
    # if superfamily is "unspecified" (any case/whitespace) AND it's a curated "Other",
    # infer a better superfamily from the family/ID (extend these rules as needed)
    superfamily_fixed = dplyr::if_else(
      subclass == "Other" & is_curated &
        stringr::str_detect(stringr::str_to_lower(superfamily_trim), "^\\s*unspecified\\s*$"),
      dplyr::case_when(
        stringr::str_detect(ID, "ERVL") ~ "LTR_ERVL",
        stringr::str_detect(ID, "ERVK") ~ "LTR_ERVK",
        stringr::str_detect(ID, "^CR1") ~ "LINE_CR1",
        TRUE ~ superfamily_trim
      ),
      superfamily_trim
    ),
    
    # derive subclass from the prefix of the (possibly fixed) superfamily
    subclass_from_sf = toupper(sub("[-_/].*$", "", superfamily_fixed)),
    subclass_from_sf = dplyr::if_else(subclass_from_sf == "PENELOPE", "PLE", subclass_from_sf),
    
    # apply the reassignment only to curated "Other" rows, and only to valid subclasses
    subclass = dplyr::if_else(
      subclass == "Other" & is_curated & subclass_from_sf %in% valid_subclasses,
      subclass_from_sf, subclass
    ),
    
    # write back the corrected superfamily: keep only the suffix after the first separator
    superfamily = dplyr::case_when(
      grepl("_", superfamily_fixed) ~ sub("^[^_]+_", "", superfamily_fixed),  # LINE_CR1 -> CR1
      grepl("/", superfamily_fixed) ~ sub("^.*/", "", superfamily_fixed),     # DNA/PIF-Harbinger -> PIF-Harbinger
      grepl("-", superfamily_fixed) ~ sub("^[^-]+-", "", superfamily_fixed),  # LTR-ERVK -> ERVK
      TRUE ~ superfamily_fixed
    ),
    
    # keep your human-friendly label in sync
    named_subclass = dplyr::case_when(
      subclass == "DNA"  ~ "DNA Transposon",
      subclass == "LTR"  ~ "LTR Retrotransposon",
      subclass == "PLE"  ~ "Penelope",
      subclass == "RC"   ~ "Rolling Circle",
      TRUE               ~ subclass
    )
  ) %>%
  dplyr::select(-is_curated, -superfamily_trim, -superfamily_fixed, -subclass_from_sf)
# for pDom
# Reclassify curated entries: if subclass == "Other" and ID ends with _PASHIS/_PASHISC/_PASDOM,
# set subclass from the prefix of superfamily (before _, -, or /).
# Also fix the special "unspecified" superfamily for GGERVL etc. using the family/ID.
divergence_eg_tes_gff_pDom <- divergence_eg_tes_gff_pDom %>%
  dplyr::mutate(
    # identify curated families by suffix
    is_curated = stringr::str_detect(ID, "(_PASHIS|_PASHISC|_PASDOM)$"),
    
    # clean superfamily text
    superfamily_trim = stringr::str_trim(superfamily),
    
    # if superfamily is "unspecified" (any case/whitespace) AND it's a curated "Other",
    # infer a better superfamily from the family/ID (extend these rules as needed)
    superfamily_fixed = dplyr::if_else(
      subclass == "Other" & is_curated &
        stringr::str_detect(stringr::str_to_lower(superfamily_trim), "^\\s*unspecified\\s*$"),
      dplyr::case_when(
        stringr::str_detect(ID, "ERVL") ~ "LTR_ERVL",
        stringr::str_detect(ID, "ERVK") ~ "LTR_ERVK",
        stringr::str_detect(ID, "^CR1") ~ "LINE_CR1",
        TRUE ~ superfamily_trim
      ),
      superfamily_trim
    ),
    
    # derive subclass from the prefix of the (possibly fixed) superfamily
    subclass_from_sf = toupper(sub("[-_/].*$", "", superfamily_fixed)),
    subclass_from_sf = dplyr::if_else(subclass_from_sf == "PENELOPE", "PLE", subclass_from_sf),
    
    # apply the reassignment only to curated "Other" rows, and only to valid subclasses
    subclass = dplyr::if_else(
      subclass == "Other" & is_curated & subclass_from_sf %in% valid_subclasses,
      subclass_from_sf, subclass
    ),
    
    # write back the corrected superfamily: keep only the suffix after the first separator
    superfamily = dplyr::case_when(
      grepl("_", superfamily_fixed) ~ sub("^[^_]+_", "", superfamily_fixed),  # LINE_CR1 -> CR1
      grepl("/", superfamily_fixed) ~ sub("^.*/", "", superfamily_fixed),     # DNA/PIF-Harbinger -> PIF-Harbinger
      grepl("-", superfamily_fixed) ~ sub("^[^-]+-", "", superfamily_fixed),  # LTR-ERVK -> ERVK
      TRUE ~ superfamily_fixed
    ),
    
    # keep your human-friendly label in sync
    named_subclass = dplyr::case_when(
      subclass == "DNA"  ~ "DNA Transposon",
      subclass == "LTR"  ~ "LTR Retrotransposon",
      subclass == "PLE"  ~ "Penelope",
      subclass == "RC"   ~ "Rolling Circle",
      TRUE               ~ subclass
    )
  ) %>%
  dplyr::select(-is_curated, -superfamily_trim, -superfamily_fixed, -subclass_from_sf)
# for pHis
# Reclassify curated entries: if subclass == "Other" and ID ends with _PASHIS/_PASHISC/_PASDOM,
# set subclass from the prefix of superfamily (before _, -, or /).
# Also fix the special "unspecified" superfamily for GGERVL etc. using the family/ID.
divergence_eg_tes_gff_pHis <- divergence_eg_tes_gff_pHis %>%
  dplyr::mutate(
    # identify curated families by suffix
    is_curated = stringr::str_detect(ID, "(_PASHIS|_PASHISC|_PASDOM)$"),
    
    # clean superfamily text
    superfamily_trim = stringr::str_trim(superfamily),
    
    # if superfamily is "unspecified" (any case/whitespace) AND it's a curated "Other",
    # infer a better superfamily from the family/ID (extend these rules as needed)
    superfamily_fixed = dplyr::if_else(
      subclass == "Other" & is_curated &
        stringr::str_detect(stringr::str_to_lower(superfamily_trim), "^\\s*unspecified\\s*$"),
      dplyr::case_when(
        stringr::str_detect(ID, "ERVL") ~ "LTR_ERVL",
        stringr::str_detect(ID, "ERVK") ~ "LTR_ERVK",
        stringr::str_detect(ID, "^CR1") ~ "LINE_CR1",
        TRUE ~ superfamily_trim
      ),
      superfamily_trim
    ),
    
    # derive subclass from the prefix of the (possibly fixed) superfamily
    subclass_from_sf = toupper(sub("[-_/].*$", "", superfamily_fixed)),
    subclass_from_sf = dplyr::if_else(subclass_from_sf == "PENELOPE", "PLE", subclass_from_sf),
    
    # apply the reassignment only to curated "Other" rows, and only to valid subclasses
    subclass = dplyr::if_else(
      subclass == "Other" & is_curated & subclass_from_sf %in% valid_subclasses,
      subclass_from_sf, subclass
    ),
    
    # write back the corrected superfamily: keep only the suffix after the first separator
    superfamily = dplyr::case_when(
      grepl("_", superfamily_fixed) ~ sub("^[^_]+_", "", superfamily_fixed),  # LINE_CR1 -> CR1
      grepl("/", superfamily_fixed) ~ sub("^.*/", "", superfamily_fixed),     # DNA/PIF-Harbinger -> PIF-Harbinger
      grepl("-", superfamily_fixed) ~ sub("^[^-]+-", "", superfamily_fixed),  # LTR-ERVK -> ERVK
      TRUE ~ superfamily_fixed
    ),
    
    # keep your human-friendly label in sync
    named_subclass = dplyr::case_when(
      subclass == "DNA"  ~ "DNA Transposon",
      subclass == "LTR"  ~ "LTR Retrotransposon",
      subclass == "PLE"  ~ "Penelope",
      subclass == "RC"   ~ "Rolling Circle",
      TRUE               ~ subclass
    )
  ) %>%
  dplyr::select(-is_curated, -superfamily_trim, -superfamily_fixed, -subclass_from_sf)


# Create summary of families
# for pIta
summary_table_pIta <- divergence_eg_tes_gff_pIta %>%
  as_tibble %>%
  group_by(ID) %>%
  dplyr::select(width, KIMURA80, subclass, superfamily, ID) %>%
  mutate(total_bp = sum(width),
         mean_width = round(mean(width), digits = 2),
         min_width = round(min(width), digits = 2),
         max_width = round(max(width), digits = 2),
         sd_width = round(sd(width), digits = 2),
         div = as.numeric(KIMURA80),
         mean_div = round(mean(div), digits = 2),
         min_div = round(min(div), digits = 2),
         max_div = round(max(div), digits = 2),
         sd_div = round(sd(div), digits = 2)
  ) %>%
  dplyr::select(-width, -div, -KIMURA80) %>%
  base::unique() %>%
  dplyr::arrange(subclass, superfamily, ID) %>%
  dplyr::rename(family = ID)
# for pDom
summary_table_pDom <- divergence_eg_tes_gff_pDom %>%
  as_tibble %>%
  group_by(ID) %>%
  dplyr::select(width, KIMURA80, subclass, superfamily, ID) %>%
  mutate(total_bp = sum(width),
         mean_width = round(mean(width), digits = 2),
         min_width = round(min(width), digits = 2),
         max_width = round(max(width), digits = 2),
         sd_width = round(sd(width), digits = 2),
         div = as.numeric(KIMURA80),
         mean_div = round(mean(div), digits = 2),
         min_div = round(min(div), digits = 2),
         max_div = round(max(div), digits = 2),
         sd_div = round(sd(div), digits = 2)
  ) %>%
  dplyr::select(-width, -div, -KIMURA80) %>%
  base::unique() %>%
  dplyr::arrange(subclass, superfamily, ID) %>%
  dplyr::rename(family = ID)
# for pHis
summary_table_pHis <- divergence_eg_tes_gff_pHis %>%
  as_tibble %>%
  group_by(ID) %>%
  dplyr::select(width, KIMURA80, subclass, superfamily, ID) %>%
  mutate(total_bp = sum(width),
         mean_width = round(mean(width), digits = 2),
         min_width = round(min(width), digits = 2),
         max_width = round(max(width), digits = 2),
         sd_width = round(sd(width), digits = 2),
         div = as.numeric(KIMURA80),
         mean_div = round(mean(div), digits = 2),
         min_div = round(min(div), digits = 2),
         max_div = round(max(div), digits = 2),
         sd_div = round(sd(div), digits = 2)
  ) %>%
  dplyr::select(-width, -div, -KIMURA80) %>%
  base::unique() %>%
  dplyr::arrange(subclass, superfamily, ID) %>%
  dplyr::rename(family = ID)

# Save summaries to file
readr::write_tsv(x = summary_table_pIta, file = paste0(out_directory_pIta, "/", species_name_pIta, "_summary_table.tsv"))
readr::write_tsv(x = summary_table_pDom, file = paste0(out_directory_pDom, "/", species_name_pDom, "_summary_table.tsv"))
readr::write_tsv(x = summary_table_pHis, file = paste0(out_directory_pHis, "/", species_name_pHis, "_summary_table.tsv"))

# Sum lengths to create data for plots (remove subclasses not in standard set and repeats which Kimura was not calculated for)
divergence_eg_tes_rounded_for_plot_pIta  <- divergence_eg_tes_gff_pIta %>%
  filter(!is.na(KIMURA80)) %>%
  mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
  dplyr::filter(KIMURA80 <= 0.5) %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
  group_by(named_subclass, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(named_subclass, KIMURA80)
divergence_eg_tes_rounded_for_plot_pDom  <- divergence_eg_tes_gff_pDom %>%
  filter(!is.na(KIMURA80)) %>%
  mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
  dplyr::filter(KIMURA80 <= 0.5) %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
  group_by(named_subclass, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(named_subclass, KIMURA80)
divergence_eg_tes_rounded_for_plot_pHis  <- divergence_eg_tes_gff_pHis %>%
  filter(!is.na(KIMURA80)) %>%
  mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
  dplyr::filter(KIMURA80 <= 0.5) %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
  group_by(named_subclass, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(named_subclass, KIMURA80)

# Set fill colours
fill_colours_pIta <- tibble(subclass = c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Other", "Unknown"),
                            named_subclass = c("DNA Transposon", "LINE", "LTR Retrotransposon", "Penelope", "Rolling Circle", "SINE", "Other", "Unknown"),
                            fill_colour = c("#E32017", "#0098D4", "#00782A", "#7156A5", "#EE7C0E", "#9B0056", "#F3A9BB", "#A0A5A9")) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pIta$subclass) %>%
  arrange(named_subclass) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pIta$subclass)
fill_colours_pDom <- tibble(subclass = c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Other", "Unknown"),
                            named_subclass = c("DNA Transposon", "LINE", "LTR Retrotransposon", "Penelope", "Rolling Circle", "SINE", "Other", "Unknown"),
                            fill_colour = c("#E32017", "#0098D4", "#00782A", "#7156A5", "#EE7C0E", "#9B0056", "#F3A9BB", "#A0A5A9")) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pDom$subclass) %>%
  arrange(named_subclass) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pDom$subclass)
fill_colours_pHis <- tibble(subclass = c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Other", "Unknown"),
                            named_subclass = c("DNA Transposon", "LINE", "LTR Retrotransposon", "Penelope", "Rolling Circle", "SINE", "Other", "Unknown"),
                            fill_colour = c("#E32017", "#0098D4", "#00782A", "#7156A5", "#EE7C0E", "#9B0056", "#F3A9BB", "#A0A5A9")) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pHis$subclass) %>%
  arrange(named_subclass) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot_pHis$subclass)

# Create repeat landscape plots
# for pIta
kimura_plot_pIta <- ggplot(divergence_eg_tes_rounded_for_plot_pIta,
                           aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
  geom_col(position = "stack", width = 0.01) +
  theme_bw() +
  labs(title = plot_title_pIta) + theme(plot.title = element_markdown(hjust = 0.5)) +
  scale_fill_manual(values = fill_colours_pIta$fill_colour, name = "TE Subclass")
subclass_kimura_plot_pIta <- kimura_plot_pIta + scale_y_continuous(expand = c(0.01,0), name = "Base pairs")

subclass_kimura_plot_pIta <- subclass_kimura_plot_pIta +
  scale_x_continuous(limits = c(0.51, -0.01),
                     expand = c(0,0), name = "Kimura 2-Parameter Distance",
                     trans = "reverse")
# for pDom
kimura_plot_pDom <- ggplot(divergence_eg_tes_rounded_for_plot_pDom,
                           aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
  geom_col(position = "stack", width = 0.01) +
  theme_bw() +
  labs(title = plot_title_pDom) + theme(plot.title = element_markdown(hjust = 0.5)) +
  scale_fill_manual(values = fill_colours_pDom$fill_colour, name = "TE Subclass")
subclass_kimura_plot_pDom <- kimura_plot_pDom + scale_y_continuous(expand = c(0.01,0), name = "Base pairs")

subclass_kimura_plot_pDom <- subclass_kimura_plot_pDom +
  scale_x_continuous(limits = c(0.51, -0.01),
                     expand = c(0,0), name = "Kimura 2-Parameter Distance",
                     trans = "reverse")
# for pHis
kimura_plot_pHis <- ggplot(divergence_eg_tes_rounded_for_plot_pHis,
                           aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
  geom_col(position = "stack", width = 0.01) +
  theme_bw() +
  labs(title = plot_title_pHis) + theme(plot.title = element_markdown(hjust = 0.5)) +
  scale_fill_manual(values = fill_colours_pHis$fill_colour, name = "TE Subclass")
subclass_kimura_plot_pHis <- kimura_plot_pHis + scale_y_continuous(expand = c(0.01,0), name = "Base pairs")

subclass_kimura_plot_pHis <- subclass_kimura_plot_pHis +
  scale_x_continuous(limits = c(0.51, -0.01),
                     expand = c(0,0), name = "Kimura 2-Parameter Distance",
                     trans = "reverse")

# save repeat landscape plots
ggsave(plot = subclass_kimura_plot_pIta, filename = paste0(out_directory_pIta, "/", species_name_pIta, "_classification_landscape.png"), device = "png", width = 12.85, height = 8.5)
ggsave(plot = subclass_kimura_plot_pDom, filename = paste0(out_directory_pDom, "/", species_name_pDom, "_classification_landscape.png"), device = "png", width = 12.85, height = 8.5)
ggsave(plot = subclass_kimura_plot_pHis, filename = paste0(out_directory_pHis, "/", species_name_pHis, "_classification_landscape.png"), device = "png", width = 12.85, height = 8.5)

# ── Shared Y-axis for classification landscape plots ───────────────────────────
# Compute the cross-species maximum stacked bar height so all three plots share
# the same y-axis scale, making comparisons directly meaningful.
class_y_max <- max(
  divergence_eg_tes_rounded_for_plot_pDom %>%
    dplyr::group_by(KIMURA80) %>% dplyr::summarise(total = sum(KIMURA_SUM), .groups = "drop") %>% dplyr::pull(total),
  divergence_eg_tes_rounded_for_plot_pIta %>%
    dplyr::group_by(KIMURA80) %>% dplyr::summarise(total = sum(KIMURA_SUM), .groups = "drop") %>% dplyr::pull(total),
  divergence_eg_tes_rounded_for_plot_pHis %>%
    dplyr::group_by(KIMURA80) %>% dplyr::summarise(total = sum(KIMURA_SUM), .groups = "drop") %>% dplyr::pull(total)
)

# Apply shared y-axis via coord_cartesian (clips view after stacking, not before)
subclass_kimura_plot_pDom <- subclass_kimura_plot_pDom + coord_cartesian(ylim = c(0, class_y_max))
subclass_kimura_plot_pIta <- subclass_kimura_plot_pIta + coord_cartesian(ylim = c(0, class_y_max))
subclass_kimura_plot_pHis <- subclass_kimura_plot_pHis + coord_cartesian(ylim = c(0, class_y_max))


# ── Shared Y-axis limits per subclass for split_class_landscape plots ─────────
# Each named_subclass (facet panel) needs the same y-axis scale across the 3
# species. Compute the max KIMURA_SUM per named_subclass per species, then take
# the cross-species maximum.
shared_y_limits_split <- bind_rows(
  divergence_eg_tes_rounded_for_plot_pIta %>% dplyr::group_by(named_subclass) %>% dplyr::summarise(y_max = max(KIMURA_SUM), .groups = "drop"),
  divergence_eg_tes_rounded_for_plot_pDom %>% dplyr::group_by(named_subclass) %>% dplyr::summarise(y_max = max(KIMURA_SUM), .groups = "drop"),
  divergence_eg_tes_rounded_for_plot_pHis %>% dplyr::group_by(named_subclass) %>% dplyr::summarise(y_max = max(KIMURA_SUM), .groups = "drop")
) %>%
  dplyr::group_by(named_subclass) %>%
  dplyr::summarise(y_max = max(y_max), .groups = "drop")

# Build a helper that adds a coord_cartesian layer per facet using ggh4x if available,
# or falls back to manually adding fixed scales via a facetted_pos_scales workaround.
# The cleanest approach without extra packages: rebuild each plot adding
# dummy data rows at the shared max so ggplot naturally expands the axis to that value,
# then use scales = "free_y" so each facet still auto-sizes — but anchored by the
# dummy rows. We add one invisible point at y = shared_max for each named_subclass
# so each facet's ceiling is consistent across species.

make_shared_scale_data <- function(plot_data) {
  # For each named_subclass present in ANY species, ensure there is a row with
  # KIMURA_SUM = the cross-species y_max so the facet axis reaches that value.
  shared_y_limits_split %>%
    dplyr::rename(KIMURA_SUM = y_max) %>%
    dplyr::mutate(KIMURA80 = NA_real_,   # won't be plotted (NA x)
                  subclass = plot_data$subclass[match(named_subclass, plot_data$named_subclass)]) %>%
    dplyr::filter(named_subclass %in% plot_data$named_subclass)
}

build_split_plot <- function(base_kimura_plot, plot_data, fill_colours) {
  anchor_rows <- make_shared_scale_data(plot_data)
  augmented_data <- bind_rows(plot_data, anchor_rows)
  
  ggplot(augmented_data, aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
    geom_col(position = "stack", width = 0.01,
             data = ~ dplyr::filter(.x, !is.na(KIMURA80))) +   # skip the anchor NA rows
    geom_blank() +                                              # anchor rows still expand the scale
    theme_bw() +
    labs(title = base_kimura_plot$labels$title) +
    theme(plot.title = element_markdown(hjust = 0.5)) +
    scale_fill_manual(values = fill_colours$fill_colour,
                      breaks = fill_colours$named_subclass,
                      name = "TE Subclass") +
    scale_y_continuous(name = "Base pairs",
                       labels = function(x) format(x, scientific = TRUE),
                       expand = c(0.01, 0)) +
    facet_grid(subclass~., scales = "free_y") +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0, 0),
                       name = "Kimura 2-Parameter Distance",
                       trans = "reverse")
}

split_subclass_kimura_plot_pIta <- build_split_plot(kimura_plot_pIta, divergence_eg_tes_rounded_for_plot_pIta, fill_colours_pIta)
split_subclass_kimura_plot_pDom <- build_split_plot(kimura_plot_pDom, divergence_eg_tes_rounded_for_plot_pDom, fill_colours_pDom)
split_subclass_kimura_plot_pHis <- build_split_plot(kimura_plot_pHis, divergence_eg_tes_rounded_for_plot_pHis, fill_colours_pHis)

# save repeat landscapes split by subclass
ggsave(plot = split_subclass_kimura_plot_pIta, filename = paste0(out_directory_pIta, "/", species_name_pIta, "_split_class_landscape.png"), device = "png", width = 12.85, height = 8.5)
ggsave(plot = split_subclass_kimura_plot_pDom, filename = paste0(out_directory_pDom, "/", species_name_pDom, "_split_class_landscape.png"), device = "png", width = 12.85, height = 8.5)
ggsave(plot = split_subclass_kimura_plot_pHis, filename = paste0(out_directory_pHis, "/", species_name_pHis, "_split_class_landscape.png"), device = "png", width = 12.85, height = 8.5)


###################################################################################
########################## Panel figure: pDom / pIta / pHis  ####################
###################################################################################

# Strip the legend from each species plot (will be added below the panel),
# and remove the x-axis title from the top two plots to avoid repetition
strip_legend <- function(p) p + theme(legend.position = "none")
strip_xlab   <- function(p) p + labs(x = NULL) + theme(axis.title.x = element_blank())

plot_pDom_panel <- strip_legend(strip_xlab(split_subclass_kimura_plot_pDom))
plot_pIta_panel <- strip_legend(strip_xlab(split_subclass_kimura_plot_pIta))
plot_pHis_panel <- strip_legend(split_subclass_kimura_plot_pHis)

# Extract the shared legend from any one of the three plots
shared_legend <- cowplot::get_legend(
  split_subclass_kimura_plot_pDom +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10)
    ) +
    guides(fill = guide_legend(nrow = 1, title = "TE Subclass"))
)

# ── Vertical panel (pDom top, pIta middle, pHis bottom) ───────────────────────
plots_column <- plot_grid(
  plot_pDom_panel,
  plot_pIta_panel,
  plot_pHis_panel,
  ncol           = 1,
  align          = "v",
  axis           = "lr",
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold"
)

panel_figure_vertical <- plot_grid(
  plots_column,
  shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.06)
)

ggsave(
  plot     = panel_figure_vertical,
  filename = paste0(out_directory, "/passer_spp_split_class_landscape_panel_vertical.png"),
  device   = "png",
  width    = 12.85,
  height   = 22,
  dpi      = 300
)

# ── Horizontal panel (pDom left, pIta centre, pHis right) ─────────────────────
# For the horizontal version all three plots keep the x-axis label but only the
# leftmost keeps the y-axis title to avoid repetition.
strip_ylab <- function(p) p + labs(y = NULL) + theme(axis.title.y = element_blank())

plot_pDom_row <- strip_legend(split_subclass_kimura_plot_pDom)
plot_pIta_row <- strip_legend(strip_ylab(split_subclass_kimura_plot_pIta))
plot_pHis_row <- strip_legend(strip_ylab(split_subclass_kimura_plot_pHis))

plots_row <- plot_grid(
  plot_pDom_row,
  plot_pIta_row,
  plot_pHis_row,
  nrow           = 1,
  align          = "h",
  axis           = "tb",
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold"
)

panel_figure_horizontal <- plot_grid(
  plots_row,
  shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.04)
)

ggsave(
  plot     = panel_figure_horizontal,
  filename = paste0(out_directory, "/passer_spp_split_class_landscape_panel_horizontal.png"),
  device   = "png",
  width    = 38,
  height   = 10,
  dpi      = 300
  
  
)

###################################################################################
#################### TE counts and base pairs by subclass ########################
###################################################################################

summarise_te_counts <- function(gff, species) {
  gff %>%
    as_tibble() %>%
    dplyr::filter(!is.na(KIMURA80)) %>%
    dplyr::mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
    dplyr::filter(KIMURA80 <= 0.5) %>%
    dplyr::group_by(subclass, named_subclass) %>%
    dplyr::summarise(
      n_TEs    = dplyr::n(),
      total_bp = sum(width),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(subclass) %>%
    dplyr::mutate(species = species) %>%
    dplyr::select(species, subclass, named_subclass, n_TEs, total_bp)
}

te_count_summary_pIta <- summarise_te_counts(divergence_eg_tes_gff_pIta, "pIta")
te_count_summary_pDom <- summarise_te_counts(divergence_eg_tes_gff_pDom, "pDom")
te_count_summary_pHis <- summarise_te_counts(divergence_eg_tes_gff_pHis, "pHis")

# Print to console
print(te_count_summary_pIta)
print(te_count_summary_pDom)
print(te_count_summary_pHis)

# Write one combined file and one per species
te_count_summary_all <- bind_rows(
  te_count_summary_pIta,
  te_count_summary_pDom,
  te_count_summary_pHis
)

readr::write_tsv(
  te_count_summary_all,
  file = paste0(out_directory, "/te_counts_by_subclass_all_species.txt")
)
readr::write_tsv(
  te_count_summary_pIta,
  file = paste0(out_directory, "/te_counts_by_subclass_pIta.txt")
)
readr::write_tsv(
  te_count_summary_pDom,
  file = paste0(out_directory, "/te_counts_by_subclass_pDom.txt")
)
readr::write_tsv(
  te_count_summary_pHis,
  file = paste0(out_directory, "/te_counts_by_subclass_pHis.txt")
)


###################################################################################
#################### Bar charts: TE counts and bp by subclass ####################
###################################################################################

# Combine all species summaries (uses the objects already created above)
te_count_summary_all <- bind_rows(
  te_count_summary_pIta,
  te_count_summary_pDom,
  te_count_summary_pHis
) %>%
  dplyr::mutate(
    species = factor(species, levels = c("pDom", "pIta", "pHis")),
    named_subclass = factor(named_subclass, levels = c(
      "DNA Transposon", "LINE", "LTR Retrotransposon",
      "Penelope", "Rolling Circle", "SINE", "Other", "Unknown"
    ))
  )

# Shared species colour palette
species_colours <- c("pDom" = "#E69F00", "pIta" = "#56B4E9", "pHis" = "#009E73")

# ── Plot 1: Total base pairs per subclass per species ─────────────────────────
bp_bar_chart <- ggplot(te_count_summary_all,
                       aes(x = species, y = total_bp, fill = species)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = species_colours, name = "Species") +
  scale_y_continuous(
    name   = "Total base pairs",
    labels = function(x) format(x, scientific = TRUE),
    expand = c(0, 0)
  ) +
  facet_wrap(~ named_subclass, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_text(face = "italic"),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ── Plot 2: Number of TEs per subclass per species ────────────────────────────
count_bar_chart <- ggplot(te_count_summary_all,
                          aes(x = species, y = n_TEs, fill = species)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = species_colours, name = "Species") +
  scale_y_continuous(
    name   = "Number of TE copies",
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    expand = c(0, 0)
  ) +
  facet_wrap(~ named_subclass, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_text(face = "italic"),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ── Panel combining both charts ───────────────────────────────────────────────
te_summary_panel <- plot_grid(
  bp_bar_chart,
  count_bar_chart,
  ncol   = 1,
  labels = c("A", "B"),
  label_size     = 14,
  label_fontface = "bold",
  align  = "v",
  axis   = "lr"
)

# Save individual plots and combined panel
ggsave(
  plot     = bp_bar_chart,
  filename = paste0(out_directory, "/te_bp_by_subclass_barchart.png"),
  device   = "png",
  width    = 12.85,
  height   = 7,
  dpi      = 300
)
ggsave(
  plot     = count_bar_chart,
  filename = paste0(out_directory, "/te_count_by_subclass_barchart.png"),
  device   = "png",
  width    = 12.85,
  height   = 7,
  dpi      = 300
)
ggsave(
  plot     = te_summary_panel,
  filename = paste0(out_directory, "/te_summary_panel_barchart.png"),
  device   = "png",
  width    = 12.85,
  height   = 14,
  dpi      = 300
)


###################################################################################
########## Bar charts: bp and TE count side-by-side per species/subclass #########
###################################################################################

# Reshape to long format so bp and n_TEs can be plotted as adjacent bars
te_count_long <- te_count_summary_all %>%
  tidyr::pivot_longer(
    cols      = c(total_bp, n_TEs),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           "total_bp" = "Base pairs",
                           "n_TEs"    = "TE copies"
    ),
    # Create a combined grouping label for x-axis ordering:
    # pDom_bp, pDom_count, pIta_bp, pIta_count, pHis_bp, pHis_count
    x_group = factor(
      paste(species, metric, sep = "\n"),
      levels = c(
        "pDom\nBase pairs", "pDom\nTE copies",
        "pIta\nBase pairs", "pIta\nTE copies",
        "pHis\nBase pairs", "pHis\nTE copies"
      )
    ),
    # Colour by species x metric combination
    bar_fill = paste(species, metric, sep = "_")
  )

# Define fills: two shades per species (solid for bp, lighter for TE count)
sidebyside_colours <- c(
  "pDom_Base pairs" = "#E69F00",
  "pDom_TE copies"  = "#F5CF80",
  "pIta_Base pairs" = "#56B4E9",
  "pIta_TE copies"  = "#AAD9F4",
  "pHis_Base pairs" = "#009E73",
  "pHis_TE copies"  = "#80CEB9"
)

# Because bp and TE counts are on very different scales, use a dual-facet
# approach: one facet row per metric, grouped by TE class
sidebyside_chart <- ggplot(te_count_long,
                           aes(x = x_group, y = value, fill = bar_fill)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = sidebyside_colours,
    breaks = c(
      "pDom_Base pairs", "pDom_TE copies",
      "pIta_Base pairs", "pIta_TE copies",
      "pHis_Base pairs", "pHis_TE copies"
    ),
    labels = c(
      "pDom – Base pairs", "pDom – TE copies",
      "pIta – Base pairs", "pIta – TE copies",
      "pHis – Base pairs", "pHis – TE copies"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = TRUE, big.mark = ","),
    expand = c(0, 0)
  ) +
  facet_grid(metric ~ named_subclass, scales = "free") +
  theme_bw() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(size = 7, face = "italic"),
    strip.text.x    = element_text(face = "bold", size = 9),
    strip.text.y    = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.text     = element_text(face = "italic"),
    panel.spacing.x = unit(0.3, "lines")
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave(
  plot     = sidebyside_chart,
  filename = paste0(out_directory, "/te_bp_and_count_sidebyside_barchart.png"),
  device   = "png",
  width    = 18,
  height   = 8,
  dpi      = 300
)

###################################################################################
################ Stacked classification landscape panel: pDom/pIta/pHis ##########
###################################################################################

# The stacked (non-split) classification landscapes are already built above as
# subclass_kimura_plot_pIta/pDom/pHis. Here we assemble them into a panel figure
# equivalent to the split_class panel.

stack_strip_legend <- function(p) p + theme(legend.position = "none")
stack_strip_xlab   <- function(p) p + labs(x = NULL) + theme(axis.title.x = element_blank())
stack_strip_ylab   <- function(p) p + labs(y = NULL) + theme(axis.title.y = element_blank())

# Shared legend extracted from pDom plot
stack_shared_legend <- cowplot::get_legend(
  subclass_kimura_plot_pDom +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10)
    ) +
    guides(fill = guide_legend(nrow = 1, title = "TE Subclass"))
)

# ── Vertical panel ─────────────────────────────────────────────────────────────
stack_plots_column <- plot_grid(
  stack_strip_legend(stack_strip_xlab(subclass_kimura_plot_pDom)),
  stack_strip_legend(stack_strip_xlab(subclass_kimura_plot_pIta)),
  stack_strip_legend(subclass_kimura_plot_pHis),
  ncol           = 1,
  align          = "v",
  axis           = "lr",
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold"
)

stack_panel_vertical <- plot_grid(
  stack_plots_column,
  stack_shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.06)
)

ggsave(
  plot     = stack_panel_vertical,
  filename = paste0(out_directory, "/passer_spp_classification_landscape_panel_vertical.png"),
  device   = "png",
  width    = 12.85,
  height   = 14,
  dpi      = 300
)

# ── Horizontal panel ───────────────────────────────────────────────────────────
stack_plots_row <- plot_grid(
  stack_strip_legend(subclass_kimura_plot_pDom),
  stack_strip_legend(stack_strip_ylab(subclass_kimura_plot_pIta)),
  stack_strip_legend(stack_strip_ylab(subclass_kimura_plot_pHis)),
  nrow           = 1,
  align          = "h",
  axis           = "tb",
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold"
)

stack_panel_horizontal <- plot_grid(
  stack_plots_row,
  stack_shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

ggsave(
  plot     = stack_panel_horizontal,
  filename = paste0(out_directory, "/passer_spp_classification_landscape_panel_horizontal.png"),
  device   = "png",
  width    = 30,
  height   = 7,
  dpi      = 300
)


###################################################################################
########################## Summary pie charts #####################################
###################################################################################
# Two types of pie chart are generated:
#   1. TE-only pie: proportions of each TE subclass relative to total annotated TE bp
#   2. Genome pie: proportions of each TE subclass relative to total genome size,
#      with the non-repeat fraction shown as black (matching the reference figure).
#
# Genome sizes must be set manually — these are the total assembly sizes in bp.
# Update these values to match your actual assembly sizes.
genome_size_pDom <- 1060000000   # <-- update with actual pDom genome size in bp
genome_size_pIta <- 1060000000   # <-- update with actual pIta genome size in bp
genome_size_pHis <- 1060000000   # <-- update with actual pHis genome size in bp

# Colour palette (TE subclasses + Non-repeat shown as black for genome pies)
pie_colours <- tibble(
  tclassif = c("DNA Transposon", "Rolling Circle", "Penelope",
                "LINE", "SINE", "LTR Retrotransposon",
                "Satellite / Simple Repeat", "Unknown", "Non-repeat"),
  colour   = c("#E32017", "#EE7C0E", "#7156A5",
                "#0098D4", "#9B0056", "#00782A", "#F3A9BB", "#A0A5A9", "#000000")
)

# Factor level order for pie slices
pie_level_order <- c("DNA Transposon", "Rolling Circle", "Penelope",
                     "LINE", "SINE", "LTR Retrotransposon",
                     "Satellite / Simple Repeat", "Unknown", "Non-repeat")

# ── Helper: build per-class coverage summary from a GFF object ─────────────────
summarise_pie_data <- function(gff_tes) {
  gff_tes %>%
    as_tibble() %>%
    dplyr::mutate(named_subclass = dplyr::if_else(
      named_subclass == "Other", "Satellite / Simple Repeat", named_subclass
    )) %>%
    dplyr::group_by(named_subclass) %>%
    dplyr::summarise(coverage_bp = sum(width), .groups = "drop") %>%
    dplyr::rename(tclassif = named_subclass)
}

# ── Function: TE-only pie (proportions within annotated TEs only) ──────────────
make_pie <- function(gff_tes, species_name) {
  pie_data <- summarise_pie_data(gff_tes) %>%
    dplyr::left_join(pie_colours, by = "tclassif") %>%
    dplyr::mutate(
      tclassif   = factor(tclassif, levels = pie_level_order),
      proportion = coverage_bp / sum(coverage_bp)
    ) %>%
    dplyr::arrange(tclassif)

  col_vec        <- pie_data$colour
  names(col_vec) <- as.character(pie_data$tclassif)

  ggplot(pie_data, aes(x = "", y = proportion, fill = tclassif)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = col_vec, drop = FALSE) +
    labs(title = paste0("*", species_name, "*"), fill = "TE Subclass") +
    theme_classic() +
    theme(
      axis.line    = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_markdown(hjust = 0.5, size = 14, face = "italic"),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold")
    )
}

# ── Function: genome proportion pie (TE subclasses + non-repeat fraction) ───────
make_genome_pie <- function(gff_tes, species_name, genome_size) {
  te_data <- summarise_pie_data(gff_tes)
  total_te_bp   <- sum(te_data$coverage_bp)
  non_repeat_bp <- genome_size - total_te_bp

  pie_data <- bind_rows(
    te_data,
    tibble(tclassif = "Non-repeat", coverage_bp = max(non_repeat_bp, 0))
  ) %>%
    dplyr::left_join(pie_colours, by = "tclassif") %>%
    dplyr::mutate(
      tclassif   = factor(tclassif, levels = pie_level_order),
      proportion = coverage_bp / sum(coverage_bp)
    ) %>%
    dplyr::arrange(tclassif)

  col_vec        <- pie_data$colour
  names(col_vec) <- as.character(pie_data$tclassif)

  ggplot(pie_data, aes(x = "", y = proportion, fill = tclassif)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = col_vec, drop = FALSE) +
    labs(title = paste0("*", species_name, "*"), fill = "TE Subclass") +
    theme_classic() +
    theme(
      axis.line    = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_markdown(hjust = 0.5, size = 14, face = "italic"),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold")
    )
}

# Build TE-only pie charts
pie_pDom <- make_pie(divergence_eg_tes_gff_pDom, "pDom")
pie_pIta <- make_pie(divergence_eg_tes_gff_pIta, "pIta")
pie_pHis <- make_pie(divergence_eg_tes_gff_pHis, "pHis")

# Build genome proportion pie charts
genome_pie_pDom <- make_genome_pie(divergence_eg_tes_gff_pDom, "pDom", genome_size_pDom)
genome_pie_pIta <- make_genome_pie(divergence_eg_tes_gff_pIta, "pIta", genome_size_pIta)
genome_pie_pHis <- make_genome_pie(divergence_eg_tes_gff_pHis, "pHis", genome_size_pHis)

# Save individual TE-only pie charts
ggsave(plot = pie_pDom, filename = paste0(out_directory, "/", species_name_pDom, "_summaryPie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)
ggsave(plot = pie_pIta, filename = paste0(out_directory, "/", species_name_pIta, "_summaryPie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)
ggsave(plot = pie_pHis, filename = paste0(out_directory, "/", species_name_pHis, "_summaryPie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)

# Save individual genome proportion pie charts
ggsave(plot = genome_pie_pDom, filename = paste0(out_directory, "/", species_name_pDom, "_genomePie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)
ggsave(plot = genome_pie_pIta, filename = paste0(out_directory, "/", species_name_pIta, "_genomePie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)
ggsave(plot = genome_pie_pHis, filename = paste0(out_directory, "/", species_name_pHis, "_genomePie.png"),
       device = "png", width = 210, height = 210, units = "mm", dpi = 300)

# ── Shared legend for all pie panels ──────────────────────────────────────────
pie_shared_legend <- cowplot::get_legend(
  genome_pie_pDom +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10)
    ) +
    guides(fill = guide_legend(nrow = 2, title = "TE Subclass"))
)

# ── TE-only pie panel ──────────────────────────────────────────────────────────
pie_row <- plot_grid(
  pie_pDom + theme(legend.position = "none"),
  pie_pIta + theme(legend.position = "none"),
  pie_pHis + theme(legend.position = "none"),
  nrow = 1, labels = c("A", "B", "C"),
  label_size = 14, label_fontface = "bold"
)
pie_panel <- plot_grid(pie_row, pie_shared_legend, ncol = 1, rel_heights = c(1, 0.15))

ggsave(plot = pie_panel, filename = paste0(out_directory, "/passer_spp_summaryPie_panel.png"),
       device = "png", width = 24, height = 10, dpi = 300)

# ── Genome proportion pie panel ────────────────────────────────────────────────
genome_pie_row <- plot_grid(
  genome_pie_pDom + theme(legend.position = "none"),
  genome_pie_pIta + theme(legend.position = "none"),
  genome_pie_pHis + theme(legend.position = "none"),
  nrow = 1, labels = c("A", "B", "C"),
  label_size = 14, label_fontface = "bold"
)
genome_pie_panel <- plot_grid(genome_pie_row, pie_shared_legend, ncol = 1, rel_heights = c(1, 0.15))

ggsave(plot = genome_pie_panel, filename = paste0(out_directory, "/passer_spp_genomePie_panel.png"),
       device = "png", width = 24, height = 10, dpi = 300)

###################################################################################
########## Combined genome pie + classification landscape panel ###################
###################################################################################
# Recreates the reference figure style: genome proportion pie inset on the left,
# stacked repeat landscape on the right, one row per species (pDom=A, pIta=B, pHis=C)
# with a shared legend below. Y-axes are fixed across species.

# Strip titles from landscapes (species label comes from the panel row label)
# and remove x-axis label from top two rows
landscape_no_title <- function(p) p + labs(title = NULL) + theme(plot.title = element_blank())
landscape_no_xlab  <- function(p) p + labs(x = NULL) + theme(axis.title.x = element_blank())
landscape_no_leg   <- function(p) p + theme(legend.position = "none")

# Build one combined row per species: pie on left (~25%), landscape on right (~75%)
make_combined_row <- function(pie_plot, landscape_plot, label) {
  plot_grid(
    pie_plot + theme(legend.position = "none", plot.title = element_blank()),
    landscape_plot + theme(legend.position = "none"),
    nrow         = 1,
    rel_widths   = c(0.28, 0.72),
    labels       = c(label, ""),
    label_size   = 14,
    label_fontface = "bold"
  )
}

row_pDom <- make_combined_row(
  genome_pie_pDom,
  landscape_no_title(landscape_no_xlab(subclass_kimura_plot_pDom)),
  "A"
)
row_pIta <- make_combined_row(
  genome_pie_pIta,
  landscape_no_title(landscape_no_xlab(subclass_kimura_plot_pIta)),
  "B"
)
row_pHis <- make_combined_row(
  genome_pie_pHis,
  landscape_no_title(subclass_kimura_plot_pHis),
  "C"
)

# Shared landscape legend (TE subclass colours)
landscape_legend <- cowplot::get_legend(
  subclass_kimura_plot_pDom +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10)
    ) +
    guides(fill = guide_legend(nrow = 1, title = "TE Subclass"))
)

combined_panel <- plot_grid(
  row_pDom,
  row_pIta,
  row_pHis,
  ncol  = 1,
  align = "v",
  axis  = "lr"
)

combined_panel_with_legend <- plot_grid(
  combined_panel,
  landscape_legend,
  ncol        = 1,
  rel_heights = c(1, 0.05)
)

ggsave(
  plot     = combined_panel_with_legend,
  filename = paste0(out_directory, "/passer_spp_genomePie_landscape_panel.png"),
  device   = "png",
  width    = 14,
  height   = 18,
  dpi      = 300
)
