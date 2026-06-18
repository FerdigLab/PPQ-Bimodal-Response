setwd("path to metabolomics file")


library(tidyverse)
library(ggrepel)

#-----------------------------
# Load and reshape data
#-----------------------------

df_wide <- read.csv("PPQ_Bump_asp-norm-areas.csv", check.names = FALSE)

df_wide <- df_wide[rowSums(is.na(df_wide) | df_wide == "") != ncol(df_wide), ]

df_wide <- df_wide %>%
  filter(!is.na(compound), compound != "")

df_long <- df_wide %>%
  pivot_longer(
    cols = -compound,
    names_to = "sample",
    values_to = "value"
  ) %>%
  filter(!str_detect(sample, "^QC"))

df_long <- df_long %>%
  mutate(
    sample_clean = str_remove(sample, "^[0-9]+-"),
    sample_clean = str_replace(sample_clean, "ATQ2$", "ATQ-2")
  ) %>%
  separate(
    sample_clean,
    into = c("line", "treatment", "replicate"),
    sep = "-",
    fill = "right",
    extra = "merge"
  ) %>%
  mutate(
    group = paste(line, treatment, sep = "-")
  )

df_for_volcano <- df_long %>%
  select(group, sample, compound, value) %>%
  pivot_wider(
    names_from = compound,
    values_from = value
  ) %>%
  select(-sample)

#-----------------------------
# Volcano plot function
#-----------------------------

make_volcano <- function(
    data,
    group1,
    group2,
    label_sig = TRUE,
    p_cutoff = 0.05,
    log2fc_cutoff = 1
) {
  
  g1 <- data %>% filter(group == group1)
  g2 <- data %>% filter(group == group2)
  
  if (nrow(g1) == 0 || nrow(g2) == 0) {
    stop("One or both groups not found: ", group1, " vs ", group2)
  }
  
  metabolite_cols <- setdiff(colnames(data), "group")
  
  g1[metabolite_cols] <- lapply(g1[metabolite_cols], as.numeric)
  g2[metabolite_cols] <- lapply(g2[metabolite_cols], as.numeric)
  
  results <- lapply(metabolite_cols, function(met) {
    
    v1 <- as.numeric(g1[[met]])
    v2 <- as.numeric(g2[[met]])
    
    v1 <- v1[is.finite(v1) & v1 > 0]
    v2 <- v2[is.finite(v2) & v2 > 0]
    
    if (length(v1) < 2 || length(v2) < 2) {
      return(data.frame(
        metabolite = met,
        log2FC = NA,
        pvalue = NA
      ))
    }
    
    v1_log <- log2(v1)
    v2_log <- log2(v2)
    
    log2FC <- mean(v1_log, na.rm = TRUE) - mean(v2_log, na.rm = TRUE)
    
    pval <- tryCatch(
      t.test(v1_log, v2_log, var.equal = FALSE)$p.value,
      error = function(e) NA
    )
    
    data.frame(
      metabolite = met,
      log2FC = log2FC,
      pvalue = pval
    )
  }) %>%
    bind_rows()
  
  results <- results %>%
    mutate(
      neglog10p = -log10(pvalue),
      significant = !is.na(pvalue) &
        pvalue < p_cutoff &
        abs(log2FC) >= log2fc_cutoff,
      point_color = ifelse(significant, "significant", "not_significant")
    )
  
  p <- ggplot(results, aes(x = log2FC, y = neglog10p, color = point_color)) +
    geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "not_significant" = "black",
        "significant" = "darkgreen" #yellow for ATQ, red for KH004-D9, blue for KH004-E9, and darkgreen for FG0305
      )
    ) +
    geom_vline(xintercept = 0, linewidth = 1.2) +
    geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      linetype = "dashed",
      linewidth = 0.7
    ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      linewidth = 0.7
    ) +
    coord_cartesian(xlim = c(-3, 3), ylim = c(0, 3)) +
    labs(
      title = paste(group1, "vs", group2),
      x = "log2(2uM PPQ/No Drug)", #Change based on condition
      y = "-log10(p-value)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.line.y = element_blank(),
      legend.position = "none"
    )
  
  if (label_sig) {
    sig_labels <- results %>% filter(significant)
    
    if (nrow(sig_labels) > 0) {
      p <- p +
        geom_text_repel(
          data = sig_labels,
          aes(label = metabolite),
          size = 3,
          max.overlaps = 50,
          color = "black"
        )
    }
  }
  
  return(list(results = results, plot = p))
}

#-----------------------------
# Run comparison
#-----------------------------

volc2 <- make_volcano(
  df_for_volcano,
  group1 = "FG0305-2uM", #conditions are either ND, 200nM, 500nM, or 2uM for either D9, E9, or FG0305
  group2 = "FG0305-ND",
  label_sig = TRUE,
  p_cutoff = 0.05,
  log2fc_cutoff = 1
)

print(volc2$plot)

write.csv(
  volc2$results,
  "FG0305_2uM_vs_FG0305_ND_volcano_results.csv",
  row.names = FALSE
)
