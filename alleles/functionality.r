

# Nonfunctional allele plot
superpop_palette <- c(
  AFR = "#E15759",
  AMR = "#F28E2B",
  EAS = "#B07AA1",
  EUR = "#76B7B2",
  SAS = "#9C755F"
)

allele_frequency_long <- v_final_freq_df %>%
  filter(gene_name %in% c("TRAV26-1", "TRDV2", "TRBV12-4")) %>%
  pivot_longer(
    cols         = starts_with("freq_"),
    names_to     = "ancestry",
    names_prefix = "freq_",
    values_to    = "allele_frequency"
  ) %>%
  mutate(
    ancestry = toupper(ancestry)
  ) %>%
  filter(ancestry %in% names(superpop_palette)) %>%
  mutate(
    ancestry = factor(ancestry, levels = names(superpop_palette))
  )


allele_frequency_plot <- ggplot(
  allele_frequency_long,
  aes(x = ancestry, y = allele_frequency, fill = ancestry)
) +
  geom_col(color = 'black') +
  scale_fill_manual(values = superpop_palette) +
  facet_wrap(~ gene_name, nrow = 1) +
  labs(
    x    = "Ancestry",
    y    = "Frequency of Nonfunctional Alleles"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.spacing = unit(2, "lines"),
    legend.position = "none",
    axis.text.x    = element_text(angle = 45, hjust = 1),
          strip.text        = element_text(size = 22)  

  )


# Functional per person
color_map <- c(
  V = "#4E79A7",
  D = "#F28E2B",
  J = "#59A14F"
)

v_genes_data <- functional_summary_data %>%
  filter(gene_group %in% c("TRAV", "TRBV"))

v_gene_histogram <- ggplot(v_genes_data, aes(x = n_functional_alleles)) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 1,
    fill     = color_map["V"],
    color    = "black"
  ) +
  facet_wrap(~ gene_group, ncol = 2, scales = "free_x") +
  labs(
    x = "Number of Distinct Functional Alleles",
    y = "Density"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.spacing = unit(2, "cm"),
    strip.text   = element_text(size = 24)
  ) +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 5) 
  )

print(v_gene_histogram)
ggsave(
  filename = file.path(output_directory, "num_functional_alleles_v_genes.pdf"),
  plot     = v_gene_histogram,
  width    = 9,
  height   = 6
)

j_genes_data <- functional_summary_data %>%
  filter(gene_group %in% c("TRAJ", "TRBJ"))

j_gene_histogram <- ggplot(j_genes_data, aes(x = n_functional_alleles)) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 1,
    fill     = color_map["J"],
    color    = "black"
  ) +
  facet_wrap(~ gene_group, ncol = 2, scales = "free_x") +
  labs(
    x = "Number of Distinct Functional Alleles",
    y = "Density"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.spacing = unit(2, "cm"),
    strip.text   = element_text(size = 24)
  )

print(j_gene_histogram)
ggsave(
  filename = file.path(output_directory, "num_functional_alleles_j_genes.pdf"),
  plot     = j_gene_histogram,
  width    = 9,
  height   = 6
)


plot_functional_alleles_by_ancestry <- function(
  functional_summary_data,
  output_directory = "./"
) {
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  
  ## ── ancestry colours ────────────────────────────────────────────────────
  superpop_palette <- c(
    AFR = "#E15759",
    AMR = "#F28E2B",
    EAS = "#B07AA1",
    EUR = "#76B7B2",
    SAS = "#9C755F"
  )
  
  plot_data <- functional_summary_data %>%
    filter(gene_group %in% gene_groups_of_interest) %>%
    mutate(
      gene_group = factor(gene_group, levels = gene_groups_of_interest),
      ancestry   = factor(ancestry,   levels = names(superpop_palette))
    )
  dodge_width <- 0.8
  full_boxplot <- ggplot(
    plot_data,
    aes(x = gene_group, y = n_functional_alleles,
        color = ancestry, fill = ancestry)
  ) +
    geom_boxplot(
      width         = dodge_width * 0.8,
      alpha         = 0.3,
      outlier.shape = NA,
      position      = position_dodge(width = dodge_width)
    ) +
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = dodge_width * 0.3,
        dodge.width  = dodge_width
      ),
      size = 2
    ) +
    scale_color_manual(values = superpop_palette, name = "Ancestry") +
    scale_fill_manual (values = superpop_palette, name = "Ancestry") +
    labs(x = "Gene Group", y = "Number of Distinct Functional Alleles") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(full_boxplot)
  ggsave(
    filename = file.path(output_directory, "functional_alleles_by_ancestry.pdf"),
    plot     = full_boxplot,
    width    = 10,
    height   = 6
  )

}

functional_alleles_plot_data <- plot_functional_alleles_by_ancestry(functional_summary_data)


# Shapley decomposition
v_gene_palette <- c(color1 = "#172432", color2 = "#4E79A7", color3 = "#b8c9db")
j_gene_palette <- c(color1 = "#172432", color2 = "#4E79A7", color3 = "#b8c9db")

gene_classes_in_order <- c("TRAV","TRBV","TRAJ","TRBJ")

predictor_map <- c(
  n_heterozygous_genes_functional = "Heterozygous functional alleles",
  total_nonfunctional_alleles     = "Nonfunctional alleles",
  total_genes_present             = "Genomic deletions"
)

compute_r2 <- function(data_frame, predictors) {
  if (length(predictors) == 0) {
    return(0)
  }
  model_formula <- as.formula(
    paste("number_distinct_functional_alleles ~", paste(predictors, collapse = " + "))
  )
  lm_fit <- lm(model_formula, data = data_frame)
  summary(lm_fit)$r.squared
}

shapley_three <- function(data_frame, a, b, c) {
  # Singletons
  R2_a <- compute_r2(data_frame, c(a))
  R2_b <- compute_r2(data_frame, c(b))
  R2_c <- compute_r2(data_frame, c(c))
  # Pairs
  R2_ab <- compute_r2(data_frame, c(a,b))
  R2_ac <- compute_r2(data_frame, c(a,c))
  R2_bc <- compute_r2(data_frame, c(b,c))
  # Full
  R2_abc <- compute_r2(data_frame, c(a,b,c))

  phi_a <- (1/3)*R2_a + (1/6)*(R2_ab - R2_b + R2_ac - R2_c) + (1/3)*(R2_abc - R2_bc)
  phi_b <- (1/3)*R2_b + (1/6)*(R2_ab - R2_a + R2_bc - R2_c) + (1/3)*(R2_abc - R2_ac)
  phi_c <- (1/3)*R2_c + (1/6)*(R2_ac - R2_a + R2_bc - R2_b) + (1/3)*(R2_abc - R2_ab)

  list(phi = c(a = phi_a, b = phi_b, c = phi_c), R2 = R2_abc)
}

shapley_results <- list()

for (gc in gene_classes_in_order) {
  class_subset <- all_class_summary %>%
    dplyr::select(number_distinct_functional_alleles,
                  n_heterozygous_genes_functional,
                  total_nonfunctional_alleles,
                  total_genes_present) %>%
    tidyr::drop_na()

  res <- shapley_three(
    data_frame = class_subset,
    a = "n_heterozygous_genes_functional",
    b = "total_nonfunctional_alleles",
    c = "total_genes_present"
  )


  tidy_res <- tibble(
    GeneClass = gc,
    Predictor_raw = c("n_heterozygous_genes_functional",
                      "total_nonfunctional_alleles",
                      "total_genes_present"),
    ShareR2 = as.numeric(shares),
    ModelR2 = model_r2
  ) %>%
    mutate(Predictor = dplyr::recode(Predictor_raw, !!!predictor_map)) %>%
    dplyr::select(GeneClass, Predictor, ShareR2, ModelR2)

  shapley_results[[gc]] <- tidy_res
}

shapley_table <- dplyr::bind_rows(shapley_results)

shapley_percent <- shapley_table %>%
  mutate(Percent = ifelse(!is.na(ModelR2) & ModelR2 > 0, ShareR2 / ModelR2 * 100, NA_real_),
         GeneClass = factor(GeneClass, levels = gene_classes_in_order)) 

shapley_plot <- ggplot(shapley_percent,
                       aes(x = GeneClass, y = Percent, fill = fill_color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity(
    name = "Component",
    breaks = c(v_gene_palette["color1"], v_gene_palette["color2"], v_gene_palette["color3"]),
    labels = c("Heterozygous functional alleles", "Genomic deletions", "Nonfunctional alleles"),
    guide  = "legend"
  ) +
  labs(x = "Gene Group", y = "% Variance Explained") +
  theme_minimal(base_size = 23) +
  scale_y_continuous(labels = percent_format(scale = 1))


ggsave("shapley_variance_explained.pdf",
       shapley_plot,
       width = 12, height = 7.5)