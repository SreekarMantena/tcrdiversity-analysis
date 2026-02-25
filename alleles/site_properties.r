# Code for visualizing IMGT site based properties

# Define IMGT regions
imgt_regions <- tibble::tribble(
  ~region, ~start, ~end,
  "FR1",   1,      26,
  "CDR1",  27,     38,
  "FR2",   39,     55,
  "CDR2",  56,     65,
  "FR3",   66,     104,
  "CDR3",  105,    116
)

atchley_factors <- c("Hydrophobicity", "Structure", "Size", "Refractivity", "Charge")
all_positions   <- seq(min(imgt_regions$start), max(imgt_regions$end))


bootstrap_ci <- function(x, nboot = 500, conf = 0.95) {
  boots <- replicate(nboot,
    mean(sample(x, length(x), replace = TRUE), na.rm = TRUE)
  )
  alpha <- (1 - conf) / 2
  quantile(boots, probs = c(alpha, 1 - alpha), na.rm = TRUE)
}

for (file_name in annotation_files) {
  annotation_data      <- read.csv(file.path(input_dir, file_name), stringsAsFactors = FALSE)
  functional_annotations <- annotation_data %>%
    distinct(assigned_allele_onefield_name, .keep_all = TRUE)

  raw_sequences        <- functional_annotations$`v.gene_aligned_aa`
  normalized_sequences <- gsub("\.", "-", raw_sequences)
  max_length           <- max(nchar(normalized_sequences))
  padded_sequences     <- vapply(normalized_sequences, function(seq) {
    chars <- strsplit(seq, "")[[1]]
    length(chars) <- max_length
    chars[is.na(chars)] <- "-"
    paste0(chars, collapse = "")
  }, character(1))

  amino_acid_matrix       <- do.call(rbind, strsplit(padded_sequences, split = ""))
  cols_with_data          <- apply(amino_acid_matrix, 2, function(col) any(col != "-"))
  pruned_amino_acid_matrix <- amino_acid_matrix[, cols_with_data, drop = FALSE]
  positions_kept          <- which(cols_with_data)

  position_variances_df <- lapply(seq_len(ncol(pruned_amino_acid_matrix)), function(j) {
    residues     <- pruned_amino_acid_matrix[, j]
    value_table  <- atchley_table[match(residues, atchley_table$symbol),
      c("hydrophobicity", "structure", "size", "refractivity", "charge")
    ]
    apply(value_table, 2, var, na.rm = TRUE)
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    mutate(position = positions_kept)

  raw_long_df <- position_variances_df %>%
    pivot_longer(
      cols = -position,
      names_to  = "factor",
      values_to = "variance"
    ) %>%
    mutate(
      variance = replace_na(variance, 0),
      factor   = str_to_title(factor)
    )

  full_grid_df <- expand.grid(
    position = all_positions,
    factor   = atchley_factors,
    stringsAsFactors = FALSE
  )
  full_long_df <- full_grid_df %>%
    left_join(raw_long_df, by = c("position", "factor")) %>%
    mutate(
      region = map_chr(position, function(pos) {
        idx <- which(pos >= imgt_regionsstart & pos <= imgt_regions
end)
        imgt_regions$region[idx]
      })
    )
                                                
  zscored_df <- full_long_df %>%
    group_by(factor) %>%
    mutate(variance_z = (variance - min(variance, na.rm = TRUE)) /
                              (max(variance, na.rm = TRUE) - min(variance, na.rm = TRUE))) %>%
                                
    ungroup() %>%
    mutate(
      factor = factor(factor, levels = atchley_factors)
    )

  max_position <- max(all_positions)
  num_factors  <- length(atchley_factors)
                                   
   tick_positions <- c(1, seq(25, max(imgt_regions$end), by = 25))
                                
  heatmap_plot <- ggplot(zscored_df, aes(x = position, y = factor, fill = variance_z)) +
    geom_tile() +
    scale_fill_viridis_c(
      name    = "Normalized \n variance",
      na.value = "lightgrey"
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(min(all_positions), max_position),
      breaks = tick_positions
    ) +
    scale_y_discrete(expand = expansion(add = c(0, 1.5))) +
    labs(x = "IMGT position", y = "Atchley factor") +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x          = element_text(angle = 45, hjust = 1),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor.x   = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor.y   = element_blank(),
      legend.position = 'right'
    ) +
    # region bars above heatmap
    geom_segment(
      data       = imgt_regions,
      aes(x = start, xend = end,
          y = num_factors + 0.53,
          yend = num_factors + 0.53),
      inherit.aes = FALSE,
      size        = 0.8
    ) +
    geom_text(
      data       = imgt_regions,
      aes(x = (start + end)/2,
          y = num_factors + 0.70,
          label = region),
      inherit.aes = FALSE, size = 4.5,
      fontface    = "bold"
    )  +
    theme(
      legend.position      = "none"
    ) +
    # tiles for positions > 100, filled white
    geom_tile(
      data = filter(zscored_df, position > 109),
      aes(x = position, y = factor),
      fill = "white"
    ) 
                                   
  print(heatmap_plot)
  ggsave(
    filename = file.path(
      output_dir,
      paste0(tools::file_path_sans_ext(file_name), "allfunctional_heatmap.pdf")
    ),
    plot   = heatmap_plot,
    width  = 11,
    height = 7
  )
}




library(ape)

build_nj_tree_from_alignment <- function(alignment_fasta, output_pdf) {

  message("Reading alignment: ", alignment_fasta)
  aligned_dnabin <- ape::read.dna(alignment_fasta, format = "fasta") 
  sequence_count <- nrow(aligned_dnabin)

  tn93_dist <- ape::dist.dna(
    aligned_dnabin,
    model = "TN93",
    pairwise.deletion = TRUE,
    as.matrix = FALSE
  )

  nj_tree <- ape::njs(tn93_dist)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(1, 1, 1, 1))
  ape::plot.phylo(nj_tree, tip.color = "#4e79a7", type = 'radial', cex = 0.7, lwd = 1.3, no.margin = FALSE, 
                 open.angle = 0,
      rotate.tree = 0)

  message("Saving tree to PDF: ", output_pdf)
  tree_radius <- max(ape::node.depth.edgelength(nj_tree))
label_offset <- 0.02 * tree_radius    
plot_radius  <- 1.15 * tree_radius    

pdf(output_pdf, width = 9, height = 9, useDingbats = FALSE)
par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = NA)  
    
    ape::plot.phylo(
  nj_tree,
  tip.color    = "#4E79A7", type = 'radial',
  cex          = 0.6,         
  lwd          = 1.2,
  no.margin    = TRUE)

  dev.off()

  invisible(nj_tree)
}