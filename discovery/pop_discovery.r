
# Dendrograms
build_nj_tree_from_alignment <- function(alignment_fasta, output_pdf) {
  # Dependencies
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required. Install with install.packages('ape').")
  }

  # --- 1) Read aligned sequences (FASTA) as DNAbin -------------------------
  if (!file.exists(alignment_fasta)) {
    stop("Input alignment file not found: ", alignment_fasta)
  }
  message("Reading alignment: ", alignment_fasta)
  aligned_dnabin <- ape::read.dna(alignment_fasta, format = "fasta")  # assumes pre-aligned FASTA
  sequence_count <- nrow(aligned_dnabin)
  # if (is.null(sequence_count) || sequence_count < 3) {
  #   stop("Need at least 3 aligned sequences to build an NJ tree.")
  # }

  # --- 2) TN93 distance matrix --------------------------------------------
  message("Computing TN93 distance matrix (pairwise deletion for gaps) …")
  tn93_dist <- ape::dist.dna(
    aligned_dnabin,
    model = "TN93",
    pairwise.deletion = TRUE,
    as.matrix = FALSE
  )

  # --- 3) Neighbor-joining tree -------------------------------------------
  message("Building neighbor-joining tree …")
  nj_tree <- ape::njs(tn93_dist)

  # --- 4) Plot to screen ---------------------------------------------------
  message("Plotting tree to screen …")
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(1, 1, 1, 1))
  ape::plot.phylo(nj_tree, tip.color = "#4e79a7", type = 'radial', cex = 0.7, lwd = 1.3, no.margin = FALSE, #, type = 'fan'
                 open.angle = 0,
      rotate.tree = 0)

  # --- 5) Save to PDF ------------------------------------------------------
  message("Saving tree to PDF: ", output_pdf)
  tree_radius <- max(ape::node.depth.edgelength(nj_tree))
label_offset <- 0.02 * tree_radius    # push labels ~2% beyond the tips
plot_radius  <- 1.15 * tree_radius    # add 15% padding so labels aren't clipped

pdf(output_pdf, width = 9, height = 9, useDingbats = FALSE)
par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = NA)  # tiny margins; allow draw into margins

# ape::plot.phylo(
#   nj_tree,
#   type         = "fan",
#   tip.color    = "#4E79A7",
#   cex          = 0.6,          # smaller labels
#   lwd          = 1.2,
#   label.offset = label_offset, # push labels outward
#   no.margin    = TRUE,
#   x.lim        = c(-plot_radius, plot_radius),
#   y.lim        = c(-plot_radius, plot_radius),
#   open.angle   = 10            # small gap can further reduce crowding/clipping
# )
    
    ape::plot.phylo(
  nj_tree,
  tip.color    = "#4E79A7", type = 'radial',
  cex          = 0.6,          # smaller labels
  lwd          = 1.2,
  no.margin    = TRUE)


    
    # plot.phylo(nj_tree, type = "fan",
    #          tip.color = tip_col,
    #          cex       = text_cex,
    #          lwd       = 1.5,
    #          no.margin = FALSE)
  dev.off()

  invisible(nj_tree)
}




# variance plotting
atchley_table <- read.csv(
  "achley_factors.csv",
  stringsAsFactors = FALSE
)

# Atchley factors (title‑case) and full IMGT position range
atchley_factors <- c("Hydrophobicity", "Structure", "Size", "Refractivity", "Charge")
all_positions   <- seq(min(imgt_regions$start), max(imgt_regions$end))

# Process each annotation file
for (file_name in annotation_files) {
  # 1. Read & filter
  annotation_data      <- read.csv(file.path(input_dir, file_name), stringsAsFactors = FALSE)
  functional_annotations <- annotation_data %>%
    filter(functional == "Functional") %>%
    distinct(assigned_allele_onefield_name, .keep_all = TRUE)

  # 2. Extract, normalize, & pad sequences
  raw_sequences        <- functional_annotations$`v.gene_aligned_aa`
  normalized_sequences <- gsub("\\.", "-", raw_sequences)
  max_length           <- max(nchar(normalized_sequences))
  padded_sequences     <- vapply(normalized_sequences, function(seq) {
    chars <- strsplit(seq, "")[[1]]
    length(chars) <- max_length
    chars[is.na(chars)] <- "-"
    paste0(chars, collapse = "")
  }, character(1))

  # 3. Build & prune AA matrix
  amino_acid_matrix       <- do.call(rbind, strsplit(padded_sequences, split = ""))
  cols_with_data          <- apply(amino_acid_matrix, 2, function(col) any(col != "-"))
  pruned_amino_acid_matrix <- amino_acid_matrix[, cols_with_data, drop = FALSE]
  positions_kept          <- which(cols_with_data)

  # 4. Compute per‑position variance for each factor
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

  # 5. Pivot to long, replace invariant NAs with zero, title‑case factor
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

  # 6. Build full grid & join (so missing positions get NA)
  full_grid_df <- expand.grid(
    position = all_positions,
    factor   = atchley_factors,
    stringsAsFactors = FALSE
  )
  full_long_df <- full_grid_df %>%
    left_join(raw_long_df, by = c("position", "factor")) %>%
    # assign region based on IMGT ranges
    mutate(
      region = map_chr(position, function(pos) {
        idx <- which(pos >= imgt_regions$start & pos <= imgt_regions$end)
        imgt_regions$region[idx]
      })
    )
                                                                      
  # 7. Z‑score within each factor
  zscored_df <- full_long_df %>%
    group_by(factor) %>%
    mutate(variance_z = (variance - min(variance, na.rm = TRUE)) /
                              (max(variance, na.rm = TRUE) - min(variance, na.rm = TRUE))) %>%
                                
    ungroup() %>%
    mutate(
      factor = factor(factor, levels = atchley_factors)
    )

  # 8. Plot: heatmap with dark grey for NAs, no vertical grid lines
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
      inherit.aes = FALSE, size = 7 #, fontface    = "bold"
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
      paste0(tools::file_path_sans_ext(file_name), "_heatmap.pdf")
    ),
    plot   = heatmap_plot,
    width  = 15,
    height = 7
  )


}