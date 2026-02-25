library(dplyr)
library(stringr)

# Identify V genes whose alleles differ by non-CDR3 residues
amino_acid_positions_to_keep <- 104L
nucleotides_to_keep <- amino_acid_positions_to_keep * 3L

vgenealleles_with_no_cdr3 <- vgenealleles %>%
  mutate(
    seq_gapped_without_cdr3 = if_else(
      is.na(seq_gapped),
      NA_character_,
      str_sub(seq_gapped, 1, nucleotides_to_keep)
    ),
    seq_without_cdr3 = if_else(
      is.na(seq_gapped_without_cdr3),
      NA_character_,
      str_replace_all(seq_gapped_without_cdr3, "\\.", "")
    )
  )

gene_level_distinctness <- vgenealleles_with_no_cdr3 %>%
  filter(!is.na(gene_name), !is.na(assigned_allele_onefield_name), !is.na(seq_without_cdr3)) %>%
  distinct(gene_name, assigned_allele_onefield_name, seq_without_cdr3) %>%
  group_by(gene_name) %>%
  summarise(
    n_distinct_alleles = n(),
    n_distinct_sequences = n_distinct(seq_without_cdr3),
    all_alleles_have_unique_seq_without_cdr3 = (n_distinct_sequences == n_distinct_alleles),
    .groups = "drop"
  )

percent_genes_all_unique <- mean(
  gene_level_distinctness$all_alleles_have_unique_seq_without_cdr3
) * 100



extract_gene_name_columns <- function(tcr_dataframe) {
  output_dataframe <- tcr_dataframe
  best_hit_columns <- grep("_best_hit$", names(output_dataframe), value = TRUE)

  for (allele_column in best_hit_columns) {
    gene_column <- sub("vallele", "vgene", sub("jallele", "jgene", allele_column))
    output_dataframe[[gene_column]] <- sub("\\*.*$", "", output_dataframe[[allele_column]])
  }

  output_dataframe
}

make_clonotype_df <- function(input_df, dataset_name) {

  message("Before filtering:\n",
          "  Rows:  ", nrow(input_df), "\n",
          "  cells: ", dplyr::n_distinct(input_df$singlecell_identifier), "\n",
          "  TRA:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRA")$singlecell_identifier), "\n",
          "  TRB:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRB")$singlecell_identifier), "\n")

  input_df <- input_df %>%
    dplyr::filter(vgene_best_hit %in% validvgenes)

  message("after valid v gene:\n",
          "  Rows:  ", nrow(input_df), "\n",
          "  cells: ", dplyr::n_distinct(input_df$singlecell_identifier), "\n",
          "  TRA:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRA")$singlecell_identifier), "\n",
          "  TRB:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRB")$singlecell_identifier), "\n")

  
  input_df <- input_df %>% dplyr::filter(vallele_perfect_hit == TRUE)
  
  message("after allele filter:\n",
          "  Rows:  ", nrow(input_df), "\n",
          "  cells: ", dplyr::n_distinct(input_df$singlecell_identifier), "\n",
          "  TRA:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRA")$singlecell_identifier), "\n",
          "  TRB:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRB")$singlecell_identifier), "\n")

  input_df <- input_df %>%
    dplyr::add_count(donor_identifier, vallele_best_hit, name = "n_v_allele_in_donor") %>%
    dplyr::filter(n_v_allele_in_donor >= 3) %>%
    dplyr::select(-n_v_allele_in_donor)

  message("after filtering 3 cells per donor per V-allele:\n",
          "  Rows:  ", nrow(input_df), "\n",
          "  cells: ", dplyr::n_distinct(input_df$singlecell_identifier), "\n",
          "  TRA:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRA")$singlecell_identifier), "\n",
          "  TRB:   ", dplyr::n_distinct(dplyr::filter(input_df, chain_type == "TRB")$singlecell_identifier), "\n")

  cell_chain_counts <- input_df %>%
    dplyr::group_by(singlecell_identifier) %>%
    dplyr::summarise(
      num_TRA = sum(chain_type == "TRA"),
      num_TRB = sum(chain_type == "TRB"),
      .groups = "drop"
    )

  valid_cells <- cell_chain_counts %>%
    dplyr::filter(num_TRA == 1, num_TRB == 1) %>%
    dplyr::pull(singlecell_identifier)

  input_df <- input_df %>%
    dplyr::filter(singlecell_identifier %in% valid_cells)

  message("after filtering exactly 1 productive TRA and 1 productive TRB:\n",
          "  Rows:  ", nrow(input_df), "\n",
          "  cells: ", dplyr::n_distinct(input_df$singlecell_identifier), "\n")

  repertoire_data <- input_df %>%
    dplyr::mutate(
      cdr3   = sub(".*\\|", "", metadata),
      cdr3nt = sub(".*\\|([^|]+)\\|[^|]+$", "\\1", metadata)
    )

  chain_specific_columns <- c(
    "cdr3", "cdr3nt",
    "vallele_best_hit", "dallele_best_hit", "jallele_best_hit",
    "vgene_best_hit", "dgene_best_hit", "jgene_best_hit"
  )

  cells_one_per_cell <- repertoire_data %>%
    tidyr::pivot_wider(
      id_cols     = c(singlecell_identifier, cluster_identifier, donor_identifier, hascovid),
      names_from  = chain_type,
      values_from = dplyr::all_of(chain_specific_columns),
      names_glue  = "{chain_type}_{.value}"
    )

  clonotype_table <- cells_one_per_cell %>%
    dplyr::mutate(
      clonotype_key = paste(
        TRA_vallele_best_hit, TRA_jallele_best_hit, TRA_cdr3nt,
        TRB_vallele_best_hit, TRB_jallele_best_hit, TRB_cdr3nt,
        sep = "|"
      )
    ) %>%
    dplyr::group_by(donor_identifier, clonotype_key) %>%
    dplyr::slice_sample(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-clonotype_key)

  message("after sampling 1 cell per clonotype (within donor):\n",
          "  Rows:  ", nrow(clonotype_table), "\n",
          "  cells: ", dplyr::n_distinct(clonotype_table$singlecell_identifier), "\n")

  clonotype_table
}




# Plot donor detection rate vs AF in EAS
compute_donor_level_allele_frequencies <- function(tcr_table,
                                                  gene_type_label,
                                                  vgene_column,
                                                  vallele_column,
                                                  donor_id_column = "donor_identifier") {

  donor_gene_allele_presence <- tcr_table %>%
    transmute(
      donor_id = .data[[donor_id_column]],
      gene_type = gene_type_label,
      gene_name = .data[[vgene_column]],
      assigned_allele_onefield_name = .data[[vallele_column]]
    ) %>%
    filter(!is.na(donor_id), !is.na(gene_name), !is.na(assigned_allele_onefield_name)) %>%
    dplyr::distinct(donor_id, gene_type, gene_name, assigned_allele_onefield_name)

  total_donors_with_gene <- donor_gene_allele_presence %>%
    distinct(donor_id, gene_type, gene_name) %>%
    dplyr::count(gene_type, gene_name, name = "total_num_donors_with_gene")

  donors_with_allele <- donor_gene_allele_presence %>%
    dplyr::count(gene_type, gene_name, assigned_allele_onefield_name, name = "num_donors_with_allele")

  donors_with_allele %>%
    inner_join(total_donors_with_gene, by = c("gene_type", "gene_name")) %>%
    dplyr::mutate(allele_frequency = num_donors_with_allele / total_num_donors_with_gene) %>%
    dplyr::arrange(gene_type, gene_name, desc(allele_frequency), assigned_allele_onefield_name)
}

donor_level_allele_frequencies <- bind_rows(
  compute_donor_level_allele_frequencies(
    tcr_table      = edahiro_data,
    gene_type_label = "TRA",
    vgene_column    = "TRA_vgene_best_hit",
    vallele_column  = "TRA_vallele_best_hit"
  ),
  compute_donor_level_allele_frequencies(
    tcr_table      = edahiro_data,
    gene_type_label = "TRB",
    vgene_column    = "TRB_vgene_best_hit",
    vallele_column  = "TRB_vallele_best_hit"
  )
)

merged_for_plot <- donor_level_allele_frequencies %>%
  inner_join(
    allele_frequencies %>%
      transmute(
        assigned_allele_onefield_name,
        allele_haplotypes_eas = as.integer(allele_haplotypes_eas),
        total_eas_haplotypes  = as.integer(total_eas_haplotypes),
        freq_eas              = as.numeric(freq_eas)
      ),
    by = "assigned_allele_onefield_name"
  )

merged_for_plot_with_ci <- merged_for_plot %>%
  rowwise() %>%
  mutate(
    donor_ci = list(binom.test(num_donors_with_allele, total_num_donors_with_gene)$conf.int),
    donor_ci_low  = list(binom.test(num_donors_with_allele, total_num_donors_with_gene)$conf.int)[[1]][1],
    donor_ci_high = list(binom.test(num_donors_with_allele, total_num_donors_with_gene)$conf.int)[[1]][2],
    donor_width = donor_ci_high - donor_ci_low,

    eas_ci = list(binom.test(allele_haplotypes_eas, total_eas_haplotypes)$conf.int),
    eas_ci_low  = list(binom.test(allele_haplotypes_eas, total_eas_haplotypes)$conf.int)[[1]][1],
    eas_ci_high = list(binom.test(allele_haplotypes_eas, total_eas_haplotypes)$conf.int)[[1]][2]
  ) %>%
  ungroup() %>%
  select(-donor_ci, -eas_ci)

allele_freq_comparison_plot <- ggplot(
 merged_for_plot_with_ci,
  aes(x = freq_eas, y = allele_frequency)
) +
  geom_errorbarh(aes(xmin = eas_ci_low, xmax = eas_ci_high), height = 0) +
  geom_errorbar(aes(ymin = donor_ci_low, ymax = donor_ci_high), width = 0) +
  geom_point() +
  theme_minimal(base_size = 20) +
  labs(
    x = "Allele frequency in EAS ancestry",
    y = "Donor detection rate"
  ) 