# Allele and allele properties visualization

color_map <- c(
  V = "#4E79A7",
  D = "#F28E2B",
  J = "#59A14F"
)

superpop_palette <- c(
  AFR = "#E15759",
  AMR = "#F28E2B",
  EAS = "#B07AA1",
  EUR = "#76B7B2",
  SAS = "#9C755F"
)




generate_all_allele_plots <- function(
  frequency_dir,
  coords_csv,
  broken_genes,
  output_dir = "./"
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)       # for unit.pmax()

  tableau10_palette <- palette.colors(10, "Classic Tableau")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  gene_coordinates <- read.csv(coords_csv, stringsAsFactors = FALSE) %>%
    select(gene_name, t2t_start) 

  segment_prefix_map <- c(V = "V")

  for (segment_label in names(segment_prefix_map)) {

    freq_pattern <- paste0("^", segment_prefix_map[[segment_label]], ".*_onefield_frequency\.csv$")
    frequency_file <- list.files(frequency_dir,
                                 pattern   = freq_pattern,
                                 full.names = TRUE)

    if (length(frequency_file) != 1) {
      warning("Skipping ", segment_label,
              ": expected 1 file, found ", length(frequency_file))
      next
    }

    allele_frequency_long_df <- read.csv(frequency_file,
                                         stringsAsFactors = FALSE) %>%
      pivot_longer(
        cols      = c(overall_frequency, freq_afr, freq_amr,
                      freq_eas, freq_eur, freq_sas),
        names_to  = "ancestry_code",
        values_to = "frequency"
      ) %>%
      mutate(
        ancestry_label = case_when(
          ancestry_code == "overall_frequency" ~ "All ancestries",
          ancestry_code == "freq_afr"          ~ "AFR",
          ancestry_code == "freq_amr"          ~ "AMR",
          ancestry_code == "freq_eas"          ~ "EAS",
          ancestry_code == "freq_eur"          ~ "EUR",
          ancestry_code == "freq_sas"          ~ "SAS"
        )
      ) %>%
      left_join(gene_coordinates, by = "gene_name") %>%
      arrange(t2t_start) %>%
      mutate(gene_name = factor(gene_name, levels = unique(gene_name)))

    
      allele_ranked_df <- allele_frequency_long_df %>%
        filter(ancestry_label == ancestry_name) %>%
        group_by(gene_name) %>%
        arrange(desc(frequency), .by_group = TRUE) %>%
        mutate(
          allele_rank    = row_number(),
          rank_factor    = factor(allele_rank, levels = seq_len(max(allele_rank))),
          stacking_order = allele_rank
        ) %>%
        ungroup() 
        

      count_bar_plot <- ggplot(common_allele_count_df,
                               aes(x = gene_name, y = common_allele_count)) +
        geom_col(width = 0.90, fill = "grey") +
        scale_y_continuous(
          expand = c(0, 0),
          breaks = scales::pretty_breaks(n = 3)
        ) +
        labs(y = "# Common Alleles") +
        theme_minimal() +
        theme(
          axis.title.x     = element_blank(),
          axis.title.y     = element_blank(), 
          axis.text.x      = element_blank(),
          axis.text.y      = element_text(size = 12),
          axis.ticks.x     = element_blank(),
          axis.ticks.y     = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          plot.margin      = margin(t = 10, r = 1, b = -2, l = 1, unit = "pt")
        )

      # Stacked allele-frequency plot with y-axis and ticks
      frequency_stack_plot <- ggplot(allele_ranked_df,
                                     aes(x = gene_name,
                                         y = frequency,
                                         fill = rank_factor,
                                         order = stacking_order)) +
        geom_col(width = 0.90,
                 colour = "white",
                 position = "stack",
                 alpha    = 0.8,
                 size     = 0.04) +
        scale_fill_manual(
          values = tableau10_palette[seq_len(max(allele_ranked_df$allele_rank))],
          guide  = "none"
        ) +
        scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.25)
        ) +
        labs(x = "V gene", y = "Allele frequency") +
        theme_minimal() +
        theme(
          axis.title.x        = element_text(size = 15),
          axis.title.y        = element_text(size = 15),
          axis.text.x         = element_text(
            angle  = 90,
            hjust  = 1,
            vjust  = 0.5,
            margin = margin(t = 10, unit = "pt"),
            size   = 7
          ),
          axis.text.y         = element_text(size = 13),
           axis.ticks.x        = element_blank(),
          axis.ticks.y        = element_line(color = "black"),
          panel.grid         = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin         = margin(t = -2, r = 1, b = 1, l = 1, unit = "pt")
        )
        

      frequency_stack_plot <- ggplot(allele_ranked_df,
                                     aes(x = gene_name,
                                         y = frequency,
                                         fill = rank_factor,
                                         order = stacking_order)) +
        geom_col(width = 0.90,
                 colour = "white",
                 position = "stack",
                 alpha    = 0.8,
                 size     = 0.04) +
        scale_fill_manual(
          values = tableau10_palette[seq_len(max(allele_ranked_df$allele_rank))],
          guide  = "none"
        ) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = "V gene", y = "Allele frequency") +
        theme_minimal() +
        theme(
          axis.title.x        = element_text(size = 23),
          axis.title.y        = element_text(size = 23), 
          axis.text.x         = element_text(
            angle  = 90,
            hjust  = 1,
            vjust  = 0.5,
            margin = margin(t = -20, unit = "pt"),
            size   = 7
          ),
          axis.text.y         = element_text(size = 13),
          axis.ticks.x        = element_blank(),
          axis.ticks.y        = element_line(color = "black"),
          panel.grid         = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin         = margin(t = -2, r = 1, b = 1, l = 1, unit = "pt")
       ) +
  scale_x_discrete(
    labels = function(gene_labels) {
      ifelse(
        grepl("TRAV", gene_labels),
        gsub("DV", "/DV", gene_labels),
        gene_labels
      )
    }
  ) 


      count_bar_grob <- ggplotGrob(count_bar_plot)
      frequency_grob <- ggplotGrob(frequency_stack_plot)
      unified_widths <- unit.pmax(count_bar_grob
widths)
      count_bar_grob
widths[] <- unified_widths

      arranged_grob <- arrangeGrob(
        count_bar_grob,
        frequency_grob,
        heights = c(0.4, 4)
      )

      grid.newpage()
      grid.draw(arranged_grob)

      ggsave(
        filename = file.path(
          output_dir,
          paste0(segment_label, "_",
                 gsub(" ", "_", ancestry_name),
                 "_diversity_with_counts.pdf")
        ),
        plot   = arranged_grob,
        width  = 20,
        height = 10
      )
    }
  
}




v_gene_palette <- c(Known = "#b8c9eb", Novel = "#4e79a7")
j_gene_palette <- c(Novel = "#59a14f", Known = "#acd0a7")
d_gene_palette <- c(Known = "#ffcb9e", Novel = "#ff7f0e")  

process_frequency_file <- function(file_path, locus, frequency_type) {
  frequency_df <- read_csv(file_path, show_col_types = FALSE)
  
  allele_name_col <- names(frequency_df)[str_detect(names(frequency_df),
                                                    "^assigned_allele.*name$")][1]
  
    print(file_path)
    frequency_df %>% filter(common == TRUE) %>%
    mutate(
      allele_status = if_else(str_detect(.data[[allele_name_col]], "N$"),
                              "Novel", "Known"),
      allele_status = factor(allele_status, levels = c("Known", "Novel"))
    ) %>%
    mutate(
      gene_groups = map(gene_name, ~{
        groups <- character()
        if (locus == "V") {
          if (str_detect(.x, "AV")) groups <- c(groups, "TRAV")
          if (str_detect(.x, "BV")) groups <- c(groups, "TRBV")
          if (str_detect(.x, "GV")) groups <- c(groups, "TRGV")
          if (str_detect(.x, "DV")) groups <- c(groups, "TRDV")
        } else if (locus == "J") {
          if (str_detect(.x, "AJ")) groups <- c(groups, "TRAJ")
          if (str_detect(.x, "BJ")) groups <- c(groups, "TRBJ")
          if (str_detect(.x, "GJ")) groups <- c(groups, "TRGJ")
          if (str_detect(.x, "DJ")) groups <- c(groups, "TRDJ")
        } else if (locus == "D") {
          if (str_detect(.x, "BD")) groups <- c(groups, "TRBD")
          if (str_detect(.x, "DD")) groups <- c(groups, "TRDD")
        }
        groups
      })
    ) %>%
    unnest(gene_groups) %>%                 
    rename(gene_group = gene_groups) %>%
    count(gene_group, allele_status, name = "allele_count") %>%
    mutate(frequency_type = frequency_type, locus = locus)
}


allele_summary <- pmap_dfr(
  file_specs,
  function(locus, frequency_type, file_name) {
    file_path <- file.path(frequency_dir, file_name)
    process_frequency_file(file_path, locus, frequency_type)
  }
)

plot_and_save <- function(summary_df, palette, title_text, out_path) {
  p <- ggplot(summary_df,
              aes(x = gene_group, y = allele_count, fill = allele_status)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = palette, name = "") +
    labs(x = "Gene Group", y = title_text) +
    theme_minimal(base_size = 25) +
    theme(
      legend.text  = element_text(size = 23),      
      legend.title = element_text(size = 23),        
      legend.key.size = unit(1.5, "cm")             
    )
  
  print(p)                        
}

file_specs %>%
  pmap(function(locus, frequency_type, file_name) {
    subset_df <- allele_summary %>%
      filter(locus == !!locus, frequency_type == !!frequency_type)
    
    palette <- if (locus == "V") {
      v_gene_palette
    } else if (locus == "J") {
      j_gene_palette
    } else {
      d_gene_palette    
    }
    
    if (frequency_type == 'onefield') {
      s <- '1-field'
    } else {
      s <- '2-field'
    }
      
    out_file <- file.path(
      output_dir,
      paste0(tolower(locus), "_gene_", frequency_type, "_common_novel_bar.pdf")
    )
    
    plot_title <- paste0("Number of Common ", s, " ", locus, " alleles")
    
    plot_and_save(subset_df, palette, plot_title, out_file)
  })


plot_shannon_entropy_by_ancestry <- function(
  frequency_dir,
  output_dir = "./"
) { 
  library(tidyverse)                        

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## 1) locate V/D/J one-field tables
  onefield_files <- list.files(
    frequency_dir,
    pattern    = "^[VDJ]_.*_onefield_frequency\.csv$",
    full.names = TRUE
  )

  ## 2) read & reshape each table
  entropy_tbl <- map_dfr(onefield_files, function(csv_path) {
    df <- read_csv(csv_path, show_col_types = FALSE)
    
    df %>% 
      pivot_longer(
        cols      = starts_with("freq_"),
        names_to  = "ancestry",
        values_to = "allele_freq",
        values_drop_na = TRUE
      ) %>%
      mutate(
        ancestry   = toupper(sub("^freq_", "", ancestry)),   # eur → EUR etc.
        gene_group = substr(gene_name, 1, 4)
      ) %>%
      ## keep ancestries that have colours defined
      filter(ancestry %in% names(superpop_palette)) %>%
      ## Shannon entropy per gene × ancestry
      group_by(gene_name, ancestry, gene_group) %>%
      summarise(
        shannon_entropy = -sum(allele_freq * log2(allele_freq), na.rm = TRUE),
        .groups         = "drop"
      )
  })
  
  ## restrict to the groups you were plotting before
  gene_groups_of_interest <- c("TRAV", "TRAJ", "TRBV", "TRBJ", "TRGV", "TRGJ")
  plot_df <- entropy_tbl %>%
    filter(gene_group %in% gene_groups_of_interest) %>%
    mutate(
      gene_group = factor(gene_group, levels = gene_groups_of_interest),
      ancestry   = factor(ancestry,     levels = names(superpop_palette))
    )
  
  ## 3) build the plot
  dodge_width <- 0.6
  p <- ggplot(plot_df,
              aes(x = gene_group,
                  y = shannon_entropy,
                  colour = ancestry,
                  fill   = ancestry)) +
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
    scale_colour_manual(values = superpop_palette, name = "Ancestry") +
    scale_fill_manual  (values = superpop_palette, name = "Ancestry") +
    labs(
      x = "Gene Group",
      y = "Shannon Entropy"
    ) +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  ggsave(
    filename = file.path(output_dir, "shannon_entropy_by_ancestry_onefield.pdf"),
    plot     = p,
    width    = 12,
    height   = 8.1,
    dpi      = 300
  )
  
  print(p)
  invisible(plot_df)
}