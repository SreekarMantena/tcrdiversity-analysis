library(tidyverse)
library(lme4)

control_bobyqa <- glmerControl(
  optimizer = "nlminbwrap",
  optCtrl   = list(iter.max = 5e8, eval.max = 5e8)
)

perform_genewise_mixed_model <- function(
  vmodeling_dataframe,
  outcome_column_name,
  valid_genes_by_chain,
  gene_type
) {
  # ensure gene_type is valid
  if (!gene_type %in% c("V", "J")) {
    stop("`gene_type` must be either 'V' or 'J'")
  }
  
  # choose genes to test
  genes_to_test <- if (gene_type == "V") {
    c(valid_genes_by_chain$TRAV, valid_genes_by_chain$TRBV)
  } else {
    c(valid_genes_by_chain$TRAJ, valid_genes_by_chain$TRBJ)
  }
  
  results_list <- vector("list", length(genes_to_test))
    
    nlopt_options <- list(
      algorithm = "NLOPT_LN_BOBYQA",
      xtol_rel   = 1e-6,
      maxeval    = 1e5
    )
  
  for (i in seq_along(genes_to_test)) {
    gene_name <- genes_to_test[i]
    
    # infer chain letter (A or B) from gene prefix
    chain_letter <- substr(gene_name, 3, 3)  # "A" for TRA*, "B" for TRB*
    
    # dynamic column names
    if (gene_type == "V") {
      allele_col <- if (chain_letter == "A") "TRA_vallele_best_hit" else "TRB_vallele_best_hit"
      gene_col   <- if (chain_letter == "A") "TRA_vgene_best_hit"     else "TRB_vgene_best_hit"
    } else {
      allele_col <- if (chain_letter == "A") "TRA_jallele_best_hit" else "TRB_jallele_best_hit"
      gene_col   <- if (chain_letter == "A") "TRA_jgene_best_hit"   else "TRB_jgene_best_hit"
    }
    
    # subset to rows carrying this gene
    subset_for_gene <- vmodeling_dataframe %>%
      filter(.data[[gene_col]] == gene_name)
    
    n_obs     <- nrow(subset_for_gene)
    n_alleles <- n_distinct(subset_for_gene[[allele_col]])
        
    # build formulas
    null_formula   <- as.formula(paste0(outcome_column_name, " ~ hascovid + dataset + Age + Sex + (1 | donor_identifier)"))
    allele_formula <- as.formula(paste0(outcome_column_name, " ~  hascovid + dataset + Age + Sex + ",
                                        allele_col, " + (1 | donor_identifier)"))
    
    # fit models
    message('real models')
    print(table(subset_for_gene$dataset))
    print(table(subset_for_gene[[allele_col]]))
    null_model   <- glmer(null_formula,   data = subset_for_gene, family = binomial, control = control_bobyqa) 
    test_model   <- glmer(allele_formula, data = subset_for_gene, family = binomial, control = control_bobyqa)

    message(warnings())
    # model comparison (likelihood‐ratio chi‐square test)
    anova_table         <- anova(null_model, test_model)
    observed_chisq      <- anova_table[2, "Chisq"]
    observed_p_value    <- anova_table[2, "Pr(>Chisq)"]
    
    # collect results
    results_list[[i]] <- tibble(
      gene              = gene_name,
      chain             = paste0("TR", chain_letter),
      gene_type         = gene_type,
      allele_column     = allele_col,
      n_observations    = n_obs,
      n_unique_alleles  = n_alleles,
      observed_chisq    = observed_chisq,
      observed_p_value  = observed_p_value
    )
  }
  
  # bind all rows into one data frame
  bind_rows(results_list)
}


valid_genes_by_chain = readRDS('valid_genes_by_chain.rds')
comb_data = readRDS('datasets6.rds')
      
results = perform_genewise_mixed_model(vmodeling_df, 'tmem_outcome', valid_genes_by_chain, 'V', 0)
warnings()
saveRDS(results, 'tmem_vgene_v5_alldatasets_results.rds')  
      
      
results = perform_genewise_mixed_model(vmodeling_df, 'cd8_outcome', valid_genes_by_chain, 'V', 0)
warnings()
saveRDS(results, 'cd8_vgene_v5_alldatasets_results.rds')  
