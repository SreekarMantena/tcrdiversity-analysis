
"""
Code to perform scans of selection

Data from the AoU CDRv8 cohort was used for these analyses. The AoU consortium has restrictions regarding the release of data and code from their portal.
The steps to perform variant QC and generate the callset used in this paper are detailed in the methods section of the manuscript.

"""

# Fst scan
#!/usr/bin/env bash
set -euo pipefail

#── INPUT & OUTPUT PATHS ───────────────────────────────────────────────────────
INPUT_DIR="merged_vcfs_filtered_and_phased_by_ancestry"
OUTDIR_BASE="fst_vcftools"

#── MAKE BASE OUTPUT DIR ───────────────────────────────────────────────────────
mkdir -p "${OUTDIR_BASE}"

echo "[$(timestamp)] Starting FST computations for loci: ${LOCI[*]} and ancestries (ordered): ${ANCESTRIES[*]}"

# Loop over loci
for locus in "${LOCI[@]}"; do
  locus_outdir="${OUTDIR_BASE}/${locus}"
  mkdir -p "${locus_outdir}"

  # Loop over ordered ancestry pairs (p1 != p2)
  for pop1 in "${ANCESTRIES[@]}"; do
    for pop2 in "${ANCESTRIES[@]}"; do
      if [[ "${pop1}" == "${pop2}" ]]; then
        continue
      fi

      # Input VCFs (expecting files like TRA_phased_eas.vcf.gz)
      bgz_pop1="${INPUT_DIR}/${locus}_phased_${pop1}.vcf.gz"
      bgz_pop2="${INPUT_DIR}/${locus}_phased_${pop2}.vcf.gz"

      # Check inputs exist
      if [[ ! -f "${bgz_pop1}" ]]; then
        echo "[$(timestamp)] [SKIP] Missing file: ${bgz_pop1}"
        continue
      fi
      if [[ ! -f "${bgz_pop2}" ]]; then
        echo "[$(timestamp)] [SKIP] Missing file: ${bgz_pop2}"
        continue
      fi

      pair_prefix="${locus,,}_${pop1}_vs_${pop2}"
      pair_outdir="${locus_outdir}/${pair_prefix}"
      mkdir -p "${pair_outdir}"

      merged_vcf="${pair_outdir}/${pair_prefix}_merged.vcf.gz"
      pop1_samples="${pair_outdir}/${pop1}_sample_list.txt"
      pop2_samples="${pair_outdir}/${pop2}_sample_list.txt"
      fst_prefix="${pair_outdir}/${pair_prefix}"

      #── Merge the two bgzipped VCFs ──────────────────────────────────────────
      echo "[$(timestamp)] Merging ${locus}: ${pop1} + ${pop2} -> ${merged_vcf}"
      bcftools merge \
        "${bgz_pop1}" \
        "${bgz_pop2}" \
        -Oz -o "${merged_vcf}"

      tabix -f -p vcf "${merged_vcf}"

      #── Extract sample lists ─────────────────────────────────────────────────
      echo "[$(timestamp)] Extracting sample names for ${pair_prefix}..."
      bcftools query -l "${bgz_pop1}" > "${pop1_samples}"
      bcftools query -l "${bgz_pop2}" > "${pop2_samples}"

      #── Compute per-site Weir & Cockerham FST ───────────────────────────────
      echo "[$(timestamp)] Computing FST for ${pair_prefix}..."
      vcftools \
        --gzvcf "${merged_vcf}" \
        --weir-fst-pop "${pop1_samples}" \
        --weir-fst-pop "${pop2_samples}" \
        --out      "${fst_prefix}"

      echo "[$(timestamp)] Done: ${pair_prefix}"
      ls -lh "${pair_outdir}" || true
      echo "------------------------------------------------------------------"
    done
  done
done

# isafe run
isafe --input eas_qced.vcf.gz --output isafe_results  --MaxRank 300  --IgnoreGaps

# selscan run
selscan --xpehh --pmap --threads 60 --vcf eas_qced.vcf.gz --vcf-ref afr_qced.vcf.gz --out "eas_vs_afr_xpehh"
selscan --xpehh --pmap --threads 60 --vcf eas_qced.vcf.gz --vcf-ref eur_qced.vcf.gz --out "eas_vs_eur_xpehh"
