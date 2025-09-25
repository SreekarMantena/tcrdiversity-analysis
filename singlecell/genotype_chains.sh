#!/usr/bin/env bash
set -euo pipefail

# Source and destination directories
readonly SRC_DIR="allele_discovery/allele_files/"
readonly DEST_DIR="singlecell_models/igblast/ncbi-igblast-1.22.0//"

# Ensure destination exists
mkdir -p "${DEST_DIR}"

# Array of gene types and their file patterns
declare -A gene_files=(
  [V]="V_gene_unique_region_forDB_shortheader.fasta"
  [D]="D_gene_unique_region_forDB_shortheader.fasta"
  [J]="J_gene_unique_region_forDB_shortheader.fasta"
)

readonly V_DB="${DEST_DIR}/${gene_files[V]}"
readonly D_DB="${DEST_DIR}/${gene_files[D]}"
readonly J_DB="${DEST_DIR}/${gene_files[J]}"

threadsused=60

# caushi
AUX_DATA="optional_file/human_gl.aux"
OUTPUT_DIR=".//caushi_tumor/"
CLONO_OUT=".//caushi_tumor/caushi_igblast_out.out"

mkdir -p "${OUTPUT_DIR}" 

echo "Running igblastn with threads..."
bin/igblastn \
  -germline_db_V "${V_DB}" \
  -germline_db_D "${D_DB}" \
  -germline_db_J "${J_DB}" \
  -organism human \
  -query "${QUERY_FASTA}" \
  -auxiliary_data "${AUX_DATA}" \
  -ig_seqtype TCR \
  -clonotype_out "${CLONO_OUT}" \
  -out "${OUTPUT_DIR}/caushi_igblast_out.txt" \
  -outfmt 7 \
  -num_threads $threadsused

# COMBAT
AUX_DATA="optional_file/human_gl.aux"
OUTPUT_DIR=".//combat/"
CLONO_OUT=".//combat/combat_data_igblast_genotyping.out"

mkdir -p "${OUTPUT_DIR}"

echo "Running igblastn with threads... COMBAT"
bin/igblastn \
  -germline_db_V "${V_DB}" \
  -germline_db_D "${D_DB}" \
  -germline_db_J "${J_DB}" \
  -organism human \
  -query "${QUERY_FASTA}" \
  -auxiliary_data "${AUX_DATA}" \
  -ig_seqtype TCR \
  -clonotype_out "${CLONO_OUT}" \
  -out "${OUTPUT_DIR}/combat_data_igblast_genotyping.txt" \
  -outfmt 7 \
  -num_threads $threadsused

# echo "Done. Results in ${OUTPUT_DIR}/edahiro_igblast_results_v2_fullcohort.txt, clonotypes in ${CLONO_OUT}"

## TENX ANTIGEN

# Other parameters
AUX_DATA="optional_file/human_gl.aux"
OUTPUT_DIR=".//tenx_antigen_data/"
CLONO_OUT=".//tenx_antigen_data/tenx_antigen_data_igblast_genotyping.out"

mkdir -p "${OUTPUT_DIR}"

echo "Running igblastn with threads... TENX"
bin/igblastn \
  -germline_db_V "${V_DB}" \
  -germline_db_D "${D_DB}" \
  -germline_db_J "${J_DB}" \
  -organism human \
  -query "${QUERY_FASTA}" \
  -auxiliary_data "${AUX_DATA}" \
  -ig_seqtype TCR \
  -clonotype_out "${CLONO_OUT}" \
  -out "${OUTPUT_DIR}/igblast_results.txt" \
  -outfmt 7 \
  -num_threads $threadsused


# EDAHIRO

# Other parameters
AUX_DATA="optional_file/human_gl.aux"
OUTPUT_DIR=".//edahiro/"
CLONO_OUT=".//edahiro/edahiro_igblast_genotyping.out"

mkdir -p "${OUTPUT_DIR}"

echo "Running igblastn with threads... EDAHIRO"
bin/igblastn \
  -germline_db_V "${V_DB}" \
  -germline_db_D "${D_DB}" \
  -germline_db_J "${J_DB}" \
  -organism human \
  -query "${QUERY_FASTA}" \
  -auxiliary_data "${AUX_DATA}" \
  -ig_seqtype TCR \
  -clonotype_out "${CLONO_OUT}" \
  -out "${OUTPUT_DIR}/igblast_results.txt" \
  -outfmt 7 \
  -num_threads $threadsused



# TEREKHOVA

# Other parameters
AUX_DATA="optional_file/human_gl.aux"
OUTPUT_DIR=".//terekhova/"
CLONO_OUT=".//terekhova/terekhova_igblast_genotyping.out"

mkdir -p "${OUTPUT_DIR}"

echo "Running igblastn with threads... TEREKHOVA"
bin/igblastn \
  -germline_db_V "${V_DB}" \
  -germline_db_D "${D_DB}" \
  -germline_db_J "${J_DB}" \
  -organism human \
  -query "${QUERY_FASTA}" \
  -auxiliary_data "${AUX_DATA}" \
  -ig_seqtype TCR \
  -clonotype_out "${CLONO_OUT}" \
  -out "${OUTPUT_DIR}/igblast_results.txt" \
  -outfmt 7 \
  -num_threads $threadsused

