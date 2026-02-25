#!/usr/bin/env bash

#-------------------------------------------------------------------------------
# CONFIGURATION
#-------------------------------------------------------------------------------
HIFI_SAMPLE_INDEX="data_hifi_pre_release.index.csv"
HIFI_RAW_DIRECTORY="hifi_raw"
ALLELE_FASTA_DIR="alleles_read_validation"
PAF_OUTPUT_DIRECTORY="readsmapped_to_alleles_"

THREAD_COUNT=8

MINIMAP2_OPTIONS="-t${THREAD_COUNT}"

#-------------------------------------------------------------------------------
# PREPARE OUTPUT DIRECTORY & MODULES
#-------------------------------------------------------------------------------
echo "Ensuring PAF output directory exists: ${PAF_OUTPUT_DIRECTORY}"
mkdir -p "${PAF_OUTPUT_DIRECTORY}"
module load samtools

#-------------------------------------------------------------------------------
# DISCOVER ALLELE FASTA FILES & CHUNK FOR THIS TASK
#-------------------------------------------------------------------------------
mapfile -t allele_fasta_paths < <(ls "${ALLELE_FASTA_DIR}"/*.fasta | sort)
total_alleles=${#allele_fasta_paths[@]}

# number of array tasks
ARRAY_JOBS=60

# compute how many FASTAs per task (ceiling division)
chunk_size=$(( (total_alleles + ARRAY_JOBS - 1) / ARRAY_JOBS ))

task_id=${SLURM_ARRAY_TASK_ID}
start_idx=$(( (task_id - 1) * chunk_size ))
end_idx=$(( task_id * chunk_size - 1 ))

# if this task has no FASTAs assigned, exit
if (( start_idx >= total_alleles )); then
  echo "[$(date)] Task ${task_id}: no FASTA files to process, exiting."
  exit 0
fi

# clamp end_idx to last element
if (( end_idx >= total_alleles )); then
  end_idx=$(( total_alleles - 1 ))
fi

echo "[$(date)] Task ${task_id}: processing FASTAs ${start_idx}–${end_idx} of ${total_alleles}"

#-------------------------------------------------------------------------------
# LOOP OVER ASSIGNED FASTA FILES
#-------------------------------------------------------------------------------
for idx in $(seq "${start_idx}" "${end_idx}"); do
  allele_fasta="${allele_fasta_paths[$idx]}"
  sample_id=$(basename "${allele_fasta}" .fasta)
  echo "[$(date)]  → mapping sample ${sample_id} (FASTA index ${idx})"

  # gather input BAMs for this sample
  mapfile -t bam_names < <(
    awk -F',' -v sid="${sample_id}" '$1==sid{print $3}' "${HIFI_SAMPLE_INDEX}"
  )

  raw_bam_paths=()
  for bam_name in "${bam_names[@]}"; do
    bam_full=$(find "${HIFI_RAW_DIRECTORY}" -type f -name "${bam_name}" | head -n1)
    if [[ -z "${bam_full}" ]]; then
      echo "Warning: raw BAM ${bam_name} not found" >&2
    else
      raw_bam_paths+=("${bam_full}")
    fi
  done

  if (( ${#raw_bam_paths[@]} == 0 )); then
    echo "No raw BAMs for ${sample_id}; skipping."
    continue
  fi

  # run merge → fastq → minimap2 → PAF
  paf_out="${PAF_OUTPUT_DIRECTORY}/${sample_id}.paf"
  echo "    → Writing PAF to ${paf_out}"
  samtools merge -n -u - "${raw_bam_paths[@]}" \
    | samtools fastq - \
    | minimap2 \
        ${MINIMAP2_OPTIONS} "${allele_fasta}" - \
    > "${paf_out}"

  echo "[$(date)]  ← done ${sample_id}"
done

echo "[$(date)] Task ${task_id} complete."