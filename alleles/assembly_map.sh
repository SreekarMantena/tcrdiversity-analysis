#-------------------------------------------------------------------------------
# CONFIGURATION
#-------------------------------------------------------------------------------
ASSEMBLY_DIR="unzipped_assemblies"
MAPPING_TABLE="all_loci_for_mapping.tsv"
OUT_DIR="assemblies_mapped_locus"
TEMP_FASTA_DIR="$OUT_DIR/temp_fastas"

mkdir -p "$OUT_DIR" "$TEMP_FASTA_DIR"
mkdir -p "assembly_mapping"
  
module load samtools 

#-------------------------------------------------------------------------------
# READ MAPPING TABLE (skip header)
#-------------------------------------------------------------------------------
mapfile -t _rows < <( tail -n +2 "$MAPPING_TABLE" )
TOTAL_ROWS=${#_rows[@]}

JOBS=50
TASK_ID=${SLURM_ARRAY_TASK_ID}

# Determine chunk size & slice
CHUNK_SIZE=$(( (TOTAL_ROWS + JOBS - 1) / JOBS ))
START_INDEX=$(( (TASK_ID - 1) * CHUNK_SIZE ))
if (( START_INDEX >= TOTAL_ROWS )); then
  echo "Task $TASK_ID: no rows to process."
  exit 0
fi
selected_rows=( "${_rows[@]:START_INDEX:CHUNK_SIZE}" )
echo "Task $TASK_ID: processing ${#selected_rows[@]} rows (indices $START_INDEX–$((START_INDEX+CHUNK_SIZE-1)))"

#-------------------------------------------------------------------------------
# MAIN LOOP: for each selected line
#-------------------------------------------------------------------------------
for row in "${selected_rows[@]}"; do
  IFS=$'\t' read -r assembly_id locus_start locus_end contig ucsc genbank sample_id locus reference_file <<< "$row"

  echo ">> Sample: $sample_id | Contig: $contig [$locus_start–$locus_end]"

  asm_fa="$ASSEMBLY_DIR/${sample_id}.fa"
  temp_slice="$TEMP_FASTA_DIR/${sample_id}_${contig}_${locus_start}_${locus_end}.fa"

  # 1) extract the slice into OUT_DIR/temp_fastas
  samtools faidx "$asm_fa" "${contig}:${locus_start}-${locus_end}" > "$temp_slice"

  # 2) ensure region FASTA is indexed
  [[ -f "${reference_file}.fai" ]] || samtools faidx "$reference_file"

  # 3) map → sort → BAM → index (no extra filtering)
  out_bam="$OUT_DIR/${sample_id}_${contig}_${locus}.bam"
  minimap2 -t $SLURM_CPUS_PER_TASK -ax asm10 "$reference_file" "$temp_slice" \
    | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$out_bam"
  samtools index "$out_bam"

  echo "   → Wrote $(basename "$out_bam")"
done

echo "All done. BAMs in $OUT_DIR; sliced FASTAs in $TEMP_FASTA_DIR"