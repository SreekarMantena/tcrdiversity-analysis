
# polarization.py
# Script used to polarize alleles for analyses.

import subprocess, sys
from pathlib import Path
ENSEMBL_DIR = Path('./homo_sapiens_ancestor_GRCh38/').resolve()       # extracted FASTA folder
CHAINFILE   = Path('hg38ToGCA_009914755.4.over.chain.gz')
OUT_DIR     = Path('/liftover/')
OUT_DIR.mkdir(exist_ok=True)

# temp / final filenames
BED_HG38 = OUT_DIR / "chr14_window.hg38.bed"
BED_T2T  = OUT_DIR / "chr14_window.t2t.bed"
UNMAPPED = OUT_DIR / "chr14_window.unmapped"
ISAFE_TXT = OUT_DIR / "chr14_window.AA.t2t.txt"

# ------------------------------------------------------------------
# helper: stream FASTA and emit only bases inside the window
# ------------------------------------------------------------------
def fasta_window_to_bed4(fa_path: Path, bed_path: Path):
    """Write BED4 for A/C/G/T bases whose 1-based position is in [start, end]."""
    with fa_path.open("rb") as fa, bed_path.open("w") as bed:
        header = fa.readline().decode().strip()[1:]
        # header looks like ANCESTOR_for_chromosome:GRCh38:14:1:...
        chrom = "chr14"                                     # forced
        pos1  = 1                                           # 1-based
        for raw in fa:
            for base in raw.decode().rstrip():
                if base.upper() in "ACGT":
                    start0 = pos1 - 1                       # 0-based start
                    bed.write(f"{chrom}\t{start0}\t{start0+1}\t{base.upper()}\n")
                pos1 += 1
                if pos1 > HG38_END:                         # stop early
                    return

# ------------------------------------------------------------------
# 1. Create BED4 on hg38 coordinates
# ------------------------------------------------------------------
fasta_path = ENSEMBL_DIR / "homo_sapiens_ancestor_14.fa"
print("▶ Building BED4 for requested window …")
fasta_window_to_bed4(fasta_path, BED_HG38)

# ------------------------------------------------------------------
# 2. liftOver to CHM13 v2.0
# ------------------------------------------------------------------
print("▶ liftOver …")
subprocess.run([
    "liftOver", "-bedPlus=4", "-tab",
    BED_HG38, CHAINFILE, BED_T2T, UNMAPPED
], check=True)
BED_HG38.unlink()

# ------------------------------------------------------------------
# 3. Convert BED4 → iSAFE table (1-based POS) & compress / index
# ------------------------------------------------------------------
with BED_T2T.open() as b_in, ISAFE_TXT.open("w") as iso:
    for row in b_in:
        chrom, start0, _end, base = row.rstrip("\n").split("\t")
        iso.write(f"chr14\t{int(start0)+1}\t{base}\n")    # 1-based POS
BED_T2T.unlink()

print("▶ bgzip + tabix …")
subprocess.run(["bgzip", "-f", ISAFE_TXT], check=True)
subprocess.run(["tabix", "-s1", "-b2", "-e2", f"{ISAFE_TXT}.gz"], check=True)

print(f"✓ DONE → {ISAFE_TXT.name}.gz  (plus .tbi) in {OUT_DIR}")
