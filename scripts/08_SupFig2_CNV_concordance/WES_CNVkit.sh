#!/bin/bash
############################################################
# WES_CNVkit.sh
# ----------------------------------------------------------
# WES copy-number calling with CNVkit for Supplementary Fig. 2.
# Matched tumor / adjacent-normal pairs (Agilent SureSelect V5, hg38):
#   T35/N35, T36/N36, T39/N39, T44/N44
# One pooled reference is built from the four normals.
# Run on a Linux SLURM server (8 CPUs, 16 GB); CNVkit installed via
#   conda create -n cnvkit -c bioconda -c conda-forge cnvkit samtools
#
# Inputs:
#   $WES_BAM_DIR/Lung_{ID}/Lung_{ID}.recal.bam (+ .bai)
#   $REF_FASTA   : chr-prefixed GRCh38 FASTA matching the BAM contigs
#   $TARGET_BED  : SureSelect V5 covered BED (S04380110_Covered.bed)
#
# Outputs ($OUT_DIR):
#   Lung_{T}.cnr (bin-level log2; input to SupFig2_CNV_concordance.R),
#   Lung_{T}.cns, ref_pooledNormals.cnn, WES_CNV_heatmap_4samples.pdf
############################################################
set -euo pipefail

WES_BAM_DIR=${WES_BAM_DIR:-wes_bams}
REF_FASTA=${REF_FASTA:-genome.fa}
TARGET_BED=${TARGET_BED:-S04380110_Covered.bed}
OUT_DIR=${OUT_DIR:-cnvkit_out}
CNVKIT=${CNVKIT:-cnvkit.py}
SAMTOOLS=${SAMTOOLS:-samtools}
THREADS=${SLURM_CPUS_PER_TASK:-8}

TUMOR_IDS="T35 T36 T39 T44"
NORMAL_IDS="N35 N36 N39 N44"

mkdir -p "$OUT_DIR"

# 1) Clean target BED (drop browser/track lines; keep chr, start, end)
TARGET_ALL="$OUT_DIR/SureSelect_V5_targets.allchr.bed"
grep -E '^chr' "$TARGET_BED" | cut -f1-3 > "$TARGET_ALL"

# 2) Preflight: contig naming must match across BAM, FASTA and BED
[ -f "$REF_FASTA" ] || { echo "[ERROR] REF_FASTA not found: $REF_FASTA"; exit 1; }
[ -f "${REF_FASTA}.fai" ] || $SAMTOOLS faidx "$REF_FASTA"
ONE_BAM="$WES_BAM_DIR/Lung_T35/Lung_T35.recal.bam"
[ -f "$ONE_BAM" ] || { echo "[ERROR] BAM not found: $ONE_BAM"; exit 1; }
bam_chr=$($SAMTOOLS view -H "$ONE_BAM" | awk '/^@SQ/{sub("SN:","",$2); print $2}' | grep -E '^(chr)?1$' | head -1)
ref_chr=$(cut -f1 "${REF_FASTA}.fai" | grep -E '^(chr)?1$' | head -1)
bed_chr=$(head -1 "$TARGET_ALL" | cut -f1)
echo "[check] chr1 naming: BAM='$bam_chr' REF='$ref_chr' BED='$bed_chr'"
if [ "$bam_chr" != "$ref_chr" ] || [ "$bam_chr" != "$bed_chr" ]; then
  echo "[ERROR] contig naming mismatch among BAM/REF/BED"; exit 1
fi

# 3) Restrict to primary chromosomes (chr1-22, X, Y) present in the BAM
$SAMTOOLS view -H "$ONE_BAM" | awk '/^@SQ/{sub("SN:","",$2); print $2}' \
  | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)$' | sort -V > "$OUT_DIR/main_contigs.txt"
TARGET="$OUT_DIR/SureSelect_V5_targets.bed"
awk 'NR==FNR{a[$1];next} ($1 in a)' "$OUT_DIR/main_contigs.txt" "$TARGET_ALL" > "$TARGET"

# 4) Accessible regions on the same contigs
$CNVKIT access "$REF_FASTA" -o "$OUT_DIR/access.full.bed"
ACCESS="$OUT_DIR/access.hg38.bed"
awk 'NR==FNR{a[$1];next} ($1 in a)' "$OUT_DIR/main_contigs.txt" "$OUT_DIR/access.full.bed" > "$ACCESS"

# 5) BAM lists
TUMORS="";  for s in $TUMOR_IDS;  do TUMORS="$TUMORS $WES_BAM_DIR/Lung_$s/Lung_$s.recal.bam";  done
NORMALS=""; for s in $NORMAL_IDS; do NORMALS="$NORMALS $WES_BAM_DIR/Lung_$s/Lung_$s.recal.bam"; done
for f in $TUMORS $NORMALS; do
  [ -f "$f" ] || { echo "[ERROR] missing BAM: $f"; exit 1; }
done

# 6) CNVkit batch (hybrid capture; pooled-normal reference)
$CNVKIT batch $TUMORS \
  --normal $NORMALS \
  --targets "$TARGET" \
  --fasta "$REF_FASTA" \
  --access "$ACCESS" \
  --output-reference "$OUT_DIR/ref_pooledNormals.cnn" \
  --output-dir "$OUT_DIR" \
  --scatter --diagram \
  -p "$THREADS"

# 7) Genome-wide segment heatmap (QC); CNVkit drops ".recal" from output names
TUMOR_CNS=""; for s in $TUMOR_IDS; do TUMOR_CNS="$TUMOR_CNS $OUT_DIR/Lung_$s.cns"; done
$CNVKIT heatmap $TUMOR_CNS -d -o "$OUT_DIR/WES_CNV_heatmap_4samples.pdf"

echo "[done] .cnr/.cns per tumor in $OUT_DIR"
