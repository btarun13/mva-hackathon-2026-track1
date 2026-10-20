#!/usr/bin/env bash
# End-to-end Track 1 run: patient VCF -> submission CSV.
# Needs: bcftools/htslib/samtools >= 1.20, Docker (VEP 113 image), uv, and ref/ from prepare_reference.sh.
# Usage: scripts/run_pipeline.sh [HF_USERNAME]
set -euo pipefail
cd "$(dirname "$0")/.."
USER_TAG=${1:-Tarun2712}
VCF=data/WGS_EX2312012_HGWCNDSX7.vcf.gz
PY="uv run --no-project --with pandas python"
mkdir -p work submissions

# 0. Patient VCF uses bare contigs (15); the scorer and ref/ use chr15. Extract panel regions + rename.
[ -f $VCF.tbi ] && [ $VCF.tbi -nt $VCF ] || tabix -f -p vcf $VCF
sed 's/^chr//' ref/panel_regions.bed > work/panel_regions.nochr.bed
bcftools view -R work/panel_regions.nochr.bed $VCF -Ou \
  | bcftools annotate --rename-chrs ref/chrmap.txt -Oz -o work/proband.panel.chr.vcf.gz
tabix -f work/proband.panel.chr.vcf.gz

# 1-2. Gene list -> BED; normalise + subset + genotype table
$PY mva_track1.py genes  --genes genes.txt --vcf work/proband.panel.chr.vcf.gz --out work/genes.bed
$PY mva_track1.py subset --vcf work/proband.panel.chr.vcf.gz --bed work/genes.bed \
  --ref ref/GRCh38.panel.fa.gz --out work/subset.vcf.gz

# 3. VEP 113 (Docker), cache-free: GTF + FASTA + ClinVar + gnomAD v4.1 + AlphaMissense
$PY mva_track1.py annotate --vep bin/vep --vcf work/subset.vcf.gz --out work/vep.tsv \
  --gtf ref/GRCh38.113.panel.gtf.gz --fasta ref/GRCh38.panel.fa.gz --clinvar ref/clinvar.panel.vcf.gz \
  --gnomade ref/gnomad.exomes.panel.vcf.gz --gnomadg ref/gnomad.genomes.panel.vcf.gz \
  --am ref/AlphaMissense.panel.tsv.gz

# 4-6. Filter, score, compound-het grouping, EPCR, write + validate submission
OUT=submissions/${USER_TAG}_vep-rules-v1.csv
$PY mva_track1.py rank --vep work/vep.tsv --gt work/subset.gt.tsv --genes genes.txt --out $OUT
$PY mva_track1.py validate --csv $OUT
