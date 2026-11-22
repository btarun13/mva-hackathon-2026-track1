#!/usr/bin/env bash
# Tier A (cache-free): coding-region annotation of the proband VCF with VEP 113 in GTF mode, then gnomAD v4.1 frequencies
# fetched by remote tabix slicing for the HIGH/MODERATE variants only (no 25 GB cache needed).
# Inputs : work/genome/coding.raw.vcf.gz (patient VCF restricted to protein-coding exons +/-25 bp, chr-prefixed)
#          /Volumes/SAN_DRIVE/mva_refs: GRCh38.primary.chr.fa.gz, GRCh38.113.chr.gtf.gz, clinvar.chr.vcf.gz, AlphaMissense_hg38.bgz.tsv.gz
# Outputs: work/genome/coding.vcf.gz, coding.gt.tsv, vep_coding.raw.tsv, gnomad_lookup.tsv, vep_coding.tsv  (patient-derived: never commit)
set -euo pipefail
cd "$(dirname "$0")/.."
SAN=/Volumes/SAN_DRIVE; REF=$SAN/mva_refs; W=work/genome
FORK=${FORK:-4}; PY="uv run --no-project --with pandas python"

# 1. normalise (split multiallelics, left-align), set ID, genotype table
if [ ! -s $W/coding.gt.tsv ]; then
  bcftools norm -m -both -f $REF/GRCh38.primary.chr.fa.gz $W/coding.raw.vcf.gz -Ou \
    | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o $W/coding.vcf.gz
  tabix -f $W/coding.vcf.gz
  bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%GQ\t%AD]\n' $W/coding.vcf.gz > $W/coding.gt.tsv
fi
echo "normalised: $(wc -l < $W/coding.gt.tsv) variants" >&2

# 2. VEP 113 cache-free: transcripts from GTF, ClinVar + AlphaMissense as custom/plugin
if [ ! -s $W/vep_coding.raw.tsv ]; then
  docker run --rm --platform linux/amd64 -v "$PWD":/d -v $SAN:/san -w /d ensemblorg/ensembl-vep:release_113.0 vep \
    -i $W/coding.vcf.gz -o $W/vep_coding.raw.tsv --tab --force_overwrite --no_stats --assembly GRCh38 \
    --gtf /san/mva_refs/GRCh38.113.chr.gtf.gz --fasta /san/mva_refs/GRCh38.primary.chr.fa.gz \
    --pick --symbol --canonical --hgvs --numbers --biotype \
    --custom file=/san/mva_refs/clinvar.chr.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT \
    --dir_plugins /plugins --plugin AlphaMissense,file=/san/mva_refs/AlphaMissense_hg38.bgz.tsv.gz \
    --fork $FORK --buffer_size 5000 \
    --fields "Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,BIOTYPE,CANONICAL,Consequence,IMPACT,EXON,INTRON,Protein_position,Amino_acids,HGVSc,HGVSp,ClinVar,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,am_class,am_pathogenicity"
fi
echo "VEP done: $(grep -vc '^#' $W/vep_coding.raw.tsv) annotated variants" >&2

# 3. gnomAD v4.1 exomes+genomes for HIGH/MODERATE variants only, by remote tabix slicing (Google storage)
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
awk -F'\t' '!/^#/ && ($10=="HIGH"||$10=="MODERATE"){split($1,a,"_"); print a[1]"\t"a[2]-1"\t"a[2]}' $W/vep_coding.raw.tsv | sort -k1,1 -k2,2n -u > $W/hm_regions.bed
echo "gnomAD lookup regions: $(wc -l < $W/hm_regions.bed)" >&2
mkdir -p $W/gnomad
for set in exomes genomes; do
  for c in $(seq 1 22) X Y; do
    awk -v c=chr$c '$1==c' $W/hm_regions.bed > $W/gnomad/reg.$set.chr$c.bed
    [ -s $W/gnomad/reg.$set.chr$c.bed ] || { rm -f $W/gnomad/reg.$set.chr$c.bed; continue; }
    printf '%s\n' "bcftools view -G -R $W/gnomad/reg.$set.chr$c.bed $G/$set/gnomad.$set.v4.1.sites.chr$c.vcf.bgz -Ou 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/AF_grpmax\\t%INFO/grpmax\\n' > $W/gnomad/$set.chr$c.tsv.part && mv $W/gnomad/$set.chr$c.tsv.part $W/gnomad/$set.chr$c.tsv" > $W/gnomad/job.$set.chr$c.sh
  done
done
ls $W/gnomad/job.*.sh | xargs -P 6 -n 1 bash
(printf "set\tCHROM\tPOS\tREF\tALT\tAF\tAF_grpmax\tgrpmax\n"; for set in exomes genomes; do cat $W/gnomad/$set.chr*.tsv 2>/dev/null | awk -v s=$set 'BEGIN{OFS="\t"}{print s,$0}'; done) > $W/gnomad_lookup.tsv
echo "gnomAD rows fetched: $(($(wc -l < $W/gnomad_lookup.tsv)-1))" >&2

# 4. append MAX_AF / MAX_AF_POPS to the VEP table
$PY scripts/add_gnomad_maxaf.py $W/vep_coding.raw.tsv $W/gnomad_lookup.tsv $W/vep_coding.tsv
echo "wrote $W/vep_coding.tsv" >&2
