#!/usr/bin/env bash
# Build the small, panel-only annotation bundle in ref/ (~250 MB) so VEP can run without its ~25 GB cache.
# Sources: Ensembl 113 (FASTA, GTF), ClinVar GRCh38 VCF, gnomAD v4.1 exomes+genomes, AlphaMissense hg38.
# Remote VCFs are read by tabix region queries, so only the panel regions are downloaded.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p ref && cd ref

CHROMS="2 4 5 11 15"   # chromosomes carrying the genes in genes.txt
printf "2 chr2\n4 chr4\n5 chr5\n11 chr11\n15 chr15\n" > chrmap.txt

# Panel regions: gene bodies (Ensembl 113) +/- 5 kb, 0-based BED, chr-prefixed
cat > panel_regions.bed <<'EOF'
chr2	110630467	110683098	BUB1
chr4	103100348	103203463	CENPE
chr5	887848	924357	TRIP13
chr11	28015617	28113186	KIF18A
chr11	95784964	95842070	CEP57
chr11	118993050	119020811	CENATAC
chr15	40155983	40226137	BUB1B
EOF

# --- Ensembl 113 FASTA (panel chromosomes only, chr-renamed) and GTF
E=https://ftp.ensembl.org/pub/release-113
for c in $CHROMS; do
  curl -sSf $E/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$c.fa.gz | gzip -dc | sed "1s/^>.*/>chr$c/"
done | bgzip > GRCh38.panel.fa.gz
samtools faidx GRCh38.panel.fa.gz
curl -sSf $E/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz | gzip -dc \
  | awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && ($1=="2"||$1=="4"||$1=="5"||$1=="11"||$1=="15"){$1="chr"$1; print}' \
  | sort -t$'\t' -k1,1 -k4,4n -k5,5n | bgzip > GRCh38.113.panel.gtf.gz
tabix -p gff GRCh38.113.panel.gtf.gz

# --- ClinVar (upstream has no chr prefix)
REGS_NOCHR=$(awk '{printf "%s%s:%d-%d", (n++?",":""), $1,$2+1,$3}' panel_regions.bed | sed 's/chr//g')
bcftools view -r "$REGS_NOCHR" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -Ou \
  | bcftools annotate --rename-chrs chrmap.txt -Oz -o clinvar.panel.vcf.gz
tabix -f clinvar.panel.vcf.gz

# --- gnomAD v4.1 exomes + genomes (sites only, AF / AF_grpmax / nhomalt)
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
for set in exomes genomes; do
  for c in $CHROMS; do
    regs=$(awk -v c=chr$c '$1==c{printf "%s%s:%d-%d", (n++?",":""), $1,$2+1,$3}' panel_regions.bed)
    bcftools view -G -r "$regs" $G/$set/gnomad.$set.v4.1.sites.chr$c.vcf.bgz -Ou \
      | bcftools annotate -x '^INFO/AF,INFO/AF_grpmax,INFO/nhomalt,INFO/AC' -Ob -o tmp.$set.chr$c.bcf
  done
  bcftools concat $(for c in $CHROMS; do echo tmp.$set.chr$c.bcf; done) -Oz -o gnomad.$set.panel.vcf.gz
  tabix -f gnomad.$set.panel.vcf.gz && rm tmp.$set.*.bcf
done

# --- AlphaMissense (no remote index: download once, slice, delete)
curl -sSfO https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
(gzip -dc AlphaMissense_hg38.tsv.gz | head -4
 gzip -dc AlphaMissense_hg38.tsv.gz | awk -F'\t' 'NR==FNR{s[NR]=$1;b[NR]=$2;e[NR]=$3;n=NR;next}
   !/^#/{for(i=1;i<=n;i++) if($1==s[i] && $2>b[i] && $2<=e[i]){print;break}}' panel_regions.bed -) \
  | bgzip > AlphaMissense.panel.tsv.gz
tabix -f -s 1 -b 2 -e 2 -S 4 AlphaMissense.panel.tsv.gz
rm -f AlphaMissense_hg38.tsv.gz clinvar.vcf.gz.tbi gnomad.*.vcf.bgz.tbi
echo "reference bundle ready in ref/"
