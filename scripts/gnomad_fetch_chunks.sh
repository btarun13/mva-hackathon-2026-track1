#!/usr/bin/env bash
# Robust remote gnomAD v4.1 lookup: regions are fetched in chunks of 100 with retries, so a dropped HTTP stream
# only costs one chunk. Usage: gnomad_fetch_chunks.sh <gnomad_dir> "<set chrN> <set chrM> ..."
set -uo pipefail
D=$1; shift; LIST=$1; G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
cd "$D"; mkdir -p chunks
fetch_chunk() {  # $1=set $2=chr $3=chunk bed
  local set=$1 c=$2 bed=$3 out=${3%.bed}.tsv
  [ -s "$out" ] && return 0
  for try in 1 2 3 4 5 6; do
    bcftools view -G -R "$bed" $G/$set/gnomad.$set.v4.1.sites.$c.vcf.bgz -Ou 2>>chunks/err.log \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_grpmax\t%INFO/grpmax\n' > "$out.part"
    if [ "${PIPESTATUS[0]}" -eq 0 ]; then mv "$out.part" "$out"; return 0; fi
    sleep $((try*5))
  done
  echo "FAILED $set $c $bed" >> chunks/failed.log; return 1
}
export -f fetch_chunk; export G
jobs=()
set -- $LIST
while [ $# -ge 2 ]; do
  set_=$1; c=$2; shift 2
  split -l 100 -a 3 reg.$set_.$c.bed chunks/$set_.$c.
  for f in chunks/$set_.$c.[a-z][a-z][a-z]; do mv "$f" "$f.bed"; echo "$set_ $c $f.bed"; done
done > chunks/jobs.txt
echo "chunks to fetch: $(wc -l < chunks/jobs.txt)" >&2
caffeinate -i xargs -P 8 -L 1 bash -c 'fetch_chunk "$0" "$1" "$2"' < chunks/jobs.txt
# assemble per set/chr
set -- $LIST
while [ $# -ge 2 ]; do
  set_=$1; c=$2; shift 2
  n=$(ls chunks/$set_.$c.*.bed | wc -l); d=$(ls chunks/$set_.$c.*.tsv 2>/dev/null | wc -l)
  cat chunks/$set_.$c.*.tsv > $set_.$c.tsv.new 2>/dev/null && mv $set_.$c.tsv.new $set_.$c.tsv
  echo "$set_ $c: $d/$n chunks ok, rows=$(wc -l < $set_.$c.tsv)" >&2
done
echo "failed chunks: $(wc -l < chunks/failed.log 2>/dev/null || echo 0)" >&2
