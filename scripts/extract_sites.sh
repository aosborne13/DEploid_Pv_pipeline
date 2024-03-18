#!/usr/bin/env bash

outdir='DEploid_input'
bcf='DEploid_input/merged_dataset.filt.DEploid.bcf.gz'
samples='DEploid_input/samples.txt'

while read -r sample; do
    # First extraction
    tag='AD{0}'
    outfile="$sample."$(sed -e 's/{/./' -e 's/}//' <<< "$tag")".txt.gz"
    outfile="$outdir/$outfile"
    format="%CHROM\t%POS[\t%$tag]\n"
    bcftools view --output-type u --samples "$sample" "$bcf" |
    bcftools query --format "$format" --print-header |
    sed -e '1s/# //' -e '1s/\[[[:digit:]]*\]//g' -e "1s/:"${tag%\{*}"//g" |
    gzip > "$outfile" &

    # Second extraction
    tag='AD{1}'
    outfile="$sample."$(sed -e 's/{/./' -e 's/}//' <<< "$tag")".txt.gz"
    outfile="$outdir/$outfile"
    format="%CHROM\t%POS[\t%$tag]\n"
    bcftools view --output-type u --samples "$sample" "$bcf" |
    bcftools query --format "$format" --print-header |
    sed -e '1s/# //' -e '1s/\[[[:digit:]]*\]//g' -e "1s/:"${tag%\{*}"//g" |
    gzip > "$outfile" &
done < "$samples"

# Wait for all background jobs to finish
wait
