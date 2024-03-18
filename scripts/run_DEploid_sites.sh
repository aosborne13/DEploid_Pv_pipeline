#!/usr/bin/env bash

plaf='DEploid_input/merged.bi.filt.GT.miss0.4.snps.PLAF.txt'
samples='samples.txt'

while read -r sample; do
    ~/tools/DEploid/dEploid \
        -ref "DEploid_input/$sample.AD.0.txt.gz" \
        -alt "DEploid_input/$sample.AD.1.txt.gz" \
        -plaf "$plaf" \
        -exclude "run_DEploid/${sample}PotentialOutliers.txt" \
        -o "run_DEploid/$sample" \
        -k '4' \
        -noPanel
done < "$samples"
