#!/usr/bin/env bash

plaf='DEploid_input/merged.bi.filt.GT.miss0.4.snps.PLAF.txt'
panel='BEST/filtered_non-swga_AF_monoclonals.DEploid.panel.txt.gz'
samples='run_DEploid/filtered_non-swga_AF.DEploid.polyclonals.txt'

while read -r sample; do
    ~/tools/DEploid/dEploid \
        -ref "DEploid_input/$sample.AD.0.txt.gz" \
        -alt "DEploid_input/$sample.AD.1.txt.gz" \
        -plaf "$plaf" \
        -panel "$panel" \
        -exclude "run_DEploid/${sample}PotentialOutliers.txt" \
        -o "BEST/$sample" \
        -k '4' \
        -best
done < "$samples"
