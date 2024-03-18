#!/usr/bin/env bash

plaf='DEploid_input/merged_dataset.filt.snps.PLAF.txt'
samples='samples.txt'

while read -r sample; do
    /home/ashley/tools/DEploid/utilities/dataExplore.r \
        -ref "DEploid_input/$sample.AD.0.txt.gz" \
        -alt "DEploid_input/$sample.AD.1.txt.gz" \
        -plaf "$plaf" \
        -o "run_DEploid/$sample"
done < "$samples"
