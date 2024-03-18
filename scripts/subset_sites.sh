#!/usr/bin/env bash

sites='DEploid_input/merged.bi.filt.GT.miss0.4.snps.PLAF.sites.txt'

bcftools view \
	--output /DEploid_input/merged_dataset.filt.DEploid.bcf.gz \
	--output-type b \
	--threads 8 \
	--regions-file "$sites" \
	merged.bi.filt.GT.miss0.4.snps.vcf.gz
