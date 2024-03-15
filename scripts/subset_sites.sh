#!/usr/bin/env bash

sites='/DEploid_input/merged_dataset.filt.snps.PLAF.sites.txt'

bcftools view \
	--output /DEploid_input/merged_dataset.filt.DEploid.bcf.gz \
	--output-type b \
	--threads 8 \
	--regions-file "$sites" \
	/DEploid_input/merged_dataset.filt.snps.vcf.gz
