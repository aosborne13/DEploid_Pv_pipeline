#!/usr/bin/env bash

# --print-header generates the following example:
# # [1]CHROM	[2]POS	[3]sample1:GT	[4]SAMPLE2:GT	[5]SampleN:GT
#
# this sed script removes all nonessentials:
# sed -e '1s/# //' -e '1s/\[[[:digit:]]*\]//g' -e "1s/:$tag//g"
#
# final result:
# CHROM	POS	sample1	SAMPLE2	SampleN

# vcf2tab <TAG> <VCF/BCF> <TXT.GZ>
vcf2tab() {
	tag="${1%\{*}"
	format="%CHROM\t%POS[\t%$1]\n"

	bcftools query --format "$format" --print-header "$2" |
	sed -e '1s/# //' -e '1s/\[[[:digit:]]*\]//g' -e "1s/:$tag//g" |
	gzip > "$3"
}


# helper <TAG> <infile> <outdir>
helper() {
	tag="$(sed -e 's/{/./' -e 's/}//' <<< "$1")"
	outfile="${infile%*.vcf.gz}.$tag.txt.gz"
	outfile="${outfile##*/}"
	outfile="$3/$outfile"
	printf '%s\n' "Subsetting "$infile" into "$outfile""

	vcf2tab "$1" "$infile" "$outfile"
}


infile='merged.bi.filt.GT.miss0.4.snps.vcf.gz'
outdir='DEploid_input'

helper 'GT' "$infile" "$outdir" &
helper 'AD{0}' "$infile" "$outdir" &
helper 'AD{1}' "$infile" "$outdir" &

wait
