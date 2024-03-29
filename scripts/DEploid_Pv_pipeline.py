#! /usr/bin/env python
import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm
from fastq2matrix import run_cmd
from collections import defaultdict
import gzip
import os
import shutil

def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")

def main(args):

    run_cmd("mkdir DEploid_input")
    run_cmd("bcftools view -v snps %(input)s -Oz -o merged_snps.vcf.gz" % vars(args))
    run_cmd("bcftools filter -i 'FMT/DP>4' -S . merged_snps.vcf.gz -Oz -o merged.filt.snps.vcf.gz")
    run_cmd("plink --const-fid --vcf merged.filt.snps.vcf.gz --mind 0.4 --recode vcf --allow-extra-chr --out merged_plink" % vars(args))
    run_cmd("grep -P \"^#CHROM\" merged_plink.vcf | awk '{ $1=\"\"; $2=\"\";$3=\"\"; $4=\"\";$5=\"\"; $6=\"\";$7=\"\"; $8=\"\";$9=\"\"; print}' | sed 's/ /\\n/g' | tail -n+10 | sed 's/^0_//'  > merged_plink_new")
    run_cmd("bcftools view -S merged_plink_new --threads 4 -O z -o  merged.miss0.4.filt.snps.vcf.gz merged.filt.snps.vcf.gz")
    #run_cmd("bcftools view merged.miss0.4.filt.snps.vcf.gz | setGT.py --fraction 0.8 | bcftools view -O z -c 1 -o merged.filt.GT.miss0.4.snps.vcf.gz")
    run_cmd("bcftools view -m2 -M2 merged.miss0.4.filt.snps.vcf.gz --threads 4 -O z -o merged.bi.filt.GT.miss0.4.snps.vcf.gz")
    run_cmd("tabix -f merged.bi.filt.GT.miss0.4.snps.vcf.gz")

    run_cmd("bash extract_DEploid.sh")
    run_cmd("Rscript ~/tools/DEploid_Pv_pipeline/scripts/calculate_PLAF_DEploid.R DEploid_input/merged.bi.filt.GT.miss0.4.snps.GT.txt.gz DEploid_input/merged.bi.filt.GT.miss0.4.snps.AD.0.txt.gz DEploid_input/merged.bi.filt.GT.miss0.4.snps.AD.1.txt.gz")
    run_cmd("bash subset_sites.sh")
    run_cmd("bash extract_sites.sh")
    run_cmd("mkdir run_DEploid")
    run_cmd("bash dataExplore_sites.sh")
    run_cmd("bash run_DEploid_sites.sh")
    run_cmd("bash compile_props_DEploid_sites.sh")
    run_cmd("Rscript ~/tools/DEploid_Pv_pipeline/scripts/compare_clonality_classification.R")
    run_cmd("mkdir BEST")
    run_cmd("Rscript ~/tools/DEploid_Pv_pipeline/scripts/resolve_mixed_sites.R")
    run_cmd("Rscript ~/tools/DEploid_Pv_pipeline/scripts/create_BEST_reference_panel.R")
    run_cmd("bash run_DEploid_BEST_sites.sh")
    run_cmd("Rscript ~/tools/DEploid_Pv_pipeline/scripts/GT2hmmIBD_BEST_sites.R")

# use filtered VCF for DEploid
# Set up the parser
parser = argparse.ArgumentParser(description='DEploid pipeline wrapper for use on P. vivax',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input',type=str,help='Merged vcf file containing all required samples in compressed format',required=True)
#parser.add_argument('--miss', default=0.4, type=int, help='Percentage missingness allowed; 0.4 allows 40% missingness')
parser.add_argument('--version', action='version', version='1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

