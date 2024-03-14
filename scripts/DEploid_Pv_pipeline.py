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
import threading
import multiprocessing
import glob
import re

def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")

def main(args):
    
    run_cmd("mkdir individual_output")
    run_cmd("mkdir BEST")
    run_cmd("bcftools view -v snps %(vcf)s -Oz -o merged_snps.vcf.gz" % vars(args))
    #run_cmd(f'bcftools view -v snps {args.vcf} -Oz -o merged_snps.vcf.gz')
    #run_cmd("bcftools filter -i 'FMT/DP>4' -S . merged_snps.vcf.gz -Oz -o merged.filt.snps.vcf.gz")
    #run_cmd("plink --const-fid --vcf merged.filt.snps.vcf.gz --mind %(miss)s --recode vcf --allow-extra-chr --out merged_plink" % vars(args))
    #run_cmd("grep -P \"^#CHROM\" merged_plink.vcf | awk '{ $1=\"\"; $2=\"\";$3=\"\"; $4=\"\";$5=\"\"; $6=\"\";$7=\"\"; $8=\"\";$9=\"\"; print}' | sed 's/ /\\n/g' | tail -n+10 | sed 's/^0_//'  > merged_plink_new")
    #run_cmd("bcftools view -S mhaps_merged_plink_new --threads 4 -O z -o  merged.miss0.4.filt.snps.vcf.gz merged.filt.snps.vcf.gz")
    #run_cmd("bcftools view merged.miss0.4.filt.snps.vcf.gz | setGT.py --fraction 0.8 | bcftools view -O z -c 1 -o merged.filt.GT.miss0.4.snps.vcf.gz")
    #run_cmd("bcftools view -m2 -M2 merged.filt.GT.miss0.4.snps.vcf.gz --threads 4 -O z -o merged.bi.filt.GT.miss0.4.snps.vcf.gz")

# use filtered VCF for DEploid

# Set up the parser
parser = argparse.ArgumentParser(description='DEploid pipeline wrapper for use on P. vivax',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='Merged vcf file containing all required samples; Use zipped format, i.e. vcf.gz',required=True)
parser.add_argument('--miss', default=0.4,type=int,help='Percentage missingness allowed; 0.4 allows 40% missingness')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
