#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse
import pysam

parser = argparse.ArgumentParser(description="Convert insertions to duplications via bwa mem realignments")
parser.add_argument("-i", "--vcf_path", required=True, help="Input VCF")
parser.add_argument("-o", "--out_vcf", required=True, help="Output VCF (output is decompressed, .vcf)")
args = parser.parse_args()

# Open input VCF
vcf_in = pysam.VariantFile(args.vcf_path)
vcf_out = pysam.VariantFile(args.out_vcf, "w", header=vcf_in.header)

for record in vcf_in:
    repeat = record.info.get("TRFrepeat")
    if repeat is not None:
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()