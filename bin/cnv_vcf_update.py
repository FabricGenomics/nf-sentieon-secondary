#!/usr/bin/env python3

import sys
import os
import re
import subprocess
import argparse
import pysam

parser = argparse.ArgumentParser(description="Update ALT and SVTYPE fields based on sample OT values")
parser.add_argument("-i", "--vcf_in", required=True, help="Input VCF")
parser.add_argument("-o", "--vcf_out", required=True, help="Output VCF")
args = parser.parse_args()

# Open input VCF (first pass to discover required ALT IDs)
vcf_in = pysam.VariantFile(args.vcf_in)
samples = list(vcf_in.header.samples)

# Collect symbolic ALT IDs that will be written
needed_alt_ids = set()
for record in vcf_in:
    original_type1 = record.samples[samples[0]]["OT"]
    original_type2 = record.samples[samples[1]]["OT"]
    if original_type1 == "CNV":
        alt_id = original_type2
    else:
        alt_id = original_type1
    if isinstance(alt_id, bytes):
        alt_id = alt_id.decode()
    needed_alt_ids.add(str(alt_id))

# Prepare output header by rebuilding and overriding INFO/STRANDS
orig_header = vcf_in.header
header = pysam.VariantHeader()

# Copy all header lines except existing STRANDS definitions (so we can override INFO only)
for rec in orig_header.records:
    try:
        if rec.get('ID') == 'STRANDS' and rec.key in ('INFO', 'FORMAT'):
            continue
    except Exception:
        pass
    header.add_line(str(rec))

# Copy samples
for s in orig_header.samples:
    header.add_sample(s)

# Gather existing ALT IDs from rebuilt header
try:
    existing_alt_ids = set(header.alts.keys())
except Exception:
    existing_alt_ids = {rec.get('ID') for rec in header.records if getattr(rec, 'key', None) == 'ALT'}

# Descriptions for common SV types; default to generic if unknown
alt_descriptions = {
    "DEL": "Deletion",
    "DUP": "Duplication",
    "INS": "Insertion",
    "INV": "Inversion",
    "CNV": "copy-number variants",
    "BND": "Breakend",
}

for alt_id in sorted(needed_alt_ids):
    if alt_id and alt_id not in existing_alt_ids:
        desc = alt_descriptions.get(alt_id, "Symbolic allele")
        header.add_meta('ALT', items=[('ID', alt_id), ('Description', desc)])

# Ensure SVTYPE info definition exists
if 'SVTYPE' not in header.info:
    header.info.add('SVTYPE', number=1, type='String', description='Type of structural variant')
# Intentionally do NOT add STRANDS; we remove it entirely from output
header.formats.add('CN', number=1, type='Integer', description='Copy number per GT')
header.info.add('COPIES', number=1, type='Integer', description='Copy number')

# Open output with augmented header
vcf_out = pysam.VariantFile(args.vcf_out, "w", header=header)

# Second pass: reopen input and write updated records
vcf_in.close()
vcf_in = pysam.VariantFile(args.vcf_in)

if 'CN' not in vcf_in.header.formats:
    vcf_in.header.formats.add('CN', number=1, type='Integer', description='Copy number per GT')
if 'COPIES' not in vcf_in.header.info:
    vcf_in.header.info.add('COPIES', number=1, type='Integer', description='Copy number')
if 'STRANDS' not in vcf_in.header.info:
    vcf_in.header.info.add('STRANDS', number=1, type='String', 
    description='Strand orientation for SV (e.g. +-, -+, ++, --)')

for record in vcf_in:
    original_type1 = record.samples[samples[0]]["OT"]
    original_type2 = record.samples[samples[1]]["OT"]
    # Extract copy number from ID (robust): matches ..._copy_num_15 or ..._copynum_15
    num = None
    id_string = record.id or ''
    m = re.search(r'(?:^|_)copy_?num_(\d+)(?:$|_)', id_string)
    if m:
        num = int(m.group(1))
    # Uncomment for quick debugging if needed
    else:
        print(f"WARN: no NCOPIES parsed for ID={id_string}", file=sys.stderr)
    if isinstance(original_type1, bytes):
        original_type1 = original_type1.decode()
    if isinstance(original_type2, bytes):
        original_type2 = original_type2.decode()

    if original_type1 == "CNV":
        new_alt_id = original_type2
        new_svtype = original_type1
    else:
        new_alt_id = original_type1
        new_svtype = original_type2

    new_alt = f"<{new_alt_id}>"
    record.alts = (new_alt,)
    record.info["SVTYPE"] = new_svtype
    # Remove STRANDS entirely (both INFO and FORMAT) and any stray INFO/CN
    record.info.pop('STRANDS', None)
    record.info.pop('CN', None)
    record.info["COPIES"] = num
    for s in samples:
        record.samples[s]["CN"] = num
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
