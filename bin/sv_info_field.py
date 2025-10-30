#!/usr/bin/env python3

import re
import argparse

parser = argparse.ArgumentParser(description="Append SVTYPE and SVLEN to VCF INFO column.")
parser.add_argument("-i", "--input", required=True, help="Input VCF (must be gzipped, .vcf.gz)")
parser.add_argument("-o", "--output", required=True, help="Output VCF (output is decompressed, .vcf)")
args = parser.parse_args()

header_added = False

with open(args.input, 'rt') as f_in, open(args.output, "wt") as f_out:
    for line in f_in:
        if line.startswith("##"):
            # Write all existing meta-information lines
            f_out.write(line)
        elif line.startswith("#CHROM") and not header_added:
            # Insert new INFO definitions BEFORE the #CHROM line
            f_out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
            f_out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            f_out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n')
            
            # Write the #CHROM line itself
            f_out.write(line)
            header_added = True
        else:
            # This is a variant line (INFO column processing)
            fields = line.strip().split("\t")
            info = fields[7]
            pos = int(fields[1])

            # Extract ID part
            id_match = re.search(r'ID=([^;]+)', info)
            if id_match:
                id_str = id_match.group(1)
                # Extract SVTYPE (DEL, INS, etc.) from ID
                svtype_match = re.search(r'-(DEL|INS|COMPLEX)-', id_str)
                if svtype_match:
                    svtype = svtype_match.group(1)
                elif "SNV" in id_str:
                    svtype = "SNV"
                else:
                    svtype = "NA"

                # Extract SVLEN from last number after last hyphen
                svlen_match = re.search(r'-(\d+)$', id_str)
                svlen = int(svlen_match.group(1)) if svlen_match else 0

                if svtype == "DEL":
                    end = pos + svlen
                else:
                    end = pos + 1

                # Make negative for deletions
                if svtype == "DEL":
                    svlen = -svlen

                # Append SVTYPE and SVLEN to INFO
                info += f";END={end};SVTYPE={svtype};SVLEN={svlen}"
                fields[7] = info

            f_out.write("\t".join(fields) + "\n")
