#!/usr/bin/env python3
import sys, argparse, gzip

def open_maybe_gzip(path, mode="rt"):
    if path == "-" or path is None:
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def parse_gt(sample_field, fmt_keys):
    # Return GT string or None if not present
    parts = sample_field.split(":")
    if "GT" not in fmt_keys:
        return None
    idx = fmt_keys.index("GT")
    if idx >= len(parts):
        return None
    return parts[idx]

def is_missing_gt(gt):
    # Missing if GT is None, "." or contains any "."
    if gt is None:
        return True
    return gt == "." or "." in gt

def is_nonmissing_gt(gt):
    return (gt is not None) and (gt != ".") and ("." not in gt)

def main():
    ap = argparse.ArgumentParser(description="Filter VCF: keep records where A has missing GT and B,C do not.")
    ap.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz). Use '-' for stdin.")
    ap.add_argument("-o", "--output", default="-", help="Output VCF (plain text). Use '-' for stdout.")
    ap.add_argument("--sample-a", required=True, help="Sample name for A (must exist in VCF).")
    ap.add_argument("--sample-b", required=True, help="Sample name for B (must exist in VCF).")
    ap.add_argument("--sample-c", required=True, help="Sample name for C (must exist in VCF).")
    args = ap.parse_args()

    with open_maybe_gzip(args.input, "rt") as fin, open_maybe_gzip(args.output, "wt") as fout:
        sample_to_idx = {}
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                fout.write(line)
                header_cols = line.rstrip("\n").split("\t")
                samples = header_cols[9:]
                name_to_idx = {name: i for i, name in enumerate(samples)}
                for name in (args.sample_a, args.sample_b, args.sample_c):
                    if name not in name_to_idx:
                        sys.stderr.write(f"ERROR: Sample '{name}' not found in VCF header.\n")
                        sys.exit(2)
                # absolute column indices for samples in the full row
                sample_to_idx = {
                    "A": 9 + name_to_idx[args.sample_a],
                    "B": 9 + name_to_idx[args.sample_b],
                    "C": 9 + name_to_idx[args.sample_c],
                }
                continue

            if not sample_to_idx:
                sys.stderr.write("ERROR: No header line (#CHROM) found before records.\n")
                sys.exit(2)

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                # Not a valid VCF record with samples, pass through unchanged
                fout.write(line)
                continue

            fmt_keys = cols[8].split(":")
            try:
                gtA = parse_gt(cols[sample_to_idx["A"]], fmt_keys)
                gtB = parse_gt(cols[sample_to_idx["B"]], fmt_keys)
                gtC = parse_gt(cols[sample_to_idx["C"]], fmt_keys)
            except IndexError:
                # Malformed line; skip
                continue

            if is_missing_gt(gtA) and is_nonmissing_gt(gtB) and is_nonmissing_gt(gtC):
                fout.write(line)

if __name__ == "__main__":
    main()
