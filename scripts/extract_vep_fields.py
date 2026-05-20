#!/usr/bin/env python3
"""
VCF CSQ 필드에서 canonical 전사체의 SYMBOL + IMPACT 추출.
bcftools +split-vep 플러그인을 사용할 수 없을 때 fallback으로 사용.

출력 컬럼 (탭 구분):
    CHROM  POS  REF  ALT  SYMBOL  IMPACT  Consequence  GT

CSQ 형식은 VCF 헤더의 ##INFO=<ID=CSQ,...,Description="...Format: ..."> 에서 자동 파싱.
"""

import sys
import gzip
import argparse
import re


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_csq_format(header_lines):
    """##INFO=<ID=CSQ,...> 헤더에서 필드 목록 파싱."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ,"):
            m = re.search(r'Format: ([^"]+)"', line)
            if m:
                return m.group(1).strip().split("|")
    return []


def extract_canonical(csq_str, fields):
    """CSQ 문자열에서 CANONICAL=YES인 첫 번째 entry 반환 (dict)."""
    if not csq_str or csq_str == ".":
        return None

    canonical_idx  = fields.index("CANONICAL")  if "CANONICAL"  in fields else -1
    symbol_idx     = fields.index("SYMBOL")      if "SYMBOL"     in fields else -1
    impact_idx     = fields.index("IMPACT")      if "IMPACT"     in fields else -1
    conseq_idx     = fields.index("Consequence") if "Consequence" in fields else -1

    for entry in csq_str.split(","):
        vals = entry.split("|")
        if len(vals) <= max(canonical_idx, symbol_idx, impact_idx, conseq_idx):
            continue
        if canonical_idx >= 0 and vals[canonical_idx] != "YES":
            continue
        return {
            "SYMBOL":      vals[symbol_idx]  if symbol_idx  >= 0 else ".",
            "IMPACT":      vals[impact_idx]  if impact_idx  >= 0 else ".",
            "Consequence": vals[conseq_idx]  if conseq_idx  >= 0 else ".",
        }
    return None


def main():
    parser = argparse.ArgumentParser(description="VCF CSQ canonical impact extractor")
    parser.add_argument("--input",  required=True, help="VEP-annotated VCF(.gz)")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    header_lines = []
    with open_vcf(args.input) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line.rstrip())
            else:
                break

    csq_fields = parse_csq_format(header_lines)
    if not csq_fields:
        print("WARNING: CSQ format not found in VCF header; output will be empty.", file=sys.stderr)

    with open_vcf(args.input) as fh, open(args.output, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tSYMBOL\tIMPACT\tConsequence\tGT\n")

        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            chrom, pos, _, ref, alt, _, _, info = parts[:8]
            fmt_keys = parts[8].split(":")
            samples   = parts[9:]

            # CSQ 추출
            csq_str = "."
            for field in info.split(";"):
                if field.startswith("CSQ="):
                    csq_str = field[4:]
                    break

            canon = extract_canonical(csq_str, csq_fields) if csq_fields else None
            symbol  = canon["SYMBOL"]      if canon else "."
            impact  = canon["IMPACT"]      if canon else "."
            conseq  = canon["Consequence"] if canon else "."

            # GT 수집 (샘플이 여러 개이면 첫 번째만)
            gt = "."
            if samples and "GT" in fmt_keys:
                gt_idx = fmt_keys.index("GT")
                vals = samples[0].split(":")
                gt = vals[gt_idx] if gt_idx < len(vals) else "."

            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{symbol}\t{impact}\t{conseq}\t{gt}\n")


if __name__ == "__main__":
    main()
