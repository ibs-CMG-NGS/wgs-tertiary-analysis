#!/usr/bin/env python3
"""이형접합 PASS SNV의 VAF 전체 분포(히스토그램)를 샘플별로 뽑아 군간 비교.

het_hom_peak_frac(VAF>0.9 이진 문턱)은 우세도(f)와 LOH 누적량이 곱해진 값이라
둘을 못 가른다 — 대신 VAF 분포 전체 모양을 보면, 아직 0.9를 못 넘었지만
진행 중인 "subclonal" LOH(예: VAF 0.6~0.85)가 얼마나 있는지 볼 수 있어 우세도와
어느 정도 독립적인 각도를 제공한다 (완전히 분리되진 않음, 참고 지표).

sample_level_summary.py와 동일한 방식(bcftools query, PASS만, 이형접합만,
멀티알렐릭은 max VAF)으로 VAF를 뽑아 일관성을 유지.

출력 1: {output-summary} — sample, group, n_het, frac_[bin] 컬럼들 +
        frac_subclonal_mid(0.6~0.85, 참고용 핵심 요약치)
출력 2: {output-hist} — long format(sample, group, bin_start, count) — 밀도곡선 그리기용
"""
import argparse
import csv
import subprocess


BIN_WIDTH = 0.02
BINS = [round(i * BIN_WIDTH, 2) for i in range(int(1 / BIN_WIDTH))]  # 0.00, 0.02, ..., 0.98


def run_bcftools_query(vcf, exclude_chroms):
    cmd = ["bcftools", "query"]
    if exclude_chroms:
        cmd += ["-t", "^" + ",".join(exclude_chroms)]
    cmd += ["-f", '%FILTER\t[%GT]\t[%VAF]\n', vcf]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return out.stdout


def max_vaf(vaf_str):
    try:
        return max(float(v) for v in vaf_str.split(",") if v not in (".", ""))
    except (ValueError, TypeError):
        return None


def het_vafs_for_sample(vcf, exclude_chroms):
    raw = run_bcftools_query(vcf, exclude_chroms)
    vafs = []
    for line in raw.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        filt, gt, vaf_str = parts
        if filt != "PASS":
            continue
        alleles = gt.replace("|", "/").split("/")
        if len(alleles) != 2 or alleles[0] == alleles[1] or "." in alleles:
            continue
        v = max_vaf(vaf_str)
        if v is not None:
            vafs.append(v)
    return vafs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filelist", default="filelist_h2o2.csv")
    ap.add_argument("--exclude-chroms", default="chrX,chrY")
    ap.add_argument("--output-summary", required=True)
    ap.add_argument("--output-hist", required=True)
    args = ap.parse_args()
    exclude_chroms = [c for c in args.exclude_chroms.split(",") if c]

    summary_rows = []
    hist_rows = []
    with open(args.filelist) as f:
        for row in csv.DictReader(f):
            sample_id, group, vcf = row["sample_id"], row["group"], row["phased_small_variant_vcf"]
            vafs = het_vafs_for_sample(vcf, exclude_chroms)
            n = len(vafs)
            counts = [0] * len(BINS)
            for v in vafs:
                idx = min(int(v / BIN_WIDTH), len(BINS) - 1)
                counts[idx] += 1
            for b, c in zip(BINS, counts):
                hist_rows.append({"sample": sample_id, "group": group, "bin_start": b,
                                   "count": c, "density": c / n if n > 0 else 0})
            frac_subclonal_mid = sum(1 for v in vafs if 0.6 <= v < 0.85) / n if n > 0 else float("nan")
            frac_peak = sum(1 for v in vafs if v > 0.9) / n if n > 0 else float("nan")
            summary_rows.append({
                "sample": sample_id, "group": group, "n_het": n,
                "frac_subclonal_mid_0.6_0.85": round(frac_subclonal_mid, 6),
                "frac_peak_gt_0.9": round(frac_peak, 6),
            })
            print(f"{sample_id}: n_het={n}, frac[0.6,0.85)={frac_subclonal_mid:.4f}, frac>0.9={frac_peak:.4f}")

    with open(args.output_summary, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "group", "n_het", "frac_subclonal_mid_0.6_0.85", "frac_peak_gt_0.9"])
        writer.writeheader()
        writer.writerows(summary_rows)

    with open(args.output_hist, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "group", "bin_start", "count", "density"])
        writer.writeheader()
        writer.writerows(hist_rows)

    print(f"저장: {args.output_summary}, {args.output_hist}")


if __name__ == "__main__":
    main()
