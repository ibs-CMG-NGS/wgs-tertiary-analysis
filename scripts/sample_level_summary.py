#!/usr/bin/env python3
"""샘플 단위 요약통계 계산 (군간비교용, 유전자별 다중검정 대신 샘플당 1행).

- total_pass_variants / total_pass_sv: PASS 레코드 수
- het_hom_peak_frac: heterozygous(0/1) PASS SNV 중 VAF>0.9인 비율
  (이형접합으로 콜됐는데 VAF가 동형접합처럼 높다 = 세포집단 내 클론 우세/LOH 신호)
- median_het_vaf: heterozygous PASS SNV의 VAF 중앙값 (이질적 집단일수록 0.5에 가까움,
  균질/클론성 집단일수록 0.5에서 벗어나 한쪽으로 치우침)
"""
import argparse
import statistics
import subprocess
import sys


def run_bcftools_query(vcf, fmt, exclude_chroms=None):
    cmd = ["bcftools", "query"]
    if exclude_chroms:
        cmd += ["-t", "^" + ",".join(exclude_chroms)]
    cmd += ["-f", fmt, vcf]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return out.stdout


def count_pass(vcf, exclude_chroms=None):
    cmd = ["bcftools", "view", "-H", "-f", "PASS"]
    if exclude_chroms:
        cmd += ["-t", "^" + ",".join(exclude_chroms)]
    cmd += [vcf]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return sum(1 for _ in out.stdout.splitlines() if _.strip())


def max_vaf(vaf_str):
    """멀티알렐릭 콤마구분 VAF 중 최댓값 반환 (float 변환 실패 시 None)"""
    try:
        return max(float(v) for v in vaf_str.split(",") if v not in (".", ""))
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="phased_small_variant_vcf (single-sample)")
    ap.add_argument("--sv-vcf", required=True, help="phased_sv_vcf (single-sample)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--exclude-chroms", default="",
                     help="콤마구분 제외 염색체 (예: chrX,chrY) — 군간 성별 구성이 다르면 성염색체 신호가 교란요인이 됨")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    exclude_chroms = [c for c in args.exclude_chroms.split(",") if c]

    total_pass_variants = count_pass(args.vcf, exclude_chroms)
    total_pass_sv = count_pass(args.sv_vcf, exclude_chroms)

    # heterozygous PASS SNV의 VAF 분포
    raw = run_bcftools_query(args.vcf, '%FILTER\t[%GT]\t[%VAF]\n', exclude_chroms)
    het_vafs = []
    for line in raw.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        filt, gt, vaf_str = parts
        if filt != "PASS":
            continue
        # 0/1, 0|1 등 이형접합만 (단상성/동형접합 제외)
        alleles = gt.replace("|", "/").split("/")
        if len(alleles) != 2 or alleles[0] == alleles[1] or "." in alleles:
            continue
        v = max_vaf(vaf_str)
        if v is not None:
            het_vafs.append(v)

    n_het = len(het_vafs)
    hom_peak_frac = (sum(1 for v in het_vafs if v > 0.9) / n_het) if n_het > 0 else float("nan")
    median_het_vaf = statistics.median(het_vafs) if n_het > 0 else float("nan")

    with open(args.output, "w") as f:
        f.write("sample\tgroup\ttotal_pass_variants\ttotal_pass_sv\tn_het_snv\thet_hom_peak_frac\tmedian_het_vaf\n")
        f.write(f"{args.sample}\t{args.group}\t{total_pass_variants}\t{total_pass_sv}\t"
                f"{n_het}\t{hom_peak_frac:.6f}\t{median_het_vaf:.6f}\n")

    print(f"저장: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
