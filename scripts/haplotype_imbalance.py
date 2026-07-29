#!/usr/bin/env python3
"""Haplotype 방향성 allelic imbalance 분석.

phased VCF의 PS(phase set) 블록 단위로, 이형접합 사이트들의 phase-0 allele
VAF가 0.5에서 한쪽으로 일관되게 치우치는지 검정한다. 개별 SNP를 독립적으로
보는 방식(§2)으로는 못 잡는, "특정 haplotype이 통째로 선택되는" 신호를 찾는다.

방법: 블록 내 각 사이트에서 phase-0 allele VAF가 0.5보다 큰지/작은지 부호를
매기고, 그 부호가 무작위(50/50)에서 벗어나는지 이항검정(exact binomial test,
scipy 없이 math.comb로 직접 계산)으로 확인한다.
"""
import argparse
import math
import subprocess
import sys
from collections import defaultdict


def _log_binom_pmf(n, i, p):
    log_choose = math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
    return log_choose + i * math.log(p) + (n - i) * math.log(1 - p)


def binom_test_two_sided(k, n, p=0.5):
    """정확 이항검정 양측 p-value (scipy.stats.binomtest와 동일한 정의).
    log-space로 계산해 n이 커도(math.comb 오버플로우 없이) 안전하게 처리."""
    if n == 0:
        return 1.0
    log_pmf = [_log_binom_pmf(n, i, p) for i in range(n + 1)]
    log_p_k = log_pmf[k]
    eps = 1e-9
    return min(1.0, sum(math.exp(lp) for lp in log_pmf if lp <= log_p_k + eps))


def max_vaf(vaf_str):
    try:
        return max(float(v) for v in vaf_str.split(",") if v not in (".", ""))
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--min-block-size", type=int, default=10,
                     help="검정 대상으로 삼을 최소 PS 블록 내 이형접합 사이트 수")
    ap.add_argument("--pvalue", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    out = subprocess.run(
        ["bcftools", "query", "-f", '%CHROM\t%FILTER\t[%GT]\t[%PS]\t[%VAF]\n', args.vcf],
        capture_output=True, text=True, check=True,
    )

    blocks = defaultdict(list)  # (chrom, ps) -> list of phase0_frac
    for line in out.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 5:
            continue
        chrom, filt, gt, ps, vaf_str = parts
        if filt != "PASS" or ps in (".", ""):
            continue
        if gt not in ("0|1", "1|0"):
            continue
        vaf = max_vaf(vaf_str)
        if vaf is None:
            continue
        # phase-0(파이프 왼쪽) allele의 VAF: GT가 1|0이면 phase0=ALT이므로 VAF 그대로,
        # 0|1이면 phase0=REF이므로 1-VAF
        phase0_frac = vaf if gt == "1|0" else (1 - vaf)
        blocks[(chrom, ps)].append(phase0_frac)

    n_blocks_tested = 0
    n_biased_blocks = 0
    abs_biases = []
    for (chrom, ps), fracs in blocks.items():
        if len(fracs) < args.min_block_size:
            continue
        n_blocks_tested += 1
        k = sum(1 for f in fracs if f > 0.5)
        n = sum(1 for f in fracs if f != 0.5)
        if n == 0:
            continue
        pval = binom_test_two_sided(k, n)
        mean_dev = sum(f - 0.5 for f in fracs) / len(fracs)
        abs_biases.append(abs(mean_dev))
        if pval < args.pvalue:
            n_biased_blocks += 1

    frac_biased = (n_biased_blocks / n_blocks_tested) if n_blocks_tested > 0 else float("nan")
    mean_abs_bias = (sum(abs_biases) / len(abs_biases)) if abs_biases else float("nan")

    with open(args.output, "w") as f:
        f.write("sample\tgroup\tn_ps_blocks_tested\tn_biased_blocks\tfrac_biased\tmean_abs_bias\n")
        f.write(f"{args.sample}\t{args.group}\t{n_blocks_tested}\t{n_biased_blocks}\t"
                f"{frac_biased:.6f}\t{mean_abs_bias:.6f}\n")

    print(f"저장: {args.output} (블록 {n_blocks_tested}개 중 편향 {n_biased_blocks}개)", file=sys.stderr)


if __name__ == "__main__":
    main()
