#!/usr/bin/env python3
"""분자(read) 단위 haplotype 일관성 확인 (1차: HP 태그 미배정 비율 프록시).

HiFi 롱리드는 whatshap/hiphase 등으로 phasing된 뒤 BAM에 HP:i:1/HP:i:2 태그가
붙는다. 진짜 깔끔한 2-클론(diploid) 모델이면 대부분의 read가 HP=1 또는 HP=2로
명확히 배정돼야 한다. 서브클론이 2개보다 많이 섞여 있으면 phasing 알고리즘이
국소적으로 흔들려 "미배정(HP 태그 없음)" read 비율이 올라갈 수 있다는 가설을
간단한 프록시로 검사한다.

전체 게놈을 다 스캔하면 12개 샘플(각 ~40GB BAM) 기준 너무 느려서, 상염색체
여러 곳에서 대표 윈도우를 뽑아 표본 조사한다.
"""
import argparse
import subprocess
import sys

# 상염색체 대표 윈도우 (5Mb씩, chrX/Y 제외 — 성별 교란 방지)
WINDOWS = [
    "chr1:20000000-25000000",
    "chr3:30000000-35000000",
    "chr5:50000000-55000000",
    "chr7:40000000-45000000",
    "chr9:60000000-65000000",
    "chr11:70000000-75000000",
    "chr14:30000000-35000000",
    "chr17:20000000-25000000",
]


def hp_tag_stats(bam, region):
    out = subprocess.run(["samtools", "view", bam, region], capture_output=True, text=True, check=True)
    total = hp1 = hp2 = 0
    for line in out.stdout.splitlines():
        total += 1
        if "\tHP:i:1" in line:
            hp1 += 1
        elif "\tHP:i:2" in line:
            hp2 += 1
    return total, hp1, hp2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    total_all = hp1_all = hp2_all = 0
    for region in WINDOWS:
        t, h1, h2 = hp_tag_stats(args.bam, region)
        total_all += t
        hp1_all += h1
        hp2_all += h2

    unassigned = total_all - hp1_all - hp2_all
    unassigned_frac = unassigned / total_all if total_all > 0 else float("nan")

    with open(args.output, "w") as f:
        f.write("sample\tgroup\ttotal_reads_sampled\thp1_reads\thp2_reads\tunassigned_reads\tunassigned_frac\n")
        f.write(f"{args.sample}\t{args.group}\t{total_all}\t{hp1_all}\t{hp2_all}\t{unassigned}\t{unassigned_frac:.6f}\n")

    print(f"저장: {args.output} (미배정 비율={unassigned_frac:.4f})", file=sys.stderr)


if __name__ == "__main__":
    main()
