#!/usr/bin/env python3
"""bcftools_roh_bed(WDL 파이프라인 원본 산출물, 이 저장소 cohort 분석엔 아직
안 쓰이던 자료)을 모아 군간 ROH(runs of homozygosity) 부담을 비교.

het_hom_peak_frac(VAF>0.9 문턱)과 달리 bcftools roh는 HMM 기반으로 동형접합
구간을 추정하므로 약간 다른 민감도 특성을 가짐 — 완전히 독립적인 지표는 아니지만
교차검증용으로 사용.

입력: filelist_h2o2.csv (phased_small_variant_vcf 경로에서 bcftools_roh_bed
경로를 유추: .../out/phased_small_variant_vcf/... -> .../out/bcftools_roh_bed/{sample}.GRCm39.bcftools_roh.bed)
출력: sample, group, n_roh_blocks, total_roh_bp, mean_roh_bp, median_roh_bp
"""
import argparse
import csv
import os
import statistics


def derive_roh_bed_path(phased_vcf_path, sample_id):
    # .../out/phased_small_variant_vcf/{sample}....vcf.gz -> .../out/bcftools_roh_bed/{sample}.GRCm39.bcftools_roh.bed
    out_dir = os.path.dirname(os.path.dirname(phased_vcf_path))  # .../out
    return os.path.join(out_dir, "bcftools_roh_bed", f"{sample_id}.GRCm39.bcftools_roh.bed")


def parse_roh_bed(path, exclude_chroms):
    """반환: (블록별 길이 리스트, {chrom: [(start,end), ...]} — gap 계산용)"""
    blocks = []
    by_chrom = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom in exclude_chroms:
                continue
            blocks.append(end - start)
            by_chrom.setdefault(chrom, []).append((start, end))
    return blocks, by_chrom


def load_chrom_sizes(fai_path, exclude_chroms):
    sizes = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.split("\t")
            chrom, length = parts[0], int(parts[1])
            if chrom not in exclude_chroms:
                sizes[chrom] = length
    return sizes


def non_roh_bp(by_chrom, chrom_sizes):
    """각 염색체에서 ROH 블록들 '사이/바깥' 구간(=잔여 이형접합 가능 구간) 총 길이.
    inbred 마우스는 게놈 대부분이 원래 동형접합이라 total_roh_bp가 거의 게놈
    크기만큼 나와 포화되므로(샘플 간 변별력 낮음), 그 여집합인 non_roh_bp가
    실제 잔여 이형접합/재조합 신호를 담은 더 민감한 지표가 된다."""
    total_non_roh = 0
    for chrom, size in chrom_sizes.items():
        blocks = sorted(by_chrom.get(chrom, []))
        covered = 0
        prev_end = 0
        for start, end in blocks:
            start = max(start, prev_end)
            if end > start:
                covered += end - start
                prev_end = max(prev_end, end)
        total_non_roh += max(size - covered, 0)
    return total_non_roh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filelist", default="filelist_h2o2.csv")
    ap.add_argument("--ref-fai", default="/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta.fai")
    ap.add_argument("--exclude-chroms", default="chrX,chrY")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    exclude_chroms = set(c for c in args.exclude_chroms.split(",") if c)
    chrom_sizes = load_chrom_sizes(args.ref_fai, exclude_chroms)

    rows = []
    with open(args.filelist) as f:
        for row in csv.DictReader(f):
            sample_id, group = row["sample_id"], row["group"]
            roh_bed = derive_roh_bed_path(row["phased_small_variant_vcf"], sample_id)
            if not os.path.exists(roh_bed):
                print(f"WARNING: 파일 없음, 스킵: {roh_bed}")
                continue
            blocks, by_chrom = parse_roh_bed(roh_bed, exclude_chroms)
            n = len(blocks)
            total_bp = sum(blocks)
            mean_bp = total_bp / n if n > 0 else float("nan")
            median_bp = statistics.median(blocks) if n > 0 else float("nan")
            non_roh = non_roh_bp(by_chrom, chrom_sizes)
            rows.append({
                "sample": sample_id, "group": group,
                "n_roh_blocks": n, "total_roh_bp": total_bp,
                "mean_roh_bp": round(mean_bp, 1), "median_roh_bp": median_bp,
                "non_roh_bp": non_roh,
            })
            print(f"{sample_id}: {n} blocks, {total_bp:,} bp ROH, {non_roh:,} bp non-ROH(잔여 이형접합 가능구간)")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "group", "n_roh_blocks", "total_roh_bp", "mean_roh_bp", "median_roh_bp", "non_roh_bp"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"저장: {args.output}")


if __name__ == "__main__":
    main()
