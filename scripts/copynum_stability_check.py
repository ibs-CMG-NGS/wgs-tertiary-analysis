#!/usr/bin/env python3
"""다중카피 유전자군 로커스의 웰 내 3-replicate 카피수 안정성 검사.

§1(SV burden)에서 반복적으로 문제됐던 다중카피 유전자군(Vmn1r/Mup/Pramel/Ear 등)
로커스에서, Sawfish가 원본 humanwgs_singleton WDL 실행 시 이미 계산해둔
depth 기반 카피수(sv_copynum_bedgraph)를 이용해 웰 내 3개 독립 생물학적
replicate 사이에 카피수가 얼마나 안정적인지 확인한다.

안정적(항상 같은 정수)이면 germline(개체 고유) 다형성, replicate마다 달라지면
배양·처치 중 발생하는 클론 동역학(우세 계통 교체)을 시사.
"""
import csv
import glob
import os
import sys

sys.path.insert(0, "/data_4tb/shared/wgs-tertiary-analysis/bin")
from make_h2o2_filelist import SAMPLE_RUN_DIRS

# §1에서 joint SV burden 유의 유전자들의 위치를 추적해 확인한 문제 로커스
LOCI = [
    ("chr7", 20124000, 20704000, "Vmn1r-a"),
    ("chr7", 20398000, 21210000, "Vmn1r-b"),
    ("chr7", 19762000, 20123000, "Vmn1r-c"),
    ("chr14", 43936000, 44406000, "Ear"),
    ("chr4", 143586000, 143953000, "Pramel"),
    ("chr4", 60256000, 60805000, "Mup"),
    ("chr6", 134876000, 135193000, "misc-DEL"),
]

GROUPS = {
    "A": ["H2O2-A01-ctrl1", "H2O2-2-A01", "H2O2-3-A01"],
    "B": ["H2O2-B01", "H2O2-2-B01", "H2O2-3-B01"],
    "C": ["H2O2-C01", "H2O2-2-C01", "H2O2-3-C01"],
    "D": ["H2O2-D01", "H2O2-2-D01", "H2O2-3-D01"],
}


def get_cn_at_locus(bedgraph_path, chrom, start, end):
    """구간과 겹치는 bedgraph 영역들 중 겹침 길이로 가중한 최빈 카피수."""
    overlaps = {}
    with open(bedgraph_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4 or parts[0] != chrom:
                continue
            s, e, cn = int(parts[1]), int(parts[2]), int(parts[3])
            ov_start, ov_end = max(s, start), min(e, end)
            if ov_start < ov_end:
                overlaps[cn] = overlaps.get(cn, 0) + (ov_end - ov_start)
    return max(overlaps, key=overlaps.get) if overlaps else None


def main():
    cn_by_sample = {}
    for sample_id, (group, run_dir) in SAMPLE_RUN_DIRS.items():
        bg = glob.glob(os.path.join(run_dir, "out", "sv_copynum_bedgraph", "*.bedgraph"))[0]
        cn_by_sample[sample_id] = {
            name: get_cn_at_locus(bg, c, s, e) for c, s, e, name in LOCI
        }

    raw_out = "/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results/cohort/copynum_locus_raw.csv"
    with open(raw_out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["locus", "chrom", "start", "end", "sample", "group", "copy_number"])
        for chrom, start, end, name in LOCI:
            for sample_id, (group, _) in SAMPLE_RUN_DIRS.items():
                writer.writerow([name, chrom, start, end, sample_id, group, cn_by_sample[sample_id][name]])
    print(f"저장: {raw_out}")

    rank_rows = []
    for _, _, _, name in LOCI:
        ranges = {}
        for g, samples in GROUPS.items():
            vals = [cn_by_sample[s][name] for s in samples]
            ranges[g] = max(vals) - min(vals)
        vals_sorted = sorted(ranges.values())
        ranks = {}
        for g in ranges:
            tied = [i + 1 for i, v in enumerate(vals_sorted) if v == ranges[g]]
            ranks[g] = sum(tied) / len(tied)
        rank_rows.append((name, ranges, ranks))

    mean_rank = {g: sum(r[2][g] for r in rank_rows) / len(rank_rows) for g in "ABCD"}

    detail_out = "/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results/cohort/copynum_instability_detail.csv"
    with open(detail_out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["locus", "group", "cn_replicate_1", "cn_replicate_2", "cn_replicate_3", "range", "rank"])
        for name, ranges, ranks in rank_rows:
            for g in "ABCD":
                vals = [cn_by_sample[s][name] for s in GROUPS[g]]
                writer.writerow([name, g, vals[0], vals[1], vals[2], ranges[g], ranks[g]])
    print(f"저장: {detail_out}")

    summary_out = "/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results/cohort/copynum_instability_rank.csv"
    with open(summary_out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["group", "mean_instability_rank"])
        for g in "ABCD":
            writer.writerow([g, round(mean_rank[g], 2)])
    print(f"저장: {summary_out}")
    print(mean_rank)


if __name__ == "__main__":
    main()
