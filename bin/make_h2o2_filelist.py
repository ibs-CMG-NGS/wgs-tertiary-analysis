#!/usr/bin/env python3
"""H2O2 12샘플(3개 드라이브 분산)용 filelist.csv 생성.

prepare_inputs.py의 find_wdl_output_dir/find_required_files를 재사용하되,
batch_results_dir 단일 가정 대신 샘플별 실행 디렉토리를 직접 지정한다.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_inputs import find_required_files, validate_files

SAMPLE_RUN_DIRS = {
    "H2O2-A01-ctrl1": ("A", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-A01-ctrl1/20260517_163205_humanwgs_singleton"),
    "H2O2-B01":       ("B", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-B01/20260518_163817_humanwgs_singleton"),
    "H2O2-C01":       ("C", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-C01/20260519_205414_humanwgs_singleton"),
    "H2O2-D01":       ("D", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-D01/20260520_184847_humanwgs_singleton"),
    "H2O2-2-A01":     ("A", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-A01/20260523_170354_humanwgs_singleton"),
    "H2O2-2-B01":     ("B", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-B01/20260524_171651_humanwgs_singleton"),
    "H2O2-2-C01":     ("C", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-C01/20260525_182243_humanwgs_singleton"),
    "H2O2-2-D01":     ("D", "/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-D01/20260526_181743_humanwgs_singleton"),
    "H2O2-3-A01":     ("A", "/mnt/hdd1_1tb/h2o2-wgs/H2O2-3-A01/20260710_102952_humanwgs_singleton"),
    "H2O2-3-B01":     ("B", "/mnt/hdd1_1tb/h2o2-wgs/H2O2-3-B01/20260711_070817_humanwgs_singleton"),
    "H2O2-3-C01":     ("C", "/mnt/hdd2_1tb/h2o2-wgs/H2O2-3-C01/20260712_023952_humanwgs_singleton"),
    "H2O2-3-D01":     ("D", "/mnt/hdd2_1tb/h2o2-wgs/H2O2-3-D01/20260713_134139_humanwgs_singleton"),
}

OUT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "filelist_h2o2.csv")


def main():
    rows = []
    errors = []
    for sample_id, (group, run_dir) in SAMPLE_RUN_DIRS.items():
        out_dir = os.path.join(run_dir, "out")
        if not os.path.isdir(out_dir):
            errors.append(f"{sample_id}: out 디렉토리 없음 -> {out_dir}")
            continue
        info = find_required_files(out_dir, sample_id)
        info["group"] = group
        missing = validate_files(info)
        status = "OK" if not missing else f"MISSING:{missing}"
        print(f"{sample_id:20s} [{group}]  {status}")
        if missing:
            errors.append(f"{sample_id}: 필수 파일 누락 - {missing}")
        rows.append(info)

    fieldnames = [
        "sample_id", "group",
        "phased_small_variant_vcf", "phased_small_variant_vcf_index",
        "phased_sv_vcf", "phased_sv_vcf_index",
        "cpg_combined_bed", "cpg_combined_bed_index", "cpg_combined_bw",
        "phased_trgt_vcf", "phased_trgt_vcf_index",
        "cpg_hap1_bed", "cpg_hap2_bed",
    ]
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (r.get(k) or "") for k in fieldnames})

    print(f"\n저장: {OUT_CSV} ({len(rows)}개 샘플)")
    if errors:
        print("\n=== 경고/오류 ===")
        for e in errors:
            print(f"  - {e}")


if __name__ == "__main__":
    main()
