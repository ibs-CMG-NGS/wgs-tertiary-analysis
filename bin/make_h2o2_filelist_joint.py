#!/usr/bin/env python3
"""H2O2 12샘플용 joint-called VCF 기반 filelist.csv 생성.

기존 make_h2o2_filelist.py는 개별(single-sample) humanwgs_singleton 출력을 사용했다.
이 스크립트는 phased_small_variant_vcf/phased_sv_vcf만 joint calling(GLnexus+Sawfish,
웰 내 3개 시퀀싱 배치 합침) 결과로 교체하고, methylation(cpg_combined_bed 등)과 TRGT VCF는
joint calling 대상이 아니므로 기존 개별 파일을 그대로 재사용한다.
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_h2o2_filelist import SAMPLE_RUN_DIRS
from prepare_inputs import find_required_files, validate_files

WELL_JOINT_DIRS = {
    "A01": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-A01/20260714_054055_joint",
    "B01": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-B01/20260714_062818_joint",
    "C01": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-C01/20260714_071318_joint",
    "D01": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-D01/20260714_075510_joint",
}

OUT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "filelist_h2o2_joint.csv")


def find_joint_vcf(joint_dir, subdir_name, sample_id):
    """out/{subdir_name}/{0,1,2}/ 안에서 sample_id로 시작하는 VCF를 찾는다 (인덱스 순서에 의존 안 함)."""
    pattern = os.path.join(joint_dir, "out", subdir_name, "*", f"{sample_id}.*.vcf.gz")
    matches = glob.glob(pattern)
    if not matches:
        return None, None
    vcf = matches[0]
    idx_pattern = os.path.join(joint_dir, "out", subdir_name + "_indices", "*", f"{sample_id}.*.vcf.gz.tbi")
    idx_matches = glob.glob(idx_pattern)
    return vcf, (idx_matches[0] if idx_matches else None)


def main():
    rows = []
    errors = []
    for sample_id, (group, run_dir) in SAMPLE_RUN_DIRS.items():
        # 기존 개별 실행 결과에서 methylation/TRGT 등 joint 대상이 아닌 파일은 그대로 재사용
        out_dir = os.path.join(run_dir, "out")
        if not os.path.isdir(out_dir):
            errors.append(f"{sample_id}: out 디렉토리 없음 -> {out_dir}")
            continue
        info = find_required_files(out_dir, sample_id)
        info["group"] = group

        # well 판별: 샘플ID에서 "A01"/"B01"/"C01"/"D01" 추출
        well = next((w for w in WELL_JOINT_DIRS if w in sample_id), None)
        if well is None:
            errors.append(f"{sample_id}: well을 판별할 수 없음")
            rows.append(info)
            continue
        joint_dir = WELL_JOINT_DIRS[well]

        sv_vcf, sv_idx = find_joint_vcf(joint_dir, "split_joint_small_variant_vcfs", sample_id)
        struct_vcf, struct_idx = find_joint_vcf(joint_dir, "split_joint_structural_variant_vcfs", sample_id)

        if sv_vcf is None:
            errors.append(f"{sample_id}: joint small variant VCF 없음 ({joint_dir})")
        else:
            info["phased_small_variant_vcf"] = sv_vcf
            info["phased_small_variant_vcf_index"] = sv_idx

        if struct_vcf is None:
            errors.append(f"{sample_id}: joint SV VCF 없음 ({joint_dir})")
        else:
            info["phased_sv_vcf"] = struct_vcf
            info["phased_sv_vcf_index"] = struct_idx

        missing = validate_files(info)
        status = "OK" if not missing else f"MISSING:{missing}"
        print(f"{sample_id:20s} [{group}] well={well}  {status}")
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
