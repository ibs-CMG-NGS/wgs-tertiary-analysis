#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDL 파이프라인 결과에서 3차 분석 입력 파일 자동 탐색 및 매핑
"""

import os
import sys
import glob
import argparse
from pathlib import Path
import csv
import yaml


def find_wdl_output_dir(batch_results_dir, sample_id):
    """
    batch_results 디렉토리에서 특정 샘플의 WDL 출력 디렉토리 찾기
    
    구조: batch_results/{sample_id}/{timestamp}_humanwgs_singleton/out/
    """
    sample_dir = os.path.join(batch_results_dir, sample_id)
    
    if not os.path.exists(sample_dir):
        raise FileNotFoundError(f"샘플 디렉토리를 찾을 수 없습니다: {sample_dir}")
    
    # timestamp_humanwgs_singleton 형태의 디렉토리 찾기
    wdl_dirs = glob.glob(os.path.join(sample_dir, "*_humanwgs_singleton"))
    
    if not wdl_dirs:
        raise FileNotFoundError(f"WDL 실행 디렉토리를 찾을 수 없습니다: {sample_dir}/*_humanwgs_singleton")
    
    # 가장 최근 것 선택 (여러 개 있을 경우)
    wdl_dir = sorted(wdl_dirs)[-1]
    out_dir = os.path.join(wdl_dir, "out")
    
    if not os.path.exists(out_dir):
        raise FileNotFoundError(f"'out' 디렉토리를 찾을 수 없습니다: {out_dir}")
    
    return out_dir


def _find_vcf(out_dir, subdir_name):
    """out_dir 하위 subdir_name 폴더에서 VCF 파일과 인덱스를 반환한다."""
    d = os.path.join(out_dir, subdir_name)
    if not (os.path.exists(d) and os.path.isdir(d)):
        return None, None
    vcf_files = glob.glob(os.path.join(d, "*.vcf.gz"))
    if not vcf_files:
        return None, None
    vcf = vcf_files[0]
    idx = next((vcf + ext for ext in (".tbi", ".csi") if os.path.exists(vcf + ext)), None)
    return vcf, idx


def _find_bed(out_dir, subdir_name):
    """out_dir 하위 subdir_name 폴더에서 BED(.gz) 파일과 인덱스를 반환한다."""
    d = os.path.join(out_dir, subdir_name)
    if not (os.path.exists(d) and os.path.isdir(d)):
        return None, None
    bed_files = glob.glob(os.path.join(d, "*.bed.gz")) or glob.glob(os.path.join(d, "*.bed"))
    if not bed_files:
        return None, None
    bed = bed_files[0]
    idx = bed + ".tbi" if os.path.exists(bed + ".tbi") else None
    return bed, idx


def find_required_files(out_dir, sample_id):
    """
    out 디렉토리에서 3차 분석에 필요한 파일들 찾기

    필수 파일:
    - phased_small_variant_vcf, phased_sv_vcf, cpg_combined_bed, cpg_combined_bw

    선택 파일 (없어도 경고만):
    - phased_trgt_vcf      : TRGT 반복서열 분석용
    - cpg_hap1_bed         : Allele-specific methylation 분석용
    - cpg_hap2_bed         : Allele-specific methylation 분석용
    """
    files = {'sample_id': sample_id}

    # ── 필수 파일 ──────────────────────────────────────────────────────────────
    vcf, idx = _find_vcf(out_dir, "phased_small_variant_vcf")
    files['phased_small_variant_vcf'] = vcf
    files['phased_small_variant_vcf_index'] = idx

    sv, sv_idx = _find_vcf(out_dir, "phased_sv_vcf")
    files['phased_sv_vcf'] = sv
    files['phased_sv_vcf_index'] = sv_idx

    bed, bed_idx = _find_bed(out_dir, "cpg_combined_bed")
    files['cpg_combined_bed'] = bed
    files['cpg_combined_bed_index'] = bed_idx

    bw_dir = os.path.join(out_dir, "cpg_combined_bw")
    bw_files = (glob.glob(os.path.join(bw_dir, "*.bw")) +
                glob.glob(os.path.join(bw_dir, "*.bigWig"))) if os.path.isdir(bw_dir) else []
    files['cpg_combined_bw'] = bw_files[0] if bw_files else None

    # ── 선택 파일 ──────────────────────────────────────────────────────────────
    trgt, trgt_idx = _find_vcf(out_dir, "phased_trgt_vcf")
    files['phased_trgt_vcf'] = trgt or ''
    files['phased_trgt_vcf_index'] = trgt_idx or ''
    if trgt:
        print(f"     └─ TRGT VCF 발견: {os.path.basename(trgt)}")

    hap1, _ = _find_bed(out_dir, "cpg_hap1_bed")
    files['cpg_hap1_bed'] = hap1 or ''
    if hap1:
        print(f"     └─ Hap1 BED 발견: {os.path.basename(hap1)}")

    hap2, _ = _find_bed(out_dir, "cpg_hap2_bed")
    files['cpg_hap2_bed'] = hap2 or ''
    if hap2:
        print(f"     └─ Hap2 BED 발견: {os.path.basename(hap2)}")

    return files


def validate_files(file_info):
    """필수 파일이 모두 존재하는지 확인"""
    required = ['phased_small_variant_vcf', 'phased_sv_vcf', 'cpg_combined_bed']
    missing = []
    
    for req in required:
        if file_info[req] is None:
            missing.append(req)
    
    return missing


def create_filelist_csv(batch_results_dir, samples, output_csv):
    """
    모든 샘플에 대한 파일 목록 CSV 생성
    
    CSV 형식:
    sample_id,group,phased_small_variant_vcf,phased_sv_vcf,cpg_combined_bed,cpg_combined_bw
    """
    all_files = []
    errors = []
    
    print(f"\n{'='*80}")
    print(f"WDL 결과 파일 탐색 시작")
    print(f"{'='*80}\n")
    print(f"Batch Results 디렉토리: {batch_results_dir}\n")
    
    for group, sample_list in samples.items():
        print(f"[{group.upper()}] 그룹 처리 중...")
        
        for sample_id in sample_list:
            print(f"  └─ {sample_id} 탐색 중...", end=" ")
            
            try:
                # WDL 출력 디렉토리 찾기
                out_dir = find_wdl_output_dir(batch_results_dir, sample_id)
                
                # 필수 파일 찾기
                file_info = find_required_files(out_dir, sample_id)
                file_info['group'] = group
                
                # 검증
                missing = validate_files(file_info)
                if missing:
                    error_msg = f"샘플 {sample_id}: 필수 파일 누락 - {', '.join(missing)}"
                    errors.append(error_msg)
                    print(f"⚠️  경고 (일부 파일 누락)")
                else:
                    print(f"✓ 완료")
                
                all_files.append(file_info)
                
            except Exception as e:
                error_msg = f"샘플 {sample_id}: {str(e)}"
                errors.append(error_msg)
                print(f"✗ 실패")
                print(f"     오류: {str(e)}")
    
    # CSV 저장
    if all_files:
        print(f"\n파일 목록 저장 중: {output_csv}")
        
        with open(output_csv, 'w', newline='') as f:
            fieldnames = [
                'sample_id', 'group',
                'phased_small_variant_vcf', 'phased_small_variant_vcf_index',
                'phased_sv_vcf', 'phased_sv_vcf_index',
                'cpg_combined_bed', 'cpg_combined_bed_index',
                'cpg_combined_bw',
                # 선택적 파일
                'phased_trgt_vcf', 'phased_trgt_vcf_index',
                'cpg_hap1_bed', 'cpg_hap2_bed',
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_files)
        
        print(f"✓ 완료: {len(all_files)}개 샘플의 파일 정보 저장\n")
    
    # 에러 요약
    if errors:
        print(f"\n{'='*80}")
        print(f"⚠️  경고 및 오류 ({len(errors)}건)")
        print(f"{'='*80}\n")
        for error in errors:
            print(f"  - {error}")
        print()
    
    # 최종 요약
    print(f"{'='*80}")
    print(f"처리 요약")
    print(f"{'='*80}")
    print(f"  총 샘플 수: {len(all_files)}")
    print(f"  성공: {len(all_files) - len(errors)}")
    print(f"  경고/오류: {len(errors)}")
    print(f"  출력 파일: {output_csv}")
    print(f"{'='*80}\n")
    
    return len(errors) == 0


def load_config(config_file):
    """config.yaml 로드"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description='WDL 파이프라인 결과에서 3차 분석 입력 파일 자동 탐색',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  # config.yaml을 사용하여 자동 탐색
  python bin/prepare_inputs.py -c config/config.yaml

  # 수동으로 경로 지정
  python bin/prepare_inputs.py -b /data_4tb/hifi-human-wgs-wdl-custom/batch_results \\
      -s KTY9537,LDK6217 -g control -s 15_JANG_SJ,16_BYUN_HY -g experimental

출력:
  filelist.csv - 발견된 모든 파일의 경로 목록
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        type=str,
        default='config/config.yaml',
        help='설정 파일 경로 (기본값: config/config.yaml)'
    )
    
    parser.add_argument(
        '-b', '--batch-results',
        type=str,
        help='WDL batch_results 디렉토리 경로 (config 파일보다 우선)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='filelist.csv',
        help='출력 CSV 파일 경로 (기본값: filelist.csv)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='실제 파일을 생성하지 않고 탐색만 수행'
    )
    
    args = parser.parse_args()
    
    # 설정 로드
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        print(f"오류: 설정 파일을 찾을 수 없습니다: {args.config}")
        print("먼저 'cp config/config.example.yaml config/config.yaml'을 실행하세요.")
        sys.exit(1)
    except Exception as e:
        print(f"오류: 설정 파일 로드 실패: {e}")
        sys.exit(1)
    
    # batch_results 경로 결정
    if args.batch_results:
        batch_results_dir = args.batch_results
    else:
        # config에서 wdl_out_dir의 상위 디렉토리를 batch_results로 가정
        # 또는 config에 batch_results_dir 필드가 있으면 그것 사용
        if 'batch_results_dir' in config['paths']:
            batch_results_dir = config['paths']['batch_results_dir']
        else:
            print("오류: batch_results 디렉토리를 지정해야 합니다.")
            print("  방법 1: --batch-results 옵션 사용")
            print("  방법 2: config.yaml에 paths.batch_results_dir 추가")
            sys.exit(1)
    
    # 디렉토리 존재 확인
    if not os.path.exists(batch_results_dir):
        print(f"오류: batch_results 디렉토리를 찾을 수 없습니다: {batch_results_dir}")
        sys.exit(1)
    
    # 샘플 정보
    samples = config['samples']
    
    if args.dry_run:
        print("\n[DRY-RUN 모드] 파일을 생성하지 않습니다.\n")
        args.output = None
    
    # 파일 탐색 및 CSV 생성
    success = create_filelist_csv(batch_results_dir, samples, args.output)
    
    if not success:
        print("⚠️  일부 샘플에서 오류가 발생했습니다. 위의 오류 메시지를 확인하세요.")
        sys.exit(1)
    
    print("✓ 모든 샘플 처리 완료!\n")
    
    if not args.dry_run:
        print(f"다음 단계:")
        print(f"  1. {args.output} 파일을 확인하세요")
        print(f"  2. Snakefile을 실행하여 3차 분석을 시작하세요:")
        print(f"     ./bin/run_pipeline.sh --cores 16\n")


if __name__ == '__main__':
    main()
