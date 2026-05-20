#!/usr/bin/env python3
"""
TRGT VCF에서 비정상적으로 긴 반복서열 탐지

TRGT(Tandem Repeat Genotyper)가 생성한 phased VCF를 분석하여
카탈로그 기대 범위를 초과한 repeat를 'expanded'로 플래그한다.

입력 FORMAT 필드:
  AL   - 실측 allele length (bp)
  ALLR - allele length range (카탈로그 기대 범위, 예: "10-25,10-25")
  SD   - spanning read depth
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


def parse_allr(allr_str):
    """
    ALLR 문자열에서 각 allele의 (min, max) 기대 범위를 파싱한다.
    예: "10-25,10-25"  → [(10, 25), (10, 25)]
         "."           → None (정보 없음)
    """
    if not allr_str or allr_str == '.':
        return None
    results = []
    for part in allr_str.split(','):
        m = re.match(r'(\d+)-(\d+)', part.strip())
        if m:
            results.append((int(m.group(1)), int(m.group(2))))
        else:
            results.append(None)
    return results if results else None


def parse_al(al_str):
    """
    AL 문자열에서 allele length 목록을 파싱한다.
    예: "15,18" → [15, 18]
         "."    → []
    """
    if not al_str or al_str == '.':
        return []
    return [int(x) for x in al_str.split(',') if x.strip().isdigit()]


def classify_allele(al_len, expected_range, expansion_buffer):
    """
    단일 allele가 정상(normal) / 확장(expanded)인지 판정한다.
    expected_range가 None이면 'unknown' 반환.
    """
    if expected_range is None:
        return 'unknown'
    lo, hi = expected_range
    if al_len > hi + expansion_buffer:
        return 'expanded'
    return 'normal'


def query_trgt_vcf(vcf_path, min_sd):
    """
    bcftools query로 TRGT VCF에서 필요한 필드를 추출한다.
    PASS 필터 + SD >= min_sd 조건 적용.

    반환: list of dict
    """
    cmd = [
        'bcftools', 'view', '-f', 'PASS', str(vcf_path),
        '|',
        'bcftools', 'query',
        '-f', r'%CHROM\t%POS\t%INFO/TRID\t%INFO/MOTIFS\t%INFO/STRUC\t[%AL\t%ALLR\t%SD]\n',
    ]
    # bcftools view | bcftools query 파이프라인은 subprocess shell=True로 실행
    cmd_str = (
        f"bcftools view -f PASS {vcf_path} "
        f"| bcftools query "
        f"-f '%CHROM\\t%POS\\t%INFO/TRID\\t%INFO/MOTIFS\\t%INFO/STRUC"
        f"\\t[%AL\\t%ALLR\\t%SD]\\n'"
    )
    result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"bcftools 오류:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    records = []
    for line in result.stdout.splitlines():
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        chrom, pos, trid, motifs, struc = parts[0], parts[1], parts[2], parts[3], parts[4]
        al_str, allr_str, sd_str = parts[5], parts[6], parts[7]

        try:
            sd = int(sd_str)
        except ValueError:
            sd = 0

        if sd < min_sd:
            continue

        records.append({
            'chrom': chrom,
            'pos': pos,
            'trid': trid,
            'motifs': motifs,
            'struc': struc,
            'al_str': al_str,
            'allr_str': allr_str,
            'sd': sd,
        })

    return records


def analyze_records(records, expansion_buffer):
    """
    각 레코드에 대해 allele별 분류를 수행하고 결과 dict 목록을 반환한다.
    """
    all_rows = []
    for r in records:
        al_list   = parse_al(r['al_str'])
        allr_list = parse_allr(r['allr_str'])

        allele1_len = al_list[0] if len(al_list) > 0 else None
        allele2_len = al_list[1] if len(al_list) > 1 else None

        range1 = allr_list[0] if allr_list and len(allr_list) > 0 else None
        range2 = allr_list[1] if allr_list and len(allr_list) > 1 else None

        status1 = classify_allele(allele1_len, range1, expansion_buffer) if allele1_len is not None else 'unknown'
        status2 = classify_allele(allele2_len, range2, expansion_buffer) if allele2_len is not None else 'unknown'

        # 어느 allele 하나라도 expanded면 expanded 판정
        if 'expanded' in (status1, status2):
            overall = 'expanded'
        elif 'unknown' in (status1, status2):
            overall = 'unknown'
        else:
            overall = 'normal'

        expected_str = ''
        if range1:
            expected_str += f"{range1[0]}-{range1[1]}"
        if range2:
            expected_str += f",{range2[0]}-{range2[1]}"

        all_rows.append({
            'chrom':         r['chrom'],
            'pos':           r['pos'],
            'trid':          r['trid'],
            'motif':         r['motifs'],
            'allele1_len':   allele1_len if allele1_len is not None else '',
            'allele2_len':   allele2_len if allele2_len is not None else '',
            'expected_range': expected_str,
            'spanning_depth': r['sd'],
            'status':        overall,
        })

    return all_rows


def write_tsv(rows, filepath, filter_expanded=False):
    fieldnames = ['chrom', 'pos', 'trid', 'motif',
                  'allele1_len', 'allele2_len', 'expected_range',
                  'spanning_depth', 'status']
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            if filter_expanded and row['status'] != 'expanded':
                continue
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description='TRGT VCF에서 반복서열 이상값(expanded) 탐지'
    )
    parser.add_argument('--input',  required=True, help='phased TRGT VCF.gz 경로')
    parser.add_argument('--output-outliers', required=True,
                        help='이상값(expanded)만 출력할 TSV 경로')
    parser.add_argument('--output-summary', required=True,
                        help='전체 repeat 요약 TSV 경로 (PASS + SD 필터 통과)')
    parser.add_argument('--min-sd', type=int, default=5,
                        help='최소 spanning read depth (기본값: 5)')
    parser.add_argument('--expansion-buffer', type=int, default=0,
                        help='기대 범위 초과 허용 bp (기본값: 0)')
    args = parser.parse_args()

    vcf_path = Path(args.input)
    if not vcf_path.exists():
        print(f"오류: 파일을 찾을 수 없습니다: {vcf_path}", file=sys.stderr)
        sys.exit(1)

    Path(args.output_outliers).parent.mkdir(parents=True, exist_ok=True)

    print(f"TRGT VCF 분석 시작: {vcf_path}", file=sys.stderr)
    print(f"  최소 spanning depth: {args.min_sd}", file=sys.stderr)
    print(f"  expansion buffer: {args.expansion_buffer} bp", file=sys.stderr)

    records = query_trgt_vcf(vcf_path, args.min_sd)
    print(f"  PASS + SD≥{args.min_sd} 통과 레코드: {len(records)}", file=sys.stderr)

    all_rows = analyze_records(records, args.expansion_buffer)

    expanded_rows = [r for r in all_rows if r['status'] == 'expanded']
    print(f"  Expanded repeat 수: {len(expanded_rows)}", file=sys.stderr)

    write_tsv(expanded_rows, args.output_outliers, filter_expanded=False)
    write_tsv(all_rows,      args.output_summary,  filter_expanded=False)

    print(f"완료: outliers → {args.output_outliers}", file=sys.stderr)
    print(f"      summary  → {args.output_summary}",  file=sys.stderr)


if __name__ == '__main__':
    main()
