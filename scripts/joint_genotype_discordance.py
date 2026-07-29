#!/usr/bin/env python3
"""웰(그룹) 내 3개 시퀀싱배치 replicate 간 joint genotype 불일치율 계산.

split_joint_small_variant_vcfs/{0,1,2}/ 의 개별 샘플 VCF 3개를 bcftools merge로 다시 합쳐
같은 사이트에 대해 3개 replicate의 GT가 일치하는지 비교한다.
(각 split VCF는 소스가 다른 필터링을 거쳐 row 수가 다를 수 있어, merge로 union set을 만들고
세 샘플 모두 유효한 GT를 가진 사이트만 비교 대상으로 삼는다.)

불일치율 = (3개 중 GT가 다른 사이트 수) / (3개 모두 유효 GT이고 그 중 최소 1개는 non-ref인 사이트 수)
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile


def normalize_gt(gt):
    """위상(phase) 무시하고 정렬된 allele 튜플로 정규화. 미콜('.')이 있으면 None."""
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or "." in alleles:
        return None
    return tuple(sorted(alleles))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joint-dir", required=True, help="joint-H2O2-XXX/<timestamp>_joint 디렉토리")
    ap.add_argument("--well", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    vcf_glob = os.path.join(args.joint_dir, "out", "split_joint_small_variant_vcfs", "*", "*.vcf.gz")
    vcfs = sorted(glob.glob(vcf_glob))
    if len(vcfs) != 3:
        print(f"경고: {args.well} joint VCF가 3개가 아님 ({len(vcfs)}개 발견): {vcfs}", file=sys.stderr)

    for v in vcfs:
        subprocess.run(["bcftools", "index", "-t", "-f", v], check=True,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tmp:
        merged = tmp.name
    subprocess.run(["bcftools", "merge"] + vcfs + ["-O", "z", "-o", merged], check=True,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    out = subprocess.run(["bcftools", "query", "-f", "[%GT\\t]\\n", merged],
                          capture_output=True, text=True, check=True)
    os.remove(merged)
    for ext in (".tbi",):
        for v in vcfs:
            idx = v + ext
            if os.path.exists(idx):
                pass  # 인덱스는 원본 VCF 옆에 그대로 둠 (재사용 가능)

    n_compared = 0
    n_discordant = 0
    for line in out.stdout.splitlines():
        gts = [g for g in line.strip().split("\t") if g]
        if len(gts) != 3:
            continue
        norm = [normalize_gt(g) for g in gts]
        if any(n is None for n in norm):
            continue  # 3개 중 하나라도 미콜이면 비교 제외
        if all(n == ("0", "0") for n in norm):
            continue  # 3개 다 REF/REF면 비교 대상 아님 (변이 없음)
        n_compared += 1
        if len(set(norm)) > 1:
            n_discordant += 1

    rate = (n_discordant / n_compared) if n_compared > 0 else float("nan")

    with open(args.output, "w") as f:
        f.write("well\tgroup\tn_sites_compared\tn_discordant\tdiscordance_rate\n")
        f.write(f"{args.well}\t{args.group}\t{n_compared}\t{n_discordant}\t{rate:.6f}\n")

    print(f"저장: {args.output} (compared={n_compared}, discordant={n_discordant}, rate={rate:.4f})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
