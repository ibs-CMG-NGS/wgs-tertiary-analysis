#!/usr/bin/env python3
"""MGI GO annotation(GAF)에서 산화스트레스/DNA손상복구 관련 유전자셋 추출.

주의(단순화): GO term의 상위/하위 계층(is_a, part_of)을 전개하지 않고, 아래 명시된
GO ID에 직접 주석된 유전자만 모은다. 즉 더 구체적인 하위 term에만 annotated된 유전자는
누락될 수 있음 — 완전한 GO closure가 아니라 실용적 근사임을 문서에 명시.
"""
import argparse
import gzip
import json
from collections import defaultdict

# 기존 2세트 (Snakefile 기본 rule이 참조 — 동작 불변 유지)
GENESETS = {
    "oxidative_stress_response": [
        "GO:0006979",  # response to oxidative stress
        "GO:0034599",  # cellular response to oxidative stress
        "GO:0000302",  # response to reactive oxygen species
        "GO:0034614",  # cellular response to reactive oxygen species
    ],
    "dna_repair": [
        "GO:0006281",  # DNA repair
        "GO:0006974",  # DNA damage response
    ],
}

# 큐레이션 확대 세트 (H2O2 성상세포 겨냥, 다카피 아티팩트 GO 의도적 제외)
# 계획: H2O2_curated_geneset_burden_plan_2026-07-27.md
CURATED_GENESETS = {
    "oxidative_stress_response": ["GO:0006979", "GO:0034599", "GO:0000302", "GO:0034614"],
    "glutathione_redox":         ["GO:0006749", "GO:0045454", "GO:0016209"],
    "dna_repair_core":           ["GO:0006281", "GO:0006974"],
    "dna_repair_ber":            ["GO:0006284"],
    "dna_repair_ner":            ["GO:0006289"],
    "dna_repair_mmr":            ["GO:0006298"],
    "dna_repair_dsb_hr":         ["GO:0000724"],
    "dna_repair_dsb_nhej":       ["GO:0006303"],
    "dna_damage_checkpoint":     ["GO:0031570", "GO:0000077"],
    "apoptosis":                 ["GO:0006915", "GO:0043065"],
    "intrinsic_apoptosis":       ["GO:0097193", "GO:0008630"],
    "ferroptosis":               ["GO:0097707"],
    "autophagy":                 ["GO:0006914", "GO:0016236"],
    "cellular_senescence":       ["GO:0090398", "GO:0007050"],
    "cell_cycle":                ["GO:0007049"],
    "er_stress_upr":             ["GO:0034976", "GO:0030968", "GO:0006986"],
    "mitochondrion":             ["GO:0006119", "GO:0007005"],
    "inflammatory_response":     ["GO:0006954"],
    "gliogenesis_astrocyte":     ["GO:0042063", "GO:0014002", "GO:0048708"],
    "hypoxia_response":          ["GO:0001666", "GO:0036293"],
    "calcium_signaling":         ["GO:0006816", "GO:0070588"],
}

PRESETS = {"original": GENESETS, "curated": CURATED_GENESETS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaf", required=True, help="mgi.gaf.gz 경로")
    ap.add_argument("--output", required=True, help="출력 JSON 경로")
    ap.add_argument("--preset", choices=list(PRESETS), default="original",
                    help="세트 컬렉션 선택 (기본 original=기존 2세트, curated=확대 세트)")
    args = ap.parse_args()

    genesets = PRESETS[args.preset]

    go_to_set = {}
    for set_name, go_ids in genesets.items():
        for go_id in go_ids:
            go_to_set[go_id] = set_name

    genes_by_set = defaultdict(set)
    with gzip.open(args.gaf, "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            symbol, go_id = parts[2], parts[4]
            set_name = go_to_set.get(go_id)
            if set_name:
                genes_by_set[set_name].add(symbol)

    out = {name: sorted(genes) for name, genes in genes_by_set.items()}
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    for name, genes in out.items():
        print(f"{name}: {len(genes)}개 유전자")
    print(f"저장: {args.output}")


if __name__ == "__main__":
    main()
