#!/bin/bash
# ================================================================================
# PacBio HiFi WGS 3차 분석 파이프라인 - 빠른 시작 스크립트
# ================================================================================

set -e  # 오류 발생 시 중단

# 기본 설정 파일
CONFIG_FILE="${SNAKEMAKE_CONFIG:-config/config.yaml}"

# 도움말
show_help() {
    cat << EOF
PacBio HiFi WGS 3차 분석 파이프라인 - 설정 확인 스크립트

사용법:
    $0 [--config FILE]

옵션:
    --config FILE    설정 파일 경로 (기본: config/config.yaml)
    -h, --help       이 도움말 표시

환경 변수:
    SNAKEMAKE_CONFIG    설정 파일 경로 (기본: config/config.yaml)

예제:
    # 기본 설정 파일 검증
    $0

    # 커스텀 설정 파일 검증
    $0 --config config/config_hifisolve.yaml

    # 환경 변수로 설정
    SNAKEMAKE_CONFIG=config/config_hifisolve.yaml $0
EOF
}

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "알 수 없는 옵션: $1"
            show_help
            exit 1
            ;;
    esac
done

echo "================================================================================"
echo "PacBio HiFi WGS 3차 분석 파이프라인 - 설정 확인"
echo "================================================================================"
echo "설정 파일: $CONFIG_FILE"
echo ""

# 색상 코드
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 함수: 명령어 존재 확인
check_command() {
    if command -v $1 &> /dev/null; then
        echo -e "${GREEN}✓${NC} $1 발견됨"
        return 0
    else
        echo -e "${RED}✗${NC} $1 없음"
        return 1
    fi
}

# 함수: 파일 존재 확인
check_file() {
    if [ -f "$1" ]; then
        echo -e "${GREEN}✓${NC} $1"
        return 0
    else
        echo -e "${RED}✗${NC} $1 (없음)"
        return 1
    fi
}

echo ""
echo "1. 필수 소프트웨어 확인..."
echo "-------------------------------------------"

ALL_OK=true

check_command snakemake  || ALL_OK=false
check_command python     || ALL_OK=false
check_command bcftools   || ALL_OK=false
check_command bgzip      || ALL_OK=false
check_command tabix      || ALL_OK=false

# bedtools (SV consequence / ASM에 필요)
check_command bedtools || echo -e "${YELLOW}⚠${NC} bedtools 없음 (conda install -c bioconda bedtools 권장)"

# Docker 또는 Singularity 중 하나는 필요
if command -v docker &> /dev/null; then
    echo -e "${GREEN}✓${NC} docker 발견됨"
elif command -v singularity &> /dev/null; then
    echo -e "${GREEN}✓${NC} singularity 발견됨"
else
    echo -e "${RED}✗${NC} docker 또는 singularity 필요"
    ALL_OK=false
fi

echo ""
echo "2. 파이프라인 파일 확인..."
echo "-------------------------------------------"

check_file "$CONFIG_FILE"             || ALL_OK=false
check_file "Snakefile"                || ALL_OK=false
check_file "scripts/dmr_analysis.R"  || ALL_OK=false
check_file "scripts/asm_analysis.R"  || ALL_OK=false
check_file "scripts/trgt_outlier.py" || ALL_OK=false
check_file "scripts/svpack"          || ALL_OK=false

echo ""
echo "3. config.yaml 검증..."
echo "-------------------------------------------"

# Python을 사용한 YAML 파싱 검증
python3 << PYEOF
import yaml
import sys

CONFIG_FILE = "${CONFIG_FILE}"

try:
    with open(CONFIG_FILE, 'r') as f:
        config = yaml.safe_load(f)
    
    # 필수 키 확인
    required_keys = ['paths', 'samples', 'parameters']
    for key in required_keys:
        if key not in config:
            print(f"\033[0;31m✗\033[0m config.yaml에 '{key}' 섹션 누락")
            sys.exit(1)
    
    # paths 확인 (batch_results_dir 또는 wdl_out_dir 중 하나 필요)
    if 'batch_results_dir' not in config['paths'] and 'wdl_out_dir' not in config['paths']:
        print("\033[0;31m✗\033[0m paths.batch_results_dir 또는 paths.wdl_out_dir 누락")
        sys.exit(1)

    # gff3_file 확인 (SV consequence 주석에 필요)
    import os
    gff3 = config['paths'].get('gff3_file', '')
    if gff3 and os.path.exists(gff3):
        print(f"\033[0;32m✓\033[0m gff3_file 발견됨: {os.path.basename(gff3)}")
    elif gff3:
        print(f"\033[1;33m⚠\033[0m gff3_file 경로가 존재하지 않음: {gff3}")
        print(f"  SV consequence 주석 (annotate_sv_consequence) 실행 불가")
    else:
        print(f"\033[1;33m⚠\033[0m paths.gff3_file 미설정 — SV consequence 주석 비활성화")

    # samples 확인
    if not config['samples'].get('control') or not config['samples'].get('experimental'):
        print("\033[0;31m✗\033[0m 샘플 그룹 (control/experimental) 누락")
        sys.exit(1)

    print("\033[0;32m✓\033[0m config.yaml 유효성 검사 통과")

    # 샘플 정보 출력
    n_control = len(config['samples']['control'])
    n_exp = len(config['samples']['experimental'])
    print(f"  - Control 샘플: {n_control}개")
    print(f"  - Experimental 샘플: {n_exp}개")
    wdl_path = config['paths'].get('batch_results_dir') or config['paths'].get('wdl_out_dir', '(미설정)')
    print(f"  - WDL 결과 경로: {wdl_path}")

    # 선택적 분석 모듈 안내
    print("\n  [선택적 분석 모듈 - filelist.csv에 해당 컬럼이 있을 때 자동 실행]")
    print("  - TRGT 반복서열 분석  : phased_trgt_vcf 컬럼 필요")
    print("  - ASM 메틸화 분석     : cpg_hap1_bed, cpg_hap2_bed 컬럼 필요")
    
except yaml.YAMLError as e:
    print(f"\033[0;31m✗\033[0m {CONFIG_FILE} 파싱 오류: {e}")
    sys.exit(1)
except FileNotFoundError:
    print(f"\033[0;31m✗\033[0m {CONFIG_FILE} 파일을 찾을 수 없습니다")
    sys.exit(1)
PYEOF

if [ $? -ne 0 ]; then
    ALL_OK=false
fi

echo ""
echo "4. 입력 파일 확인..."
echo "-------------------------------------------"

# config 파일에서 wdl_out_dir 추출
WDL_OUT_DIR=$(python3 -c "import yaml; config=yaml.safe_load(open('${CONFIG_FILE}')); print(config['paths']['wdl_out_dir'])")

if [ -d "$WDL_OUT_DIR" ]; then
    echo -e "${GREEN}✓${NC} WDL 출력 디렉토리: $WDL_OUT_DIR"
    
    # 샘플 파일 확인
    SAMPLE_FILES=$(ls -1 "$WDL_OUT_DIR"/phased_*.vcf.gz 2>/dev/null | wc -l)
    if [ $SAMPLE_FILES -gt 0 ]; then
        echo -e "${GREEN}✓${NC} phased VCF 파일: ${SAMPLE_FILES}개 발견"
    else
        echo -e "${YELLOW}⚠${NC} phased VCF 파일을 찾을 수 없습니다"
        echo "  경로: $WDL_OUT_DIR/phased_*.vcf.gz"
    fi
else
    echo -e "${RED}✗${NC} WDL 출력 디렉토리를 찾을 수 없음: $WDL_OUT_DIR"
    ALL_OK=false
fi

echo ""
echo "================================================================================"

if [ "$ALL_OK" = true ]; then
    echo -e "${GREEN}모든 확인 통과!${NC}"
    echo ""
    echo "다음 명령어로 파이프라인을 실행할 수 있습니다:"
    echo ""
    echo "  # Dry-run (테스트)"
    echo "  snakemake --dry-run --cores 1"
    echo ""
    echo "  # 실제 실행 (Docker 사용)"
    echo "  snakemake --use-singularity --cores 16"
    echo ""
    echo "  # 클러스터 실행 (SLURM)"
    echo "  snakemake --use-singularity --cluster 'sbatch -p normal -n {threads}' --jobs 10"
    echo ""
    exit 0
else
    echo -e "${RED}일부 확인 실패!${NC}"
    echo ""
    echo "위의 오류를 수정한 후 다시 시도하세요."
    echo "자세한 내용은 README.md를 참조하세요."
    echo ""
    exit 1
fi
