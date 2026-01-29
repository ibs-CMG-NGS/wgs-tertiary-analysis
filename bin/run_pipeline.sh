#!/bin/bash
# ================================================================================
# PacBio HiFi WGS 3차 분석 파이프라인 - 실행 스크립트
# ================================================================================

set -e

# 기본 설정
CORES=${CORES:-16}
USE_SINGULARITY=${USE_SINGULARITY:-true}
CONFIG_FILE=${SNAKEMAKE_CONFIG:-"config/config.yaml"}

# 색상 코드
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# 도움말
show_help() {
    cat << EOF
PacBio HiFi WGS 3차 분석 파이프라인 실행 스크립트

사용법:
    $0 [옵션]

옵션:
    -h, --help              이 도움말 표시
    -d, --dry-run           Dry-run 모드 (실제 실행 안 함)
    -c, --cores N           사용할 코어 수 (기본: 16)
    --config FILE           설정 파일 경로 (기본: config/config.yaml)
    -t, --target RULE       특정 규칙만 실행
    --unlock                워크플로우 잠금 해제
    --dag                   DAG 그래프 생성 (dag.pdf)
    --report                HTML 리포트 생성

예제:
    # Dry-run
    $0 --dry-run

    # 커스텀 설정 파일 사용
    $0 --config config/config_hifisolve.yaml --cores 8

    # 전체 파이프라인 실행 (8 코어)
    $0 --cores 8

    # Small variant 필터링만 실행
    $0 --target filter_small_variants

    # DAG 시각화
    $0 --dag

환경 변수:
    CORES                   사용할 코어 수 (기본: 16)
    USE_SINGULARITY         Singularity 사용 여부 (기본: true)
    SNAKEMAKE_CONFIG        설정 파일 경로 (기본: config/config.yaml)

EOF
}

# 인자 파싱
DRY_RUN=false
TARGET=""
UNLOCK=false
DAG=false
REPORT=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            export SNAKEMAKE_CONFIG="$2"
            shift 2
            ;;
        -t|--target)
            TARGET="$2"
            shift 2
            ;;
        --unlock)
            UNLOCK=true
            shift
            ;;
        --dag)
            DAG=true
            shift
            ;;
        --report)
            REPORT=true
            shift
            ;;
        *)
            echo "알 수 없는 옵션: $1"
            show_help
            exit 1
            ;;
    esac
done

# 설정 파일 확인
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${YELLOW}경고: 설정 파일을 찾을 수 없습니다: $CONFIG_FILE${NC}"
    echo "config/config.example.yaml을 복사하여 설정 파일을 만드세요:"
    echo "  cp config/config.example.yaml $CONFIG_FILE"
    exit 1
fi

# 로그 디렉토리 생성
mkdir -p logs

echo -e "${BLUE}================================================================================${NC}"
echo -e "${BLUE}PacBio HiFi WGS 3차 분석 파이프라인${NC}"
echo -e "${BLUE}================================================================================${NC}"
echo -e "설정 파일: ${GREEN}$CONFIG_FILE${NC}"
echo ""

# 잠금 해제
if [ "$UNLOCK" = true ]; then
    echo -e "${YELLOW}워크플로우 잠금 해제 중...${NC}"
    snakemake --unlock
    echo -e "${GREEN}완료!${NC}"
    exit 0
fi

# DAG 생성
if [ "$DAG" = true ]; then
    echo -e "${YELLOW}DAG 그래프 생성 중...${NC}"
    snakemake --dag | dot -Tpdf > dag.pdf
    echo -e "${GREEN}dag.pdf 생성 완료!${NC}"
    
    snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
    echo -e "${GREEN}rulegraph.pdf 생성 완료!${NC}"
    exit 0
fi

# Snakemake 명령어 구성
SNAKEMAKE_CMD="snakemake"

# Singularity 사용
if [ "$USE_SINGULARITY" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --use-singularity"
fi

# 코어 수
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cores $CORES"

# Dry-run
if [ "$DRY_RUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run --printshellcmds"
    echo -e "${YELLOW}Dry-run 모드${NC}"
fi

# 특정 타겟
if [ -n "$TARGET" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD $TARGET"
    echo -e "${YELLOW}타겟: $TARGET${NC}"
fi

# 추가 옵션
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --printshellcmds --keep-going --rerun-incomplete"

# 실행
echo ""
echo -e "${GREEN}실행 명령어:${NC}"
echo "$SNAKEMAKE_CMD"
echo ""

if [ "$DRY_RUN" = false ]; then
    echo -e "${YELLOW}파이프라인 실행 시작...${NC}"
    echo ""
fi

# Snakemake 실행
$SNAKEMAKE_CMD

# 리포트 생성
if [ "$REPORT" = true ] && [ "$DRY_RUN" = false ]; then
    echo ""
    echo -e "${YELLOW}HTML 리포트 생성 중...${NC}"
    snakemake --report report.html
    echo -e "${GREEN}report.html 생성 완료!${NC}"
fi

if [ "$DRY_RUN" = false ]; then
    echo ""
    echo -e "${BLUE}================================================================================${NC}"
    echo -e "${GREEN}파이프라인 실행 완료!${NC}"
    echo -e "${BLUE}================================================================================${NC}"
    echo ""
    echo "결과 파일 위치:"
    echo "  - Small Variants: 3rd_analysis_results/small_variants/"
    echo "  - Structural Variants: 3rd_analysis_results/structural_variants/"
    echo "  - DMR 분석: 3rd_analysis_results/dmr_analysis/"
    echo "  - 메틸화 데이터: 3rd_analysis_results/methylation/"
    echo ""
fi
