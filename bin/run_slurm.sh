#!/bin/bash
# ================================================================================
# SLURM 클러스터 환경에서 파이프라인 실행
# ================================================================================

#SBATCH --job-name=pacbio-tertiary
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal

set -e

# 로그 디렉토리 생성
mkdir -p logs

echo "================================================================================"
echo "PacBio HiFi WGS 3차 분석 파이프라인 - SLURM 실행"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Hostname: $(hostname)"
echo "Start time: $(date)"
echo "================================================================================"

# 환경 모듈 로드 (필요시 수정)
# module load python/3.10
# module load singularity

# Conda 환경 활성화 (옵션)
# source activate pacbio-tertiary

# Snakemake 실행 (SLURM 클러스터 모드)
snakemake \
    --use-singularity \
    --jobs 50 \
    --cluster "sbatch \
        --partition={cluster.partition} \
        --cpus-per-task={threads} \
        --mem={resources.mem_mb}M \
        --time={cluster.time} \
        --output=logs/slurm_{rule}_{wildcards}.out \
        --error=logs/slurm_{rule}_{wildcards}.err" \
    --cluster-config config/cluster_config.yaml \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

echo ""
echo "================================================================================"
echo "파이프라인 실행 완료"
echo "End time: $(date)"
echo "================================================================================"
