# 📁 프로젝트 구조

```
wgs-tertiary-pipeline/
│
├── 📄 Snakefile                    # 메인 워크플로우 정의
├── 📄 README.md                    # 상세 사용 설명서
├── 📄 .gitignore                   # Git 제외 파일 목록
├── 🐍 environment.yaml             # Conda 환경 정의
│
├── 📁 config/                      # 설정 파일 디렉토리
│   ├── config.yaml                 # 실제 설정 (사용자 수정) [gitignore]
│   ├── config.example.yaml         # 설정 템플릿 (버전 관리)
│   └── cluster_config.yaml         # SLURM 리소스 설정
│
├── 📁 bin/                         # 실행 스크립트 디렉토리
│   ├── run_pipeline.sh             # 메인 실행 스크립트
│   ├── setup_check.sh              # 환경 검증 스크립트
│   └── run_slurm.sh                # SLURM 클러스터 실행 스크립트
│
├── 📁 scripts/                     # 분석 스크립트 디렉토리
│   └── dmr_analysis.R              # DMR 분석 R 스크립트 (DSS)
│
└── 📁 docs/                        # 문서 디렉토리
    ├── QUICKSTART.md               # 빠른 시작 가이드
    └── SUMMARY.md                  # 프로젝트 전체 요약
```

---

## 🗂️ 디렉토리 설명

### 📁 config/
설정 파일을 모아둔 디렉토리입니다.

- **config.yaml**: 실제 사용하는 설정 파일 (개인 정보 포함, `.gitignore`에 추가됨)
- **config.example.yaml**: 설정 템플릿 (버전 관리됨, 사용자는 이를 복사하여 `config.yaml` 생성)
- **cluster_config.yaml**: SLURM 클러스터 환경에서 각 규칙별 리소스 할당 정의

**사용 방법:**
```bash
cp config/config.example.yaml config/config.yaml
nano config/config.yaml  # 경로 및 샘플 정보 수정
```

### 📁 bin/
실행 가능한 스크립트들을 모아둔 디렉토리입니다.

- **run_pipeline.sh**: 로컬/워크스테이션에서 파이프라인 실행
- **setup_check.sh**: 환경 및 설정 검증 (실행 전 확인용)
- **run_slurm.sh**: SLURM 클러스터에서 파이프라인 실행

**사용 방법:**
```bash
# 환경 확인
bash bin/setup_check.sh

# 로컬 실행
./bin/run_pipeline.sh --cores 16

# 클러스터 실행
sbatch bin/run_slurm.sh
```

### 📁 scripts/
분석에 사용되는 스크립트들을 모아둔 디렉토리입니다.

- **dmr_analysis.R**: DSS 패키지를 사용한 차이 메틸화 영역(DMR) 분석 R 스크립트

이 스크립트는 Snakemake가 자동으로 호출하므로 직접 실행할 필요는 없습니다.

### 📁 docs/
문서 파일들을 모아둔 디렉토리입니다.

- **QUICKSTART.md**: 5분 안에 시작하는 가이드 (초보자용)
- **SUMMARY.md**: 프로젝트 전체 기술 요약 (개발자/연구자용)

**추천 읽기 순서:**
1. 처음 사용자 → `QUICKSTART.md`
2. 상세 정보 필요 → `README.md`
3. 전체 아키텍처 이해 → `docs/SUMMARY.md`

---

## 🔗 파일 간 참조 관계

### Snakefile
```python
configfile: "config/config.yaml"  # 설정 파일 로드
```

### bin/run_slurm.sh
```bash
--cluster-config config/cluster_config.yaml  # 클러스터 리소스 설정
```

### bin/setup_check.sh
```bash
check_file "config/config.yaml"  # 설정 파일 존재 확인
```

---

## 🚀 빠른 시작

### 1단계: 설정 파일 준비
```bash
cp config/config.example.yaml config/config.yaml
nano config/config.yaml
```

### 2단계: 환경 확인
```bash
bash bin/setup_check.sh
```

### 3단계: 실행
```bash
./bin/run_pipeline.sh --cores 16
```

---

## 📝 버전 관리 전략

### Git에 포함되는 파일
- ✅ `config/config.example.yaml` (템플릿)
- ✅ `config/cluster_config.yaml` (클러스터 설정)
- ✅ `bin/*.sh` (실행 스크립트)
- ✅ `scripts/*.R` (분석 스크립트)
- ✅ `docs/*.md` (문서)
- ✅ `Snakefile` (워크플로우)
- ✅ `README.md` (설명서)

### Git에서 제외되는 파일 (`.gitignore`)
- ❌ `config/config.yaml` (개인 설정)
- ❌ `3rd_analysis_results/` (결과 파일)
- ❌ `logs/` (로그 파일)
- ❌ `.snakemake/` (메타데이터)
- ❌ `*.vcf.gz`, `*.bam`, `*.bed` (데이터 파일)

---

## 🔄 권장 워크플로우

```
1. git clone <repository>
   └─> 프로젝트 클론
   
2. cp config/config.example.yaml config/config.yaml
   └─> 설정 파일 생성
   
3. nano config/config.yaml
   └─> 경로 및 샘플 정보 수정
   
4. bash bin/setup_check.sh
   └─> 환경 검증
   
5. ./bin/run_pipeline.sh --dry-run
   └─> Dry-run 테스트
   
6. ./bin/run_pipeline.sh --cores 16
   └─> 실제 실행
   
7. Check results in 3rd_analysis_results/
   └─> 결과 확인
```

---

## 📚 추가 문서

- **메인 가이드**: [README.md](README.md)
- **빠른 시작**: [docs/QUICKSTART.md](docs/QUICKSTART.md)
- **기술 요약**: [docs/SUMMARY.md](docs/SUMMARY.md)

---

**업데이트**: 2026-01-29
