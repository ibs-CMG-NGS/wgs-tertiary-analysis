# Tertiary Analysis Parameters Guide

This guide provides detailed explanations of all configurable parameters in the WGS tertiary analysis pipeline, their impact on results, and recommendations for different analysis scenarios.

## Table of Contents

1. [Small Variant Filtering (slivar)](#1-small-variant-filtering-slivar)
2. [Structural Variant Processing (svpack)](#2-structural-variant-processing-svpack)
3. [Variant Annotation (VEP)](#3-variant-annotation-vep)
4. [Differential Methylation Analysis (DSS)](#4-differential-methylation-analysis-dss)
5. [Parameter Selection Guidelines](#5-parameter-selection-guidelines)
6. [Common Use Cases](#6-common-use-cases)

---

## 1. Small Variant Filtering (slivar)

### Overview
Slivar filters small variants (SNVs/indels) based on quality metrics, population frequency, and predicted impact. These parameters directly affect the sensitivity and specificity of variant calling.

### Parameters

#### `min_gq` (Default: 20)
```yaml
parameters:
  slivar:
    min_gq: 20
```

**Description**: Minimum Genotype Quality score for a variant to pass filtering.

**Range**: 0-99 (Phred-scaled probability that genotype is wrong)
- GQ 10 = 10% error probability
- GQ 20 = 1% error probability  
- GQ 30 = 0.1% error probability

**Impact on Results**:
- **Lower values (10-15)**: 
  - ✅ Higher sensitivity - captures more true variants
  - ❌ More false positives - includes lower confidence calls
  - 📊 Use for: Discovery studies, rare disease with high suspicion
  
- **Higher values (25-30)**:
  - ✅ Higher specificity - fewer false positives
  - ❌ Lower sensitivity - may miss real variants in difficult regions
  - 📊 Use for: Clinical reporting, high-confidence variant sets

**Recommendation**: 
- Research: GQ ≥ 15
- Clinical: GQ ≥ 20
- Ultra-high confidence: GQ ≥ 30

---

#### `max_gnomad_af` (Default: 0.01)
```yaml
parameters:
  slivar:
    max_gnomad_af: 0.01
```

**Description**: Maximum allele frequency in gnomAD database for a variant to be retained.

**Range**: 0.0-1.0 (proportion of population carrying the variant)
- 0.01 = 1% of population
- 0.001 = 0.1% of population
- 0.0001 = 0.01% of population

**Impact on Results**:
- **Higher values (0.05-0.10)**:
  - ✅ Retains common variants with potential effects
  - ❌ Large number of variants to review
  - 📊 Use for: Pharmacogenomics, complex trait analysis
  
- **Lower values (0.001-0.0001)**:
  - ✅ Focuses on rare/novel variants
  - ❌ May miss disease variants with higher frequency in specific populations
  - 📊 Use for: Mendelian disease, rare disease diagnosis

**Population Considerations**:
- Asian populations: Some pathogenic variants have AF > 0.01 in East Asian gnomAD
- Founder populations: Adjust based on disease allele frequency in specific groups

**Recommendation**:
- Rare disease: AF ≤ 0.001
- Research/discovery: AF ≤ 0.01
- Common variant studies: AF ≤ 0.05

---

#### `impact_severity` (Default: "HIGH,MODERATE")
```yaml
parameters:
  slivar:
    impact_severity: "HIGH,MODERATE"
```

**Description**: Variant Effect Predictor (VEP) consequence severity levels to retain.

**Available Levels**:
1. **HIGH**: Loss-of-function variants
   - Stop gained, frameshift, splice donor/acceptor
   - ~0.1% of variants per genome
   
2. **MODERATE**: Missense and other protein-altering
   - Missense, inframe indels, splice region
   - ~5-10% of variants per genome
   
3. **LOW**: Synonymous and other minor effects
   - Synonymous, stop retained
   - ~40-50% of variants per genome
   
4. **MODIFIER**: Non-coding and intergenic
   - Intron, upstream, downstream
   - ~45-50% of variants per genome

**Impact on Results**:
- **"HIGH" only**:
  - ✅ Smallest, most actionable variant set
  - ❌ Misses missense pathogenic variants (e.g., BRCA1 p.C61G)
  - 📊 Use for: Initial screening, high-penetrance disorders
  
- **"HIGH,MODERATE"** (recommended):
  - ✅ Balanced sensitivity/specificity
  - ✅ Captures most clinically relevant coding variants
  - 📊 Use for: Most clinical/research applications
  
- **"HIGH,MODERATE,LOW"**:
  - ✅ Includes synonymous variants (may affect splicing)
  - ❌ Significantly increases variant count
  - 📊 Use for: RNA splicing disorders, comprehensive analysis

**Recommendation**: "HIGH,MODERATE" for 95% of use cases

---

## 2. Structural Variant Processing (svpack)

### Overview
Svpack filters and processes structural variants (SVs: deletions, duplications, inversions, translocations). SV calling is inherently noisier than SNV calling, requiring careful parameter tuning.

### Parameters

#### `min_sv_size` (Default: 50)
```yaml
parameters:
  svpack:
    min_sv_size: 50
```

**Description**: Minimum size (in base pairs) for structural variants to be retained.

**Range**: 30-10,000 bp

**Impact on Results**:
- **Smaller values (30-50 bp)**:
  - ✅ Captures small deletions/duplications in coding regions
  - ❌ Higher false positive rate (alignment artifacts)
  - 📊 Use for: Exon-level CNV detection, repeat expansions
  
- **Larger values (100-1000 bp)**:
  - ✅ Higher confidence calls
  - ❌ Misses small pathogenic SVs (e.g., BRCA1 exon deletions)
  - 📊 Use for: Large CNV analysis, low-coverage data

**Biology Context**:
- Exon sizes: typically 100-300 bp → use ≤50 bp to detect exon deletions
- Gene sizes: 1-100 kb → use ≥50 bp for gene-level CNVs
- Repeat elements: 100-6000 bp → adjust based on target repeat

**Recommendation**:
- High-coverage PacBio (≥30x): 50 bp
- Low-coverage data (15-20x): 100 bp
- Large SV focus: 1000 bp

---

#### `sv_types` (Default: "DEL,DUP,INV,BND")
```yaml
parameters:
  svpack:
    sv_types: "DEL,DUP,INV,BND"
```

**Description**: Types of structural variants to include in analysis.

**Available Types**:
1. **DEL** (Deletion): Loss of DNA segment
   - Haploinsufficiency, loss-of-function
   - ~10,000-20,000 per genome
   
2. **DUP** (Duplication): Gain of DNA segment
   - Gene dosage increase, fusion genes
   - ~5,000-10,000 per genome
   
3. **INV** (Inversion): DNA segment flipped orientation
   - Can disrupt genes at breakpoints
   - ~500-1,000 per genome
   
4. **BND** (Breakend/Translocation): Inter-chromosomal rearrangement
   - Fusion genes (e.g., BCR-ABL)
   - ~50-200 per genome
   
5. **INS** (Insertion): Addition of DNA sequence
   - Mobile element insertions
   - ~1,000-5,000 per genome (often filtered)

**Impact on Results**:
- **All types** (recommended):
  - ✅ Comprehensive SV detection
  - ❌ More variants to review
  
- **DEL,DUP only**:
  - ✅ Focus on dosage-sensitive disorders
  - ❌ Misses fusion genes from translocations
  - 📊 Use for: CNV-focused analysis
  
- **BND only**:
  - ✅ Cancer gene fusion detection
  - 📊 Use for: Oncology applications

**Recommendation**: Include all types, filter downstream by gene overlap

---

#### `max_sv_size` (Default: 1000000)
```yaml
parameters:
  svpack:
    max_sv_size: 1000000
```

**Description**: Maximum size (in base pairs) for structural variants to be retained.

**Range**: 10,000-5,000,000 bp (1 Mb = 1,000,000 bp)

**Impact on Results**:
- **Smaller values (100 kb)**:
  - ✅ Excludes whole-chromosome aneuploidies
  - ❌ Misses large pathogenic CNVs (e.g., Williams syndrome)
  - 📊 Use for: Gene-level SV analysis
  
- **Larger values (1-5 Mb)**:
  - ✅ Detects microdeletion/microduplication syndromes
  - ✅ Identifies chromosomal abnormalities
  - 📊 Use for: Constitutional disorder diagnosis

**Common Syndromes**:
- DiGeorge syndrome (22q11.2): 1.5-3 Mb deletion → need max ≥ 3 Mb
- Williams syndrome: 1.5 Mb deletion → need max ≥ 2 Mb
- Prader-Willi/Angelman: 5-7 Mb deletion → need max ≥ 10 Mb

**Recommendation**: 
- Constitutional analysis: 5,000,000 bp (5 Mb)
- Somatic/cancer: 1,000,000 bp (1 Mb)

---

## 3. Variant Annotation (VEP)

### Overview
VEP annotates variants with functional predictions, population frequencies, and clinical databases. These parameters control which annotation sources are used and how predictions are made.

### Parameters

#### `assembly` (Default: "GRCh38")
```yaml
parameters:
  vep:
    assembly: "GRCh38"
```

**Description**: Reference genome assembly version for annotation.

**Options**: "GRCh38" or "GRCh37"

**Impact on Results**:
- **GRCh38** (recommended):
  - ✅ Latest genome build with corrections
  - ✅ Better representation of alternate loci and patch scaffolds
  - ❌ Some older databases not fully updated
  
- **GRCh37/hg19**:
  - ✅ More legacy clinical databases available
  - ❌ Known errors in genome sequence
  - 📊 Use only if required for compatibility

**Important**: Assembly MUST match your VCF input files. Mismatches cause incorrect annotations.

**Recommendation**: Always use GRCh38 for new projects

---

#### `cache_version` (Default: 110)
```yaml
parameters:
  vep:
    cache_version: 110
```

**Description**: Ensembl VEP cache version for transcript and gene annotations.

**Range**: 95-112 (as of 2024, updates quarterly)

**Impact on Results**:
- **Newer versions**:
  - ✅ Updated gene models and transcripts
  - ✅ New clinical annotations
  - ❌ May change transcript selection (MANE transcripts updated)
  
- **Older versions**:
  - ✅ Consistency with historical analyses
  - ❌ Missing recent gene discoveries

**Recommendation**: Use latest version (110+) for new analyses, but document version for reproducibility

---

#### `plugins` (Default: "CADD,REVEL,dbNSFP")
```yaml
parameters:
  vep:
    plugins: "CADD,REVEL,dbNSFP"
```

**Description**: VEP plugins for additional variant annotation and pathogenicity prediction.

**Available Plugins**:

##### **CADD** (Combined Annotation Dependent Depletion)
- **Function**: Integrative deleteriousness score
- **Range**: 0-99 (higher = more deleterious)
- **Interpretation**:
  - CADD ≥ 10: Top 10% most deleterious
  - CADD ≥ 20: Top 1% most deleterious
  - CADD ≥ 30: Top 0.1% most deleterious
- **Use**: General variant prioritization, non-missense variants

##### **REVEL** (Rare Exome Variant Ensemble Learner)
- **Function**: Missense variant pathogenicity predictor
- **Range**: 0-1 (higher = more pathogenic)
- **Interpretation**:
  - REVEL ≥ 0.5: Likely pathogenic (sensitivity 85%)
  - REVEL ≥ 0.75: High confidence pathogenic
- **Use**: Missense variant interpretation in disease genes

##### **dbNSFP** (Database of Non-Synonymous Functional Predictions)
- **Function**: Meta-database with multiple prediction algorithms
- **Includes**: 
  - SIFT, PolyPhen-2, MutationTaster
  - Conservation scores (phyloP, phastCons)
  - Functional annotations
- **Use**: Comprehensive missense variant assessment

**Impact on Results**:
- **All plugins** (recommended):
  - ✅ Multiple independent pathogenicity predictions
  - ✅ Consensus scoring improves accuracy
  - ❌ Longer runtime (adds ~30% processing time)
  
- **CADD only**:
  - ✅ Fast, single score
  - ❌ May miss variants with discordant predictions
  
- **No plugins**:
  - ✅ Fastest VEP runtime
  - ❌ No pathogenicity predictions, manual interpretation required

**Recommendation**: Use all three plugins for clinical analysis, CADD only for large-scale screening

---

#### `pick_order` (Default: "mane_select,canonical,biotype,length")
```yaml
parameters:
  vep:
    pick_order: "mane_select,canonical,biotype,length"
```

**Description**: Priority order for selecting one transcript per gene when multiple exist.

**Available Criteria**:
1. **mane_select**: MANE Select transcript (gold standard clinical transcript)
2. **canonical**: Ensembl canonical transcript (longest CDS)
3. **biotype**: Protein coding > other biotypes
4. **length**: Longest transcript
5. **tsl**: Transcript Support Level (experimental evidence)

**Impact on Results**:
- **mane_select first** (recommended):
  - ✅ Consistent with clinical reporting guidelines (ACMG/ClinGen)
  - ✅ Matches HGVS nomenclature in ClinVar
  - 📊 Use for: Clinical variant interpretation
  
- **canonical first**:
  - ✅ Stable over time
  - ❌ May differ from clinical databases
  - 📊 Use for: Research, consistency with older studies

**Recommendation**: Always prioritize MANE Select for clinical applications

---

## 4. Differential Methylation Analysis (DSS)

### Overview
DSS identifies differentially methylated regions (DMRs) between control and experimental groups. These parameters control the statistical stringency and biological significance thresholds.

### Parameters

#### `p_threshold` (Default: 0.001)
```yaml
parameters:
  dmr:
    p_threshold: 0.001
```

**Description**: Statistical significance threshold (p-value) for calling DMRs.

**Range**: 0.0001-0.05

**Impact on Results**:
- **Lower p-values (0.0001-0.001)**:
  - ✅ High confidence DMRs, low false discovery rate
  - ❌ May miss subtle but real methylation changes
  - 📊 Use for: Clinical biomarkers, follow-up validation
  
- **Higher p-values (0.01-0.05)**:
  - ✅ Higher sensitivity for exploratory analysis
  - ❌ More false positives requiring validation
  - 📊 Use for: Discovery studies, large sample sizes

**Multiple Testing Context**:
- Genome has ~28 million CpGs
- DMR analysis tests thousands of regions
- Recommend p < 0.001 to control family-wise error rate

**Recommendation**:
- Discovery: p ≤ 0.01
- Clinical: p ≤ 0.001
- High-confidence: p ≤ 0.0001

---

#### `delta_threshold` (Default: 0.1)
```yaml
parameters:
  dmr:
    delta_threshold: 0.1
```

**Description**: Minimum methylation difference (beta value) between groups to call a DMR.

**Range**: 0.05-0.5 (proportion of methylation)
- 0.1 = 10% methylation difference
- 0.2 = 20% methylation difference
- 0.5 = 50% methylation difference

**Impact on Results**:
- **Lower thresholds (0.05-0.1)**:
  - ✅ Detects subtle methylation changes
  - ❌ May include biologically insignificant changes
  - 📊 Use for: Tissue differentiation, aging studies
  
- **Higher thresholds (0.2-0.5)**:
  - ✅ Focuses on major methylation alterations
  - ❌ Misses subtle but cumulative effects
  - 📊 Use for: Cancer vs normal, imprinting disorders

**Biological Context**:
- Promoter methylation: 20-30% change often affects transcription
- Enhancer methylation: 10-15% change can alter activity
- Imprinted loci: 50% difference (monoallelic vs biallelic)

**Recommendation**:
- General analysis: Δ ≥ 0.1 (10%)
- Cancer studies: Δ ≥ 0.2 (20%)
- Imprinting/X-inactivation: Δ ≥ 0.3 (30%)

---

#### `min_length` (Default: 50)
```yaml
parameters:
  dmr:
    min_length: 50
```

**Description**: Minimum length (in base pairs) for a region to be called a DMR.

**Range**: 25-500 bp

**Impact on Results**:
- **Shorter regions (25-50 bp)**:
  - ✅ Detects small regulatory elements (TF binding sites)
  - ❌ More susceptible to technical noise
  - 📊 Use for: High-coverage data (≥30x)
  
- **Longer regions (100-500 bp)**:
  - ✅ Higher confidence, spans multiple CpGs
  - ❌ May miss focal methylation changes
  - 📊 Use for: Low-coverage data, enhancer regions

**CpG Density Context**:
- CpG islands: ~1 CpG per 10 bp → 50 bp = ~5 CpGs
- CpG shores: ~1 CpG per 30 bp → 50 bp = ~2 CpGs
- Recommend ≥3 CpGs per DMR for reliability

**Recommendation**: 
- High-coverage (≥30x): 50 bp
- Standard coverage (20-30x): 100 bp
- Low-coverage (≤20x): 200 bp

---

#### `min_cpg` (Default: 3)
```yaml
parameters:
  dmr:
    min_cpg: 3
```

**Description**: Minimum number of CpG sites in a region to be called a DMR.

**Range**: 2-10 CpGs

**Impact on Results**:
- **Fewer CpGs (2-3)**:
  - ✅ More DMRs called, higher sensitivity
  - ❌ Single-site noise can drive false positives
  - 📊 Use for: Sparse CpG regions (gene bodies)
  
- **More CpGs (5-10)**:
  - ✅ Robust to technical variation
  - ❌ Biased toward CpG islands, misses shores/shelves
  - 📊 Use for: CpG island-focused analysis

**Statistical Power**:
- 3 CpGs: Minimum for statistical testing
- 5 CpGs: Good balance of sensitivity/specificity
- 10 CpGs: High confidence but reduced coverage

**Recommendation**: 3-5 CpGs for most applications

---

## 5. Parameter Selection Guidelines

### Decision Tree

```
START: What is your analysis goal?
│
├─ Clinical Diagnosis
│  ├─ Rare Disease
│  │  ├─ slivar: min_gq=20, max_gnomad_af=0.001, impact="HIGH,MODERATE"
│  │  ├─ svpack: min_sv_size=50, all types, max_sv_size=5000000
│  │  └─ vep: CADD+REVEL+dbNSFP, mane_select first
│  │
│  └─ Cancer/Somatic
│     ├─ slivar: min_gq=25, max_gnomad_af=0.0001, impact="HIGH,MODERATE"
│     ├─ svpack: min_sv_size=100, BND+DEL+DUP, max_sv_size=1000000
│     └─ dmr: p=0.001, delta=0.2, min_cpg=5
│
├─ Research/Discovery
│  ├─ Exploratory
│  │  ├─ slivar: min_gq=15, max_gnomad_af=0.01, impact="HIGH,MODERATE,LOW"
│  │  ├─ svpack: min_sv_size=50, all types
│  │  └─ dmr: p=0.01, delta=0.1, min_cpg=3
│  │
│  └─ Validation/Follow-up
│     ├─ slivar: min_gq=20, max_gnomad_af=0.001, impact="HIGH,MODERATE"
│     └─ dmr: p=0.001, delta=0.1, min_cpg=5
│
└─ Method Development
   ├─ High Sensitivity
   │  └─ slivar: min_gq=10, max_gnomad_af=0.05, all impacts
   │
   └─ High Specificity
      └─ slivar: min_gq=30, max_gnomad_af=0.0001, impact="HIGH"
```

### Sample Size Considerations

| Sample Size | DMR p-threshold | Delta Threshold | Min CpG |
|-------------|-----------------|-----------------|---------|
| 3-5 per group | 0.01 | 0.2 | 5 |
| 6-10 per group | 0.005 | 0.15 | 3-5 |
| 11-20 per group | 0.001 | 0.1 | 3 |
| 20+ per group | 0.0001 | 0.05-0.1 | 3 |

### Coverage-Based Recommendations

| Coverage | min_gq | min_sv_size | min_length | Min CpG |
|----------|--------|-------------|------------|---------|
| 15-20x | 15 | 100 bp | 200 bp | 5 |
| 20-30x | 20 | 50 bp | 100 bp | 3-5 |
| 30-40x | 20 | 50 bp | 50 bp | 3 |
| 40x+ | 25 | 50 bp | 50 bp | 3 |

---

## 6. Common Use Cases

### Use Case 1: Rare Mendelian Disease Diagnosis

**Clinical Scenario**: 5-year-old patient with intellectual disability, dysmorphic features, and developmental delay. Exome sequencing negative.

**Recommended Parameters**:
```yaml
parameters:
  slivar:
    min_gq: 20
    max_gnomad_af: 0.001  # Rare variants only
    impact_severity: "HIGH,MODERATE"
  
  svpack:
    min_sv_size: 50  # Detect exon-level deletions
    sv_types: "DEL,DUP,INV,BND"
    max_sv_size: 5000000  # Include microdeletion syndromes
  
  vep:
    assembly: "GRCh38"
    cache_version: 110
    plugins: "CADD,REVEL,dbNSFP"  # Multiple pathogenicity scores
    pick_order: "mane_select,canonical,biotype,length"
```

**Rationale**: 
- Low AF threshold captures rare pathogenic variants
- Comprehensive SV detection for CNV syndromes
- Clinical-grade variant annotation with MANE transcripts

**Expected Output**: 
- 50-200 high/moderate impact rare variants
- 10-50 rare SVs overlapping genes
- Focus review on known disease genes

---

### Use Case 2: Cancer Genomics (Tumor vs Normal)

**Clinical Scenario**: Paired tumor-normal whole genome sequencing for precision oncology.

**Recommended Parameters**:
```yaml
parameters:
  slivar:
    min_gq: 25  # Higher stringency for somatic variants
    max_gnomad_af: 0.0001  # Exclude common germline polymorphisms
    impact_severity: "HIGH,MODERATE"
  
  svpack:
    min_sv_size: 100
    sv_types: "BND,DEL,DUP"  # Focus on gene fusions and CNVs
    max_sv_size: 1000000
  
  dmr:
    p_threshold: 0.001
    delta_threshold: 0.2  # Cancer shows large methylation changes
    min_length: 100
    min_cpg: 5
```

**Rationale**:
- Higher GQ threshold reduces false positives in tumor samples
- BND detection critical for fusion gene discovery
- Large delta threshold for cancer-specific hypermethylation

**Expected Output**:
- 5-50 somatic driver mutations
- 2-20 gene fusions/large SVs
- 100-500 DMRs in promoters/enhancers

---

### Use Case 3: Population Genomics Study

**Research Scenario**: Identify genetic variants associated with complex trait in 100 individuals.

**Recommended Parameters**:
```yaml
parameters:
  slivar:
    min_gq: 15  # More lenient for discovery
    max_gnomad_af: 0.05  # Include common variants
    impact_severity: "HIGH,MODERATE,LOW"
  
  svpack:
    min_sv_size: 50
    sv_types: "DEL,DUP,INV,BND"
    max_sv_size: 1000000
  
  dmr:
    p_threshold: 0.01  # Discovery threshold
    delta_threshold: 0.1
    min_length: 50
    min_cpg: 3
```

**Rationale**:
- Higher AF threshold captures common trait-associated variants
- Broad impact severity for comprehensive variant catalog
- Relaxed DMR thresholds for exploratory analysis

**Expected Output**:
- 500-2000 variants per sample for association testing
- 100-300 SVs per sample
- 1000-5000 DMRs for epigenome-wide association

---

### Use Case 4: Epigenetic Biomarker Discovery

**Research Scenario**: Identify DNA methylation signatures distinguishing disease subtypes (n=10 per group).

**Recommended Parameters**:
```yaml
parameters:
  dmr:
    p_threshold: 0.005  # Balance sensitivity and FDR
    delta_threshold: 0.15  # Moderate effect size
    min_length: 100  # Span multiple CpGs
    min_cpg: 5  # Robust to technical variation
```

**Rationale**:
- Moderate p-value with medium sample size
- 15% methylation difference has biological relevance
- Minimum 5 CpGs for reliable biomarker

**Expected Output**:
- 200-1000 candidate DMRs
- Prioritize promoter/enhancer DMRs
- Validate top 20-50 DMRs in independent cohort

---

### Use Case 5: Pharmacogenomics Panel

**Clinical Scenario**: Pre-emptive screening for drug metabolism variants.

**Recommended Parameters**:
```yaml
parameters:
  slivar:
    min_gq: 20
    max_gnomad_af: 0.10  # Include common PGx alleles
    impact_severity: "HIGH,MODERATE,LOW"  # Include synonymous PGx variants
  
  vep:
    assembly: "GRCh38"
    plugins: "CADD,REVEL,dbNSFP"
    pick_order: "mane_select,canonical,biotype,length"
```

**Rationale**:
- Higher AF threshold captures common metabolizer alleles (e.g., CYP2D6*4 AF=0.2)
- Include synonymous variants (some affect splicing)
- Focus on known pharmacogenes (CYP, UGT, SLCO families)

**Expected Output**:
- 10-30 actionable PGx variants
- Genotype-to-phenotype translation (e.g., *1/*4 = intermediate metabolizer)
- Clinical decision support for drug dosing

---

## Parameter Interaction Effects

### Interaction 1: GQ threshold × Coverage
Higher coverage allows higher GQ thresholds without losing sensitivity:
- 15x coverage: min_gq ≤ 15
- 30x coverage: min_gq = 20 (optimal)
- 60x coverage: min_gq = 25 (very high confidence)

### Interaction 2: AF threshold × Population ancestry
Adjust gnomAD AF based on patient population:
- European ancestry: Use global AF
- East Asian: Check gnomAD EAS subset (some pathogenic variants AF > 0.01)
- African: Use AFR subset (highest genetic diversity)
- Admixed: Use most conservative (lowest) AF across populations

### Interaction 3: DMR p-value × delta threshold
Stricter p-value allows relaxed delta, and vice versa:
- High stringency: p ≤ 0.0001, delta ≥ 0.05
- Balanced: p ≤ 0.001, delta ≥ 0.1
- Discovery: p ≤ 0.01, delta ≥ 0.2

### Interaction 4: SV size × SV type
Different SV types have different size distributions:
- DEL: 50 bp - 5 Mb (most common: 100-500 bp)
- DUP: 100 bp - 10 Mb (most common: 1-50 kb)
- INV: 1 kb - 10 Mb (most common: 10-100 kb)
- BND: N/A (breakpoint, not size-dependent)

Adjust min/max_sv_size based on sv_types selected.

---

## Troubleshooting

### Problem: Too many variants (>10,000 small variants)
**Solutions**:
1. Increase `min_gq` to 25-30
2. Decrease `max_gnomad_af` to 0.001-0.0001
3. Use only `impact_severity: "HIGH"`
4. Filter post-hoc by gene panels (e.g., ACMG genes)

### Problem: Too few variants (<10 small variants)
**Solutions**:
1. Decrease `min_gq` to 15
2. Increase `max_gnomad_af` to 0.01-0.05
3. Add "LOW" to impact_severity
4. Check VCF input quality (may have pre-filtering)

### Problem: No DMRs detected
**Solutions**:
1. Increase `p_threshold` to 0.01-0.05
2. Decrease `delta_threshold` to 0.05-0.1
3. Decrease `min_cpg` to 2-3
4. Check sample sizes (n ≥ 3 per group recommended)
5. Verify methylation data quality (coverage ≥ 10x per CpG)

### Problem: Too many DMRs (>5,000)
**Solutions**:
1. Decrease `p_threshold` to 0.0001
2. Increase `delta_threshold` to 0.2-0.3
3. Increase `min_cpg` to 5-10
4. Check for batch effects or outlier samples

### Problem: SV calling misses known pathogenic CNV
**Solutions**:
1. Decrease `min_sv_size` to 30-50 bp (for small deletions)
2. Increase `max_sv_size` to 5-10 Mb (for large syndromes)
3. Verify SV type is in `sv_types` (check if INV or BND)
4. Check input VCF: some WDL pipelines pre-filter SVs

---

## References and Further Reading

1. **Slivar Documentation**: https://github.com/brentp/slivar
2. **gnomAD Allele Frequency Guidance**: https://gnomad.broadinstitute.org/help/faf
3. **VEP Consequence Terms**: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
4. **ACMG/AMP Variant Interpretation Guidelines**: PMID 25741868
5. **MANE Transcripts**: https://www.ncbi.nlm.nih.gov/refseq/MANE/
6. **DSS User Guide**: https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html
7. **CADD Score Interpretation**: https://cadd.gs.washington.edu/info
8. **Structural Variant Standards**: PMID 31690885

---

## Appendix: Quick Reference Tables

### Table A1: Parameter Presets by Analysis Type

| Analysis Type | Preset Name | Key Parameters |
|---------------|-------------|----------------|
| Clinical Rare Disease | `strict_rare` | GQ≥20, AF≤0.001, HIGH+MODERATE, SV 50bp-5Mb |
| Clinical Cancer | `somatic_strict` | GQ≥25, AF≤0.0001, DMR delta≥0.2 |
| Research Discovery | `discovery` | GQ≥15, AF≤0.01, DMR p≤0.01 |
| Pharmacogenomics | `pgx` | GQ≥20, AF≤0.10, all impacts |
| Population Study | `population` | GQ≥15, AF≤0.05, all SV types |

### Table A2: Expected Runtime by Parameter Stringency

| Stringency | Small Variants | SVs | VEP | DMR | Total |
|------------|----------------|-----|-----|-----|-------|
| Lenient (GQ≥10, AF≤0.05) | ~500K variants | ~50K SVs | 4-6h | 2-4h | 6-10h |
| Moderate (GQ≥20, AF≤0.01) | ~50K variants | ~10K SVs | 1-2h | 1-2h | 2-4h |
| Strict (GQ≥30, AF≤0.001) | ~5K variants | ~1K SVs | 15-30m | 30m-1h | 1-2h |

*Runtime estimates for 30x WGS on 16 cores*

### Table A3: Disk Space Requirements

| Analysis Component | Intermediate Files | Final Output | Total |
|-------------------|-------------------|--------------|-------|
| Small Variants | 5-10 GB | 100-500 MB | ~10 GB |
| Structural Variants | 2-5 GB | 50-200 MB | ~5 GB |
| VEP Annotation | 10-20 GB | 500 MB-2 GB | ~20 GB |
| Methylation DMR | 5-10 GB | 200 MB-1 GB | ~10 GB |
| **Total Pipeline** | **25-45 GB** | **1-4 GB** | **~50 GB** |

*Per-sample estimates. Multiply by sample count for cohort studies.*

---

## Version History

- **v1.0** (2026-01-29): Initial comprehensive parameter guide
  - All tertiary analysis parameters documented
  - Use case examples added
  - Troubleshooting section included

---

**Document Maintainer**: WGS Tertiary Analysis Pipeline  
**Last Updated**: 2026-01-29  
**For Questions**: See main [README.md](../README.md) for contact information
