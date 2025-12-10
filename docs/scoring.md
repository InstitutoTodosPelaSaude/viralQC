# Scoring Logic

ViralQC implements its **own quality scoring system** that is more flexible than Nextclade's standard system.

## Why a Custom System?

While Nextclade uses three categories (good, mediocre, bad) based on fixed scores, ViralQC offers **four levels** (A, B, C, D):

- **A** - Sequences meeting all quality criteria
- **B** - Sequences suitable for most analyses
- **C** - Use with caution
- **D** - Use not recommended

## Quality Metrics

ViralQC calculates scores for **6 metrics**:

### 1. missingDataQuality - Missing Data Quality

Based on **genome coverage**:

| Coverage | Score |
|----------|-------|
| ≥ 90% | **A** |
| ≥ 75% | **B** |
| ≥ 50% | **C** |
| < 50% | **D** |

### 2. privateMutationsQuality - Private Mutations Quality

Based on the **number of private mutations** relative to a configurable threshold:

| Private Mutations | Score |
|-------------------|-------|
| ≤ threshold | **A** |
| ≤ threshold × 1.05 | **B** |
| ≤ threshold × 1.10 | **C** |
| > threshold × 1.10 | **D** |

**Default threshold**: 10 mutations (configurable per virus via `private_mutation_total_threshold`)

### 3. mixedSitesQuality - Mixed Sites Quality

| Mixed Sites | Score |
|-------------|-------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 4. snpClustersQuality - SNP Clusters Quality

| SNPs Clustered | Score |
|----------------|-------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 5. frameShiftsQuality - Frameshifts Quality

| Frameshifts | Score |
|-------------|-------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 6. stopCodonsQuality - Stop Codons Quality

| Stop Codons | Score |
|-------------|-------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

## genomeQuality Calculation

The overall genome score combines all 6 metrics:

1. **Conversion**: A=4, B=3, C=2, D=1 points
2. **Normalization**: `score = (total / 24) × 24`
3. **Classification**:

| Score | genomeQuality |
|-------|---------------|
| = 24 | **A** |
| ≥ 18 | **B** |
| ≥ 12 | **C** |
| < 12 | **D** |

## Target Region Extraction

ViralQC extracts sequences based on quality hierarchy:

```{mermaid}
graph TD
    A[Analyze Sequence] --> B{genomeQuality = A or B?}
    B -->|Yes| C[Extract COMPLETE GENOME]
    B -->|No| D{targetRegionsQuality = A or B?}
    D -->|Yes| E[Extract TARGET REGIONS]
    D -->|No| F{targetGeneQuality = A or B?}
    F -->|Yes| G[Extract TARGET GENE]
    F -->|No| H[DO NOT EXTRACT]
```
