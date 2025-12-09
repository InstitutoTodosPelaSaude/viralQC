# Lógica de Pontuação

O ViralQC implementa um **sistema próprio de pontuação** mais flexível que o Nextclade.

## Por Que Um Sistema Próprio?

Enquanto o Nextclade usa três categorias (good, mediocre, bad), o ViralQC oferece **quatro níveis** (A, B, C, D):

- **A** - Perfeito: Sequências que atendem todos os critérios
- **B** - Bom: Sequências adequadas para maioria das análises
- **C** - Aceitável: Usar com cautela
- **D** - Ruim: Uso não recomendado

## Métricas de Qualidade

O ViralQC calcula pontuações para **6 métricas**:

### 1. missingDataQuality - Qualidade de Dados Ausentes

Baseada na **cobertura do genoma**:

| Cobertura | Pontuação |
|-----------|-----------|
| ≥ 90% | **A** |
| ≥ 75% | **B** |
| ≥ 50% | **C** |
| < 50% | **D** |

### 2. privateMutationsQuality - Qualidade de Mutações Privadas

Baseada no **número de mutações privadas** em relação a um limite configurável:

| Mutações Privadas | Pontuação |
|-------------------|-----------|
| ≤ threshold | **A** |
| ≤ threshold × 1.05 | **B** |
| ≤ threshold × 1.10 | **C** |
| > threshold × 1.10 | **D** |

**Threshold padrão**: 10 mutações (configurável via `private_mutation_total_threshold`)

### 3. mixedSitesQuality - Qualidade de Sítios Mistos

| Sítios Mistos | Pontuação |
|---------------|-----------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 4. snpClustersQuality - Qualidade de Clusters de SNPs

| SNPs Clustered | Pontuação |
|----------------|-----------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 5. frameShiftsQuality - Qualidade de Frameshifts

| Frameshifts | Pontuação |
|-------------|-----------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

### 6. stopCodonsQuality - Qualidade de Códons de Parada

| Stop Codons | Pontuação |
|-------------|-----------|
| 0 | **A** |
| 1 | **B** |
| 2 | **C** |
| ≥ 3 | **D** |

## Cálculo do genomeQuality

A pontuação geral combina as 6 métricas:

1. **Conversão**: A=4, B=3, C=2, D=1 pontos
2. **Normalização**: `score = (total / 24) × 24`
3. **Classificação**:

| Score | genomeQuality |
|-------|---------------|
| = 24 | **A** |
| ≥ 18 | **B** |
| ≥ 12 | **C** |
| < 12 | **D** |

## Extração de Regiões-Alvo

O ViralQC extrai sequências baseado na hierarquia de qualidade:

```{mermaid}
graph TD
    A[Analisar Sequência] --> B{genomeQuality = A ou B?}
    B -->|Sim| C[Extrair GENOMA COMPLETO]
    B -->|Não| D{targetRegionsQuality = A ou B?}
    D -->|Sim| E[Extrair REGIÕES-ALVO]
    D -->|Não| F{targetGeneQuality = A ou B?}
    F -->|Sim| G[Extrair GENE-ALVO]
    F -->|Não| H[NÃO EXTRAIR]
```
