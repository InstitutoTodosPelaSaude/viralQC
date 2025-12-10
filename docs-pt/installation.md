# Instalação

## Pré-requisitos

O ViralQC requer Python 3.8 ou superior (mas inferior a 3.12) e várias dependências de bioinformática.

## Instalação via pip

### Passo 1: Instalar Dependências

Instale as dependências usando `micromamba` ou `conda`:

```bash
micromamba install \
  -c conda-forge \
  -c bioconda \
  "python>=3.8.0,<3.12.0" \
  "snakemake-minimal>=7.32.0,<7.33.0" \
  "blast>=2.16.0,<2.17.0" \
  "nextclade>=3.15.0,<3.16.0" \
  "seqtk>=1.5.0,<1.6.0" \
  "ncbi-datasets-cli>=18.9.0,<18.10.0" \
  "taxonkit>=0.20.0,<0.21.0"
```

**Descrição das Dependências:**

| Pacote | Versão | Descrição |
|--------|--------|-----------|
| Python | ≥3.8.0, <3.12.0 | Linguagem de programação base |
| Snakemake | ≥7.32.0, <7.33.0 | Gerenciador de workflows |
| BLAST | ≥2.16.0, <2.17.0 | Ferramenta de busca de similaridade |
| Nextclade | ≥3.15.0, <3.16.0 | Análise de clados e controle de qualidade |
| Seqtk | ≥1.5.0, <1.6.0 | Utilitário para arquivos FASTA |
| NCBI Datasets CLI | ≥18.9.0, <18.10.0 | Download de dados do NCBI |
| TaxonKit | ≥0.20.0, <0.21.0 | Manipulação de taxonomia |

### Passo 2: Instalar ViralQC

```bash
pip install viralQC
```

### Passo 3: Verificar Instalação

```bash
vqc --help
```

## Instalação a partir do código-fonte

### Passo 1: Clonar o Repositório

```bash
git clone https://github.com/InstitutoTodosPelaSaude/viralQC.git
cd viralQC
```

### Passo 2: Criar Ambiente Conda

```bash
micromamba env create -f env.yml
micromamba activate viralQC
```

### Passo 3: Instalar ViralQC

```bash
pip install -e .
```

### Passo 4: Verificar Instalação

```bash
vqc --help
```
