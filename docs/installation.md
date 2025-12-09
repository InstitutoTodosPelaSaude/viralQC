# Installation

## Prerequisites

ViralQC requires Python 3.8 or higher (but lower than 3.12) and several bioinformatics dependencies.

## Installation via pip

### Step 1: Install Dependencies

First, install the required dependencies using `micromamba` or `conda`:

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

**Dependency Descriptions:**

| Package | Version | Description |
|---------|---------|-------------|
| Python | ≥3.8.0, <3.12.0 | Base programming language |
| Snakemake | ≥7.32.0, <7.33.0 | Workflow manager for analysis orchestration |
| BLAST | ≥2.16.0, <2.17.0 | Sequence similarity search tool |
| Nextclade | ≥3.15.0, <3.16.0 | Clade analysis and quality control tool |
| Seqtk | ≥1.5.0, <1.6.0 | Utility for processing FASTA files |
| NCBI Datasets CLI | ≥18.9.0, <18.10.0 | Tool for downloading data from NCBI |
| TaxonKit | ≥0.20.0, <0.21.0 | Utility for NCBI taxonomy manipulation |

### Step 2: Install ViralQC

```bash
pip install viralQC
```

### Step 3: Verify Installation

```bash
vqc --help
```

You should see the ViralQC help message listing the available commands.

## Installation from Source Code

### Step 1: Clone the Repository

```bash
git clone https://github.com/InstitutoTodosPelaSaude/viralQC.git
cd viralQC
```

### Step 2: Create Conda Environment

```bash
micromamba env create -f env.yml
micromamba activate viralQC
```

### Step 3: Install ViralQC

```bash
pip install -e .
```

### Step 4: Verify Installation

```bash
vqc --help
```
