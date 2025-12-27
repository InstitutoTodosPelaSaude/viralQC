# ViralQC Documentation

```{admonition} 🌐 Language / Idioma
:class: tip

Esta documentação também está disponível em **Português**: [Clique aqui](https://viralqc.readthedocs.io/pt-br/latest/)
```

**ViralQC** is a Python tool and package developed for virus identification and quality control from FASTA files.

The tool uses [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/), [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), and a series of internal logics to classify viral sequences and perform quality control of complete genomes, regions, or target genes.

The tool was developed with the goal of automating the use of Nextclade and BLAST tools, focusing on integration with viral genomic databases, in order to reduce submission errors (such as submitting genomes with incorrect taxonomic information) and also to provide the user viral genome quality metrics of submitted sequences, combining Nextclade metrics provided by standardized datasets, as well as providing general metrics and annotations based on reference genomes (refseq) provided by NCBI.

## Main Features

- **Automatic virus identification** using Nextclade and BLAST
- **Quality control** of viral genomes using Nextclade
- **Target region extraction** (CDS or specific genes)
- **Analysis of multiple viruses** in a single FASTA file
- **Flexible configuration** through the `datasets.yml` file

## Documentation Contents

```{toctree}
:maxdepth: 2

installation
configuration
scoring
commands
adding-datasets
output
examples
troubleshooting
developer-guide
```

## Quick Links

- [GitHub Repository](https://github.com/InstitutoTodosPelaSaude/viralQC)
- [Issues/Bugs](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- **License:** MIT

## References

When using viralQC for academic purposes, also cite:

- **Nextclade**: Aksamentov, I., Roemer, C., Hodcroft, E. B., & Neher, R. A., (2021). Nextclade: clade assignment, mutation calling and quality control for viral genomes. JOSS, 6(67), 3773.

- **BLAST**: Altschul SF, et al. (1990). Basic local alignment search tool. J Mol Biol. 215(3):403-10.
