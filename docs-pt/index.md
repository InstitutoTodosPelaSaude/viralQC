# Documentação ViralQC

```{admonition} 🌐 Language / Idioma
:class: tip

This documentation is also available in **English**: [Click here](https://viralqc.readthedocs.io/en/latest/)
```

**ViralQC** é uma ferramenta e pacote Python desenvolvido para identificação de vírus e controle de qualidade a partir de arquivos FASTA.

A ferramenta utiliza o [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/), o [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) e uma série de lógicas internas para classificar sequências virais e realizar controle de qualidade de genomas completos, regiões ou genes-alvo.

## Principais Funcionalidades

- **Identificação automática de vírus** usando Nextclade e BLAST
- **Controle de qualidade** de genomas virais usando Nextclade
- **Extração de regiões-alvo** (CDS ou genes específicos)
- **Análise de múltiplos vírus** em um único arquivo FASTA
- **Configuração flexível** através do arquivo `datasets.yml`

## Conteúdo da Documentação

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
```

## Links Rápidos

- [Repositório GitHub](https://github.com/InstitutoTodosPelaSaude/viralQC)
- [Issues/Bugs](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- **Licença:** MIT

## Referências

Ao utilizar o viralQC para fins acadêmicos, cite também:

- **Nextclade**: Aksamentov, I., Roemer, C., Hodcroft, E. B., & Neher, R. A., (2021). Nextclade: clade assignment, mutation calling and quality control for viral genomes. JOSS, 6(67), 3773.

- **BLAST**: Altschul SF, et al. (1990). Basic local alignment search tool. J Mol Biol. 215(3):403-10.
