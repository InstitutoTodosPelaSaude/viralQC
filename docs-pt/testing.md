# Testes

O `viralqc` usa [pytest](https://docs.pytest.org) para testes unitários. Todos os testes estão no diretório `tests/` na raiz do repositório.

## Configuração

Instale o pacote e as dependências de desenvolvimento no ambiente micromamba `viralQC`:

```bash
micromamba run -n viralQC pip install -e ".[dev]"
```

## Executando os testes

A partir da raiz do repositório:

```bash
# Todos os testes (verboso, tracebacks curtos — configurado em pyproject.toml)
micromamba run -n viralQC pytest

# Apenas um módulo
micromamba run -n viralQC pytest tests/test_ncbi_submission.py

# Uma classe ou função específica
micromamba run -n viralQC pytest tests/test_ncbi_submission.py::TestFilterSequences
micromamba run -n viralQC pytest tests/test_ncbi_submission.py::TestFilterSequences::test_high_n_content_dropped
```

## Estrutura dos testes

```
tests/
├── conftest.py              # Fixtures e utilitários compartilhados
├── test_ncbi_submission.py  # Testes de viralqc/core/ncbi_submission.py
├── test_run_analysis.py     # Testes de viralqc/core/run_analysis.py
└── test_scripts_python.py   # Testes de viralqc/scripts/python/
```

### `conftest.py`

Fornece helpers e fixtures pytest compartilhados:

| Símbolo | Descrição |
|---------|-----------|
| `make_fasta(tmp_path, name, records)` | Escreve um arquivo FASTA mínimo |
| `make_results_tsv(tmp_path, rows)` | Escreve um TSV de resultados do ViralQC |
| `make_metadata_csv(tmp_path, rows)` | Escreve um CSV de metadados |
| `tmp_fasta` (fixture) | FASTA com duas sequências que passam nos filtros de qualidade |
| `tmp_results` (fixture) | TSV de resultados SARS-CoV-2 com duas linhas |
| `tmp_metadata` (fixture) | CSV de metadados de vírus padrão com duas linhas |

### `test_ncbi_submission.py`

Cobre todas as funções públicas de `viralqc/core/ncbi_submission.py`:

| Classe de teste | Funções testadas |
|---|---|
| `TestSanitizeSeqName` | `sanitize_seq_name` |
| `TestNFraction` | `_n_fraction` |
| `TestLoadResults` | `load_results` (TSV / CSV / JSON / extensão inválida) |
| `TestLoadFasta` | `load_fasta` |
| `TestFilterSequences` | `filter_sequences` (conteúdo N, comprimento, limites personalizados) |
| `TestWriteFasta` | `write_fasta` (arquivo único, lotes, rename log) |
| `TestBatchPath` | `_batch_path` |
| `TestWriteHelpers` | `write_dropped_table`, `write_rename_log` |
| `TestCopyTbl` | `copy_tbl` |
| `TestIsPlainHeaderVirus` | `is_plain_header_virus` |
| `TestSerotypeFields` | `_serotype_fields` |
| `TestLoadMetadata` | `load_metadata` |
| `TestWriteSubmissionMetadata` | `write_submission_metadata` |

### `test_run_analysis.py`

Cobre `viralqc/core/run_analysis.RunAnalysis`. O Snakemake é **mockado** com `unittest.mock.patch`, portanto nenhuma ferramenta externa é necessária:

| Classe de teste | O que é testado |
|---|---|
| `TestGetOutputFormat` | Extensões válidas / inválidas |
| `TestRun` | Snakemake é chamado, chaves de config corretas, caminhos absolutos, formato inválido lança exceção |

### `test_scripts_python.py`

Testa funções auxiliares puras de cada script. Pontos de entrada com muita I/O (`main()`, `__main__`) não são testados aqui — eles dependem de testes de integração com dados reais.

| Classe de teste | Script / função |
|---|---|
| `TestValidateFasta` | `validate_fasta.validate_fasta_file` |
| `TestParseCdsCoverage` / `TestReorderCdsCoverage` / `TestReadGffGeneOrder` / `TestProcessNextcladeTsv` | `reorder_cds` |
| `TestSanitizeName` / `TestLoadIdMapping` / `TestParseTblBlocks` | `split_tbl_by_sample` |
| `TestCreateFastaPath` / `TestMapDatasetsToLocalPaths` / `TestWriteUnmappedSequences` | `format_nextclade_sort` |
| `TestParseFastaLengths` / `TestCleanCdsName` | `jsonl_to_gff` |

## Escrevendo novos testes

1. Coloque os arquivos de teste em `tests/` com o prefixo `test_`.
2. Importe utilitários de `tests/conftest.py` quando precisar de arquivos mockados.
3. Use `tmp_path` (fixture embutida do pytest) para diretórios temporários — o pytest os limpa automaticamente.
4. Prefira testes parametrizados (`@pytest.mark.parametrize`) para casos repetitivos.
5. Mocke chamadas externas (`run_snakemake`, requisições de rede, etc.) com `unittest.mock.patch`.

```python
# Exemplo: mockar uma dependência externa
from unittest.mock import patch, MagicMock

@patch("viralqc.core.run_analysis.run_snakemake")
def test_algo(mock_snake, tmp_path):
    mock_snake.return_value = MagicMock(success=True)
    ...
```
