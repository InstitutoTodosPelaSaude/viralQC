# Guia do Desenvolvedor

Este guia fornece documentação detalhada para os scripts Python usados no pipeline do ViralQC. Destina-se a desenvolvedores que desejam entender a lógica interna, escolhas de implementação e o fluxo de execução da ferramenta.

Os scripts estão localizados em `viralqc/scripts/python/` e são orquestrados por workflows do Snakemake.

## Gerenciamento de datasets (Datasets)

Estes scripts são usados pelos comandos `get-nextclade-datasets` e `get-blast-database` para baixar e preparar dados de referência.

### get_github_dataset.py

**Propósito**: Baixa diretórios específicos de datasets de um repositório GitHub sem usar a API do GitHub (para evitar limites de taxa).

**Execução**: Chamado por `get_public_datasets.smk` ao processar vírus configurados como `github`.

**Detalhes de Implementação**:
- **Sem Uso de API**: Em vez de usar a API do GitHub, baixa o arquivo compactado do repositório via `https://codeload.github.com/.../zip/refs/heads/main`.
- **Extração Seletiva**: Faz o stream do arquivo zip e extrai apenas os arquivos que correspondem ao `dataset-path` solicitado.
- **Achatamento de Estrutura**: Lida com a remoção da pasta raiz (ex: `repo-main/`) para colocar os arquivos diretamente no diretório de destino.

### jsonl_to_gff.py

**Propósito**: Converte relatórios de anotação JSONL do NCBI Datasets para o formato GFF3, garantindo compatibilidade com o Nextclade.

**Execução**: Chamado por `get_blast_database.smk` após baixar dados do NCBI RefSeq.

**Funções Principais**:
- `clean_cds_name(cds_name)`: Sanitiza nomes de CDS removendo caracteres especiais, truncando para 20 caracteres e padronizando a formatação. Isso é crucial porque o Nextclade pode falhar com nomes de genes complexos, duplicados ou muito longos.
- `jsonl_to_gff(...)`:
  - **Validação**: Verifica se os comprimentos das CDS são múltiplos de 3. Se não, o número de acesso é marcado como inválido e excluído.
  - **Agrupamento**: Agrupa entradas de CDS divididas (ex: genes unidos) por nome para criar características gênicas únicas.
  - **Gene vs CDS**: Se os dados de CDS estiverem ausentes, tenta criar uma característica de gene usando o comprimento total do genoma (se divisível por 3).

### get_minimizer_index.py

**Propósito**: Gera um arquivo JSON de índice de minimizadores a partir de arquivos FASTA de referência, permitindo que o Nextclade mapeie sequências para datasets externos (hospedados no GitHub).

**Execução**: Chamado por `get_public_datasets.smk` para datasets hospedados no GitHub.

**Detalhes de Implementação**:
- **Origem**: Esta é uma adaptação simplificada do script `minimizer.py` do projeto Nextclade.
- **Personalização**: A função `fasta_read` foi modificada para incluir o nome do dataset nas anotações do registro de sequência. Isso garante que o índice gerado associe corretamente minimizadores de sequência com seus caminhos de datasets locais correspondentes.

## Pipeline de Análise

Estes scripts são executados durante o fluxo de trabalho principal `vqc run` para processar sequências, identificar vírus e avaliar a qualidade.

### format_nextclade_sort.py

**Propósito**: Processa a saída do `nextclade sort` para vincular datasets identificados com seus caminhos de arquivos locais e identificar sequências que não corresponderam a nenhum conjunto de dados.

**Execução**: Executado imediatamente após `nextclade sort` em `run_analysis.smk`.

**Funções Principais**:
- `map_datasets_to_local_paths(...)`: Lê a configuração YAML para construir um mapeamento entre nomes de datasets remotos (ex: do Nextclade) e caminhos de armazenamento local.
- `format_nextclade_output(...)`: Mescla resultados de classificação "local" e "externo". Adiciona uma coluna `localDataset` apontando para o diretório do vírus identificado.
- `write_unmapped_sequences(...)`: Extrai sequências que não têm conjunto de dados atribuído (`localDataset` é NaN) e grava seus nomes em `unmapped_sequences.txt` para análise subsequente do BLAST.

### blast_wrapper.py

**Propósito**: Um wrapper em torno do comando `blastn` para lidar com segurança com cabeçalhos FASTA contendo espaços.

**Execução**: Executado pela regra `blast` em `run_analysis.smk` para sequências que não foram identificadas pelo Nextclade.

**Detalhes de Implementação**:
- **Sanitização**: O BLAST pode truncar cabeçalhos no primeiro espaço, levando a incompatibilidades de ID em etapas de pós processamento. Este script verifica espaços nos cabeçalhos.
- **Fluxo de Renomeação**:
  1. Se espaços forem encontrados, gera um FASTA temporário onde as sequências são renomeadas para índices simples (1, 2, 3...).
  2. Salva um arquivo de mapeamento (`mapping.tsv`) vinculando índices aos cabeçalhos originais.
  3. Executa o BLAST com o arquivo renomeado.
  4. Restaura os cabeçalhos originais no TSV de saída do BLAST usando o arquivo de mapeamento.

### reorder_cds.py

**Propósito**: Reordena a string `cdsCoverage` na saída TSV do Nextclade para corresponder à ordem dos genes definida no arquivo GFF.

**Execução**: Executado após cada execução do `nextclade run`.

**Lógica**:
- O Nextclade emite genes em ordem alfabética e omite genes com zero de cobertura.
- Este script lê o GFF para estabelecer a ordem canônica dos genes (posição de `start`).
- Analisa o `cdsCoverage` existente (formato `Gene:Cov,...`), reordena-o e insere `Gene:0.0` para quaisquer genes ausentes. Isso garante ordenação consistente das colunas para processamento posterior.

### post_process_nextclade.py

**Propósito**: O script central de agregação que combina resultados do Nextclade, BLAST e análises genéricas em um relatório final. Calcula métricas de qualidade categóricas (notas A-D) e produz a saída final (TSV, CSV ou JSON).

**Execução**: O passo final da regra `post_process_nextclade`.

**Gerenciamento de Memória e Geradores**:

Este script é projetado para processar datasets massivos (ex: milhões de sequências) com uma pegada de memória constante para saídas CSV/TSV.

1.  **Carregamento Preguiçoso (Lazy Loading) com Geradores**:
    A função `format_dfs` é implementada como um **Gerador** Python. Em vez de retornar uma lista de todos os DataFrames (o que carregaria todos os arquivos na RAM), ela `yields` (produz) um DataFrame processado por vez.
    
    *   **Lógica**: Itera pela lista de arquivos de entrada. Para cada arquivo, lê os dados, enriquece com metadados, otimiza tipos, entrega ao consumidor e **imediatamente** deleta a referência e força a "coleta de lixo" (garbage collection).

2.  **Escrita em Fluxo (Streamed Writing)**:
    A função `write_combined_df` (e sua auxiliar `_write_csv_tsv_output`) consome este gerador. Itera sobre o gerador, escrevendo cada pedaço produzido no disco imediatamente usando `mode='a'` (append/anexar).
    
    *   **Resultado**: Em qualquer momento dado, apenas os dados de um único arquivo de entrada existem na memória.
    
    *   **Limitação JSON**: Para saída JSON (`_write_json_output`), o script *deve* acumular todos os dados para formar um array JSON válido. No entanto, ainda emprega coleta de lixo para descartar artefatos de processamento intermediários assim que são anexados à lista principal.

3.  **Coleta de Lixo Explícita**:
    A coleta de lixo automática do Python pode não ser acionada rápido o suficiente ao lidar com loops grandes e apertados de carregamento de dados. Chamadas explícitas de `del` combinadas com `gc.collect()` são colocadas estrategicamente para garantir que a memória seja liberada de volta para o SO antes de alocar o próximo pedaço.

**Funções Principais**:

*   `format_dfs(files, config_file, blast_metadata_df)`:
    O gerador primário. Determina se um arquivo de resultado pertence a um vírus conhecido (com dataset configurado no YAML) ou é uma execução genérica (nextclade executado com referências informadas pela análise de BLAST). Chama a lógica de processamento apropriada (`_process_with_virus_info` ou `_process_generic_run`) e produz o resultado.


*   `load_blast_metadata(metadata_path)`:
    Carrega os metadados do banco de dados BLAST e normaliza nomes de colunas (ex: `accession` -> `virus`) para garantir consistência com as saídas do Nextclade.

*   `optimize_dataframe_memory(df)`:
    Analisa colunas do DataFrame. Se uma coluna de string (como `virus`, `clade`, `dataset`) tiver uma baixa cardinalidade (número de valores únicos < 50% do total de linhas), converte a coluna para o tipo `category`. Isso reduz drasticamente o uso de RAM.

*   `add_qualities(df, virus_info)`:
    Aplica a lógica de pontuação de qualidade linha por linha. Invoca o auxiliar `_compute_metrics_qualities` e então usa as funções `get_*_quality` para atribuir notas.

*   `add_coverages(df, virus_info)`:
    Analisa e formata a string `cdsCoverage`. Também calcula a cobertura para regiões alvo específicas e as adiciona como novas colunas (`targetRegionsCoverage`, `targetGeneCoverage`).

*   `format_sc2_clade(df, dataset_name)`:
    Contém lógica específica para SARS-CoV-2. Como o Nextclade emite linhagens Pango em uma coluna específica (`Nextclade_pango`), esta função mapeia para a coluna `clade` padrão para consistência.

*   `create_unmapped_df(unmapped_sequences, blast_results, blast_metadata_df)`:
    Lida com sequências que falharam tanto na identificação do Nextclade quanto no BLAST. Lê o arquivo bruto `unmapped_sequences.txt` e cria um DataFrame rotulado como "Unclassified".

*   `write_combined_df(df_iterator, output_file, output_format, ...)`:
    O despachante principal. Recebe o gerador `format_dfs` (iterador) e o direciona para o backend de escrita apropriado (CSV/TSV ou JSON).


**Funções de Controle de Qualidade**:

*   `get_genome_quality(scores)`:
    Agrega pontuações de métricas individuais em uma Qualidade de Genoma final. Soma as pontuações (A=4, B=3, C=2, D=1) e normaliza para uma escala de 24 pontos para atribuir a nota final.
*   `get_target_regions_quality(...)`:
    Determina a qualidade de regiões alvo específicas. Lógica: Se o genoma inteiro for A/B, retorna vazio (implícito bom). Caso contrário, calcula a cobertura média das regiões alvo para atribuir uma nota.
*   `get_cds_cov_quality(...)`:
    Verifica cada CDS contra limites de cobertura para atribuir notas A/B/C/D por gene.
*   `get_missing_data_quality(coverage)`: Pontuado com base em limites (0.9, 0.75, 0.5).
*   `get_private_mutations_quality(total, threshold)`: Pontuado com base no desvio de um limite definido.
*   `get_qc_quality(total)`: Pontuação geral para métricas de contagem (0=A, 1=B, 2=C, >2=D).

### extract_target_regions.py

**Propósito**: Extrai as coordenadas genômicas de regiões de "boa qualidade" para uso posterior (ex: geração de consenso ou design de primers).

**Execução**: Executado pela regra `extract_target_regions` após o pós-processamento.

**Lógica de Seleção**:

A função `check_target_regions` determina a melhor região para extrair com base na qualidade:

1. **Genoma**: Se o `genomeQuality` geral for A ou B, o genoma inteiro é selecionado.
2. **Região Alvo**: Senão, se `targetRegionsQuality` for A ou B, as regiões alvo específicas são selecionadas.
3. **Gene Alvo**: Senão, se `targetGeneQuality` for A ou B, o gene alvo é selecionado.

**Mapeamento de Coordenadas**:
- Usa `get_regions` para consultar as coordenadas de início/fim da característica selecionada (gene ou genoma completo) no arquivo GFF.
- Produz um arquivo BED compatível com `seqtk subseq`.
