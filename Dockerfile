FROM mambaorg/micromamba:1.5.1

USER root
ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /app
ENV DATASETS_DIR="/app/datasets"

RUN apt-get update && apt-get install -y unzip && apt-get clean

RUN micromamba install -y -n base \
    -c conda-forge \
    -c bioconda \
    snakemake-minimal=7.32.4 \
    blast=2.16.0 \
    nextclade=3.15.0 \
    seqtk=1.5 \
    ncbi-datasets-cli=18.9.0 \
    taxonkit=0.20.0 \
    pip && \
    micromamba clean -ay

RUN pip install viralQC

RUN set -eo pipefail && \
    vqc get-nextclade-datasets --datasets-dir $DATASETS_DIR -v && \
    vqc get-blast-database --output-dir $DATASETS_DIR -v

ENTRYPOINT ["vqc"]