from viralqc import PKG_PATH
from viralqc.core.defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_OUTPUT_FILE,
    DEFAULT_OUTPUT_FORMAT,
    DEFAULT_DATASETS_LOCAL_PATH,
    DEFAULT_NEXTCLADE_SORT_MIN_SCORE,
    DEFAULT_NEXTCLADE_SORT_MIN_HITS,
    DEFAULT_BLAST_DATABASE,
    DEFAULT_BLAST_DATABASE_METADATA,
    DEFAULT_BLAST_IDENTITY_THRESHOLD,
    DEFAULT_BLAST_EVALUE,
    DEFAULT_BLAST_QCOV,
    DEFAULT_BLAST_TASK,
    DEFAULT_CORES,
)
import csv
import logging
import os

workflow_logger = logging.getLogger("viralqc.workflow")
workflow_logger.setLevel(logging.WARNING)

rule parameters:
    params:
        sequences_fasta = config["sequences_fasta"],
        output_dir = DEFAULT_OUTPUT_DIR,
        output_file = config.get("output_file", DEFAULT_OUTPUT_FILE),
        output_format = config.get("output_format", DEFAULT_OUTPUT_FORMAT),
        config_file = config.get("config_file", DEFAULT_CONFIG_FILE),
        datasets_local_path = config.get("datasets_local_path", DEFAULT_DATASETS_LOCAL_PATH),
        external_datasets_minimizers = f"{config.get('datasets_local_path', DEFAULT_DATASETS_LOCAL_PATH)}/external_datasets_minimizers.json",
        nextclade_sort_min_score = config.get("nextclade_sort_min_score", DEFAULT_NEXTCLADE_SORT_MIN_SCORE),
        nextclade_sort_min_hits = config.get("nextclade_sort_min_hits", DEFAULT_NEXTCLADE_SORT_MIN_HITS),
        blast_database = config.get("blast_database", DEFAULT_BLAST_DATABASE),
        blast_database_metadata = config.get("blast_database_metadata", DEFAULT_BLAST_DATABASE_METADATA),
        blast_identity_threshold = config.get("blast_identity_threshold", DEFAULT_BLAST_IDENTITY_THRESHOLD),
        blast_evalue = config.get("blast_evalue", DEFAULT_BLAST_EVALUE),
        blast_qcov = config.get("blast_qcov", DEFAULT_BLAST_QCOV),
        blast_task = config.get("blast_task", DEFAULT_BLAST_TASK),
        threads = config.get("threads", DEFAULT_CORES)

parameters = rules.parameters.params


rule itemize_sequences:
    message:
        "Renaming sequences to sequential IDs"
    input:
        sequences = parameters.sequences_fasta
    output:
        fasta = f"{parameters.output_dir}/sequences_sanitized.fasta",
        mapping = f"{parameters.output_dir}/identified_datasets/id_mapping.tsv"
    shell:
        """
        python {PKG_PATH}/scripts/python/itemize_sequences.py \
            --input {input.sequences} \
            --output-fasta {output.fasta} \
            --output-mapping {output.mapping}
        """

_run_counters = {}
def get_nextclade_outputs(wildcards):
    output_dir = parameters.output_dir
    if output_dir not in _run_counters:
        _run_counters[output_dir] = 0

    datasets_selected_file = checkpoints.select_datasets_from_nextclade.get(**wildcards).output.datasets_selected
    viruses = set()
    dataset_not_found = []

    with open(datasets_selected_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('localDataset', None):
                virus_name = row.get('localDataset').split('/')[-1]
                viruses.add(virus_name)
            else:
                if _run_counters[output_dir] == 0 and row['dataset'] not in dataset_not_found:
                    workflow_logger.warning(f"The '{row['dataset']}' dataset was not found locally.")
                    dataset_not_found.append(row['dataset'])

    nextclade_results = [f"{parameters.output_dir}/nextclade_results/{virus}.nextclade.tsv" for virus in viruses]
    if not nextclade_results and _run_counters[output_dir] == 0:
        workflow_logger.warning(f"Nextclade will not run for any input sequence.")
    _run_counters[output_dir] += 1
    return nextclade_results

rule all:
    input:
        viruses_identified = f"{parameters.output_dir}/identified_datasets/viruses.tsv",
        datasets_selected = f"{parameters.output_dir}/identified_datasets/datasets_selected.tsv",
        unmapped_sequences = f"{parameters.output_dir}/identified_datasets/unmapped_sequences.txt",
        nextclade_outputs = get_nextclade_outputs,
        output = f"{parameters.output_dir}/{parameters.output_file}",
        target_regions_bed = f"{parameters.output_dir}/sequences_target_regions.bed",
        target_regions_sequences = f"{parameters.output_dir}/sequences_target_regions.fasta"

rule nextclade_sort:
    message:
        "Run nextclade sort to identify datasets"
    input:
        sequences = rules.itemize_sequences.output.fasta
    params:
        output_dir = parameters.output_dir,
        min_score = parameters.nextclade_sort_min_score,
        min_hits = parameters.nextclade_sort_min_hits,
        external_datasets_minimizers = parameters.external_datasets_minimizers
    output:
        viruses_identified =  f"{parameters.output_dir}/identified_datasets/viruses.tsv",
        viruses_identified_external =  f"{parameters.output_dir}/identified_datasets/viruses.external_datasets.tsv"
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/nextclade_sort.log"
    shell:
        """
        set -euo pipefail

        mkdir -p {params.output_dir}/identified_datasets 2>{log}
        mkdir -p {params.output_dir}/logs 2>>{log}

        nextclade sort {input.sequences} \
            --output-path '{params.output_dir}/identified_datasets/{{name}}/sequences.fa' \
            --output-results-tsv {output.viruses_identified} \
            --min-score {params.min_score} \
            --min-hits {params.min_hits} \
            --jobs {threads} 2>>{log}

        # Run nextclade sort again using only sequences that were not mapped in the datasets from nextclade_data
        # This only runs if external_datasets_minimizers.json exists (i.e., github datasets are configured)
        awk -F"\\t" '{{if ($3 == "") print $2}}' \
            {output.viruses_identified} > \
            {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt 2>>{log}

        if [ -f "{params.external_datasets_minimizers}" ] && [ -s {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt ]; then
            seqtk subseq {input.sequences} {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt > \
                {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta 2>>{log}

            if [ -s {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta ]; then
                nextclade sort {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta \
                    --input-minimizer-index-json {params.external_datasets_minimizers} \
                    --output-path '{params.output_dir}/identified_datasets/{{name}}/sequences.fa' \
                    --output-results-tsv {output.viruses_identified_external} \
                    --min-score {params.min_score} \
                    --min-hits {params.min_hits} \
                    --jobs {threads} 2>>{log}
            else
                echo -e "seqName\tdataset\tscore\tnumHits" > {output.viruses_identified_external} 2>>{log}
            fi
        else
            # No external datasets configured or no unmapped sequences
            echo -e "seqName\tdataset\tscore\tnumHits" > {output.viruses_identified_external} 2>>{log}
        fi

        rm {params.output_dir}/identified_datasets/tmp_unmapped_sequences.* || true 2>>{log}
        """

checkpoint select_datasets_from_nextclade:
    message:
        "Select datasets based on nextclade sort output."
    input:
        viruses_identified = rules.nextclade_sort.output.viruses_identified,
        viruses_identified_external = rules.nextclade_sort.output.viruses_identified_external,
        config_file = parameters.config_file,
    params:
        datasets_local_path = parameters.datasets_local_path,
        output_dir = parameters.output_dir
    output:
        datasets_selected = f"{parameters.output_dir}/identified_datasets/datasets_selected.tsv",
        unmapped_sequences = f"{parameters.output_dir}/identified_datasets/unmapped_sequences.txt"
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/format_nextclade_sort.log"
    shell:
        """
        set -euo pipefail

        python {PKG_PATH}/scripts/python/format_nextclade_sort.py \
            --nextclade-output {input.viruses_identified} \
            --nextclade-external-output {input.viruses_identified_external} \
            --config-file {input.config_file} \
            --local-datasets-path {params.datasets_local_path}/ \
            --output-path {params.output_dir}/identified_datasets 2>>{log}
        """


rule blast:
    message:
        "Run BLAST for unmapped sequences"
    input:
        sequences = rules.itemize_sequences.output.fasta,
        unmapped_sequences = rules.select_datasets_from_nextclade.output.unmapped_sequences,
        blast_database = parameters.blast_database
    params:
        identity_threshold = parameters.blast_identity_threshold,
        blast_evalue = parameters.blast_evalue,
        blast_qcov = parameters.blast_qcov,
        blast_task = parameters.blast_task,
        output_dir = parameters.output_dir
    output:
        viruses_identified = f"{parameters.output_dir}/blast_results/unmapped_sequences.blast.tsv",
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/blast.log"
    shell:
        """
        set -euo pipefail

        mkdir -p {params.output_dir}/blast_results 2>>{log}
        if [ -s {input.unmapped_sequences} ]; then
            seqtk subseq {input.sequences} {input.unmapped_sequences} > {params.output_dir}/blast_results/unmapped_sequences.fasta 2>>{log}

            blastn \
                -query {params.output_dir}/blast_results/unmapped_sequences.fasta \
                -db {input.blast_database} \
                -out {output.viruses_identified} \
                -task {params.blast_task} \
                -evalue {params.blast_evalue} \
                -qcov_hsp_perc {params.blast_qcov} \
                -perc_identity {params.identity_threshold} \
                -num_threads {threads} \
                -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp" \
                -max_hsps 1 \
                -max_target_seqs 1 2>>{log}

            rm {params.output_dir}/blast_results/unmapped_sequences.fasta 2>>{log}
        else
            touch {output.viruses_identified}
        fi
        """

_virus_info_cache = {}
def get_virus_info(wildcards, field):
    output_dir = parameters.output_dir
    if output_dir not in _virus_info_cache:
        _virus_info_cache[output_dir] = {}
        datasets_selected_file = checkpoints.select_datasets_from_nextclade.get().output.datasets_selected
        with open(datasets_selected_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                localDataset = row.get('localDataset')
                virus_name = localDataset.split('/')[-1]
                if virus_name not in _virus_info_cache[output_dir]:
                    _virus_info_cache[output_dir][virus_name] = {
                        'splittedFasta': row.get('splittedFasta'),
                        'localDataset': row.get('localDataset')
                    }

    return _virus_info_cache[output_dir][wildcards.virus][field]

def get_fasta_for_virus(wildcards):
    return get_virus_info(wildcards, 'splittedFasta')

def get_dataset_for_virus(wildcards):
    return get_virus_info(wildcards, 'localDataset')



def get_nextclade_threads(wildcards):
    try:
        datasets_selected_file = checkpoints.select_datasets_from_nextclade.get().output.datasets_selected
        viruses = set()
        with open(datasets_selected_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('localDataset', None):
                    viruses.add(row.get('localDataset'))
        
        count = len(viruses)
        if count == 0:
            return 1
        return max(1, int(parameters.threads / count))
    except Exception:
        return 1


rule nextclade:
    message:
        "Run nextclade for virus {wildcards.virus}"
    input:
        fasta = get_fasta_for_virus,
        dataset = get_dataset_for_virus
    output:
        nextclade_tsv = f"{parameters.output_dir}/nextclade_results/{{virus}}.nextclade.tsv",
        nextclade_gff = f"{parameters.output_dir}/gff_files/{{virus}}.nextclade.gff"
    threads:
        get_nextclade_threads
    log:
        f"{parameters.output_dir}/logs/nextclade.{{virus}}.log"
    shell:
        """
        set -euo pipefail

        nextclade run \
            --input-dataset {input.dataset} \
            --output-tsv {output.nextclade_tsv}.tmp \
            --output-annotation-gff {output.nextclade_gff} \
            --min-seed-cover 0.05 \
            --jobs {threads} \
            {input.fasta} 2>{log}

        # Reorder cdsCoverage column based on GFF gene order
        python {PKG_PATH}/scripts/python/reorder_cds.py \
            --nextclade-tsv {output.nextclade_tsv}.tmp \
            --gff {input.dataset}/genome_annotation.gff3 \
            --output {output.nextclade_tsv} 2>>{log}

        rm {output.nextclade_tsv}.tmp
        """

checkpoint process_blast_results:
    input:
        blast_results = rules.blast.output.viruses_identified
    output:
        blast_viruses_list = f"{parameters.output_dir}/blast_results/blast_viruses.list"
    shell:
        """
        cut -f 3 {input.blast_results} | sort | uniq > {output.blast_viruses_list}
        """

def get_generic_nextclade_outputs(wildcards):
    blast_viruses_list = checkpoints.process_blast_results.get(**wildcards).output.blast_viruses_list
    blast_viruses = []
    if os.path.exists(blast_viruses_list):
        with open(blast_viruses_list, 'r') as f:
            blast_viruses = [line.strip() for line in f if line.strip()]

    return [f"{parameters.output_dir}/nextclade_results/{virus}.generic.nextclade.tsv" for virus in blast_viruses]

def get_generic_nextclade_threads(wildcards):
    try:
        blast_viruses_list = checkpoints.process_blast_results.get().output.blast_viruses_list
        count = 0
        if os.path.exists(blast_viruses_list):
            with open(blast_viruses_list, 'r') as f:
                 count = sum(1 for line in f if line.strip())
        
        if count == 0:
            return 1
        return max(1, int(parameters.threads / count))
    except Exception:
        return 1


rule run_generic_nextclade:
    message:
        "Run generic nextclade for blast identified virus {wildcards.virus}"
    input:
        sequences = rules.itemize_sequences.output.fasta,
        blast_results = rules.blast.output.viruses_identified,
        blast_database = parameters.blast_database
    params:
        output_dir = parameters.output_dir,
        datasets_dir = parameters.datasets_local_path
    output:
        nextclade_tsv = f"{parameters.output_dir}/nextclade_results/{{virus}}.generic.nextclade.tsv",
        nextclade_gff = f"{parameters.output_dir}/gff_files/{{virus}}.generic.nextclade.gff"
    threads:
        get_generic_nextclade_threads
    log:
        f"{parameters.output_dir}/logs/generic_nextclade.{{virus}}.log"
    shell:
        """
        set -euo pipefail

        # Get sequences for this virus from blast results
        grep "{wildcards.virus}" {input.blast_results} | cut -f 1 > {params.output_dir}/blast_results/{wildcards.virus}.ids 2>>{log}
        seqtk subseq {input.sequences} {params.output_dir}/blast_results/{wildcards.virus}.ids > {params.output_dir}/blast_results/{wildcards.virus}.fasta 2>>{log}

        # Get reference sequence from blast db
        echo "{wildcards.virus}" > {params.output_dir}/blast_results/{wildcards.virus}.ref.id 2>>{log}
        seqtk subseq {input.blast_database} {params.output_dir}/blast_results/{wildcards.virus}.ref.id > {params.output_dir}/blast_results/{wildcards.virus}.ref.fasta 2>>{log}

        # Check if GFF exists
        GFF_FILE="{params.datasets_dir}/blast_gff/{wildcards.virus}.gff"
        
        if [ -f "$GFF_FILE" ]; then
            echo "Valid GFF found for {wildcards.virus}. Running with annotation." >> {log}
            nextclade run \
                --input-ref {params.output_dir}/blast_results/{wildcards.virus}.ref.fasta \
                --input-annotation "$GFF_FILE" \
                --output-tsv {output.nextclade_tsv}.tmp \
                --output-annotation-gff {output.nextclade_gff} \
                --min-seed-cover 0.05 \
                --jobs {threads} \
                {params.output_dir}/blast_results/{wildcards.virus}.fasta 2>>{log}

            # Reorder cdsCoverage column based on GFF gene order
            python {PKG_PATH}/scripts/python/reorder_cds.py \
                --nextclade-tsv {output.nextclade_tsv}.tmp \
                --gff "$GFF_FILE" \
                --output {output.nextclade_tsv} 2>>{log}

            rm {output.nextclade_tsv}.tmp
        else
            echo "No GFF found for {wildcards.virus}. Running without annotation." >> {log}
            nextclade run \
                --input-ref {params.output_dir}/blast_results/{wildcards.virus}.ref.fasta \
                --output-tsv {output.nextclade_tsv} \
                --min-seed-cover 0.05 \
                --jobs {threads} \
                {params.output_dir}/blast_results/{wildcards.virus}.fasta 2>>{log}
            touch {output.nextclade_gff}
        fi

        rm {params.output_dir}/blast_results/{wildcards.virus}.ids {params.output_dir}/blast_results/{wildcards.virus}.ref.id {params.output_dir}/blast_results/{wildcards.virus}.fasta {params.output_dir}/blast_results/{wildcards.virus}.ref.fasta
        """


rule post_process_nextclade:
    message:
        "Process nextclade outputs"
    input:
        nextclade_results = get_nextclade_outputs,
        generic_nextclade_results = get_generic_nextclade_outputs,
        blast_results = rules.blast.output.viruses_identified,
        unmapped_sequences = f"{parameters.output_dir}/identified_datasets/unmapped_sequences.txt",
        config_file = parameters.config_file,
        blast_database_metadata = {parameters.blast_database_metadata},
        id_mapping = rules.itemize_sequences.output.mapping
    params:
        output_format = parameters.output_format
    output:
        output_file = f"{parameters.output_dir}/{parameters.output_file}"
    log:
        f"{parameters.output_dir}/logs/pp_nextclade.log"
    shell:
        """
        set -euo pipefail

        python {PKG_PATH}/scripts/python/post_process_nextclade.py \
            --files {input.nextclade_results} \
            --generic-files {input.generic_nextclade_results} \
            --unmapped-sequences {input.unmapped_sequences} \
            --blast-results {input.blast_results} \
            --blast-metadata {input.blast_database_metadata} \
            --id-mapping {input.id_mapping} \
            --config-file {input.config_file} \
            --output {output.output_file} \
            --output-format {params.output_format} 2>{log}
        """

rule extract_target_regions:
    message:
        "Extracts the regions marked as good"
    input:
        sequences = parameters.sequences_fasta,
        post_processed_data = rules.post_process_nextclade.output.output_file,
        id_mapping = rules.itemize_sequences.output.mapping
    params:
        output_format = parameters.output_format
    output:
        target_regions_bed = f"{parameters.output_dir}/sequences_target_regions.bed",
        target_regions_sequences = f"{parameters.output_dir}/sequences_target_regions.fasta"
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/extract_target_regions.log"
    shell:
        """
        set -euo pipefail

        python {PKG_PATH}/scripts/python/extract_target_regions.py \
            --pp-results {input.post_processed_data} \
            --output-format {params.output_format} \
            --id-mapping {input.id_mapping} \
            --output {output.target_regions_bed} 2>{log}

        # Remove range values (:start-end) that seqtk subseq includes in the header.
        seqtk subseq {input.sequences} {output.target_regions_bed} | \
            sed -e 's/\:[0-9]*-[0-9]*$//g' > {output.target_regions_sequences} 2>>{log}
        
        # Clean up intermediate files
        rm -rf {parameters.output_dir}/sequences_sanitized.fasta
        find {parameters.output_dir}/identified_datasets/ -mindepth 1 -type d -exec rm -rf {{}} +
        """