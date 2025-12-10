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

logging.basicConfig(level=logging.WARNING, format='%(levelname)s:%(message)s')

rule parameters:
    params:
        sequences_fasta = config["sequences_fasta"],
        output_dir = config.get("output_dir", DEFAULT_OUTPUT_DIR),
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



count_get_nextclade_outputs_run = 0
def get_nextclade_outputs(wildcards):
    datasets_selected_file = checkpoints.select_datasets_from_nextclade.get(**wildcards).output.datasets_selected
    viruses = set()
    dataset_not_found = []
    global count_get_nextclade_outputs_run

    with open(datasets_selected_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('localDataset', None):
                virus_name = row.get('localDataset').split('/')[-1]
                viruses.add(virus_name)
            else:
                if count_get_nextclade_outputs_run == 0 and row['dataset'] not in dataset_not_found:
                    logging.warning(f"The '{row['dataset']}' dataset was not found locally.")
                    dataset_not_found.append(row['dataset'])

    nextclade_results = [f"{parameters.output_dir}/nextclade_results/{virus}.nextclade.tsv" for virus in viruses]
    if not nextclade_results and count_get_nextclade_outputs_run == 0:
        logging.warning(f"Nextclade will not run for any input sequence.")
    count_get_nextclade_outputs_run += 1
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
        sequences = parameters.sequences_fasta,
        external_datasets_minimizers = parameters.external_datasets_minimizers
    params:
        output_dir = parameters.output_dir,
        min_score = parameters.nextclade_sort_min_score,
        min_hits = parameters.nextclade_sort_min_hits
    output:
        viruses_identified =  f"{parameters.output_dir}/identified_datasets/viruses.tsv",
        viruses_identified_external =  f"{parameters.output_dir}/identified_datasets/viruses.external_datasets.tsv"
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/nextclade_sort.log"
    shell:
        """
        mkdir -p {params.output_dir}/identified_datasets
        mkdir -p {params.output_dir}/logs

        nextclade sort {input.sequences} \
            --output-path '{params.output_dir}/identified_datasets/{{name}}/sequences.fa' \
            --output-results-tsv {output.viruses_identified} \
            --min-score {params.min_score} \
            --min-hits {params.min_hits} \
            --jobs {threads} 2>{log}

        # Run nextclade sort again using only sequences that were not mapped in the datasets from nextclade_data
        awk -F"\\t" '{{if ($3 == "") print $2}}' \
            {output.viruses_identified} > \
            {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt

        if [ -s {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt ]; then
            seqtk subseq {input.sequences} {params.output_dir}/identified_datasets/tmp_unmapped_sequences.txt > \
                {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta
        fi

        if [ -s {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta ]; then
            nextclade sort {params.output_dir}/identified_datasets/tmp_unmapped_sequences.fasta \
                --input-minimizer-index-json {input.external_datasets_minimizers} \
                --output-path '{params.output_dir}/identified_datasets/{{name}}/sequences.fa' \
                --output-results-tsv {output.viruses_identified_external} \
                --min-score {params.min_score} \
                --min-hits {params.min_hits} \
                --jobs {threads} 2>>{log}
        else
            echo -e "seqName\tdataset\tscore\tnumHits" > {output.viruses_identified_external}
        fi

        rm {params.output_dir}/identified_datasets/tmp_unmapped_sequences.*
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
    shell:
        """
        python {PKG_PATH}/scripts/python/format_nextclade_sort.py \
            --nextclade-output {input.viruses_identified} \
            --nextclade-external-output {input.viruses_identified_external} \
            --config-file {input.config_file} \
            --local-datasets-path {params.datasets_local_path}/ \
            --output-path {params.output_dir}/identified_datasets 
        """


rule blast:
    message:
        "Run BLAST for unmapped sequences"
    input:
        sequences = parameters.sequences_fasta,
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
        mkdir -p {params.output_dir}/blast_results
        if [ -s {input.unmapped_sequences} ]; then
            seqtk subseq {input.sequences} {input.unmapped_sequences} > {params.output_dir}/blast_results/unmapped_sequences.fasta
            blastn -db {input.blast_database} \
                -query {params.output_dir}/blast_results/unmapped_sequences.fasta \
                -out {output.viruses_identified} \
                -task {params.blast_task} \
                -evalue {params.blast_evalue} \
                -qcov_hsp_perc {params.blast_qcov} \
                -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp" \
                -max_hsps 1 \
                -max_target_seqs 1 \
                -perc_identity {params.identity_threshold} \
                -num_threads {threads} 2> {log}
            rm {params.output_dir}/blast_results/unmapped_sequences.fasta
        else
            touch {output.viruses_identified}
        fi
        """

virus_info = None
def get_virus_info(wildcards, field):
    global virus_info
    if virus_info is None:
        virus_info = {}
        datasets_selected_file = checkpoints.select_datasets_from_nextclade.get().output.datasets_selected
        with open(datasets_selected_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                localDataset = row.get('localDataset')
                virus_name = localDataset.split('/')[-1]
                if virus_name not in virus_info:
                    virus_info[virus_name] = {
                        'splittedFasta': row.get('splittedFasta'),
                        'localDataset': row.get('localDataset')
                    }

    return virus_info[wildcards.virus][field]

def get_fasta_for_virus(wildcards):
    return get_virus_info(wildcards, 'splittedFasta')

def get_dataset_for_virus(wildcards):
    return get_virus_info(wildcards, 'localDataset')

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
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/nextclade.{{virus}}.log"
    shell:
        """
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
            --output {output.nextclade_tsv}

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

rule run_generic_nextclade:
    message:
        "Run generic nextclade for blast identified virus {wildcards.virus}"
    input:
        sequences = parameters.sequences_fasta,
        blast_results = rules.blast.output.viruses_identified,
        blast_database = parameters.blast_database
    params:
        output_dir = parameters.output_dir,
        datasets_dir = parameters.datasets_local_path
    output:
        nextclade_tsv = f"{parameters.output_dir}/nextclade_results/{{virus}}.generic.nextclade.tsv",
        nextclade_gff = f"{parameters.output_dir}/gff_files/{{virus}}.generic.nextclade.gff"
    threads:
        parameters.threads
    log:
        f"{parameters.output_dir}/logs/generic_nextclade.{{virus}}.log"
    shell:
        """
        # Get sequences for this virus from blast results
        grep "{wildcards.virus}" {input.blast_results} | cut -f 1 > {params.output_dir}/blast_results/{wildcards.virus}.ids
        seqtk subseq {input.sequences} {params.output_dir}/blast_results/{wildcards.virus}.ids > {params.output_dir}/blast_results/{wildcards.virus}.fasta

        # Get reference sequence from blast db
        echo "{wildcards.virus}" > {params.output_dir}/blast_results/{wildcards.virus}.ref.id
        seqtk subseq {input.blast_database} {params.output_dir}/blast_results/{wildcards.virus}.ref.id > {params.output_dir}/blast_results/{wildcards.virus}.ref.fasta

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
                --output {output.nextclade_tsv}

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
        blast_database_metadata = {parameters.blast_database_metadata}
    params:
        output_format = parameters.output_format
    output:
        output_file = f"{parameters.output_dir}/{parameters.output_file}"
    log:
        f"{parameters.output_dir}/logs/pp_nextclade.log"
    shell:
        """
        python {PKG_PATH}/scripts/python/post_process_nextclade.py \
            --files {input.nextclade_results} \
            --generic-files {input.generic_nextclade_results} \
            --unmapped-sequences {input.unmapped_sequences} \
            --blast-results {input.blast_results} \
            --blast-metadata {input.blast_database_metadata} \
            --config-file {input.config_file} \
            --output {output.output_file} \
            --output-format {params.output_format} 2>{log}
        """

rule extract_target_regions:
    message:
        "Extracts the regions marked as good"
    input:
        sequences = parameters.sequences_fasta,
        post_processed_data = rules.post_process_nextclade.output.output_file
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
        python {PKG_PATH}/scripts/python/extract_target_regions.py \
            --pp-results {input.post_processed_data} \
            --output-format {params.output_format} \
            --output {output.target_regions_bed} 2>{log}

        # Remove range values (:start-end) that seqtk subseq includes in the header.
        seqtk subseq {input.sequences} {output.target_regions_bed} | \
            sed -e 's/\:[0-9]*-[0-9]*$//g' > {output.target_regions_sequences} 2>{log}
        """
