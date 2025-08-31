from viralqc import PKG_PATH
import csv

rule parameters:
    params:
        sequences_fasta = config["sequences_fasta"],
        sort_mode = config["sort_mode"],
        output_dir = config["output_dir"],
        output_file = config["output_file"],
        output_format = config["output_format"],
        config_file = config["config_file"],
        datasets_local_path = config["datasets_local_path"],
        nextclade_sort_min_score = config["nextclade_sort_min_score"],
        nextclade_sort_min_hits = config["nextclade_sort_min_hits"],
        blast_database = config["blast_database"],
        threads = config["threads"]

parameters = rules.parameters.params

# Checkpoint to ensure datasets_selected.tsv is processed
checkpoint create_datasets_selected:
    input:
        ""
    output:
        tsv = f"{parameters.output_dir}/datasets_selected.tsv"

def get_nextclade_outputs(wildcards):
    datasets_selected_file = checkpoints.create_datasets_selected.get(**wildcards).output.tsv
    viruses = set()
    with open(datasets_selected_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            virus_name = row['localDataset'].split('/')[-1]
            viruses.add(virus_name)
    return [f"{parameters.output_dir}/{virus}.nextclade.tsv" for virus in viruses]

rule all:
    input:
        viruses_identified = f"{parameters.output_dir}/viruses.tsv",
        datasets_selected = f"{parameters.output_dir}/datasets_selected.tsv",
        unmapped_sequences = f"{parameters.output_dir}/unmapped_sequences.txt",
        nextclade_outputs = get_nextclade_outputs,
        output = f"{parameters.output_dir}/{parameters.output_file}"

if parameters.sort_mode == "nextclade":
    rule nextclade_sort:
        message:
            "Run nextclade sort to identify datasets"
        input:
            sequences = parameters.sequences_fasta
        params:
            output_dir = parameters.output_dir,
            min_score = parameters.nextclade_sort_min_score,
            min_hits = parameters.nextclade_sort_min_hits
        output:
            viruses_identified =  f"{parameters.output_dir}/viruses.tsv"
        threads:
            parameters.threads
        log:
            "logs/nextclade_sort.log"
        shell:
            """
            mkdir -p {params.output_dir}

            nextclade sort {input.sequences} \
                --output-path '{params.output_dir}/{{name}}/sequences.fa' \
                --output-results-tsv {output.viruses_identified} \
                --min-score {params.min_score} \
                --min-hits {params.min_hits} \
                --jobs {threads} 2>{log}
            """

    checkpoint select_datasets_from_nextclade:
        message:
            "Select datasets based on nextclade sort output."
        input:
            viruses_identified = rules.nextclade_sort.output.viruses_identified,
            config_file = parameters.config_file,
        params:
            datasets_local_path = parameters.datasets_local_path,
            output_dir = parameters.output_dir
        output:
            datasets_selected = f"{parameters.output_dir}/datasets_selected.tsv",
            unmapped_sequences = f"{parameters.output_dir}/unmapped_sequences.txt"
        threads:
            parameters.threads
        shell:
            """
            python {PKG_PATH}/scripts/python/format_nextclade_sort.py \
                --nextclade-output {input.viruses_identified} \
                --config-file {input.config_file} \
                --local-datasets-path {params.datasets_local_path}/ \
                --output-path {params.output_dir} 
            """
else:
    rule blast:
        message:
            "Run blast to identify datasets"
        input:
            sequences = parameters.sequences_fasta,
            blast_db = parameters.blast_database
        params:
            output_dir = parameters.output_dir,
            datasets_local_path = parameters.datasets_local_path
        output:
            viruses_identified =  f"{parameters.output_dir}/viruses.tsv",
        threads:
            parameters.threads
        log:
            "logs/blast.log"
        shell:
            """
            mkdir -p {params.output_dir}

            blastn -db {input.blast_db} \
                -query {input.sequences} \
                -out {output.viruses_identified} \
                -task megablast \
                -evalue 0.001 \
                -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp" \
                -max_hsps 1 \
                -max_target_seqs 1 \
                -num_threads {threads} 2>{log}

            # Split sequences by virus, similar to nextclade sort
            for virus in $(cut -f 3 {output.viruses_identified} | sort -u);do
                grep $virus {output.viruses_identified} | cut -f 1 > {params.output_dir}/tmp.$virus.sequence_names.txt && \
                seqtk subseq {input.sequences} {params.output_dir}/tmp.$virus.sequence_names.txt > {params.output_dir}/sequences.$virus.fasta;
            done

            rm {params.output_dir}/tmp.*
            """

    checkpoint select_datasets_from_blast:
        message:
            "Select datasets based on blast output."
        input:
            sequences = parameters.sequences_fasta,
            blast_results = rules.blast.output.viruses_identified,
        params:
            datasets_local_path = parameters.datasets_local_path,
            output_dir = parameters.output_dir
        output:
            datasets_selected = f"{parameters.output_dir}/datasets_selected.tsv",
            unmapped_sequences = f"{parameters.output_dir}/unmapped_sequences.txt"
        threads:
            parameters.threads
        shell:
            """
            # Create datasets information file with header similar to nextclade sort

            echo -e "seqName\tlocalDataset\tsplittedFasta" > {output.datasets_selected}
            awk -F'\t' \
                -v datasets_local_path="{params.datasets_local_path}" \
                -v output_dir="{params.output_dir}" \
                'BEGIN{{OFS="\t"}} {{print $1, datasets_local_path "/" $3, output_dir "/sequences." $3 ".fasta"}}' \
                {input.blast_results} >> {output.datasets_selected}

            # Save sequences without matches
            grep ">" {input.sequences} | sed -e "s/>//g" > {params.output_dir}/tmp.sequence_ids.txt
            cut -f 1 {input.blast_results} > {params.output_dir}/tmp.sequences_with_match.txt
            fgrep -vw -f {params.output_dir}/tmp.sequences_with_match.txt {params.output_dir}/tmp.sequence_ids.txt > {output.unmapped_sequences} || true
            
            rm {params.output_dir}/tmp.*
            """

virus_info = None
def get_virus_info(wildcards, field):
    global virus_info
    if virus_info is None:
        virus_info = {}
        datasets_selected_file = checkpoints.create_datasets_selected.get().output.tsv
        with open(datasets_selected_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                localDataset = row['localDataset']
                virus_name = localDataset.split('/')[-1]
                if virus_name not in virus_info:
                    virus_info[virus_name] = {
                        'splittedFasta': row['splittedFasta'],
                        'localDataset': row['localDataset']
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
        nextclade_tsv = f"{parameters.output_dir}/{{virus}}.nextclade.tsv"
    threads:
        parameters.threads
    log:
        "logs/nextclade.{virus}.log"
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-tsv {output.nextclade_tsv} \
            --jobs {threads} \
            {input.fasta} 2>{log}
        """

rule post_process_nextclade:
    message:
        "Process nextclade outputs"
    input:
        nextclade_results = get_nextclade_outputs,
        unmapped_sequences = f"{parameters.output_dir}/unmapped_sequences.txt",
        config_file = parameters.config_file
    params:
        output_format = parameters.output_format
    output:
        output_file = f"{parameters.output_dir}/{parameters.output_file}"
    log:
        "logs/pp_nextclade.log"
    shell:
        """
        python {PKG_PATH}/scripts/python/post_process_nextclade.py \
            --files {input.nextclade_results} \
            --unmapped-sequences {input.unmapped_sequences} \
            --config-file {input.config_file} \
            --output {output.output_file} \
            --output-format {params.output_format} 2>{log}
        """