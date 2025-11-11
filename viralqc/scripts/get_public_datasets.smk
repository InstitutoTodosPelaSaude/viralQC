from viralqc import PKG_PATH

rule parameters:
    params:
        nextclade_viruses = [
            k for k, v in config.get("nextclade_data", {}).items()
            if k != "datasets_dir" and isinstance(v, dict)
        ],
        github_viruses = [
            k for k, v in config.get("github", {}).items()
            if k != "datasets_dir" and isinstance(v, dict)
        ],
        datasets_dir = config["datasets_dir"]

parameters = rules.parameters.params

rule all:
    input:
        # Nextclade datasets
        readmes_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/README.md", virus=parameters.nextclade_viruses),
        changelogs_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md", virus=parameters.nextclade_viruses),
        genome_annotations_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3", virus=parameters.nextclade_viruses),
        pathogen_configs_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/pathogen.json", virus=parameters.nextclade_viruses),
        references_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/reference.fasta", virus=parameters.nextclade_viruses),
        sequences_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/sequences.fasta", virus=parameters.nextclade_viruses),
        trees_nextclade = expand(f"{parameters.datasets_dir}/{{virus}}/tree.json", virus=parameters.nextclade_viruses),

        # GitHub datasets
        readmes_github = expand(f"{parameters.datasets_dir}/{{virus}}/README.md", virus=parameters.github_viruses),
        changelogs_github = expand(f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md", virus=parameters.github_viruses),
        genome_annotations_github = expand(f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3", virus=parameters.github_viruses),
        pathogen_configs_github = expand(f"{parameters.datasets_dir}/{{virus}}/pathogen.json", virus=parameters.github_viruses),
        references_github = expand(f"{parameters.datasets_dir}/{{virus}}/reference.fasta", virus=parameters.github_viruses),
        sequences_github = expand(f"{parameters.datasets_dir}/{{virus}}/sequences.fasta", virus=parameters.github_viruses),
        trees_github = expand(f"{parameters.datasets_dir}/{{virus}}/tree.json", virus=parameters.github_viruses),
        minimizer_index = f"{parameters.datasets_dir}/external_datasets_minimizers.json"

rule get_nextclade_databases:
    message:
        "Get datasets provided by nextclade"
    params:
        datasets_dir = parameters.datasets_dir,
        dataset = lambda wc: config["nextclade_data"][wc.virus]["dataset"],
        tag = lambda wc: config["nextclade_data"][wc.virus]["tag"]
    output:
        readmes = f"{parameters.datasets_dir}/{{virus}}/README.md",
        changelogs = f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md",
        genome_annotations = f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3",
        pathogen_configs = f"{parameters.datasets_dir}/{{virus}}/pathogen.json",
        references = f"{parameters.datasets_dir}/{{virus}}/reference.fasta",
        sequences = f"{parameters.datasets_dir}/{{virus}}/sequences.fasta",
        trees = f"{parameters.datasets_dir}/{{virus}}/tree.json"
    wildcard_constraints:
        virus='|'.join(parameters.nextclade_viruses)
    shell:
        """
        nextclade dataset get \
            --name "{params.dataset}" \
            --tag "{params.tag}" \
            --output-dir "{params.datasets_dir}/{wildcards.virus}"

        for f in sequences.fasta tree.json; do
            [ -f "{params.datasets_dir}/{wildcards.virus}/$f" ] || touch "{params.datasets_dir}/{wildcards.virus}/$f"
        done
        """

rule get_github_databases:
    message:
        "Get datasets hosted on github repositories"
    params:
        datasets_dir = parameters.datasets_dir,
        repository = lambda wc: config["github"][wc.virus]["repository"],
        dataset = lambda wc: config["github"][wc.virus]["dataset"],
    output:
        readmes = f"{parameters.datasets_dir}/{{virus}}/README.md",
        changelogs = f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md",
        genome_annotations = f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3",
        pathogen_configs = f"{parameters.datasets_dir}/{{virus}}/pathogen.json",
        references = f"{parameters.datasets_dir}/{{virus}}/reference.fasta",
        sequences = f"{parameters.datasets_dir}/{{virus}}/sequences.fasta",
        trees = f"{parameters.datasets_dir}/{{virus}}/tree.json"
    wildcard_constraints:
        virus='|'.join(parameters.github_viruses)
    shell:
        """
        python {PKG_PATH}/scripts/python/get_github_dataset.py \
            --repository "{params.repository}" \
            --dataset-path "{params.dataset}" \
            --output-dir "{params.datasets_dir}/{wildcards.virus}"
        """

rule create_minimzer_indexes:
    message:
        "Get minimizer indexes for github datasets"
    input:
        references = expand(
            "{datasets_dir}/{virus}/reference.fasta",
            datasets_dir=parameters.datasets_dir,
            virus=parameters.github_viruses
        )
    output:
        minimizer_index = f"{parameters.datasets_dir}/external_datasets_minimizers.json"
    shell:
        """
        python {PKG_PATH}/scripts/python/get_minimizer_index.py \
            --input-refs {input.references} \
            --output-json {output.minimizer_index}
        """