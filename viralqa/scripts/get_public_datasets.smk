rule parameters:
	params:
		viruses = [
			k for k, v in config.items()
			if k != "datasets_dir" and isinstance(v, dict)
		],
        datasets_dir = config["datasets_dir"]

parameters = rules.parameters.params

rule all:
    input:
        readmes = expand(f"{parameters.datasets_dir}/{{virus}}/README.md", virus=parameters.viruses),
        changelogs = expand(f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md", virus=parameters.viruses),
        genome_annotations = expand(f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3", virus=parameters.viruses),
        pathogen_configs = expand(f"{parameters.datasets_dir}/{{virus}}/pathogen.json", virus=parameters.viruses),
        references = expand(f"{parameters.datasets_dir}/{{virus}}/reference.fasta", virus=parameters.viruses),
        sequences = expand(f"{parameters.datasets_dir}/{{virus}}/sequences.fasta", virus=parameters.viruses),
        trees = expand(f"{parameters.datasets_dir}/{{virus}}/tree.json", virus=parameters.viruses)

rule get_nextclade_databases:
    message:
        "Get datasets provided by nextclade"
    params:
        datasets_dir = parameters.datasets_dir,
        dataset = lambda wc: config[wc.virus]["dataset"],
        tag = lambda wc: config[wc.virus]["tag"]
    output:
        readmes = f"{parameters.datasets_dir}/{{virus}}/README.md",
        changelogs = f"{parameters.datasets_dir}/{{virus}}/CHANGELOG.md",
        genome_annotations = f"{parameters.datasets_dir}/{{virus}}/genome_annotation.gff3",
        pathogen_configs = f"{parameters.datasets_dir}/{{virus}}/pathogen.json",
        references = f"{parameters.datasets_dir}/{{virus}}/reference.fasta",
        sequences = f"{parameters.datasets_dir}/{{virus}}/sequences.fasta",
        trees = f"{parameters.datasets_dir}/{{virus}}/tree.json"
    shell:
        """
        nextclade dataset get \
            --name "{params.dataset}" \
            --tag "{params.tag}" \
            --output-dir "{params.datasets_dir}/{wildcards.virus}"
        """