import os
from importlib.resources import files

VIRUSES = list(config.keys())

rule all:
	input:
		readmes = expand("datasets/{virus}/README.md", virus=VIRUSES),
		changelogs = expand("datasets/{virus}/CHANGELOG.md", virus=VIRUSES),
		genome_annotations = expand("datasets/{virus}/genome_annotation.gff3", virus=VIRUSES),
		pathogen_configs = expand("datasets/{virus}/pathogen.json", virus=VIRUSES),
		references = expand("datasets/{virus}/reference.fasta", virus=VIRUSES),
		sequences = expand("datasets/{virus}/sequences.fasta", virus=VIRUSES),
		trees = expand("datasets/{virus}/tree.json", virus=VIRUSES)

rule get_nextclade_databases:
	message:
		"""
		Get datasets provided by nextclade
		"""
	params:
		virus = lambda wildcards: wildcards.virus,
		dataset = lambda wildcards: config[wildcards.virus]["dataset"],
		tag = lambda wildcards: config[wildcards.virus]["tag"]
	output:
		readmes = "datasets/{virus}/README.md",
		changelogs = "datasets/{virus}/CHANGELOG.md",
		genome_annotations = "datasets/{virus}/genome_annotation.gff3",
		pathogen_configs = "datasets/{virus}/pathogen.json",
		references = "datasets/{virus}/reference.fasta",
		sequences = "datasets/{virus}/sequences.fasta",
		trees = "datasets/{virus}/tree.json"
	shell:
		"""
		nextclade dataset get \
			--name "{params.dataset}" \
			--tag "{params.tag}" \
			--output-dir "datasets/{params.virus}"
		"""