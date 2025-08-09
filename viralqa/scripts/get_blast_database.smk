rule parameters:
	params:
		datasets = config["datasets"],
        output_dir = config["output_dir"]

parameters = rules.parameters.params

rule all:
    input:
        blast_database = f"{parameters.output_dir}/viruses.fa",

rule makeblast_db:
    message:
        "Create BLAST database"
    params:
        datasets_dir = parameters.datasets,
        output_dir = parameters.output_dir
    output:
        sort_results =  f"{parameters.output_dir}/viruses.fa",
    shell:
        """
        mkdir -p {params.output_dir}

        for i in $(ls {params.datasets_dir}/*/reference.fasta);do
            ref_name=$(grep ">" $i | sed -e "s/>//g") && \
            virus_id=$(echo $i | sed -e "s/\/reference.fasta//g" -e "s/{params.datasets_dir}\///g") && \
            sed -e "s#$ref_name#$virus_id#g" $i >> {params.output_dir}/viruses.fa;
        done

        makeblastdb -dbtype nucl -in {params.output_dir}/viruses.fa
        """