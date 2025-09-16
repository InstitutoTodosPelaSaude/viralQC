output_dir = config["output_dir"]

rule all:
    input:
        blast_database = f"{output_dir}/blast.fasta"

rule makeblast_db:
    message:
        "Create BLAST database"
    params:
        output_dir = output_dir
    output:
        blast_database =  f"{output_dir}/blast.fasta",
    shell:
        """
        mkdir -p {params.output_dir}

        datasets download virus genome taxon 10239 --refseq --include genome --fast-zip-validation 
        unzip -n ncbi_dataset.zip
        sed -e "s/,.*//g" \
            -e "s/virus.*/virus/g" \
            -e "s/MAG: //g" \
            -e "s/UNVERIFIED: //g" \
            -e "s/ /_/" \
            -e "s/ /-/g" ncbi_dataset/data/genomic.fna > {output.blast_database}
        makeblastdb -dbtype nucl -in {output.blast_database}

        rm -rf ncbi_dataset.zip ncbi_dataset/
        """