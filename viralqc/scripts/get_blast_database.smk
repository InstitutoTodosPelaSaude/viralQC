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
        blast_metadata = f"{output_dir}/blast.tsv",
    shell:
        """
        mkdir -p {params.output_dir}

        # Download viral genomes from NCBI RefSeq
        datasets download virus genome taxon 10239 --refseq --include genome --fast-zip-validation
        unzip -n ncbi_dataset.zip
        sed -e "s/ .*//g" ncbi_dataset/data/genomic.fna > {output.blast_database}
        echo -e "accession\tsegment\tvirus_name\tvirus_tax_id\tspecies_name\tspecies_tax_id\tdatabase_version" > tmp_metadata.tsv
        dataformat tsv virus-genome \
            --inputfile ncbi_dataset/data/data_report.jsonl \
            --fields accession,segment,virus-name,virus-tax-id | \
            grep -v "Accession" >> tmp_metadata.tsv

        join_files() {{
            local file1="$1"
            local file2="$2"
            
            awk -F'\\t' '
            BEGIN {{ FS = OFS = "\t" }}
            FNR==NR {{
                map[$1]=$2"\\t"$3
                next
            }}
            FNR==1 {{
                print $0
                next
            }}
            {{
                key=$4
                if(key in map){{
                    print $1"\\t"$2"\\t"$3"\\t"$4"\\t"map[key]
                }} else {{
                    print $1"\\t"$2"\\t"$3"\\t"$4"\\tna\\tna"
                }}
            }}
            ' "$file2" "$file1"
        }}
        export -f join_files

        # Download ncbi taxdump files, required for taxonkit tool
        mkdir -p tmp_taxdump
        cd tmp_taxdump
        wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -zxvf taxdump.tar.gz

        mkdir -p $HOME/.taxonkit
        cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit
        cd ..

        # Get virus species name and tax_id
        echo -e "virus_tax_id\\tspecies_tax_id\\tspecies_name" > taxid_mapping.tsv
        cut -f 4 tmp_metadata.tsv | grep -v "virus_tax_id" | taxonkit lineage | taxonkit reformat2 -t -f "{{species}}"  | cut -f 1,3,4 >> taxid_mapping.tsv
        
        # Join the taxid mapping to the metadata
        join_files tmp_metadata.tsv taxid_mapping.tsv > tmp_metadata_with_species.tsv

        # Add database version
        awk -v version="ncbi-refseq-virus_$(date +%Y-%m-%d)" 'BEGIN{{OFS="\\t"}} NR==1{{print $0; next}} {{print $0, version}}' tmp_metadata_with_species.tsv > {output.blast_metadata}

        makeblastdb -dbtype nucl -in {output.blast_database}
        rm -rf ncbi_dataset.zip ncbi_dataset/ tmp_taxdump/ tmp_metadata.tsv tmp_metadata_with_species.tsv taxid_mapping.tsv
        """