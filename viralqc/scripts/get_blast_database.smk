from viralqc import PKG_PATH
output_dir = config["output_dir"]
release_date = config.get("release_date", None)

rule all:
    input:
        blast_database = f"{output_dir}/blast.fasta"

rule makeblast_db:
    message:
        "Create BLAST database"
    params:
        output_dir = output_dir,
        release_date = release_date if release_date else ""
    output:
        blast_database =  f"{output_dir}/blast.fasta",
        blast_metadata = f"{output_dir}/blast.tsv",
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.output_dir}/tmp_ncbi

        # Download viral genomes from NCBI RefSeq
        cd {params.output_dir}/tmp_ncbi
        datasets download virus genome taxon 10239 --refseq --include genome,annotation --fast-zip-validation
        unzip -n ncbi_dataset.zip

        # Create metadata with release-date field
        echo -e "accession\\tsegment\\tvirus_name\\tvirus_tax_id\\trelease_date\\tspecies_name\\tspecies_tax_id\\tdatabase_version" > tmp_metadata.tsv
        dataformat tsv virus-genome \\
            --inputfile ncbi_dataset/data/data_report.jsonl \\
            --fields accession,segment,virus-name,virus-tax-id,release-date | \\
            grep -v "Accession" >> tmp_metadata.tsv

        # Filter by release date if provided
        RELEASE_DATE="{params.release_date}"
        if [ -n "$RELEASE_DATE" ]; then
            echo "Filtering sequences by release date: $RELEASE_DATE"
            # Filter metadata to include only sequences with release_date <= user-specified date
            # NCBI release_date is in ISO 8601 format (e.g., 2021-06-01T00:00:00Z), extract just the date part
            head -1 tmp_metadata.tsv > tmp_metadata_filtered.tsv
            awk -F'\\t' -v date="$RELEASE_DATE" 'NR>1 {{ split($5, dt, "T"); if (dt[1] <= date) print }}' tmp_metadata.tsv >> tmp_metadata_filtered.tsv
            mv tmp_metadata_filtered.tsv tmp_metadata.tsv

            # Get list of accessions to keep
            cut -f1 tmp_metadata.tsv | tail -n +2 > accessions_to_keep.txt
            
            # Filter genomic.fna to keep only matching sequences
            seqtk subseq ncbi_dataset/data/genomic.fna accessions_to_keep.txt > genomic_filtered.fna
            mv genomic_filtered.fna ncbi_dataset/data/genomic.fna
            rm accessions_to_keep.txt
        fi

        sed -e "s/ .*//g" ncbi_dataset/data/genomic.fna > ../blast.fasta

        join_files() {{
            local file1="$1"
            local file2="$2"
            
            awk -F'\\t' '
            BEGIN {{ FS = OFS = "\\t" }}
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
                    print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"map[key]
                }} else {{
                    print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tna\\tna"
                }}
            }}
            ' "$file2" "$file1"
        }}
        export -f join_files

        # Download ncbi taxdump files, required for taxonkit tool
        mkdir -p tmp_taxdump
        cd tmp_taxdump
        rm -f taxdump.tar.gz*
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -zxvf taxdump.tar.gz

        mkdir -p $HOME/.taxonkit
        cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit
        cd ..

        # Get virus species name and tax_id
        echo -e "virus_tax_id\\tspecies_tax_id\\tspecies_name" > taxid_mapping.tsv
        cut -f 4 tmp_metadata.tsv | grep -v "virus_tax_id" | taxonkit lineage | taxonkit reformat2 -t -f "{{species}}"  | cut -f 1,3,4 >> taxid_mapping.tsv
        
        # Join the taxid mapping to the metadata
        join_files tmp_metadata.tsv taxid_mapping.tsv > tmp_metadata_with_species.tsv

        # Add database version (use release_date if provided, otherwise current date)
        if [ -n "$RELEASE_DATE" ]; then
            VERSION_DATE="$RELEASE_DATE"
        else
            VERSION_DATE=$(date +%Y-%m-%d)
        fi
        awk -v version="ncbi-refseq-virus_$VERSION_DATE" 'BEGIN{{OFS="\\t"}} NR==1{{print $0; next}} {{print $0, version}}' tmp_metadata_with_species.tsv > ../blast.tsv

        makeblastdb -dbtype nucl -in ../blast.fasta

        # Split GFF (Convert JSONL to GFF)
        python {PKG_PATH}/scripts/python/jsonl_to_gff.py \\
            --input ncbi_dataset/data/annotation_report.jsonl \\
            --output-dir ../blast_gff \\
            --fasta ncbi_dataset/data/genomic.fna

        # Clean up temporary directory
        cd ../..
        rm -rf {params.output_dir}/tmp_ncbi
        """