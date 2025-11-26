import json
import argparse
import os
import sys
import re


def parse_fasta_lengths(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequence lengths.

    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        A dictionary where keys are accession IDs and values are sequence lengths.
    """
    lengths = {}
    current_acc = None
    current_len = 0

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_acc:
                    lengths[current_acc] = current_len
                current_acc = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_acc:
            lengths[current_acc] = current_len

    return lengths


def clean_cds_name(cds_name):
    """
    Clean and truncate CDS name. This is necessary because some CDS names are
    too long and contain special characters that break nextclade execution

    Args:
        cds_name: The original CDS name.

    Returns:
        The cleaned name.
    """
    cds_name = re.sub(r"[ ,/]+", "_", cds_name).capitalize().strip("'")

    if ";" in cds_name:
        cds_name = cds_name.split(";")[0].strip()
    if ":" in cds_name:
        cds_name = cds_name.split(":")[0].strip()
    if len(cds_name) > 20:
        cds_name = cds_name[:20]

    return cds_name


def jsonl_to_gff(jsonl_file, output_dir, fasta_file):
    """
    Convert NCBI Datasets JSONL to GFF3 format.

    Args:
        jsonl_file: Path to the input JSONL file.
        output_dir: Directory where GFF files will be saved.
        fasta_file: Path to the FASTA file for sequence lengths.

    Returns:
        None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    lengths = parse_fasta_lengths(fasta_file)
    invalid_accessions = []

    with open(jsonl_file, "r") as f:
        for line in f:
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue

            accession = entry.get("accession")
            if not accession:
                continue

            genes = entry.get("genes", [])
            if not genes:
                continue

            # Validate CDS lengths before creating GFF. This is necessary because
            # nextclade requires that the lenght of CDS to be a multiple of 3
            has_invalid_cds = False
            for gene in genes:
                gene_name = gene.get("name", "unknown")
                cds_data = gene.get("cds", [])

                if cds_data:  # Only validate if there is CDS data
                    # Group CDSs by name some CDSs have the same name e.g hypothethical
                    # protein and it breaks the nextclade run
                    cds_groups = {}
                    for cds in cds_data:
                        cds_name = cds.get("name", gene_name)
                        cds_name = clean_cds_name(cds_name)

                        if cds_name not in cds_groups:
                            cds_groups[cds_name] = []
                        cds_groups[cds_name].append(cds)

                    # Validate each CDS group
                    for cds_name, cds_list in cds_groups.items():
                        all_ranges = []
                        for cds in cds_list:
                            nucleotide = cds.get("nucleotide", {})
                            all_ranges.extend(nucleotide.get("range", []))

                        if all_ranges:
                            try:
                                # Calculate total length for this CDS group
                                min_start = min(int(r.get("begin")) for r in all_ranges)
                                max_end = max(int(r.get("end")) for r in all_ranges)
                                total_length = max_end - min_start + 1

                                if total_length % 3 != 0:
                                    has_invalid_cds = True
                                    break
                            except (ValueError, TypeError):
                                # Skip groups with invalid range data
                                continue

                if has_invalid_cds:
                    break

            # If CDS validation failed, add to invalid list and skip
            if has_invalid_cds:
                invalid_accessions.append(accession)
                continue

            # Validate genes without CDS data (genome length must be divisible by 3)
            # here we use the complete gene as a CDS
            for gene in genes:
                cds_data = gene.get("cds", [])
                if not cds_data:  # Gene has no CDS information
                    genome_length = lengths.get(accession)
                    if genome_length and genome_length % 3 != 0:
                        invalid_accessions.append(accession)
                        has_invalid_cds = True
                        break

            # Skip this accession if it was marked as invalid
            if has_invalid_cds:
                continue

            gff_file_path = os.path.join(output_dir, f"{accession}.gff")
            with open(gff_file_path, "w") as gff:
                gff.write("##gff-version 3\n")

                seq_len = lengths.get(accession)
                if seq_len:
                    gff.write(f"##sequence-region\t{accession}\t1\t{seq_len}\n")

                # Track used gene IDs to ensure uniqueness
                used_gene_ids = set()

                for gene in genes:
                    gene_name = gene.get("name", "unknown")
                    gene_id = gene.get("geneId")

                    # Check if gene has CDS information
                    cds_data = gene.get("cds", [])

                    # If no CDS data is available, create a gene entry using genome length
                    if not cds_data:
                        # Clean up gene name
                        clean_gene_name = (
                            re.sub(r"[ ,/]+", "_", gene_name).capitalize().strip("'")
                        )

                        # Generate unique ID for the gene
                        unique_gene_id = clean_gene_name
                        counter = 1
                        while unique_gene_id in used_gene_ids:
                            counter += 1
                            unique_gene_id = f"{clean_gene_name}_{counter}"
                        used_gene_ids.add(unique_gene_id)

                        # Use genome length from FASTA for gene coordinates
                        genome_length = lengths.get(accession)
                        if not genome_length:
                            continue

                        # Write gene line spanning entire genome
                        unique_gene_id = re.sub(r" ", "_", unique_gene_id).capitalize()
                        gene_attributes = f"gene_name={unique_gene_id}"
                        gff.write(
                            f"{accession}\tRefSeq\tgene\t1\t{genome_length}\t.\t+\t.\t{gene_attributes}\n"
                        )
                        continue

                    # Group CDSs by name to handle split genes
                    cds_groups = {}

                    for cds in cds_data:
                        cds_name = cds.get("name", gene_name)
                        cds_name = clean_cds_name(cds_name)

                        if cds_name not in cds_groups:
                            cds_groups[cds_name] = []

                        cds_groups[cds_name].append(cds)

                    # Process each group of CDSs as a single gene entity
                    for cds_name, cds_list in cds_groups.items():
                        # Generate unique ID for the gene
                        unique_gene_id = cds_name
                        counter = 1
                        while unique_gene_id in used_gene_ids:
                            counter += 1
                            unique_gene_id = f"{cds_name}_{counter}"
                        used_gene_ids.add(unique_gene_id)

                        # Collect all ranges to find gene boundaries
                        all_ranges = []
                        for cds in cds_list:
                            nucleotide = cds.get("nucleotide", {})
                            all_ranges.extend(nucleotide.get("range", []))

                        if not all_ranges:
                            continue

                        # Calculate global start/end for the gene
                        min_start = min(int(r.get("begin")) for r in all_ranges)
                        max_end = max(int(r.get("end")) for r in all_ranges)

                        # Determine strand (assume all segments are on same strand)
                        orientation = all_ranges[0].get("orientation", "plus")
                        strand = "-" if orientation == "minus" else "+"

                        # Write Parent Gene Feature
                        unique_gene_id = re.sub(r" ", "_", unique_gene_id).capitalize()
                        gene_attributes = f"gene_name={unique_gene_id}"
                        gff.write(
                            f"{accession}\tRefSeq\tgene\t{min_start}\t{max_end}\t.\t{strand}\t.\t{gene_attributes}\n"
                        )

    # Write not_included.txt file
    not_included_path = os.path.join(output_dir, "not_included.txt")
    with open(not_included_path, "w") as f:
        for accession in invalid_accessions:
            f.write(f"{accession}\n")

    print(
        f"Processed GFF files. {len(invalid_accessions)} accessions excluded due to invalid CDS lengths.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert NCBI Datasets JSONL to GFF3.")
    parser.add_argument("--input", required=True, help="Input JSONL file")
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for GFF files"
    )
    parser.add_argument("--fasta", help="Input FASTA file for sequence lengths")
    args = parser.parse_args()
    jsonl_to_gff(args.input, args.output_dir, args.fasta)
