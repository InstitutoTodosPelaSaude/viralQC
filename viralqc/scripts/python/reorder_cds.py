import argparse
from pathlib import Path
import csv


def read_gff_gene_order(gff_file: Path) -> list[str]:
    """
    Read GFF file and extract gene names ordered by their genomic start position.

    Args:
        gff_file: Path to the GFF file

    Returns:
        List of gene names ordered by their start position in the genome
    """
    genes = []
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 9 and fields[2] == "gene":
                start_pos = int(fields[3])
                attributes = fields[8]
                for attr in attributes.split(";"):
                    if "gene_name=" in attr:
                        gene_name = attr.split("gene_name=")[1].strip()
                        genes.append((start_pos, gene_name))
                        break
    
    # Sort by start position and extract gene names
    genes.sort(key=lambda x: x[0])
    return [gene_name for _, gene_name in genes]


def parse_cds_coverage(cds_coverage_str: str) -> dict:
    """
    Parse cdsCoverage string into a dictionary.

    Args:
        cds_coverage_str: String like "2K:1,C:1,E:1,NS1:0.997"

    Returns:
        Dictionary mapping gene name to coverage value
    """
    if not cds_coverage_str or cds_coverage_str.strip() == "":
        return {}

    cds_dict = {}
    for pair in cds_coverage_str.split(","):
        if ":" in pair:
            gene, coverage = pair.split(":", 1)
            cds_dict[gene.strip()] = coverage.strip()
    return cds_dict


def reorder_cds_coverage(cds_dict: dict, gene_order: list[str]) -> str:
    """
    Reorder cdsCoverage based on gene order from GFF.
    Add missing genes with 0.0 coverage.

    Args:
        cds_dict: Dictionary of gene:coverage pairs
        gene_order: Ordered list of gene names from GFF

    Returns:
        Reordered cdsCoverage string
    """
    ordered_pairs = []
    for gene in gene_order:
        if gene in cds_dict:
            ordered_pairs.append(f"{gene}:{cds_dict[gene]}")
        else:
            ordered_pairs.append(f"{gene}:0.0")

    return ",".join(ordered_pairs)


def process_nextclade_tsv(tsv_file: Path, gff_file: Path, output_file: Path) -> None:
    """
    Process nextclade TSV file and reorder cdsCoverage column.

    Args:
        tsv_file: Path to input nextclade TSV file
        gff_file: Path to GFF annotation file
        output_file: Path to output TSV file
    """
    gene_order = read_gff_gene_order(gff_file)

    with open(tsv_file, "r", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = reader.fieldnames

        rows = []
        for row in reader:
            if gene_order and "cdsCoverage" in fieldnames:
                cds_coverage = row.get("cdsCoverage", "")
                cds_dict = parse_cds_coverage(cds_coverage)
                row["cdsCoverage"] = reorder_cds_coverage(cds_dict, gene_order)
            rows.append(row)

    with open(output_file, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reorder cdsCoverage column in nextclade TSV based on GFF gene order"
    )
    parser.add_argument(
        "--nextclade-tsv", type=Path, required=True, help="Path to nextclade TSV file"
    )
    parser.add_argument(
        "--gff", type=Path, required=True, help="Path to GFF annotation file"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Path to output TSV file"
    )

    args = parser.parse_args()

    process_nextclade_tsv(args.nextclade_tsv, args.gff, args.output)
