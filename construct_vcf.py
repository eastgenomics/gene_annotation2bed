"""
    Construct a VCF file from a BED file containing the chromosome,
    start, and end positions to be tested with VEP or visualised in IGV Reports.
    Each row of the BED file is is converted into a row in the VCF for
    start, middle, and end coordinates using the reference genome.

    TODO:
    - Add argparse to allow user to specify input and output files
    - Add argparse to allow user to specify reference genome
    - Check right coordinates are being used for start/end positions
    - Write docstrings
    - Write tests
    - Add logging
    Example cmd for quick ref:
    construct_vcf.py -fasta hs37d5.fa -b output_hg38_general_test1.bed
    Returns
    -------
    vcf_file written to test.vcf.
"""

import pandas as pd
import pysam
import argparse


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Parameters
    ----------
    None

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(description="GFF Processing Script")
    parser.add_argument(
        "-b", "--bed_file", help="BED file to be converted to VCF",
        required=True
    )
    parser.add_argument(
        "-ref", "--reference_genome", help="Reference genome (hg19/hg38)",
        required=False, choices=('hg19', 'hg38')
    )
    parser.add_argument(
        "-fasta",
        "--reference_file",
        help="Path to Reference genome fasta file for igv_reports",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help="Path to output file",
        required=True,
    )

    args = parser.parse_args()
    return args


def fetch_nucleotides(row, reference_path, variant_dict):
    """
    Fetches the nucleotide sequence from a reference database using the chromosome,
    start, and end positions in a row of a dataframe.

    Parameters
    ----------
    row (pandas.Series): A row of a pandas dataframe containing the chromosome,
        start, and end positions of the nucleotide sequence to fetch.

    reference_path (str): The path to the reference genome database.

    variant_dict (dict): A dictionary containing the reference nucleotide as the key
        and the corresponding alternate nucleotide as the value.

    Returns
    -------
    pandas.DataFrame: A dataframe containing the variant information for the start,
        middle, and end positions in the format specified by a VCF file.
    """

    # Define the chromosome number and NCBI identifier
    chr = str(row['chr']).strip('chr')
    if chr == 'X':
        ncbi_chr = 'X'
    elif chr == 'Y':
        ncbi_chr = 'Y'
    elif int(chr) < 10 or int(chr) >= 10:
        ncbi_chr = int(chr)
    else:
        # If the chromosome is not valid, raise an error
        raise ValueError(f"Error: {chr} is not a valid chromosome")

    # Define the start, middle, and end positions
    start = int(row['start']) + 1
    end = int(row['end'])
    middle = int((start + end) / 2)
    INFO_field = row['info']
    INFO_col = f"DP=268;{INFO_field}"
    FORMAT = f"GT:GQ:AD:DP:VF:NL:SB:NC:US:AQ:LQ"
    SAMPLE_col = ("1/1:0:0,114:114:1:65:-100:0.2192:"
                  "27,12,16,14,23,22,27,12,16,14,23,22:100:100"
                )

    # Fetch the start nucleotide from reference
    print(f"Fetching nucleotide sequence for {ncbi_chr}:{start}-{end}...")
    start_seq = pysam.faidx(reference_path,
                            f"{ncbi_chr}:{start}-{start}")
    middle_seq = pysam.faidx(reference_path,
                             f"{ncbi_chr}:{middle}-{middle}")
    end_seq = pysam.faidx(reference_path,
                          f"{ncbi_chr}:{end}-{end}")
    start_nuc = start_seq.splitlines()[1]
    middle_nuc = middle_seq.splitlines()[1]
    end_nuc = end_seq.splitlines()[1]
    # # Convert the ref nucleotides to the alternate nucleotides
    start_nuc_variant = variant_dict[start_nuc]
    middle_nuc_variant = variant_dict[middle_nuc]
    end_nuc_variant = variant_dict[end_nuc]

    # Define the VCF dataframe with the variant information for the start,
    # middle, and end positions
    data = {'#CHROM': [ncbi_chr, ncbi_chr, ncbi_chr],
            'POS': [start, middle, end],
            'ID': ['.', '.', '.'],
            'REF': [start_nuc, middle_nuc, end_nuc],
            'ALT': [start_nuc_variant, middle_nuc_variant, end_nuc_variant],
            'QUAL': [100, 100, 100],
            'FILTER': ['PASS', 'PASS', 'PASS'],
            'INFO': [INFO_col, INFO_col, INFO_col],
            'FORMAT': [FORMAT, FORMAT, FORMAT],
            'test-123456-1-DNA-egg6.bam': [SAMPLE_col, SAMPLE_col, SAMPLE_col]
            }

    vcf_df = pd.DataFrame(data)

    return vcf_df


def main():
    """
    Main function to run the script
    """
    args = parse_args()
    header = ['chr', 'start', 'end', 'info', 'gene']
    bed_file_path = args.bed_file
    bed_file = pd.read_csv(bed_file_path, names=header, sep='\t')

    # Define the variant dictionary to convert ref to alt
    variant_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'A'}

    # Define the output dataframe columns
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
               'FILTER', 'INFO', 'FORMAT', 'test-123456-1-DNA-egg6.bam']

    # Initialize an empty list to store the rows of the output dataframe
    reference_path = args.reference_file

    output_df = pd.DataFrame(columns=columns)
    for i, row in bed_file.iterrows():
        # Fetch the nucleotide sequence from NCBI
        # using the fetch_nucleotides function
        seq_df = fetch_nucleotides(row, reference_path, variant_dict)
        output_df = pd.concat([output_df, seq_df])

    # saving as tsv file
    output_df.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    print("Running construct_vcf.py...")
    main()
