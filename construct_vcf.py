"""
    Construct a VCF file from a BED file containing the chromosome, start, and end positions
    to be tested with VEP or visualised in IGV Reports.
    Each row of the BED file is is converted into a row in the VCF for
    start, middle, and end coordinates using the reference genome.

    TODO:
    - Add argparse to allow user to specify input and output files
    - Add argparse to allow user to specify reference genome
    - Check right coordinates are being used for start/end positions
    - Write docstrings
    - Write tests
    - Add logging

    Returns
    -------
    vcf_file written to test.vcf.
"""

import pandas as pd
import pysam


def fetch_nucleotides(row, reference_path):
    """
    This function fetches the nucleotide sequence from a reference database using the
    chromosome, start, and end positions in a row of a dataframe.

    Parameters:
    row (pandas.Series): A row of a pandas dataframe containing the chromosome,
    start, and end positions of the nucleotide sequence to fetch.

    Returns:
    pandas.DataFrame: A dataframe containing the variant information for the start,
    middle, and end positions in the format specified by a VCF file.
    """

    # Define the chromosome number and NCBI identifier
    chr = str(row['chr'])  #.strip('chr')
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
    role_in_cancer = row['role_in_cancer']
    INFO_col = f"DP=268;ROLE={role_in_cancer}"
    FORMAT = f"GT:GQ:AD:DP:VF:NL:SB:NC:US:AQ:LQ"
    SAMPLE_col = f"1/1:0:0,114:114:1:65:-100:0.2192:27,12,16,14,23,22,27,12,16,14,23,22:100:100"

    # Variant dictionary to convert ref to alt
    variant_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'A'}

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
            'QUAL': ['.', '.', '.'],
            'FILTER': ['.', '.', '.'],
            'INFO': [INFO_col, INFO_col, INFO_col],
            'FORMAT': [FORMAT, FORMAT, FORMAT],
            'SAMPLE_col': [SAMPLE_col, SAMPLE_col, SAMPLE_col]
            }

    vcf_df = pd.DataFrame(data)
    return vcf_df


def main():
    """
    _summary_
    """
    header = ['chr', 'start', 'end', 'role_in_cancer', 'gene']
    bed_file = pd.read_csv('output_hg19_overlap_test5.bed', names=header, sep='\t')

    # # Define the variant dictionary to convert ref to alt
    # variant_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'A'}

    # Define the output dataframe columns
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'test-123456-1-DNA-egg6.bam']

    # Initialize an empty list to store the rows of the output dataframe
    rows = []
    reference_path = "/home/rswilson1/Documents/MSC_dissertation/MELT_example_data/hs37d5.fa"
    # Loop through each row of the input dataframe
    counter = 0
    for i, row in bed_file.iterrows():
        # Fetch the nucleotide sequence from NCBI using the fetch_nucleotides function
        seq_df = fetch_nucleotides(row, reference_path)
        # Iterate through the rows of the seq_df to create the output dataframe rows
        for j, seq_row in seq_df.iterrows():
            # Define the output dataframe row as a dictionary
            counter += 1
            output_row = {
                '#CHROM': seq_row['#CHROM'],
                'POS': seq_row['POS'],
                'ID': counter,
                'REF': str(seq_row['REF']),
                'ALT': seq_row['ALT'],
                'QUAL': 100,
                'FILTER': 'PASS',
                'INFO': seq_row['INFO'],
                'FORMAT': seq_row['FORMAT'],
                'test-123456-1-DNA-egg6.bam': str(seq_row['SAMPLE_col'])
            }

            # Append the output dataframe row to the rows list
            rows.append(output_row)

    # Create the output dataframe from the rows list
    output_df = pd.DataFrame(rows, columns=columns)
    # Print the output dataframe
    print(output_df.head(10))

    # saving as tsv file
    output_df.to_csv('test.vcf', sep="\t", index=False)

    # test_vcf_cancer_TSG_V4.vcf -- 0-based coordinates


if __name__ == "__main__":
    print("Running construct_vcf.py...")
    main()
