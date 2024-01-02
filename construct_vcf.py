"""
    Construct a VCF file from a BED file containing the chromosome,
    start, and end positions to be tested with VEP or visualised in IGV Reports.
    Each row of the BED file is is converted into a row in the VCF for
    start, middle, and end coordinates using the reference genome.

    TODO:
    - Check right coordinates are being used for start/end positions
    - Write tests
    - Add logging

    Example cmd for quick ref:
    construct_vcf.py -fasta /path/to/hs37d5.fa -b output_hg38_general_test1.bed

    Returns
    -------
    None

    Outputs
    -------
    vcf_file written with variants for start, middle, and end positions of
    each row in the BED file.
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


class ConstructVCF():
    # Define the variant dictionary to convert ref to alt nucleotides
    variant_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'A'}

    def __init__(self, bed_file_path, reference_path, output_file):
        self.bed_file_path = bed_file_path
        self.reference_path = reference_path
        self.output_file = output_file

    def fetch_nucleotides(self, row, reference_path):
        """
        Fetches the nucleotide sequence from a reference database using the chromosome,
        start, and end positions in a row of a dataframe.

        Parameters
        ----------
        row (pandas.Series): A row of a pandas dataframe containing the chromosome,
            start, and end positions of the nucleotide sequence to fetch.

        reference_path (str): The path to the reference genome database.

        Returns
        -------
        pandas.DataFrame: A dataframe containing the variant information for the start,
            middle, and end positions in the format specified by a VCF file.
        """
        # Define the chromosome number and NCBI identifier
        chr_str = str(row['chr']).replace('chr', '')
        if chr_str.upper() == 'X':
            ncbi_chr = 'X'
        elif chr_str.upper() == 'Y':
            ncbi_chr = 'Y'
        else:
            try:
                ncbi_chr = int(chr_str)
                assert ncbi_chr <= 22
            except:
                raise ValueError(f"Error: {chr_str} is not a valid chromosome")

        # Define the start, middle, and end positions
        start = (int(row['start']) + 1)
        end = int(row['end'])
        middle = int((int(row['start']) + end) / 2)

        # Define the content for VCF columns to make compatible with VEP
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

        # Convert the ref nucleotides to the alternate nucleotides
        start_nuc_variant = ConstructVCF.variant_dict[start_nuc]
        middle_nuc_variant = ConstructVCF.variant_dict[middle_nuc]
        end_nuc_variant = ConstructVCF.variant_dict[end_nuc]

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

    def convert_bed_to_vcf(self):
        """
        Main logic to convert bed file into vcf file.
        Using the fetch_nucleotides function.
        This outputs a VCF (TSV file).
        Parameters
        ----------
        self:
            instance of class.

        Class Parameters
        ----------
        bed_file_path : str
            path to bed file.
        reference_path : str
            path to reference genome.
        output_file : str
            path to write output file.

        Returns
        -------
        output_df
            vcf file as a pandas dataframe.

        Outputs
        -------
        Vcf File
            VCF file containing the variants for the start, middle, and end
            of each gene in bed file provided.
        """

        header = ['chr', 'start', 'end', 'info', 'gene']

        bed_file = pd.read_csv(self.bed_file_path, names=header, sep='\t')

        # Define the output dataframe columns
        columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                   'FILTER', 'INFO', 'FORMAT', 'test-123456-1-DNA-egg6.bam']

        output_df = pd.DataFrame(columns=columns)

        for i, row in bed_file.iterrows():
            # Fetch the nucleotide sequence from NCBI
            # using the fetch_nucleotides function
            seq_df = self.fetch_nucleotides(row, self.reference_path)
            output_df = pd.concat([output_df, seq_df])

        # sequences = bed_file.apply(lambda x: self.fetch_nucleotides(x, self.reference_path), axis=1)
        # output_df = pd.concat([output_df, sequences])

        # saving as tsv file
        # Define the columns and their desired data types
        columns_to_convert = {
            '#CHROM': 'int64',
            'POS': 'int64',
            'ID': 'str',
            'REF': 'str',
            'ALT': 'str',
            'QUAL': 'int64',
            'FILTER': 'str',
            'INFO': 'str',
            'FORMAT': 'str',
            'test-123456-1-DNA-egg6.bam': 'str'
        }

        # Convert the data types of the specified columns
        output_df = output_df.astype(columns_to_convert)
        # reset index
        output_df.reset_index(drop=True, inplace=True)
        output_df.to_csv(self.output_file, sep="\t", index=False)
        return output_df


def main():
    """
    Main function to run the script
    """

    args = parse_args()
    bed_file_path = args.bed_file

    # Gather file paths
    reference_path = args.reference_file
    output_file = args.output_file
    vcf_object = ConstructVCF(bed_file_path, reference_path, output_file)
    vcf_object.convert_bed_to_vcf()


if __name__ == "__main__":
    print("Running construct_vcf.py...")
    main()
