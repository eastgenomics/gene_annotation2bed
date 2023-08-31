"""
This script takes a GFF3 file and an annotation file,
producing a BED file for annotation of relevant transcripts.
Example cmd:
/home/rswilson1/anaconda3/envs/Annotation_app/bin/python
/home/rswilson1/Documents/Programming_project/gene_annotation2bed.py -gff "GCF_000001405.25_GRCh37.p13_genomic.gff"
-ig cancerGeneList_test.tsv -ref "hg19" -f 5
--assembly_summary "GCF_000001405.25_GRCh37.p13_assembly_report.txt"
-o "test6"
"""

import argparse
import argcomplete
import pandas as pd
import numpy as np
import itertools
pd.options.mode.chained_assignment = None  # default='warn'
import gff2pandas as gffpd
import igv_report as igv

def parse_gff(gff_file):
    """
    Summary: Import GFF3 file and convert to pandas DataFrame.

    Parameters
    ----------
        gff_file : gff2pandas object
            GFF object which contains the df and header.

    Returns
    -------
    transcripts_df : pandas DataFrame
        DataFrame containing the all the 'NM_' prefixed
        transcripts from the GFF3 file.
    """
    transcripts_gff = gffpd.read_gff3(gff_file)
    info = transcripts_gff.stats_dic()
    gff_df = transcripts_gff.attributes_to_columns()
    # drop columns that are not needed to reduce memory footprint
    gff_df = gff_df.drop(['Gap', 'Is_circular', 'Name', 'Note', 'Parent',
                          'Target', 'anticodon', 'assembly_bases_aln',
                          'assembly_bases_seq', 'bit_score', 'blast_aligner',
                          'blast_score', 'bound_moiety', 'chromosome',
                          'codons', 'common_component', 'consensus_splices',
                          'country', 'description', 'direction', 'e_value',
                          'end_range', 'exception', 'exon_identity',
                          'exon_number', 'experiment', 'feat_class',
                          'filter_score', 'for_remapping', 'function',
                          'gap_count', 'gene_biotype', 'gene_synonym',
                          'genome', 'hsp_percent_coverage', 'identity',
                          'idty', 'inference', 'inversion_merge_aligner',
                          'isolation-source', 'lxr_locAcc_currStat_120',
                          'lxr_locAcc_currStat_35', 'map', 'matchable_bases',
                          'matched_bases', 'matches', 'merge_aligner',
                          'mobile_element_type', 'mol_type',
                          'not_for_annotation', 'note', 'num_ident',
                          'num_mismatch', 'number', 'partial', 'pct_coverage',
                          'pct_coverage_hiqual', 'pct_identity_gap',
                          'pct_identity_gapopen_only', 'pct_identity_ungap',
                          'product', 'product_coverage', 'protein_id',
                          'pseudo', 'rank', 'recombination_class',
                          'regulatory_class', 'rpt_family', 'rpt_type',
                          'rpt_unit_range', 'rpt_unit_seq', 'satellite',
                          'splices', 'standard_name', 'start_range', 'tag',
                          'tissue-type', 'transl_except', 'transl_table',
                          'weighted_identity'], axis=1)

    # Apply extract_hgnc_id function to create 'hgnc_id' column
    gff_df['hgnc_id'] = gff_df['Dbxref'].apply(extract_hgnc_id)

    # set dtype for each column to reduce memory footprint
    dtype_mapping = {
        'ID': 'str',
        'transcript_id': 'str',
        'hgnc_id': 'Int64',
    }

    gff_df = gff_df.astype(dtype_mapping)

    # Filter GFF DataFrame to select entries with 'NM' type
    transcripts_df = gff_df[gff_df['transcript_id'].str.startswith('NM_')]

    return transcripts_df


# Function to extract HGNC ID
def extract_hgnc_id(dbxref_str):
    """
    Wrapper function to extract HGNC ID from a string of dbxrefs.

    Parameters
    ----------
    dbxref_str : str
        various ids separated by commas

    Returns
    -------
    int
        HGNC ID as an integer i.e. 427 for HGNC:427.
    """
    if dbxref_str is None:
        return None
    parts = dbxref_str.split(',')
    for part in parts:
        if 'HGNC:' in part:
            return int(part.split(':')[-1])
        elif 'hgnc:' in part:
            return int(part.split(':')[-1])
        elif 'HGNC_' in part:
            return int(part.split('_')[-1])
    return None


def read_assembly_mapping(assembly_file):
    """
    Reads in the associated assembly file and returns a dictionary mapping
    to find chromosome for each refseq accession.

    Parameters
    ----------
    assembly_file : tsv
        found at: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/

    Returns
    -------
    dictionary
        mapping of refseq accession to chromosome
    """
    accession_to_chromosome = {}
    assembly_df = pd.read_csv(assembly_file, sep='\t', comment='#', header=None)
    assembly_df = assembly_df.dropna()  # Drop rows with missing values
    for _, row in assembly_df.iterrows():
        chromosome = row[2]
        accession = row[6]
        if chromosome.startswith('na'):
            continue
        accession_to_chromosome[accession] = chromosome
    return accession_to_chromosome


def map_accession_to_chromosome(accession, accession_to_chromosome):
    """
    Simple mapping function to find chromosome for a given refseq accession.
    Calls the accession_to_chromosome dictionary and extracts the chromosome.

    Parameters
    ----------
    accession: str
    str of the refseq accession

    accession_to_chromosome: dictionary
    dictionary mapping of refseq accession to chromosome

    Returns
    -------
    str value for the corresponding chromosome for the accession key.
    Or if not present in the dictionary, returns "Unknown - {accession}"
    """
    return accession_to_chromosome.get(accession, f"Unknown - {accession}")


def parse_pickle(pickle_file):
    """
    Parses a pickle file and returns a DataFrame of transcripts.

    Parameters
    ----------
    pickle_file : pkl
        pickle file of a GFF DataFrame once parsed
        with columns from attributes_to_columns

    Returns
    -------
    transcripts_df: dataframe
        dataframe of transcripts with columns for attributes.
        Contains only transcripts with NM_ prefix.
    """
    gff_df = pd.read_pickle(f"./{pickle_file}")
    transcripts_df = gff_df[gff_df['transcript_id'].fillna('').str.startswith('NM_')]
    return transcripts_df


def merge_overlapping(bed_df):
    """
    Function to merge overlapping regions in a bed file by annotation.

    Parameters
    ----------
    bed_df : dataframe
        bed file with columns: seq_id, start_flank,
            end_flank, hgnc_id, annotation, gene, chromosome

    Returns
    -------
    merged_df: dataframe
        dataframe of merged rows with columns: chromosome, start, end, annotation
    """
    # Sort by chromosome, start, and end
    # This makes sure that overlapping regions are next to each other.

    bed_df = bed_df.sort_values(by=["annotation", "chromosome", "start_flank", "end_flank"])
    # Sort by first annotation then chromosome, start, and end.
    merged_rows = []

    current_row = bed_df.iloc[0]
    for _, row in bed_df.iterrows():
        if row['annotation'] != current_row['annotation']:
            merged_rows.append(current_row) # Append the merged row
            current_row = row # Start a new potential merged row
            # Only rows with same annotation are merged
        if row["chromosome"] != current_row["chromosome"]:
            merged_rows.append(current_row)
            current_row = row
            # Only rows with same chromosome are merged.
        if row["start_flank"] <= current_row["end_flank"]:
            current_row["end_flank"] = max(current_row["end_flank"], row["end_flank"])
            # Extend the end if overlapping
        else:
            merged_rows.append(current_row)
            current_row = row

    merged_rows.append(current_row)  # Append the last merged row
    merged_df = pd.DataFrame(merged_rows)
    return merged_df


def config_igv_report(args):
    """
    Function to call igv report script with the correct parameters.
    Generates an IGV html report using generic handling.

    Parameters
    ----------
    args : argeparse object
        argeparse object with the following attributes:
        reference_genome, output_file_suffix, gff_file/pickle_file,
        annotation_file/transcript_file, assembly_file, and flanking.

    Returns
    -------
    None
    """
    # assign vars.
    bed_file = f"output_{args.reference_genome}_{args.output_file_suffix}.maf"
    genome = args.reference_genome
    info_columns = []
    title = f"{args.output_file_suffix}_report"
    output_file = f"{title}.html"
    print("Creating IGV report...")
    print(f"Bed file: {bed_file}, Genome: {genome}, Info columns: {info_columns}, Title: {title}, Output: {output_file}")
    igv.create_igv_report(bed_file, genome, info_columns, title, output_file)
    return "IGV report created successfully!"


def main():
    parser = argparse.ArgumentParser(description="GFF Processing Script")
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('-gff', "--gff_file", help="Path to GFF file")
    group1.add_argument('-pkl', "--pickle", help="Import gff as pickle file")

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('-ig', "--annotation_file", help="Path to the annotation file (TSV)")
    group2.add_argument('-it', "--transcript_file", help="Path to transcript annotation file")

    parser.add_argument('-o', "--output_file_suffix", help="Output file suffix", required=True)
    parser.add_argument('-ref', "--reference_genome", help="Reference genome (GrCh37/38)", required=True)
    parser.add_argument('-f', "--flanking", type=int, help="Flanking size", required=True)
    parser.add_argument('--assembly_summary', help="Path to assembly summary file", required=True)
    # parser.add_argument('--report_name', help="Name for report")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # read in pickle file if provided
    if args.pickle is not None:
        transcripts_df = parse_pickle(args.pickle)
        print("Parsed pickle file")
    else:
        # Parse gff file
        transcripts_df = parse_gff(args.gff_file)

    # Read the annotation file into a pandas DataFrame
    if args.annotation_file is not None:
        annotation_df = pd.read_csv(args.annotation_file, sep="\t")
    # Read the transcript annotation file
    elif args.transcript_file is not None:
        annotation_df = pd.read_csv(args.transcript_file, sep="\t")
    else:
        SystemExit("Please provide an annotation file")

    # Merge NM entries with matching HGNC IDs
    print("Merging annotation and gff dataframes")
    merged_df = transcripts_df.merge(annotation_df, left_on="hgnc_id",
                                     right_on="hgnc_id", how="inner")
    # Create BED file with flanking regions
    print("Creating BED file")
    print("Adding flanking regions")
    merged_df["start_flank"] = merged_df["start"] - args.flanking
    merged_df["end_flank"] = merged_df["end"] + args.flanking
    bed_columns = ["seq_id", "start_flank", "end_flank", "hgnc_id", "annotation", "gene"]
    bed_df = merged_df[bed_columns]

    # Extract chromosome from seqid and create the 'chromosome' column
    accession_to_chromosome = read_assembly_mapping(args.assembly_summary)
    # Add a new column 'chromosome' by mapping accession to chromosome identifier
    bed_df.loc[:, "chromosome"] = bed_df["seq_id"].apply(lambda x: map_accession_to_chromosome(x, accession_to_chromosome))
    print(f"Summary of BED file df before collapsing \n {bed_df.head()}")

    # # Merge overlapping entries
    collapsed_df = merge_overlapping(bed_df)
    print(f"Summary of BED file df after collapsing \n {collapsed_df.head()}")
    # Reorder the columns to match the BED format
    cols = ['chromosome', 'start_flank', 'end_flank', 'annotation', 'gene']
    collapsed_df = collapsed_df[cols]
    # Rename columns
    new_column_names = {
        "start_flank": "start",
        "end_flank": "end"
    }
    collapsed_df.rename(columns=new_column_names, inplace=True)

    # Write the collapsed data to an output file
    output_file_name = f"output_{args.reference_genome}" \
                       f"_{args.output_file_suffix}.maf"
    collapsed_df.to_csv(output_file_name, sep="\t",
                        header=False, index=False)

    # Create an IGV report
    report = config_igv_report(args)

if __name__ == "__main__":
    main()
