"""
This script takes a GFF3 file and an annotation file,
producing a BED file for annotation of relevant transcripts.
Example cmd:
/home/rswilson1/anaconda3/envs/Annotation_app/bin/python
/home/rswilson1/Documents/Programming_project/gene_annotation2bed.py -gff "GCF_000001405.25_GRCh37.p13_genomic.gff"
-ig cancerGeneList_test.tsv -ref "Grch37" -f 5
--assembly_summary "GCF_000001405.25_GRCh37.p13_assembly_report.txt"
-o "test6"
"""

import argparse
import argcomplete
import pandas as pd
import numpy as np
import itertools
# import pybedtools # collapsed = bed_file.merge(c=4, o='collapse') or bed_file.merge(c=4, o='distinct')
#import gffpandas.gffpandas as gffpd
pd.options.mode.chained_assignment = None  # default='warn'
import gff2pandas as gffpd
import igv_report as igv

def parse_gff(gff_file):
    """
    Summary: Import GFF3 file and convert to pandas DataFrame.
    """
    transcripts_gff = gffpd.read_gff3(gff_file)
    info = transcripts_gff.stats_dic()
    print(info)
    gff_df = transcripts_gff.attributes_to_columns()
    print(gff_df.head())
    print(gff_df.columns)
    print(gff_df.memory_usage(deep=True).sum() / 1024 / 1024)
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
    print(gff_df.dtypes)

    # Filter GFF DataFrame to select entries with 'NM' type
    transcripts_df = gff_df[gff_df['transcript_id'].str.startswith('NM_')]
    # transcripts_df = transcripts_df.reset_index(drop=True)
    print(transcripts_df.head())
    print(transcripts_df.memory_usage(deep=True).sum() / 1024 / 1024)
    return transcripts_df


def attributes_to_columns(attribute_df) -> pd.DataFrame:
    """Saving each attribute-tag to a single column.
    Attribute column will be split by the tags in the single columns.
    For this method only a pandas DataFrame and not a Gff3DataFrame
    will be returned. Therefore, this data frame can not be saved as
    gff3 file.
    :return: pandas dataframe, whereby the attribute column of the gff3
            file are splitted into the different attribute tags
    :rtype: pandas DataFrame
    """
    df_attributes = attribute_df.loc[:, "seqid":"attributes"]
    attribute_df["attr_dic"] = attribute_df.attributes.apply(
        lambda attributes: dict(
            [
                key_value_pair.split(sep="=", maxsplit=1)
                for key_value_pair in attributes.split(";")
            ]
        )
    )
    attribute_df["attr_dic_keys"] = attribute_df["attr_dic"].apply(
        lambda attr_dic: list(attr_dic.keys())
    )
    merged_attribute_list = list(
        itertools.chain.from_iterable(attribute_df["attr_dic_keys"])
    )
    nonredundant_list = sorted(list(set(merged_attribute_list)))
    for atr in nonredundant_list:
        df_attributes[atr] = attribute_df["attr_dic"].apply(
            lambda attr_dic: attr_dic.get(atr)
        )
    df_attributes.drop(columns=["attributes"], inplace=True)
    return df_attributes


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
    """
    return accession_to_chromosome.get(accession, f"Unknown - {accession}")


def parse_pickle(pickle_file):
    gff_df = pd.read_pickle(f"./{pickle_file}")
    transcripts_df = gff_df[gff_df['transcript_id'].fillna('').str.startswith('NM_')]
    return transcripts_df


def merge_overlapping(bed_df):
    """
    Function to merge overlapping regions in a bed file.

    Parameters
    ----------
    bed_df : _type_
        bed file with columns: chromosome, start, end

    Returns
    -------
    _type_
        _description_
    """
    bed_df = bed_df.sort_values(by=["chromosome", "start_flank", "end_flank"])  # Sort by chromosome, start, and end
    merged_rows = []

    current_row = bed_df.iloc[0]
    for _, row in bed_df.iterrows():
        if row["chromosome"] != current_row["chromosome"]:
            merged_rows.append(current_row)
            current_row = row
        if row["start_flank"] <= current_row["end_flank"]:
            current_row["end_flank"] = max(current_row["end_flank"], row["end_flank"])  # Extend the end if overlapping
        else:
            merged_rows.append(current_row)  # Append the merged row
            current_row = row  # Start a new potential merged row

    merged_rows.append(current_row)  # Append the last merged row

    return pd.DataFrame(merged_rows)


def config_igv_report(args):
    # assign vars.
    bed_file = f"output_{args.reference_genome}_{args.output_file_suffix}.bed"
    genome = args.reference_genome #"hg38"
    info_columns = [] #["TCGA", "GTEx", "variant_name"]
    title = f"TEST{args.output_file_suffix}" #Sample A
    output = "example_igv_report.html"
    print("Creating IGV report...")
    print(f"Bed file: {bed_file}, Genome: {genome}, Info columns: {info_columns}, Title: {title}, Output: {output}")
    igv.create_igv_report(bed_file, genome, info_columns, title, output)
    return "IGV report created successfully!"


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


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
        print(annotation_df.dtypes)
    # Read the transcript annotation file
    elif args.transcript_file is not None:
        annotation_df = pd.read_csv(args.transcript_file, sep="\t")
    else:
        SystemExit("Please provide an annotation file")

    # Merge NM entries with matching HGNC IDs
    print("Merging annotation and gff dataframes")
    #merged_df = pd.concat([transcripts_df, annotation_df], axis=1, join="inner")
    merged_df = transcripts_df.merge(annotation_df, left_on="hgnc_id",
                                     right_on="hgnc_id", how="inner")
    print(merged_df.head())
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
    cols = ['chromosome', 'start_flank', 'end_flank', 'annotation']
    collapsed_df = collapsed_df[cols]
    # Rename columns
    new_column_names = {
        "start_flank": "start",
        "end_flank": "end"
    }
    collapsed_df.rename(columns=new_column_names, inplace=True)

    # Write the collapsed data to an output file
    output_file_name = f"output_{args.reference_genome}" \
                       f"_{args.output_file_suffix}.bed"
    print(collapsed_df.head())
    collapsed_df.to_csv(output_file_name, sep="\t",
                        header=False, index=False)

    # Create an IGV report
    report = config_igv_report(args)

if __name__ == "__main__":
    main()
