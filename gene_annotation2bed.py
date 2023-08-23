"""
This script takes a GFF3 file and an annotation file,
producing a BED file for annotation of relevant transcripts.

"""

import argparse
import argcomplete
import pandas as pd
import numpy as np
import itertools
# import pybedtools # collapsed = bed_file.merge(c=4, o='collapse') or bed_file.merge(c=4, o='distinct')
import gffpandas.gffpandas as gffpd


def import_gff(gff_file):
    """
    Summary: Import GFF3 file and convert to pandas DataFrame.
    """
    transcripts_gff = gffpd.read_gff3(gff_file)
    #gff_df = transcripts_gff.df
    gff_df = transcripts_gff.attributes_to_columns()
    print(gff_df.head())
    print(gff_df.columns)
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
        'hgnc_id': 'Int64'
    }

    gff_df = gff_df.astype(dtype_mapping)

    # Filter GFF DataFrame to select entries with 'NM' type
    transcripts_df = gff_df[gff_df['transcript_id'].str.startswith('NM_')]
    # transcripts_df = transcripts_df.reset_index(drop=True)
    print(transcripts_df.head())
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


def parse_gff3(gff3_file):
    """
    Summary: Parse a GFF3 file into a pandas DataFrame.

    Parameters
    ----------
    gff3_file : _type_
        _description_
    """
    # Process GFF file and create DataFrame with relevant columns

    gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    col_dtypes = {"seqid": "str", "source": "str", "type": "str", "start": "Int64", "end": "Int64", "score": "str",
                  "phase": "str", "strand": "str", "attributes": "str"}
    gff_df = pd.read_csv(gff3_file, sep="\t", comment="#", header=None,
                         names=gff_columns, dtype=col_dtypes)
    gff_df = attributes_to_columns(gff_df)

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

    # set dtype for each column to reduce memory footprint
    dtype_mapping = {
        'ID': 'str',
        'transcript_id': 'str',
        'hgnc_id': 'Int64'
    }

    gff_df = gff_df.astype(dtype_mapping)

    # Apply extract_hgnc_id function to create 'hgnc_id' column
    transcripts_df['hgnc_id'] = transcripts_df['Dbxref'].apply(extract_hgnc_id)

    # Filter GFF DataFrame to select entries with 'NM' type
    transcripts_df = gff_df[gff_df['transcript_id'].str.startswith('NM_')]
    # transcripts_df = transcripts_df.reset_index(drop=True)
    print(transcripts_df.head())
    return transcripts_df


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
        # if chromosome.startswith('na'):
        #     continue
        accession_to_chromosome[accession] = chromosome
    print(accession_to_chromosome)
    return accession_to_chromosome


def map_accession_to_chromosome(accession, accession_to_chromosome):
    """
    Simple mapping function to find chromosome for a given refseq accession.
    """
    return accession_to_chromosome.get(accession, f"Unknown - {accession}")


def parse_pickle(pickle_file):
    gff_df = pd.read_pickle(f"./{pickle_file}")
    transcripts_df = gff_df[gff_df['transcript_id'].fillna('').str.startswith('NM_')]
    print(transcripts_df.columns)
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


def main():
    parser = argparse.ArgumentParser(description="GFF Processing Script")
    parser.add_argument('-ig', "--annotation_file", help="Path to the annotation file (TSV)")
    parser.add_argument('-o', "--output_file_suffix", help="Output file suffix")
    parser.add_argument('-it', "--transcript_file", help="Path to transcript annotation file")
    parser.add_argument('-gff', "--gff_file", help="Path to GFF file")
    parser.add_argument('-ref', "--reference_genome", help="Reference genome (GrCh37/38)")
    parser.add_argument('-f', "--flanking", type=int, help="Flanking size")
    parser.add_argument('-pkl', "--pickle", help="Import gff as pickle file")
    parser.add_argument('--assembly_summary', help="Path to assembly summary file")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # read in pickle file if provided
    if args.pickle is not None:
        transcripts_df = parse_pickle(args.pickle)
        print("Parsed pickle file")
    else:
        # Parse gff file
        transcripts_df = import_gff(args.gff_file)

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
    print(merged_df.head())
    # Create BED file with flanking regions
    print("Creating BED file")
    print("Adding flanking regions")
    merged_df["start_flank"] = merged_df["start"] - args.flanking
    merged_df["end_flank"] = merged_df["end"] + args.flanking
    bed_columns = ["seq_id", "start_flank", "end_flank", "hgnc_id", "annotation"]
    bed_df = merged_df[bed_columns]

    # Extract chromosome from seqid and create the 'chromosome' column
    accession_to_chromosome = read_assembly_mapping(args.assembly_summary)
    # Add a new column 'chromosome' by mapping accession to chromosome identifier
    print(bed_df.head())
    bed_df["chromosome"] = bed_df["seq_id"].apply(lambda x: map_accession_to_chromosome(x, accession_to_chromosome))
    print(bed_df.head())

    # # Merge overlapping entries
    collapsed_df = merge_overlapping(bed_df)

    # Reorder the columns to match the BED format
    cols = collapsed_df.columns.tolist()
    cols = ['chromosome', 'start_flank', 'end_flank', 'annotation']
    collapsed_df = collapsed_df[cols]

    # Write the collapsed data to an output file
    output_file_name = f"output_{args.reference_genome}" \
                       f"_{args.output_file_suffix}.bed"
    collapsed_df.to_csv(output_file_name, sep="\t", index=False)

if __name__ == "__main__":
    main()
