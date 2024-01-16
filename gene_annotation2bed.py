"""
This script takes a GFF3 file and an annotation file,
producing a BED file for annotation of relevant transcripts.
Example cmd:
...

Current working cmd:
/home/rswilson1/anaconda3/envs/Annotation_app/bin/python \
/home/rswilson1/Documents/Programming_project/gene_annotation2bed/bin/gene_annotation2bed.py \
-pkl ./tests/test_data/refseq_gff_preprocessed.pkl \
-ig data/mixed_dataset.tsv \
-ref_igv ./tests/test_data/hs37d5.fa -ref hg38 -f 5 \
--assembly_summary data/GCF_000001405.25_GRCh37.p13_assembly_report.txt \
-o "test_X"
"""

import argparse

import argcomplete
import numpy as np
import pandas as pd
import re

from utils import gff2pandas as gffpd
from scripts import igv_report as igv
from utils import helper


pd.options.mode.chained_assignment = None  # default='warn'


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
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument("-gff", "--gff_file", help="Path to GFF file")
    group1.add_argument("-pkl", "--pickle", help="Import gff as pickle file")

    parser.add_argument(
        "-ig", "--annotation_file", help="Path to the annotation file (TSV)",
        required=True
    )

    parser.add_argument(
        "-o", "--output_file_suffix", help="Output file suffix", required=True
    )
    parser.add_argument(
        "-build", "--genome_build", help="Human reference genome (hg19/hg38)",
        required=True, choices=('hg19', 'hg38')
    )
    parser.add_argument(
        "-ref_igv",
        "--reference_file_for_igv",
        help="Path to Reference genome fasta file for igv_reports",
    )
    parser.add_argument(
        "-f", "--flanking", type=int, help="Flanking size", required=True
    )
    parser.add_argument(
        "--assembly_summary", help="Path to assembly summary file", required=True
    )
    parser.add_argument(
        "-gs", "--symbols_present", help="Flag to indicate gene symbols present",
        required=False, action='store_true'
    ) # is this needed?
    parser.add_argument(
        '-dump', '--hgnc_dump_path', required=False,
        help='Path to HGNC TSV file with HGNC information. \
              Required if gene symbols are present.'
    )
    # parser.add_argument('--report_name', help="Name for report")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    return args


def parse_gff(gff_file):
    """
    Import GFF3 file and convert to pandas DataFrame.

    The GFF3 file is imported into a dataframe and then all the attributes
    in the attributes column are split into seperate columns.
    It then drops many of the additional fields from the attributes column
    which are not needed to reduce memory footprint.
    The dataframe is then filtered to only include entries which have the
    'transcript_id' start with 'NM_'.

    Parameters
    ----------
        gff_file : gff2pandas object
            GFF object which contains the df and header.

    Returns
    -------
    transcripts_df : pandas DataFrame
        DataFrame containing the all the 'NM_' prefixed
        transcripts from the GFF3 file.

    Transformation from initial dataframe (gff df) to final dataframe:
    +--------------+------------+------------+-------+-----------+-------+
    |    seq_id    |   source   |    type    | start |    end    | score |
    +--------------+------------+------------+-------+-----------+-------+
    | NC_000001.10 | RefSeq     | region     |     1 | 249250621 |     . |
    | NC_000001.10 | BestRefSeq | pseudogene | 11874 |     14409 |     . |
    | NC_000001.10 | BestRefSeq | transcript | 11874 |     14409 |     . |
    +--------------+------------+------------+-------+-----------+-------+
    +--------+-------+---------------------------------------------------+
    | strand | phase |                   attributes                      |
    +--------+-------+---------------------------------------------------+
    | +      |     . | attributes string...                              |
    | +      |     . | attributes string...                              |
    | +      |     . | attributes string...                              |
    +--------+-------+---------------------------------------------------+

                                    |
                                    |
                                    |
                                    V

    Transcripts dataframe:
    +--------------+------------+------+-------+-------+-------+--------+-------+
    |    seq_id    |   source   | type | start |  end  | score | strand | phase |
    +--------------+------------+------+-------+-------+-------+--------+-------+
    | NC_000001.10 | BestRefSeq | mRNA | 65419 | 71585 |     . | +      |     . |
    | NC_000001.10 | BestRefSeq | exon | 65419 | 65433 |     . | +      |     . |
    | NC_000001.10 | BestRefSeq | exon | 65520 | 65573 |     . | +      |     . |
    +--------------+------------+------+-------+-------+-------+--------+-------+
    +-----------------------------------------------------+----------------------+
    |                       Dbxref                        |          ID          |
    +-----------------------------------------------------+----------------------+
    | GeneID:79501,Genbank:NM_001005484.2,HGNC:HGNC:14825 | rna-NM_001005484.2   |
    | GeneID:79501,Genbank:NM_001005484.2,HGNC:HGNC:14825 | exon-NM_001005484.2-1|
    | GeneID:79501,Genbank:NM_001005484.2,HGNC:HGNC:14825 | exon-NM_001005484.2-2|
    +-----------------------------------------------------+----------------------+
    +---------------------------------------------------+
    |                   attributes                      |
    +---------------------------------------------------+
    | attributes string...                              |
    | attributes string...                              |
    | attributes string...                              |
    +---------------------------------------------------+
    +-------+-------+----------------+---------+
    | gbkey | gene  | transcript_id  | hgnc_id |
    +-------+-------+----------------+---------+
    | mRNA  | OR4F5 | NM_001005484.2 |   14825 |
    | mRNA  | OR4F5 | NM_001005484.2 |   14825 |
    | mRNA  | OR4F5 | NM_001005484.2 |   14825 |
    +-------+-------+----------------+---------+

    Produced using https://ozh.github.io/ascii-tables/.
    """
    transcripts_gff = gffpd.read_gff3(gff_file)
    gff_df = transcripts_gff.attributes_to_columns()
    # drop columns that are not needed to reduce memory footprint
    gff_df = gff_df.drop(
        [
            "Gap", "Is_circular", "Name", "Note", "Parent", "Target", "anticodon",
            "assembly_bases_aln", "assembly_bases_seq", "bit_score", "blast_aligner",
            "blast_score", "bound_moiety", "chromosome", "codons", "common_component",
            "consensus_splices", "country", "description", "direction", "e_value",
            "end_range", "exception", "exon_identity", "exon_number", "experiment",
            "feat_class", "filter_score", "for_remapping", "function", "gap_count",
            "gene_biotype", "gene_synonym", "genome", "hsp_percent_coverage",
            "identity", "idty", "inference", "inversion_merge_aligner",
            "isolation-source", "lxr_locAcc_currStat_120", "lxr_locAcc_currStat_35",
            "map", "matchable_bases", "matched_bases", "matches", "merge_aligner",
            "mobile_element_type", "mol_type", "not_for_annotation", "note",
            "num_ident", "num_mismatch", "number", "partial", "pct_coverage",
            "pct_coverage_hiqual", "pct_identity_gap", "pct_identity_gapopen_only",
            "pct_identity_ungap", "product", "product_coverage", "protein_id",
            "pseudo", "rank", "recombination_class", "regulatory_class",
            "rpt_family", "rpt_type", "rpt_unit_range", "rpt_unit_seq",
            "satellite", "splices", "standard_name", "start_range", "tag",
            "tissue-type", "transl_except", "transl_table", "weighted_identity",
        ],
        axis=1,
    )

    # Apply extract_hgnc_id function to create 'hgnc_id' column
    gff_df["hgnc_id"] = gff_df["Dbxref"].apply(extract_hgnc_id)

    # set dtype for each column to reduce memory footprint
    dtype_mapping = {
        "ID": "category",
        "transcript_id": "str",
        "hgnc_id": "Int64",
    }

    gff_df = gff_df.astype(dtype_mapping)

    # Filter GFF DataFrame to select entries with 'NM' type
    gff_transcripts_df = gff_df[gff_df["transcript_id"].str.startswith("NM_")]
    return gff_transcripts_df


def replace_chromosome_prefix_suffix(chromosome):
    """
    replace chr/chromosome in chromosome column.

    Parameters
    ----------
    chromosome : str
        string from chromosome column. i.e. chr1, chromosome1, Chr1.

    Returns
    -------
    replaced string
        string with instances of chr/chromosome replaced with empty string
    """
    return re.sub(r"(?i)(chr|omosome)?", "", chromosome)


def convert_coordinates(coordinates_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert coordinates dataframe to BED format.

    Parameters
    ----------
    coordinates_df : pd.DataFrame (0-based)
        coordinate format provided by the annotation file.
        ID                  annotation
        chr1:11874-14409    promoter_of_interest

    Returns
    -------
    pd.DataFrame
        Bed format dataframe with columns: chromosome, start,
            end, annotation, gene.

        +------------------+----------------------+
        |        ID        |      annotation      |
        +------------------+----------------------+
        | chr1:11874-14409 | promoter_of_interest |
        +------------------+----------------------+

                           |
                           |
                           V

    +------------+-------+-------+----------------------+
    | chromosome | start |  end  |      annotation      |
    +------------+-------+-------+----------------------+
    | chr1       | 11874 | 14409 | promoter_of_interest |
    +------------+-------+-------+----------------------+
    """
    # If the "Coordinates" column is empty, return an empty dataframe:
    if coordinates_df.empty:
        # Define the columns and their corresponding data types

        # Create an empty DataFrame with specified columns and data types
        empty_df = pd.DataFrame(
            columns=["chromosome", "start", "end", "annotation", "gene"])
        print("No Coordinates found in the annotation file.")
        return empty_df

    # create columns
    coordinates_df[["chr", "start", "end", "gene"]] = ""

    try:
        # Split the "Coordinates" column by ':' and '-'
        coordinates_df[["chromosome", "start", "end"]] = coordinates_df[
            "Coordinates"
        ].str.split("[:-]", expand=True)
        coordinates_df["chromosome"] = coordinates_df["chromosome"].apply(
            replace_chromosome_prefix_suffix)
        coordinates_df = coordinates_df[
            ["chromosome", "start", "end", "annotation", "gene"]
        ]

    except Exception as err:
        print(f"Error: {err}")
        print("Please check the format of the coordinates in the annotation file.")
        empty_df = pd.DataFrame(
            columns=["chromosome", "start", "end", "annotation", "gene"])
        return empty_df

    try:
        coordinates_df["chromosome"] = coordinates_df["chromosome"].astype(
            'str')
        coordinates_df["start"] = coordinates_df["start"].astype('Int64')
        coordinates_df["end"] = coordinates_df["end"].astype('Int64')
        coordinates_df["annotation"] = coordinates_df["annotation"].astype(
            'str')
    except ValueError as e:
        print(f"Error: {e}")

    return coordinates_df


def parse_annotation_tsv(path: str,
                         gff_transcripts_df: pd.DataFrame,
                         hgnc_dump=None):
    """
    Parse an annotation TSV file and separate it into dataframes for HGNC IDs,
    Transcript IDs, and Coordinates, then merge them with a GFF dataframe.

    Parameters
    ----------
    path : str
        The file path to the TSV annotation file.
    gff_transcripts_df : pd.DataFrame
        A dataframe containing GFF information including transcript IDs.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple containing two dataframes:
        1. The merged dataframe for HGNC IDs and transcripts. (hgnc_merged_df)
        2. The coordinated dataframe for coordinates to be appended
           to a BED file later (coordinates_df).
        3. The gene symbol dataframe for gene symbols to be appended.
    """
    df = pd.read_csv(path, sep="\t", dtype={
                     'ID': 'string', 'annotation': 'string'})
    # Create masks for HGNC, Transcript, and Coordinates dataframes
    assert 'ID' in df.columns, 'The annotation file does not contain an "ID" column'

    hgnc_mask = df["ID"].str.startswith("HGNC:") | df["ID"].str.isnumeric()
    # Use regex to match transcript IDs/chromosome coordinates
    pattern_nm = r'^NM'
    transcript_mask = df["ID"].str.contains(pattern_nm, case=True)
    pattern_chr = r'^(chr|chromosome|Chr|Chromosome)'
    coordinates_mask = df["ID"].str.contains(pattern_chr, case=False)
    pattern_gene_symbol = r'\b[A-Z][A-Z0-9]+\b'  # r'^[A-Z]+[A-Z0-9]+'
    gene_symbol_mask = df["ID"].str.contains(pattern_gene_symbol, case=True)
    # Use masks to filter the original dataframe
    not_separated_rows = df[
        ~(hgnc_mask | transcript_mask | coordinates_mask | gene_symbol_mask)
    ]

    # for hgcnID and transcriptID don't exist
    if not_separated_rows.empty:
        print("All rows were separated successfully")
    else:
        print(f"These rows were not separated: \n {not_separated_rows}")

    # Create dataframes for HGNC IDs, Transcript IDs, and Coordinates
    hgnc_df = df[hgnc_mask]
    transcript_df = df[transcript_mask]
    coordinates_df = df[coordinates_mask]
    gene_symbol_df = df[gene_symbol_mask]
    # set dtype for each column
    dtype_mapping_hgnc = {
        "ID": "Int64",
        "annotation": "category",
    }

    dtype_mapping_transcript = {
        "ID": "str",
        "annotation": "category",
    }
    dtype_mapping_gff = {
        "hgnc_id": "Int64",
    }
    dtype_mapping_symbols = {
        "ID": "str",
        "annotation": "category",
    }

    hgnc_df = hgnc_df.astype(dtype_mapping_hgnc)
    transcript_df = transcript_df.astype(dtype_mapping_transcript)
    gff_transcripts_df = gff_transcripts_df.astype(dtype_mapping_gff)
    gene_symbol_df = gene_symbol_df.astype(dtype_mapping_symbols)
    # Rename columns for clarity
    hgnc_df = hgnc_df.rename(columns={"ID": "hgnc_id"})
    transcript_df = transcript_df.rename(columns={"ID": "transcript_id"})
    coordinates_df = coordinates_df.rename(columns={"ID": "Coordinates"})

    # Remove everything after '.' in the "transcript_id" column
    gff_transcripts_df["transcript_id"] = (
        gff_transcripts_df["transcript_id"].str.split(".").str[0]
    )
    transcript_df["transcript_id"] = (
        transcript_df["transcript_id"].str.split(".").str[0]
    )

    # Convert gene symbols to hgnc_ids for concatenation
    if hgnc_dump is not None:
        if not gene_symbol_df.empty:
            corrected_gene_symbols = convert_symbols_to_hgnc_ids(
                gene_symbol_df, hgnc_dump
            )
            gene_symbol_df_corrected = corrected_gene_symbols[
                ['hgnc_id', 'annotation']
            ].copy()
            # Stripping "HGNC:" from the 'HGNC ID' column using apply and lambda
            gene_symbol_df_corrected['hgnc_id'] = gene_symbol_df_corrected[
                'hgnc_id'
            ].apply(lambda x: x.split(":")[1])
            dtype_mapping_symbols = {
                "hgnc_id": "Int64",
                "annotation": "category",
            }
            gene_symbol_df_corrected = gene_symbol_df_corrected.astype(
                dtype_mapping_symbols
                )

        else:
            print("No gene symbols found in the annotation file.")
            gene_symbol_df_corrected = gene_symbol_df
    elif hgnc_dump is None:
        if not gene_symbol_df.empty:
            raise RuntimeError(
                "No HGNC dump file provided."
                "Gene symbols present will not be converted to HGNC IDs."
            )
        else:
            gene_symbol_df_corrected = gene_symbol_df

    return hgnc_df, transcript_df, coordinates_df, gene_symbol_df_corrected


def merge_dataframes(hgnc_df: pd.DataFrame,
                     transcript_df: pd.DataFrame,
                     coordinates_df: pd.DataFrame,
                     gff_df: pd.DataFrame,
                     gene_symbol_df: pd.DataFrame):
    """
    merge_dataframes function
    Extract the corresponding transcripts from the GFF dataframe using HGNC_ID.
    Then Merge based on the HGNC_ID field into final dataframes
    with just coordinates and annotation.

    Parameters
    ----------
    hgnc_df : pd.DataFrame
        A dataframe containing HGNC IDs with annotation.
    transcript_df : pd.DataFrame
        A dataframe containing transcript information with annotation.
    coordinates_df : pd.DataFrame
        A dataframe containing coordinates with annotation.
    gff_df : pd.DataFrame
        A dataframe containing GFF information including transcript IDs
        and HGNC_Ids and coordinate information for producing final bed.
    gene_symbol_df : pd.DataFrame
        A dataframe containing HGNC_IDs and annotation
        for provided gene symbols.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple containing two dataframes:
        1. final_merged_df
            The merged dataframe for HGNC IDs, gene symbols and transcripts.
        2. coordinates_df
            The dataframe with the coordinates to be appended
            to a BED file later.
    """
    # Check if any of the dataframes are empty
    # If so, create an empty dataframe with the correct columns
    # for the merge to work.
    if gene_symbol_df.empty:
        print("No Gene Symbols found in the annotation file.")
        gene_symbol_df = pd.DataFrame(columns=["hgnc_id", "annotation"])
    if hgnc_df.empty:
        print("No HGNC IDs found in the annotation file.")
        hgnc_df = pd.DataFrame(columns=["hgnc_id", "annotation"])
    if transcript_df.empty:
        print("No Transcript IDs found in the annotation file.")
        transcript_df = pd.DataFrame(columns=["transcript_id", "annotation"])
    # Remove everything after '.' in the "transcript_id" column for gff resource
    gff_df["transcript_id"] = (
        gff_df["transcript_id"].str.split(".").str[0]
    )
    # Merge the HGNC and Transcript dataframes with gff dataframe based on the 'ID' column
    merged_hgnc_df = gff_df.merge(
        hgnc_df, on="hgnc_id", how="inner")
    merged_transcript_df = gff_df.merge(
        transcript_df, on="transcript_id", how="inner"
    )
    merged_symbol_df = gff_df.merge(
        gene_symbol_df, on="hgnc_id", how="inner")

    # Sanity Check
    # Find any rows dropped during the merge
    dropped_hgnc_rows = hgnc_df[~hgnc_df["hgnc_id"].isin(
        merged_hgnc_df["hgnc_id"])]
    dropped_transcript_rows = transcript_df[
        ~transcript_df["transcript_id"].isin(
            merged_transcript_df["transcript_id"])
    ]
    dropped_symbol_rows = gene_symbol_df[~gene_symbol_df["hgnc_id"].isin(
        merged_symbol_df["hgnc_id"])]
    if not dropped_hgnc_rows.empty:
        print(f"Summary of dropped HGNC rows: \n {dropped_hgnc_rows}")
    else:
        if not merged_hgnc_df.empty:
            print("All HGNC rows were merged successfully")
    if not dropped_transcript_rows.empty:
        print(
            f"Summary of dropped Transcript rows: \n {dropped_transcript_rows}")
    else:
        if not merged_transcript_df.empty:
            print("All Transcript rows were merged successfully")
    if not dropped_symbol_rows.empty:
        print(f"Summary of dropped Symbol rows: \n {dropped_symbol_rows}")
    else:
        if not merged_symbol_df.empty:
            print("All Symbol rows were merged successfully")

    # Concatenate the merged dataframes to get a single merged dataframe.
    final_merged_df = pd.concat([merged_hgnc_df, merged_transcript_df])
    final_merged_df = pd.concat([final_merged_df, merged_symbol_df])

    # Coordinates dataframe split into columns
    coordinates_df = convert_coordinates(coordinates_df)

    return final_merged_df, coordinates_df


def extract_hgnc_id(dbxref_str: str):
    """
    Wrapper function to extract HGNC ID from a string of dbxrefs.

    Parameters
    ----------
    dbxref_str : str
        various ids separated by commas

    Returns
    -------
    int | None
        HGNC ID as an integer i.e. 427 for HGNC:427.
        Returns None if no HGNC ID found.

    Raises
    ------
    ValueError
        If more than one HGNC ID is found in the input string.
    """
    if not dbxref_str:
        return None
    parts = dbxref_str.split(",")
    hgnc_ids = []
    for part in parts:
        match = re.search(r"hgnc[:_][0-9]+", part, re.IGNORECASE)
        if match:
            hgnc_id = int(match.group().replace("_", ":").split(":")[-1])
            hgnc_ids.append(hgnc_id)
    try:
        if len(hgnc_ids) > 1:
            raise ValueError("Multiple HGNC IDs found: " +
                             ", ".join(map(str, hgnc_ids)))
        elif hgnc_ids:
            return hgnc_ids[0]
        else:
            return None
    except ValueError as e:
        print(f"Error: {e}")
        return hgnc_ids[0]


def convert_row(row, hgnc_dump):
    """
    Convert gene symbols in the ID column of the DataFrame to HGNC ids.
    Parameters
    ----------
    row : pd.DataFrame
        Row of a DataFrame with 'ID' column containing gene symbols.
    hgnc_dump : pd.DataFrame
        HGNC DataFrame with 'Approved', 'Previous symbols',
        'Alias symbols', 'HGNC ID', and 'Status' columns.

    Returns
    -------
    pd.Series
        Series with 'ID', 'annotation', 'hgnc_id', and 'Status' columns.
    """
    gene_symbol = row['ID']
    hgnc_info = helper.find_hgnc_id(gene_symbol, hgnc_dump)
    return pd.Series({
        'ID': gene_symbol,
        'annotation': row['annotation'],  # Preserve the 'Annotation' column
        'hgnc_id': hgnc_info[1],
        'Status': hgnc_info[2]
    })


def convert_symbols_to_hgnc_ids(df, hgnc_dump):
    """Convert gene symbols in the ID column of the DataFrame to HGNC ids.

    Args:
        df (pd.DataFrame): DataFrame with 'ID' column containing gene symbols.
        hgnc_dump (pd.DataFrame): HGNC DataFrame.

    Returns:
        pd.DataFrame: DataFrame with 'ID', 'annotation', 'HGNC ID', and 'Status' columns.
    """
    print("Converting gene symbols to HGNC IDs...")
    result_df = df.apply(lambda row: convert_row(row, hgnc_dump), axis=1)

    return result_df


def read_assembly_mapping(assembly_file: str):
    """
    Reads in the associated assembly file and returns a dictionary mapping
    to find chromosome for each refseq accession.

    Parameters
    ----------
    assembly_file : str (file path to tsv)
        found at: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/

    Returns
    -------
    dictionary
        mapping of refseq accession to chromosome
    """
    accession_to_chromosome = {}
    assembly_df = pd.read_csv(assembly_file, sep="\t",
                              comment="#", header=None)
    assembly_df = assembly_df.dropna()  # Drop rows with missing values
    # filter out na from chromosome column and turn accession and chromosome columns to dict
    assembly_df = assembly_df[~assembly_df[2].str.startswith("na")]
    accession_to_chromosome = dict(zip(assembly_df[6], assembly_df[2]))

    return accession_to_chromosome


def map_accession_to_chromosome(accession: str, accession_to_chromosome: dict):
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


def parse_pickle(pickle_file: str):
    """
    Parses a pickle file and returns a DataFrame of transcripts.

    Parameters
    ----------
    pickle_file : str (path to Pickle file)
        pickle file of a GFF DataFrame once parsed
        with columns from attributes_to_columns (1-based)

    Returns
    -------
    transcripts_df: dataframe
        dataframe of transcripts with columns for attributes.
        Contains only transcripts with NM_ prefix.
    """
    gff_df = pd.read_pickle(pickle_file)
    transcripts_df = gff_df[gff_df["transcript_id"].fillna(
        "").str.startswith("NM_")]
    return transcripts_df


def merge_overlapping(bed_df: pd.DataFrame):
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
        dataframe of merged rows with columns: chromosome, start,
        end, annotation
    """
    # Sort by chromosome, start, and end
    # This makes sure that overlapping regions are next to each other.

    bed_df = bed_df.sort_values(
        by=["annotation", "chromosome", "start_flank", "end_flank"]
    )
    # Sort by first annotation then chromosome, start, and end.
    merged_rows = []

    current_row = bed_df.iloc[0]
    for _, row in bed_df.iterrows():
        if row["annotation"] != current_row["annotation"]:
            merged_rows.append(current_row)  # Append the merged row
            current_row = row  # Start a new potential merged row
            # Only rows with same annotation are merged
        if row["chromosome"] != current_row["chromosome"]:
            merged_rows.append(current_row)
            current_row = row
            # Only rows with same chromosome are merged.
        if row["start_flank"] <= current_row["end_flank"]:
            current_row["end_flank"] = max(
                current_row["end_flank"], row["end_flank"])
            # Extend the end if overlapping
        else:
            merged_rows.append(current_row)
            current_row = row

    merged_rows.append(current_row)  # Append the last merged row
    merged_df = pd.DataFrame(merged_rows)
    return merged_df


def config_igv_report(args: argparse.Namespace):
    """
    Function to call igv report script with the correct parameters.
    Generates an IGV html report using generic handling.

    Parameters
    ----------
    args : argeparse object
        argeparse object with the following attributes:
        genome_build, output_file_suffix, gff_file/pickle_file,
        annotation_file/transcript_file, assembly_file, and flanking.

    Returns
    -------
    None
    """
    # assign vars.
    maf_file = f"output_{args.genome_build}_{args.output_file_suffix}.maf"
    bed_file = f"output_{args.genome_build}_{args.output_file_suffix}.bed"
    genome = args.genome_build
    fasta_ref = args.reference_file_for_igv
    info_columns = []
    title = f"{args.output_file_suffix}_report"
    output_file = f"{title}.html"
    print("Creating IGV report...")

    print(
        f"Bed file: {bed_file}\nGenome: {genome}\n"
        f"Info columns: {info_columns}\nTitle: {title}\nOutput: {output_file}"
    )

    igv.create_igv_report(
        bed_file, maf_file, genome, fasta_ref, info_columns, title, output_file
    )

    print("IGV report created successfully!")


def subtract_and_replace(position, flanking_int):
    """
    Define a function to apply the subtraction and replace with 0 if negative

    Parameters
    ----------
    position : int
        position to subtract from.
    flanking_int : int
        integer value to subtract from each value in the list.

    Returns
    -------
    result : list
        list of values with flanking subtracted, minimum value = 0.
    """
    return max(1, position - flanking_int)


def write_bed(annotation_df: pd.DataFrame,
              coordinates_df: pd.DataFrame,
              args: argparse.Namespace) -> None:
    """
    Combines dataframes, extracts chromosome for HGNC_ids,
    and writes to MAF & BED file for IGV visualisation and VEP annotation.

    Parameters
    ----------
    annotation_df : pd.DataFrame
        A dataframe containing annotation information.
    coordinates_df : pd.DataFrame
        A dataframe containing coordinates and annotation info..
    args : Namespace
        A namespace containing command-line arguments and options.

    Outputs
    -------
    bed file: (file) bed file containing the relevant transcripts
        for annotation for visualisation in igv.
    """
    # Create BED file with flanking regions
    print("Creating BED file")
    print("Adding flanking regions")

    # Apply the function to the specified column
    annotation_df["start_flank"] = annotation_df["start"].apply(
        subtract_and_replace, flanking_int=args.flanking
    )
    annotation_df["end_flank"] = annotation_df["end"].apply(
        subtract_and_replace, flanking_int=args.flanking
    )

    bed_columns = [
        "seq_id",
        "start_flank",
        "end_flank",
        "hgnc_id",
        "annotation",
        "gene",
    ]
    bed_df = annotation_df[bed_columns]
    bed_df = bed_df.reindex()
    print(f"Summary of BED file df before collapsing \n {bed_df.head()}")
    # Extract chromosome from seqid and create the 'chromosome' column
    accession_to_chromosome = read_assembly_mapping(args.assembly_summary)
    # Add a new column 'chromosome' by mapping accession to chromosome identifier
    bed_df.loc[:, "chromosome"] = bed_df["seq_id"].apply(
        lambda x: map_accession_to_chromosome(x, accession_to_chromosome)
    )
    print(f"Summary of BED file df before collapsing \n {bed_df.head()}")

    # Merge overlapping entries
    collapsed_df = merge_overlapping(bed_df).reset_index(drop=True)
    print(f"Summary of BED file df after collapsing \n {collapsed_df.head()}")
    # Reorder the columns to match the BED format
    cols = ["chromosome", "start_flank", "end_flank", "annotation", "gene"]
    collapsed_df = collapsed_df[cols]
    # Rename columns
    new_column_names = {"start_flank": "start", "end_flank": "end"}
    collapsed_df.rename(columns=new_column_names, inplace=True)
    collapsed_df = pd.concat(
        [collapsed_df, coordinates_df], axis=0, ignore_index=True)
    # Write the collapsed data to an output file
    output_file_name_maf = (
        f"output_{args.genome_build}_{args.output_file_suffix}.maf"
    )
    output_file_name_bed = (
        f"output_{args.genome_build}_{args.output_file_suffix}.bed"
    )
    collapsed_df.to_csv(output_file_name_maf, sep="\t",
                        header=True, index=False)
    collapsed_df.to_csv(output_file_name_bed, sep="\t",
                        header=False, index=False)


def main():
    """
    Main logic for script
    Collects arguments.
    Based on this imports the correct inputs and parses them.
    Creates a BED file for annotation of relevant transcripts.
    Creates an IGV report.
    """
    args = parse_args()

    # read in pickle file if provided
    if args.pickle:
        transcripts_df = parse_pickle(args.pickle)
        print("Parsed pickle file")
    else:
        # Parse gff file
        transcripts_gff_df = parse_gff(args.gff_file)

    if args.symbols_present & (args.hgnc_dump_path is not None):
        hgnc_dump = helper.parse_tsv(args.hgnc_dump_path)
    elif args.symbols_present & (not args.hgnc_dump_path):
        raise RuntimeError(
            "No HGNC dump file provided."
            "Gene symbols present will not be converted to HGNC IDs."
        )
    elif (not args.symbols_present) & args.hgnc_dump_path:
        raise RuntimeError(
            "HGNC dump file provided but -gs flag not set."
            "Re-run with -gs flag to convert gene symbols to HGNC IDs."
        )
    else:
        print("No gene symbols flag, returning None for hgnc_dump")
        hgnc_dump = None
    # Read the annotation file into a pandas DataFrame

    hgnc_df, transcript_df, coordinates_df, corrected_gene_symbols = parse_annotation_tsv(
        args.annotation_file, transcripts_gff_df, hgnc_dump
    )
    # Merge the dataframes
    annotation_df, coordinates_df = merge_dataframes(
        hgnc_df, transcript_df, coordinates_df, transcripts_gff_df, corrected_gene_symbols
    )
    # Merge NM entries with matching HGNC IDs

    print("Merging annotation and gff dataframes")
    write_bed(annotation_df, coordinates_df, args)
    # Create an IGV report
    config_igv_report(args)


if __name__ == "__main__":
    main()
