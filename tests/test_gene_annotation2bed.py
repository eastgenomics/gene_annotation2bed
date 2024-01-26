"""
Unit tests for gene_2annotation.py
Tests for gene_annotation2bed.py functions
- convert_coordinates
- extract_hgnc_id

Run: python -m pytest -v tests/test_gene_2annotation_script.py

"""

import os
import sys
from io import StringIO

import argparse
import unittest
import pandas as pd
import pytest
import numpy as np

# set up the path to the module
sys.path.append('../gene_annotation2bed')

from gene_annotation2bed import (convert_coordinates, extract_hgnc_id,
                                 parse_annotation_tsv, parse_pickle,
                                 merge_dataframes, write_bed,
                                 subtract_and_replace, merge_overlapping)

TEST_DATA_DIR = (
    os.path.join(os.path.dirname(__file__), 'test_data')
)


class TestExtractHGNCID(unittest.TestCase):
    def test_extract_hgnc_id_found(self):
        """
        Testing the extraction of HGNC ID from the attributes string.
        Using the extract_hgnc_id function.
        """
        attributes_str = "Dbxref=GeneID:123,HGNC:456"
        result = extract_hgnc_id(attributes_str)
        self.assertEqual(result, 456)

    def test_extract_hgnc_id_not_found(self):
        """
        Testing the extraction of HGNC ID when no HGNC ID is found.
        """
        attributes_str = "Dbxref=GeneID:123"
        result = extract_hgnc_id(attributes_str)
        self.assertIsNone(result)

    def test_extract_hgnc_id_multiple_entries(self):
        """
        Checking the passing of HGNC IDs when multiple are found.
        Currently this passes and selects the first.
        Only one HGNC ID should be found per gene.
        Should raise an error but take first HGNC_id.
        """
        attributes_str = "Dbxref=GeneID:123,HGNC:456,HGNC:789"
        result = extract_hgnc_id(attributes_str)
        self.assertEqual(result, 456)


class TestConvertCoordinates(unittest.TestCase):

    def test_convert_coordinates(self):
        """
        Testing the extraction of HGNC ID from the attributes string.
        Using the convert_coordinates function.
        """
        input_data = {
            "Coordinates": ["chr1:11874-14409", "chr2:20000-25000"],
            "annotation": ["promoter_of_interest", "enhancer"]
        }
        input_df = pd.DataFrame(input_data)
        result_df = convert_coordinates(input_df)
        expected_data = {
            "chromosome": ["1", "2"],
            "start": [11874, 20000],
            "end": [14409, 25000],
            "annotation": ["promoter_of_interest", "enhancer"],
            "gene": ["", ""]
        }
        expected_df = pd.DataFrame(expected_data)
        expected_df = expected_df.astype({'start': "Int64", 'end': "Int64"})
        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_convert_coordinates_empty(self):
        """
        Test when the input DataFrame is empty
        """
        input_df = pd.DataFrame()
        result_df = convert_coordinates(input_df)
        self.assertTrue(result_df.empty)

    def test_convert_coordinates_extra_columns(self):
        """
        Test when the input DataFrame has extra columns
        """
        input_data = {
            "Coordinates": ["chr1:11874-14409", "chr2:20000-25000"],
            "annotation": ["promoter_of_interest", "enhancer"],
            "extra_column": ["extra1", "extra2"]
        }
        input_df = pd.DataFrame(input_data)
        result_df = convert_coordinates(input_df)
        expected_data = {
            "chromosome": ["1", "2"],
            "start": [11874, 20000],
            "end": [14409, 25000],
            "annotation": ["promoter_of_interest", "enhancer"],
            "gene": ["", ""]
        }
        expected_df = pd.DataFrame(expected_data)
        expected_df = expected_df.astype({'start': "Int64", 'end': "Int64"})

        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_convert_coordinates_different_spellings_of_chr(self):
        """
        Test different spellings of "chr" (case-insensitive)
        """
        input_data = {
            "Coordinates": ["Chr1:11874-14409", "Chromosome2:20000-25000"],
            "annotation": ["promoter_of_interest", "enhancer"]
        }
        input_df = pd.DataFrame(input_data)
        result_df = convert_coordinates(input_df)
        expected_data = {
            "chromosome": ["1", "2"],
            "start": [11874, 20000],
            "end": [14409, 25000],
            "annotation": ["promoter_of_interest", "enhancer"],
            "gene": ["", ""]
        }
        expected_df = pd.DataFrame(expected_data)
        expected_df = expected_df.astype({'start': "Int64", 'end': "Int64"})

        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_convert_coordinates_large_numbers(self):
        """
        Test with large numbers
        """
        input_data = {
            "Coordinates": ["chr1:999999999-1000000000"],
            "annotation": ["large_region"]
        }
        input_df = pd.DataFrame(input_data)
        result_df = convert_coordinates(input_df)
        expected_data = {
            "chromosome": ["1"],
            "start": [999999999],
            "end": [1000000000],
            "annotation": ["large_region"],
            "gene": [""]
        }
        expected_df = pd.DataFrame(expected_data)
        expected_df = expected_df.astype({'start': "Int64", 'end': "Int64"})

        pd.testing.assert_frame_equal(result_df, expected_df)

class TestParseAnnotationTsv(unittest.TestCase):
    """
    Tests for checking correct parsing of the annotation resource.
    """

    def setUp(self):
        """
        Set up the test data. Load the preprocessed gff file.
        If not present then exit.
        # TODO: Add a test for the parsing of the gff file.
        # run script to produce gff if not present?
        """
        try:
            self.gff_transcripts_df = parse_pickle(
                f"{TEST_DATA_DIR}/refseq_gff_preprocessed.pkl"
            )
        except FileNotFoundError:
            print("File not found! Ensure the preprocessed gff file is present.")
            sys.exit(1)

    def test_parsing_prints(self):
        """
        Test parsing of transcripts from the annotation file.
        - Checks prints for correct output.
        """
        # Set-up to capture stdout
        capturedOutput = StringIO()
        sys.stdout = capturedOutput

        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
            path, self.gff_transcripts_df)

        # Check print output for coorect that all rows were separated.
        sys.stdout = sys.__stdout__
        print_output = capturedOutput.getvalue().split("\n")
        expected_output = "All rows were separated successfully"

        self.assertEqual(print_output[0], expected_output)

    def test_empty_annotation_file(self):
        """
        Test parsing of transcripts from an empty annotation file.
        Should raise error.
        """
        with self.assertRaises(RuntimeError):
            path = f"{TEST_DATA_DIR}/empty.tsv"
            hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
                path, self.gff_transcripts_df)


    def test_empty_file(self):
        """
        test no pandas import as empty file.
        """
        expected_output = (
            "The annotation file should be a tab-separated file with two columns: ",
            "'ID' and 'annotation'"
        )
        path = f"{TEST_DATA_DIR}/emptyfile.tsv"
        with self.assertRaises(pd.errors.EmptyDataError):
            path = f"{TEST_DATA_DIR}/empty_file.tsv"
            hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
                path, self.gff_transcripts_df)



class TestMerge_Dataframes(unittest.TestCase):
    """
    Tests for merging of gff dataframe with transcripts
    and the HGNC ID dataframe, Gene Symbol dataframe and Coordinates dataframe.
    """

    def setUp(self):
        """
        Set up the test data. Load the preprocessed gff file.
        If not present then exit.
        # TODO: Add a test for the parsing of the gff file.
        # run script to produce pkl file for gff if not present.
        """
        try:
            self.gff_transcripts_df = parse_pickle(
                f"{TEST_DATA_DIR}/refseq_gff_preprocessed.pkl"
            )
        except FileNotFoundError:
            print("File not found! Ensure the preprocessed gff file is present.")
            sys.exit(1)

    def test_parsing_transcripts_prints(self):
        """
        Test parsing of transcripts from the annotation file.
            - Checks prints for correct output.
        """
        # Set-up to capture stdout
        capturedOutput = StringIO()
        sys.stdout = capturedOutput

        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = merge_dataframes(
            hgnc_df, transcript_df, coordinates_df, self.gff_transcripts_df
        )
        # Check print output for coorect that all rows were separated.
        sys.stdout = sys.__stdout__
        print_output = capturedOutput.getvalue().split("\n")

        expected_output_list = [
            'All rows were separated successfully',
            'No HGNC IDs found in the annotation file.', # missing transcript line?
            'No Coordinates found in the annotation file.',
            ''
        ]

        for i, (actual_output, expected_output) in enumerate(zip(print_output, expected_output_list), 1):
            with self.subTest(f"Test case {i} - {actual_output}"):
                self.assertEqual(actual_output, expected_output)

    def test_parsing_transcripts_output(self):
        """
        Test parsing of transcripts from the annotation file.
            - Checks right output df is created with correct numbers of columns.
        """
        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = merge_dataframes(
            hgnc_df, transcript_df, coordinates_df, self.gff_transcripts_df
        )
        # Check df output
        assert annotation_df.empty is False
        transcripts_list = ["NM_000124", "NM_000059", "NM_000546"]

        for i in transcripts_list:
            with self.subTest(transcript=i):
                transcript = annotation_df.loc[annotation_df['transcript_id'] == i]
                assert transcript.empty is False, f"Transcript {i} not found"
                assert_failed_msg = (
                    f"Transcript {i} has an incorrect number of columns:",
                    f"{transcript.shape[1]}"
                )
                assert transcript.shape[1] == 16, assert_failed_msg

        with self.subTest(dataframe="coordinates_df"):
            assert coordinates_df.empty is True, "Coordinates DataFrame should be empty"

    def test_parsing_coordinates(self):
        """
        Test parsing of raw coordinates from the annotation file.
        """
        path = f"{TEST_DATA_DIR}/coordinates_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = merge_dataframes(
            hgnc_df, transcript_df, coordinates_df, self.gff_transcripts_df
        )

        expected_data = {
            'chromosome': ['1', '2', '1', '19', '17'],
            'start': [5000000, 5000, 5000, 1, 1],
            'end': [248956422, 10000, 10000, 100000, 100000],
            'annotation': ['Non-Oncogene', 'Oncogene', 'Oncogene', 'Non-Oncogene', 'Oncogene'],
            'gene': ['', '', '', '', '']
        }
        expected_df = pd.DataFrame(expected_data)
        expected_df = expected_df.astype({
            'start': "Int64", 'end': "Int64",
            'chromosome': "str", 'annotation': "str"
        })

        pd.testing.assert_frame_equal(coordinates_df, expected_df)
        self.assertEqual(annotation_df.empty, True)

    # def test_switched_coordinates_order(self):
    #     """
    #     Test if handles switched start/end or incorrect start/end coordinates
    #     """
    #     pass


class TestWriteBed(unittest.TestCase):
    """
    Tests for writing the bed file function.
    """
    def setUp(self) -> None:
        self.assembly_file = "data/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
        return super().setUp()

    def test_write_bed_creates_file(self):
        """
        Test if the write_bed function creates the correct file.
        """
        # Create sample dataframes and args
        annotation_df = pd.read_csv("example_final_merged_df.csv")
        coordinates_df = pd.read_csv("example_coordinates_df.csv")
        args = argparse.Namespace(
            flanking=10,
            assembly_summary=self.assembly_file,
            genome_build="hg38",
            output_file_suffix="test"
        )

        # Call the function
        write_bed(annotation_df, coordinates_df, args)

        # Check if the output files are created
        assert os.path.exists("output_hg38_test.bed")
        assert os.path.exists("output_hg38_test.maf")

        # Clean up - delete the generated files
        os.remove("output_hg38_test.bed")
        os.remove("output_hg38_test.maf")


    def test_write_bed_edge_cases(self):
        """
        Test empty files for write_bed function.
        Check error raised and files not created.
        """
        with self.assertRaises(RuntimeError):
            # Test with empty annotation df
            annotation_df = pd.DataFrame()
            coordinates_df = pd.read_csv("example_coordinates_df.csv")
            args = argparse.Namespace()

            # Call the function and check if it raises an error
            write_bed(annotation_df, coordinates_df, args)

        # Check if the output files are not created
        assert not os.path.exists("output_.bed")
        assert not os.path.exists("output_.maf")


    # Additional tests for normal function with example annotation df and coordinates df
    def test_write_bed_normal_function(self):
        """
        Test the write_bed function with example annotation and coordinates dataframes.
        """
        # Load example dataframes from files
        annotation_df = pd.read_csv("example_final_merged_df.csv")
        coordinates_df = pd.read_csv("example_coordinates_df.csv")
        args = argparse.Namespace(
            flanking=10,
            assembly_summary=self.assembly_file,
            genome_build="hg38",
            output_file_suffix="test"
        )

        # Call the function
        write_bed(annotation_df, coordinates_df, args)

        # Check if the output files are created
        assert os.path.exists("output_hg38_test.bed")
        assert os.path.exists("output_hg38_test.maf")

        # Clean up - delete the generated files
        os.remove("output_hg38_test.bed")
        os.remove("output_hg38_test.maf")


class TestSubtractAndReplace(unittest.TestCase):
    """
    Tests for subtract_and_replace function.
    """
    def test_subtract_and_replace_positive(self):
        """
        Test subtract_and_replace function with positive integers.
        """
        assert subtract_and_replace(10, 5) == 5
        assert subtract_and_replace(100, 50) == 50


    def test_subtract_and_replace_zero_and_negative(self):
        """
        Test subtract_and_replace function with zero and negative integers.
        """
        assert subtract_and_replace(0, 5) == 1  # Minimum value is 1 when input is 0
        assert subtract_and_replace(3, 5) == 1  # Minimum value is 1 when input is less than 5


    def test_subtract_and_replace_large_integers(self):
        """
        Test subtract_and_replace function with very large integers.
        """
        assert subtract_and_replace(100000000000, 50000000000) == 50000000000
        assert subtract_and_replace(999999999999, 500000000000) == 499999999999

class TestMergeOverlapping(unittest.TestCase):
    def test_merge_overlapping_empty_dataframe(self):
        """
        Test merge_overlapping function with an empty dataframe.
        """
        empty_df = pd.DataFrame(columns=["seq_id", "start_flank", "end_flank", "hgnc_id", "annotation", "gene", "chromosome"])
        with self.assertRaises(RuntimeError):
            merged_df = merge_overlapping(empty_df)


    def test_merge_overlapping_no_overlap(self):
        """
        Test merge_overlapping function with no overlapping regions.
        """
        bed_df = pd.DataFrame({
            "seq_id": ["seq1", "seq2", "seq3"],
            "start_flank": [100, 200, 300],
            "end_flank": [150, 250, 350],
            "hgnc_id": ["hgnc1", "hgnc2", "hgnc3"],
            "annotation": ["anno1", "anno2", "anno3"],
            "gene": ["gene1", "gene2", "gene3"],
            "chromosome": ["chr1", "chr2", "chr3"]
        })
        merged_df = merge_overlapping(bed_df)
        assert merged_df.equals(bed_df)

    def test_merge_overlapping_with_overlap(self):
        """
        Test merge_overlapping function with overlapping regions.
        """
        bed_df = pd.DataFrame({
            "seq_id": ["seq1", "seq1", "seq1", "seq2", "seq2"],
            "start_flank": [100, 120, 140, 200, 220],
            "end_flank": [150, 170, 180, 250, 270],
            "hgnc_id": ["hgnc1", "hgnc1", "hgnc1", "hgnc4", "hgnc4"],
            "annotation": ["anno1", "anno1", "anno1", "anno2", "anno2"],
            "gene": ["gene1", "gene1", "gene1", "gene2", "gene2"],
            "chromosome": ["chr1", "chr1", "chr1", "chr2", "chr2"]
        })
        merged_df = merge_overlapping(bed_df)

        # Define the expected merged dataframe
        expected_df = pd.DataFrame({
            "seq_id": ["seq1", "seq2"],
            "start_flank": [100, 200],
            "end_flank": [180, 270],
            "hgnc_id": ["hgnc1", "hgnc4"],
            "annotation": ["anno1", "anno2"],
            "gene": ["gene1", "gene2"],
            "chromosome": ["chr1", "chr2"]
        })

        assert merged_df.equals(expected_df)

if __name__ == '__main__':
    unittest.main()
