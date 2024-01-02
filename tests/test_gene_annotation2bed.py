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

import unittest
import pandas as pd
import pytest
import numpy as np

# set up the path to the module
sys.path.append('../gene_annotation2bed')
from gene_annotation2bed import convert_coordinates, extract_hgnc_id, parse_annotation_tsv, parse_pickle

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

    def test_parsing_transcripts(self):
        """
        Test parsing of transcripts from the annotation file.
        """
        # Set-up to capture stdout
        capturedOutput = StringIO()
        sys.stdout = capturedOutput

        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_merged_df, coordinates_df = parse_annotation_tsv(path, self.gff_transcripts_df)

        # Check print output for coorect that all rows were separated.
        sys.stdout = sys.__stdout__
        print_output = capturedOutput.getvalue().split("\n")
        expected_output_list = [
            "All rows were separated successfully",
            "All HGNC rows were merged successfully",
            "All Transcript rows were merged successfully",
            "No Coordinates found in the annotation file."
        ]

        for i, (actual_output, expected_output) in enumerate(zip(print_output, expected_output_list), 1):
            with self.subTest(f"Test case {i} - {actual_output}"):
                self.assertEqual(actual_output, expected_output)

        # Check df output
        assert hgnc_merged_df.empty is False
        transcripts_list = ["NM_000124", "NM_000059", "NM_000546"]
        for i in transcripts_list:
            transcript = hgnc_merged_df.loc[hgnc_merged_df['transcript_id'] == i]
            assert transcript.empty == False  # transcript found
            assert transcript.shape[1] == 16  # correct number of columns
        assert coordinates_df.empty == True


    def test_parsing_coordinates(self):
        """
        Test parsing of raw coordinates from the annotation file.
        """
        filename = f"{TEST_DATA_DIR}/coordinates_anno_test.tsv"
        hgnc_merged_df, coordinates_df = parse_annotation_tsv(filename, self.gff_transcripts_df)

        # chr1:5000000-248956422	Non-Oncogene
        # Chromosome1:5000-10000	Oncogene
        # Chr19:1-100000	Non-Oncogene
        # chr17:1-100000	Oncogene
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
        self.assertEqual(hgnc_merged_df.empty, True)


    # def test_switched_coordinates_order(self):
    #     """
    #     Test if handles switched start/end or incorrect start/end coordinates
    #     """
    #     pass


if __name__ == '__main__':
    unittest.main()
