"""
Unit tests for gene_2annotation.py
Tests for gene_annotation2bed.py functions
- convert_coordinates
- extract_hgnc_id

Run: python -m pytest -v tests/test_gene_2annotation_script.py
"""
import sys
from io import StringIO

import unittest
import pandas as pd
import numpy as np

# set up the path to the module
sys.path.append('../gene_annotation2bed')
from gene_annotation2bed import convert_coordinates, extract_hgnc_id


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
        # Test regular input
        input_data = {
            "Coordinates": ["chr1:11874-14409", "chr2:20000-25000"],
            "annotation": ["promoter_of_interest", "enhancer"]
        }
        input_df = pd.DataFrame(input_data)
        print(input_df['Coordinates'])
        result_df = convert_coordinates(input_df)
        print(result_df)
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
        # Test when input DataFrame is empty
        input_df = pd.DataFrame()
        result_df = convert_coordinates(input_df)
        print(result_df)
        self.assertTrue(result_df.empty)

    def test_convert_coordinates_extra_columns(self):
        # Test when input DataFrame has extra columns
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
        # Test different spellings of "chr" (case-insensitive)
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
        # Test with large numbers
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

if __name__ == '__main__':
    unittest.main()
