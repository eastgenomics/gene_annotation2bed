"""
Unit tests for gene_2annotation.py
Tests for gene_annotation2bed.py functions
- convert_coordinates
- extract_hgnc_id
"""
import sys
from io import StringIO

import unittest
import pandas as pd
import numpy as np

# set up the path to the module
sys.path.append('../gene_annotation2bed')
from gene_annotation2bed import convert_coordinates, extract_hgnc_id


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
        pd.testing.assert_frame_equal(result_df, expected_df, check_dtype=False)

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
        pd.testing.assert_frame_equal(result_df, expected_df)

if __name__ == '__main__':
    unittest.main()
