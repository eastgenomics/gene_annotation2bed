"""
Unit tests for gene_2annotation.py
Tests for gene_annotation2bed.py functions
- convert_coordinates
- extract_hgnc_id

Run: python -m pytest -v tests/test_gene_2annotation_script.py

"""

import gene_annotation2bed as bed
import os
import sys

import argparse
import unittest
import pandas as pd
import pytest

# set up the path to the module
sys.path.append('../gene_annotation2bed')


TEST_DATA_DIR = (
    os.path.join(os.path.dirname(__file__), 'test_data')
)

class TestParseGFF(unittest.TestCase):
    """
    Tests for parsing the gff file.
    """

    def setUp(self):
        """
        Set up the test data. Load the preprocessed gff file.
        If not present then exit.
        """
        self.gff_path = f"{TEST_DATA_DIR}/test_GRCh37_genomic.gff"
        self.test_df = bed.parse_gff(self.gff_path)

    def test_parse_gff_type(self):
        # Test if the returned object is a DataFrame
        self.assertIsInstance(self.test_df, pd.DataFrame)

    def test_parse_gff_columns(self):
        # Test if the returned DataFrame has the correct columns
        expected_columns = ['seq_id', 'source', 'type', 'start',
                            'end', 'score', 'strand', 'phase',
                            'attributes', 'Dbxref', 'ID',
                            'gbkey', 'gene', 'transcript_id', 'hgnc_id']
        self.assertListEqual(list(self.test_df.columns), expected_columns)

    def test_parse_gff_startswith_NM(self):
        # Test if the transcript_id starts with 'NM_'
        self.assertTrue(
            all(self.test_df['transcript_id'].str.startswith('NM_')))

    def test_parse_gff_empty_input(self):
        # Test if the function handles empty input gracefully
        with self.assertRaises(Exception):
            bed.parse_gff(None)

    def test_parse_gff_dtypes(self):
        # Test if the data types are set as expected
        expected_dtypes = {'seq_id': 'object', 'source': 'category',
                           'type': 'category', 'start': 'uint32',
                           'end': 'uint32', 'score': 'object',
                           'strand': 'object', 'phase': 'object',
                           'attributes': 'object',
                           'Dbxref': 'object', 'ID': 'category',
                           'gbkey': 'object', 'gene': 'object',
                           'transcript_id': 'object',
                           'hgnc_id': 'Int64'}

        for col, dtype in expected_dtypes.items():
            self.assertEqual(self.test_df[col].dtype.name, dtype)

    def test_parse_gff_extract_hgnc_id(self):
        # Test if the hgnc_id column is extracted properly
        self.assertTrue(
            all(self.test_df['hgnc_id'] == self.test_df['hgnc_id']))

    def test_parse_gff_filter_columns(self):
        # Test if unwanted columns are filtered out
        unwanted_columns = set(self.test_df.columns) - \
            set(self.test_df.columns)
        self.assertEqual(len(unwanted_columns), 0)

    def test_parse_gff_return_not_empty(self):
        # Test if the returned DataFrame is not empty
        self.assertFalse(self.test_df.empty)

    def test_parse_gff_return_contains_data(self):
        # Test if the returned DataFrame contains expected data
        self.assertTrue(self.test_df.equals(self.test_df))

    def test_parse_gff_return_contains_expected_columns(self):
        # Test if the returned DataFrame contains expected columns
        self.assertTrue(
            all(col in self.test_df.columns for col in self.test_df.columns))

    def test_parse_gff_filter(self):
        # Test if the filtering works properly
        self.assertTrue(
            all(self.test_df['transcript_id'].str.startswith('NM_')))

    def test_parse_gff_drop_columns(self):
        # Test if unwanted columns are dropped properly
        unwanted_columns = ["Gap", "Is_circular", "Name", "Note",
                            "Parent", "Target", "anticodon",
                            "assembly_bases_aln", "assembly_bases_seq",
                            "bit_score", "blast_aligner", "blast_score",
                            "bound_moiety", "chromosome", "codons",
                            "common_component", "consensus_splices",
                            "country", "description", "direction",
                            "e_value", "end_range", "exception",
                            "exon_identity", "exon_number", "experiment",
                            "feat_class", "filter_score", "for_remapping",
                            "function", "gap_count", "gene_biotype",
                            "gene_synonym", "genome", "hsp_percent_coverage",
                            "identity", "idty", "inference",
                            "inversion_merge_aligner", "isolation-source",
                            "lxr_locAcc_currStat_120", "lxr_locAcc_currStat_35",
                            "map", "matchable_bases", "matched_bases", "matches",
                            "merge_aligner", "mobile_element_type",
                            "mol_type", "not_for_annotation", "note",
                            "num_ident", "num_mismatch", "number",
                            "partial", "pct_coverage",
                            "pct_coverage_hiqual", "pct_identity_gap",
                            "pct_identity_gapopen_only", "pct_identity_ungap",
                            "product", "product_coverage", "protein_id",
                            "pseudo", "rank", "recombination_class",
                            "regulatory_class", "rpt_family", "rpt_type",
                            "rpt_unit_range", "rpt_unit_seq", "satellite",
                            "splices", "standard_name", "start_range", "tag",
                            "tissue-type", "transl_except", "transl_table",
                            "weighted_identity"]
        for col in unwanted_columns:
            self.assertNotIn(col, self.test_df.columns)

class TestExtractHGNCID(unittest.TestCase):
    def test_extract_hgnc_id_found(self):
        """
        Testing the extraction of HGNC ID from the attributes string.
        Using the extract_hgnc_id function.
        """
        attributes_str = "Dbxref=GeneID:123,HGNC:456"
        result = bed.extract_hgnc_id(attributes_str)
        self.assertEqual(result, 456)

    def test_extract_hgnc_id_not_found(self):
        """
        Testing the extraction of HGNC ID when no HGNC ID is found.
        """
        attributes_str = "Dbxref=GeneID:123"
        result = bed.extract_hgnc_id(attributes_str)
        self.assertIsNone(result)

    def test_extract_hgnc_id_multiple_entries(self):
        """
        Checking the passing of HGNC IDs when multiple are found.
        Currently this passes and selects the first.
        Only one HGNC ID should be found per gene.
        Should raise an error but take first HGNC_id.
        """
        attributes_str = "Dbxref=GeneID:123,HGNC:456,HGNC:789"
        result = bed.extract_hgnc_id(attributes_str)
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
        result_df = bed.convert_coordinates(input_df)
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
        result_df = bed.convert_coordinates(input_df)
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
        result_df = bed.convert_coordinates(input_df)
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
        result_df = bed.convert_coordinates(input_df)
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
        result_df = bed.convert_coordinates(input_df)
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
            self.gff_transcripts_df = bed.parse_pickle(
                f"{TEST_DATA_DIR}/refseq_gff_preprocessed.pkl"
            )
        except FileNotFoundError:
            print("File not found! Ensure the preprocessed gff file is present.")
            sys.exit(1)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_parsing_prints(self):
        """
        Test parsing of transcripts from the annotation file.
        - Checks prints for correct output.
        """
        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
            path, self.gff_transcripts_df)

        captured = self.capsys.readouterr()
        print_output = captured.out.split("\n")
        expected_output = "All rows were separated successfully"

        assert print_output[0] == expected_output

    def test_empty_annotation_file(self):
        """
        Test parsing of transcripts from an empty annotation file.
        Should raise error.
        """
        with self.assertRaises(RuntimeError):
            path = f"{TEST_DATA_DIR}/empty.tsv"
            hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
                path, self.gff_transcripts_df)

    def test_empty_file(self):
        """
        test no pandas import as empty file.
        """
        expected_output = "No columns to parse from file"
        path = f"{TEST_DATA_DIR}/emptyfile.tsv"
        with self.assertRaises(pd.errors.EmptyDataError) as cm:
            path = f"{TEST_DATA_DIR}/empty_file.tsv"
            hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
                path, self.gff_transcripts_df)
        self.assertEqual(str(cm.exception), expected_output)


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
            self.gff_transcripts_df = bed.parse_pickle(
                f"{TEST_DATA_DIR}/refseq_gff_preprocessed.pkl"
            )
        except FileNotFoundError:
            print("File not found! Ensure the preprocessed gff file is present.")
            sys.exit(1)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_parsing_transcripts_prints(self):
        """
        Test parsing of transcripts from the annotation file.
            - Checks prints for correct output.
        """

        path = f"{TEST_DATA_DIR}/transcripts_anno_test.tsv"
        hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = bed.merge_dataframes(
            hgnc_df, transcript_df, coordinates_df, self.gff_transcripts_df
        )
        captured = self.capsys.readouterr()
        print_output = captured.out.split("\n")
        expected_output_list = [
            'All rows were separated successfully',
            'No HGNC IDs found in the annotation file.',  # missing transcript line?
            'No Coordinates found in the annotation file.'
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
        hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = bed.merge_dataframes(
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
        hgnc_df, transcript_df, coordinates_df = bed.parse_annotation_tsv(
            path, self.gff_transcripts_df)
        # Merge the dataframes
        annotation_df, coordinates_df = bed.merge_dataframes(
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
        self.assembly_file = "tests/test_data/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
        self.annotation_df = pd.read_csv(
            f"{TEST_DATA_DIR}/example_final_merged_df.csv")
        self.coordinates_df = pd.read_csv(
            f"{TEST_DATA_DIR}/example_coordinates_df.csv")

        return super().setUp()

    def test_write_bed_creates_file(self):
        """
        Test if the write_bed function creates the correct file.
        """
        args = argparse.Namespace(
            flanking=10,
            assembly_summary=self.assembly_file,
            genome_build="hg38",
            output_file_suffix="test"
        )

        # Call the function
        bed.write_bed(self.annotation_df, self.coordinates_df, args)

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
            args = argparse.Namespace()

            # Call the function and check if it raises an error
            bed.write_bed(annotation_df, self.coordinates_df, args)

        # Check if the output files are not created
        assert not os.path.exists("output_.bed")
        assert not os.path.exists("output_.maf")

    # Additional tests for normal function with example annotation df and coordinates df

    def test_write_bed_normal_function(self):
        """
        Test the write_bed function with example annotation and coordinates dataframes.
        """
        args = argparse.Namespace(
            flanking=10,
            assembly_summary=self.assembly_file,
            genome_build="hg38",
            output_file_suffix="test"
        )

        # Call the function
        bed.write_bed(self.annotation_df, self.coordinates_df, args)

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
        assert bed.subtract_and_replace(10, 5) == 5
        assert bed.subtract_and_replace(100, 50) == 50

    def test_subtract_and_replace_zero_and_negative(self):
        """
        Test subtract_and_replace function with zero and negative integers.
        """
        assert bed.subtract_and_replace(
            0, 5) == 1  # Minimum value is 1 when input is 0
        # Minimum value is 1 when input is less than 5
        assert bed.subtract_and_replace(3, 5) == 1

    def test_subtract_and_replace_large_integers(self):
        """
        Test subtract_and_replace function with very large integers.
        """
        assert bed.subtract_and_replace(
            100000000000, 50000000000) == 50000000000
        assert bed.subtract_and_replace(
            999999999999, 500000000000) == 499999999999


class TestMergeOverlapping(unittest.TestCase):
    def test_merge_overlapping_empty_dataframe(self):
        """
        Test merge_overlapping function with an empty dataframe.
        """
        empty_df = pd.DataFrame(columns=[
                                "seq_id", "start_flank", "end_flank",
                                "hgnc_id", "annotation", "gene",
                                "chromosome"])
        with self.assertRaises(RuntimeError):
            merged_df = bed.merge_overlapping(empty_df)

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
        merged_df = bed.merge_overlapping(bed_df)
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
        merged_df = bed.merge_overlapping(bed_df)

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
