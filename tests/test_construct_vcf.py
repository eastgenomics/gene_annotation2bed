import unittest
import pandas as pd
import pytest
import os
from unittest.mock import patch
from construct_vcf import ConstructVCF, parse_args

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


class TestConstructVCF(unittest.TestCase):
    def setUp(self):
        # Define paths for sample data
        self.bed_file_path = f"{TEST_DATA_DIR}/example_bed_hg38.bed"
        self.reference_path = f"{TEST_DATA_DIR}/hs37d5.fa"
        self.output_file = "output_new_test.vcf"
        self.bed_df = pd.read_csv(f"{TEST_DATA_DIR}/example_bed_hg38.bed",
                                  sep="\t", header=None, names=["chr", "start", "end", "info", "gene"])

    def test_parse_args(self):
        # Test command line argument parsing
        with patch("sys.argv", ["script.py", "-b", "tests/test_data/example_bed_hg38.bed", "-fasta", "tests/test_data/hs37d5.fa", "-o", "output_new_test.vcf"]):
            args = parse_args()
            self.assertEqual(
                args.bed_file, "tests/test_data/example_bed_hg38.bed")
            self.assertEqual(args.reference_file, "tests/test_data/hs37d5.fa")
            self.assertEqual(args.output_file, "output_new_test.vcf")

    def test_fetch_nucleotides_normal(self):
        # Test fetch_nucleotides method

        # Mock the pysam.faidx function
        with patch("pysam.faidx") as mock_faidx:
            # Mock the return values
            # mock_faidx.return_value = ">\nA\n" # if you want to return a single value
            mock_faidx.side_effect = [">\nA\n", ">\nT\n", ">\nC\n", ">\nG\n"]

            # Initialize ConstructVCF object
            vcf_obj = ConstructVCF(
                self.bed_file_path, self.reference_path, self.output_file)

            # Mock a sample row
            row = pd.Series({'chr': 'chr14', 'start': 1,
                            'end': 10, 'info': 'test', 'gene': 'BRCA1'})

            # Call the method
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)
            print(result_df)
            print(result_df.iloc[0])
            # Assertions
            assert mock_faidx.call_args_list == [
                ((self.reference_path, "14:2-2"), {}),
                ((self.reference_path, "14:5-5"), {}),
                ((self.reference_path, "14:10-10"), {}),
            ]
            self.assertEqual(result_df.iloc[0]['REF'], 'A')
            self.assertEqual(result_df.iloc[0]['ALT'], 'T')
            self.assertEqual(result_df.iloc[1]['REF'], 'T')
            self.assertEqual(result_df.iloc[1]['ALT'], 'A')
            self.assertEqual(result_df.iloc[2]['REF'], 'C')
            self.assertEqual(result_df.iloc[2]['ALT'], 'G')

    def test_fetch_nucleotides_normal_X(self):
        # Test fetch_nucleotides method

        # Mock the pysam.faidx function
        with patch("pysam.faidx") as mock_faidx:
            # Mock the return values
            # mock_faidx.return_value = ">\nA\n" # if you want to return a single value
            mock_faidx.side_effect = [">\nA\n", ">\nT\n", ">\nC\n", ">\nG\n"]

            # Initialize ConstructVCF object
            vcf_obj = ConstructVCF(
                self.bed_file_path, self.reference_path, self.output_file)

            # Mock a sample row
            row = pd.Series({'chr': 'chrX', 'start': 1,
                            'end': 10, 'info': 'test', 'gene': 'Gene_A'})

            # Call the method
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)
            print(result_df)
            print(result_df.iloc[0])
            # Assertions
            assert mock_faidx.call_args_list == [
                ((self.reference_path, "X:2-2"), {}),
                ((self.reference_path, "X:5-5"), {}),
                ((self.reference_path, "X:10-10"), {}),
            ]
            self.assertEqual(result_df.iloc[0]['REF'], 'A')
            self.assertEqual(result_df.iloc[0]['ALT'], 'T')
            self.assertEqual(result_df.iloc[1]['REF'], 'T')
            self.assertEqual(result_df.iloc[1]['ALT'], 'A')
            self.assertEqual(result_df.iloc[2]['REF'], 'C')
            self.assertEqual(result_df.iloc[2]['ALT'], 'G')

    def test_fetch_nucleotides_invalid_chr_character(self):
        # Tests if the correct error is raised by
        # fetch_nucleotides when an invalid chromosome is passed

        # Initialize ConstructVCF object
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        # Mock a sample row
        chromosome = 'chrA'
        row = pd.Series({'chr': f'{chromosome}', 'start': 1,
                        'end': 10, 'info': 'test', 'gene': 'BRCA1'})
        # Check error is raised:
        with pytest.raises(ValueError) as excinfo:
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)
        assert str(
            excinfo.value) == f"Error: {chromosome.replace('chr', '')} is not a valid chromosome"

    def test_fetch_nucleotides_invalid_chr_number(self):
        # Tests if the correct error is raised by
        # fetch_nucleotides when an invalid chromosome is passed

        # Initialize ConstructVCF object
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        # Mock a sample row
        chromosome = 'chr100'
        row = pd.Series({'chr': f'{chromosome}', 'start': 1,
                        'end': 10, 'info': 'test', 'gene': 'BRCA1'})
        # Check error is raised:
        with pytest.raises(ValueError) as excinfo:
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)
        assert str(
            excinfo.value) == f"Error: {chromosome.replace('chr', '')} is not a valid chromosome, greater than 22"

    def test_fetch_nucleotides_invalid_case(self):
        # Tests if the correct error is raised by
        # fetch_nucleotides when an invalid chromosome is passed

        # Initialize ConstructVCF object
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        # Mock a sample row
        chromosome = 'chrx'
        row = pd.Series({'chr': f'{chromosome}', 'start': 1,
                        'end': 10, 'info': 'test', 'gene': 'BRCA1'})
        # Check error is raised:
        with pytest.raises(ValueError) as excinfo:
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)
        assert str(
            excinfo.value) == f"Error: {chromosome.replace('chr', '')} is not a valid chromosome, use uppercase X or Y"

    def test_convert_bed_to_vcf(self):
        # Test convert_bed_to_vcf method
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        vcf_df = vcf_obj.convert_bed_to_vcf()
        expected_df = pd.read_csv(
            f"{TEST_DATA_DIR}/expected_output.vcf", sep="\t")
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
        expected_df = expected_df.astype(columns_to_convert)
        pd.testing.assert_frame_equal(expected_df, vcf_df)


if __name__ == '__main__':
    unittest.main()
