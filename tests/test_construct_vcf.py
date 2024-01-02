import unittest
import pandas as pd
import pytest
import os
from unittest.mock import patch
from construct_vcf import ConstructVCF, parse_args

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


class TestConstructVCF(unittest.TestCase):
    def setUp(self):
        """
        Set-up for the test cases. Import the sample data.
        """
        self.bed_file_path = f"{TEST_DATA_DIR}/example_bed_hg38.bed"
        self.reference_path = f"{TEST_DATA_DIR}/hs37d5.fa"
        self.output_file = "output_new_test.vcf"
        self.bed_df = pd.read_csv(
            f"{TEST_DATA_DIR}/example_bed_hg38.bed",
            sep="\t", header=None,
            names=["chr", "start", "end", "info", "gene"]
        )

    def test_parse_args(self):
        """
        Test command line argument parsing
        """
        with patch("sys.argv", ["script.py", "-b", "tests/test_data/example_bed_hg38.bed", "-fasta", "tests/test_data/hs37d5.fa", "-o", "output_new_test.vcf"]):
            args = parse_args()
            self.assertEqual(
                args.bed_file, "tests/test_data/example_bed_hg38.bed")
            self.assertEqual(args.reference_file, "tests/test_data/hs37d5.fa")
            self.assertEqual(args.output_file, "output_new_test.vcf")

    def test_fetch_nucleotides_valid_chr(self):
        """
        Test the fetch_nucleotides method for expected valid chromosomes.

        Uses patch to avoid reference file dependency. The test case mocks the
        pysam.faidx function to simulate fetching nucleotide sequences for
        specific genomic coordinates.

        Test Parameters:
        - Reference path: Path to the reference file.
        - Bed file path: Path to the BED file.
        - Output file: Path to the output file.

        Mocked Data:
        - Mocked nucleotide sequences returned by pysam.faidx for specific
          coordinates on chromosome 14.

        Assertions:
        - Verifies the correct function calls to pysam.faidx for specified genomic
          coordinates.
        - Compares the results of the fetch_nucleotides method with expected
          nucleotide sequences.
        """

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

            # Assertions
            assert mock_faidx.call_args_list == [
                ((self.reference_path, "14:2-2"), {}),
                ((self.reference_path, "14:5-5"), {}),
                ((self.reference_path, "14:10-10"), {}),
            ]
            expected = [
                (result_df.iloc[0]['REF'], 'A'),
                (result_df.iloc[0]['ALT'], 'T'),
                (result_df.iloc[1]['REF'], 'T'),
                (result_df.iloc[1]['ALT'], 'A'),
                (result_df.iloc[2]['REF'], 'C'),
                (result_df.iloc[2]['ALT'], 'G')
            ]

            assert all([x[0] == x[1] for x in expected])

    def test_fetch_nucleotides_valid_sex_chr(self):
        """
        Test fetch_nucleotides method for valid X chromosome with mocked data.

        Uses patch to simulate pysam.faidx and avoid reference file dependency.
        Verifies correct handling of nucleotide sequences for specified coordinates
        on chromosome X. Mocks a ConstructVCF object, a sample row with chromosome X,
        and asserts the expected function calls and results.

        Mocked Data:
        - Mocked nucleotide sequences returned by pysam.faidx for specific
          coordinates on chromosome X.

        Assertions:
        - Verifies the correct function calls to pysam.faidx for specified genomic
          coordinates on chromosome X.
        - Compares the results of the fetch_nucleotides method with expected
          nucleotide sequences.
        """

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
            # Subtests for each coordinate

            # Assertions
            assert mock_faidx.call_args_list == [
                ((self.reference_path, "X:2-2"), {}),
                ((self.reference_path, "X:5-5"), {}),
                ((self.reference_path, "X:10-10"), {}),
            ]
            expected = [
                (result_df.iloc[0]['REF'], 'A'),
                (result_df.iloc[0]['ALT'], 'T'),
                (result_df.iloc[1]['REF'], 'T'),
                (result_df.iloc[1]['ALT'], 'A'),
                (result_df.iloc[2]['REF'], 'C'),
                (result_df.iloc[2]['ALT'], 'G')
            ]

            assert all([x[0] == x[1] for x in expected])

    def test_fetch_nucleotides_invalid_chr_character(self):
        """
        Tests if the correct error is raised by fetch_nucleotides
        when an invalid chromosome character is parsed, i.e. not chrX or chrY.
        """

        # Initialize ConstructVCF object
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        # Mock a sample row
        chromosome = 'chrA'
        row = pd.Series({'chr': f'{chromosome}', 'start': 1,
                        'end': 10, 'info': 'test', 'gene': 'BRCA1'})
        # Check error is raised:
        expected_error = f"Error: {chromosome.replace('chr', '')} is not a valid chromosome"

        with pytest.raises(ValueError, match=expected_error):
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)

    def test_fetch_nucleotides_invalid_chr_number(self):
        """
        Tests if the correct error is raised by fetch_nucleotides
        when an invalid chromosome integer is parsed i.e. Not Chr1-22.
        """

        # Initialize ConstructVCF object
        vcf_obj = ConstructVCF(
            self.bed_file_path, self.reference_path, self.output_file)
        # Mock a sample row
        chromosome = 'chr100'
        row = pd.Series({'chr': f'{chromosome}', 'start': 1,
                        'end': 10, 'info': 'test', 'gene': 'BRCA1'})
        # Check error is raised:
        expected_error = f"Error: {chromosome.replace('chr', '')} is not a valid chromosome"

        with pytest.raises(ValueError, match=expected_error):
            result_df = vcf_obj.fetch_nucleotides(row, self.reference_path)

    def test_convert_bed_to_vcf(self):
        """
        Tests if the convert_bed_to_vcf method returns the expected output
        """
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
