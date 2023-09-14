import pandas as pd
import unittest
import sys
from io import StringIO

# set up the path to the module
sys.path.append('../gene_annotation2bed')
import gff2pandas as gff2pd

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 60)
pd.set_option('display.width', 150)


class Test_parsing_gff(unittest.TestCase):
    def setUp(self):
        """
        Set-up for the test cases.
        """
        self.gff_file = "tests/test_data/test3.gff"
        self.gff = gff2pd.read_gff3(self.gff_file)
        self.gff_df = self.gff.df
        self.gff_header = self.gff.header
        self.gff_new_df = self.gff.attributes_to_columns()

    def test_gff_parsing(self):
        """
        This tests:
        The correct instances are produced by the set-up.
        The correct columns are produced.
        The correct values are produced in the columns.
        """
        self.assertIsInstance(self.gff_df, pd.DataFrame)
        self.assertIsInstance(self.gff_header, str)
        self.assertEqual(self.gff_df.shape, (3, 9))
        self.assertEqual(self.gff_df.columns.tolist(), [
            "seq_id",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ])
        self.assertEqual(self.gff_df["seq_id"].tolist(), [
                         "NC_000001.10", "NC_000001.10", "NC_000001.10"])
        self.assertEqual(self.gff_df["source"].tolist(), [
                         "BestRefSeq", "BestRefSeq", "BestRefSeq"])
        self.assertEqual(self.gff_df["type"].tolist(), [
                         "pseudogene", "transcript", "exon"])
        self.assertEqual(self.gff_df["start"].tolist(), [14362, 14362, 29321])
        self.assertEqual(self.gff_df["end"].tolist(), [29370, 29370, 29370])
        self.assertEqual(self.gff_df["score"].tolist(), [".", ".", "."])
        self.assertEqual(self.gff_df["strand"].tolist(), ["-", "-", "-"])
        self.assertEqual(self.gff_df["phase"].tolist(), [".", ".", "."])
        self.assertEqual(self.gff_df["attributes"].tolist(), [
            "ID=gene-WASH7P;Dbxref=GeneID:653635,HGNC:HGNC:38034;Name=WASH7P;description=WASP family homolog 7%2C pseudogene;gbkey=Gene;gene=WASH7P;gene_biotype=transcribed_pseudogene;gene_synonym=FAM39F,WASH5P;pseudo=true",
            "ID=rna-NR_024540.1;Parent=gene-WASH7P;Dbxref=GeneID:653635,Genbank:NR_024540.1,HGNC:HGNC:38034;Name=NR_024540.1;gbkey=misc_RNA;gene=WASH7P;product=WASP family homolog 7%2C pseudogene;pseudo=true;transcript_id=NR_024540.1",
            "ID=exon-NR_024540.1-1;Parent=rna-NR_024540.1;Dbxref=GeneID:653635,Genbank:NR_024540.1,HGNC:HGNC:38034;gbkey=misc_RNA;gene=WASH7P;product=WASP family homolog 7%2C pseudogene;pseudo=true;transcript_id=NR_024540.1",
        ])

    def test_attr_to_cols(self):
        """
        test the gff attributes column to pandas cols.
        currently this is just checking the instance is a dataframe.
        """
        self.assertIsInstance(self.gff_new_df, pd.DataFrame)
        # Add further test here to check it
        # splits the attributes column into correctly.



class TestSplitAttributes(unittest.TestCase):
    """
    Test Class for testing the _split_atts function from gff2pd.py.
    """
    def setUp(self):
        """
        Set-up for the test cases.
        """
        self._split_atts = gff2pd._split_atts
    
    def test_empty_input(self):
        """Test splitting empty input."""
        atts = ""
        result = self._split_atts(atts)
        self.assertEqual(result, {})

    def test_single_attribute(self):
        """Test splitting a single attribute."""
        atts = "ID=12345"
        result = self._split_atts(atts)
        self.assertEqual(result, {"ID": "12345"})

    def test_multiple_attributes(self):
        """Test splitting multiple attributes."""
        atts = "ID=12345;Name=GeneX;Type=Protein"
        result = self._split_atts(atts)
        self.assertEqual(result, {"ID": "12345", "Name": "GeneX", "Type": "Protein"})

    def test_attributes_with_spaces(self):
        """Test splitting attributes with spaces in values."""
        atts = "ID=12345;Name=Gene X;Type=Protein"
        result = self._split_atts(atts)
        self.assertEqual(result, {"ID": "12345", "Name": "Gene X", "Type": "Protein"})

    def test_attributes_with_special_characters(self):
        """Test splitting attributes with special characters in values."""
        atts = "ID=12345;Note=This is a test!;Score=9.8"
        result = self._split_atts(atts)
        self.assertEqual(result, {"ID": "12345", "Note": "This is a test!", "Score": "9.8"})

    def test_dbxref_attribute(self):
        """Test splitting a Dbxref attribute."""
        atts = "ID=gene-WASH7P;Dbxref=GeneID:653635,HGNC:HGNC:38034;Name=WASH7P;description=WASP family homolog 7%2C pseudogene;gbkey=Gene;gene=WASH7P;gene_biotype=transcribed_pseudogene;gene_synonym=FAM39F,WASH5P;pseudo=true"
        result = self._split_atts(atts)
        expected = {
            "ID": "gene-WASH7P",
            "Dbxref": "GeneID:653635,HGNC:HGNC:38034",
            "Name": "WASH7P",
            "description": "WASP family homolog 7%2C pseudogene",
            "gbkey": "Gene",
            "gene": "WASH7P",
            "gene_biotype": "transcribed_pseudogene",
            "gene_synonym": "FAM39F,WASH5P",
            "pseudo": "true",
        }
        self.assertEqual(result, expected)

class TestAttributesToColumns(unittest.TestCase):

    def setUp(self):
        """
        Set-up for the test cases.
        """
        self.gff_file = "tests/test_data/test3.gff"
        self.gff = gff2pd.read_gff3(self.gff_file)
        self.gff_df = self.gff.df
        self.gff_header = self.gff.header
        self.gff_new_df = self.gff.attributes_to_columns()

    def test_attributes_to_columns(self):
        # Check if the resulting DataFrame has the expected columns
        expected_columns = ['seq_id', 'source', 'type', 'start', 'end',
                            'score', 'strand', 'phase', 'attributes',
                            'Dbxref', 'ID', 'Name', 'Parent', 'description',
                            'gbkey', 'gene', 'gene_biotype', 'gene_synonym',
                            'product', 'pseudo', 'transcript_id']
        self.assertListEqual(self.gff_new_df.columns.tolist(), expected_columns)

        # Check if the values in the resulting DataFrame match the provided GFF3 data
        self.assertEqual(self.gff_new_df.loc[0, "ID"], "gene-WASH7P")
        self.assertEqual(self.gff_new_df.loc[0, "Name"], "WASH7P")
        self.assertEqual(self.gff_new_df.loc[1, "ID"], "rna-NR_024540.1")
        self.assertEqual(self.gff_new_df.loc[1, "Name"], "NR_024540.1")


    def test_attributes_to_columns_empty_attributes(self):
        # Test when attributes are empty
        gff_file_empty_attr = "tests/test_data/test_empty_attributes.gff"
        gff_empty_attr = gff2pd.read_gff3(gff_file_empty_attr)
        gff_empty_attr_new_df = gff_empty_attr.attributes_to_columns()
        # Define the expected columns based on your GFF3 structure
        expected_columns = ["seq_id", "source", "type", "start","end",
                            "score", "strand", "phase", "attributes"]

        # Check that the number of columns in the DataFrame matches the expected number
        self.assertEqual(len(gff_empty_attr_new_df.columns), len(expected_columns))

        # Optionally, you can check if the column names match the expected ones
        self.assertListEqual(gff_empty_attr_new_df.columns.tolist(), expected_columns)

    def test_attributes_to_columns_missing_attributes(self):
        # Test when attributes are missing in some rows
        gff_file_missing_attr = "tests/test_data/test_missing_attributes.gff"
        gff_missing_attr = gff2pd.read_gff3(gff_file_missing_attr)
        gff_missing_attr_new_df = gff_missing_attr.attributes_to_columns()

        # Check if the resulting DataFrame has the expected columns
        expected_columns = ['seq_id', 'source', 'type', 'start', 'end', 'score',
                            'strand', 'phase', 'attributes', 'ID', 'Name', 'gene',
                            'transcript_id']
        self.assertListEqual(gff_missing_attr_new_df.columns.tolist(), expected_columns)

        # Check if the values in the resulting DataFrame match the provided GFF3 data
        self.assertEqual(gff_missing_attr_new_df.loc[0, "ID"], "gene-WASH7P")
        self.assertEqual(gff_missing_attr_new_df.loc[0, "Name"], "WASH7P")
        self.assertEqual(gff_missing_attr_new_df.loc[1, "ID"], "rna-NR_024540.1")
        self.assertEqual(gff_missing_attr_new_df.loc[1, "gene"], "WASH7P")
        self.assertEqual(gff_missing_attr_new_df.loc[1, "transcript_id"], "NR_024540.1")



if __name__ == "__main__":
    unittest.main()
