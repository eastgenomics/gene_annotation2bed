import pandas as pd
import unittest
import sys
from io import StringIO

# set up the path to the module
sys.path.append('../gene_annotation2bed')
import gff2pandas

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 60)
pd.set_option('display.width', 150)


class Test_parsing_gff(unittest.TestCase):
    def setUp(self):
        """
        Set-up for the test cases.
        """
        self.gff_file = "tests/test_data/test3.gff"
        self.gff = gff2pandas.read_gff3(self.gff_file)
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


if __name__ == "__main__":
    unittest.main()
