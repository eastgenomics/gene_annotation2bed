import sys
sys.path.append('../gene_annotation2bed')
import unittest
import gff2pandas
import pandas as pd

class Test_parsing_gff(unittest.TestCase):
    def setUp(self):
        self.gff_file = "tests/test_data/test3.gff"
        self.gff = gff2pandas.read_gff3(self.gff_file)
        self.gff_df = self.gff.df
        self.gff_header = self.gff.header
        self.gff_new_df = self.gff.attributes_to_columns()

    def test_gff_parsing(self):
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
        self.assertEqual(self.gff_df["seq_id"].tolist(), ["NC_000001.10", "NC_000001.10", "NC_000001.10"])
        self.assertEqual(self.gff_df["source"].tolist(), ["BestRefSeq", "BestRefSeq", "BestRefSeq"])
        self.assertEqual(self.gff_df["type"].tolist(), ["pseudogene", "transcript", "exon"])
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
        """
        self.assertIsInstance(self.gff_new_df, pd.DataFrame)
        print(self.gff_new_df)

if __name__ == "__main__":
    unittest.main()
