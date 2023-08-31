import unittest
import sys
sys.path.append('../gene_annotation2bed')
from gene_annotation2bed import extract_hgnc_id


class TestExtractHGNCID(unittest.TestCase):
    def test_extract_hgnc_id_found(self):
        attributes_str = "Dbxref=GeneID:123,HGNC:456"
        result = extract_hgnc_id(attributes_str)
        self.assertEqual(result, 456)

    def test_extract_hgnc_id_not_found(self):
        attributes_str = "Dbxref=GeneID:123"
        result = extract_hgnc_id(attributes_str)
        self.assertIsNone(result)

    def test_extract_hgnc_id_multiple_entries(self):
        """
        Currently this passes and selects the first.
        TODO: review if this is the desired behaviour.
        """
        attributes_str = "Dbxref=GeneID:123,HGNC:456,HGNC:789"
        result = extract_hgnc_id(attributes_str)
        self.assertEqual(result, 456)


if __name__ == '__main__':
    unittest.main()
