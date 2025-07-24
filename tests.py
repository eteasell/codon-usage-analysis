import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rscu.rscu import RSCU, extract_window

class TestCalculateRSCU(unittest.TestCase):
    
    def test_calculate_not_multiple_of_3(self):
        rscu = RSCU("ATGGT")
        self.assertEqual(rscu.rscu_values, {}, 'The dictionary should be empty when invalid sequence')
        
    def test_calculate_funky_codons(self):
        rscu = RSCU("ATGBBB")
        self.assertEqual(rscu.rscu_values, {}, 'The dictionary should be empty when invalid sequence')
        
    def test_calculate_valid_mini_sequence(self):
        rscu = RSCU("ATGGCCGCT")
        self.assertEqual(len(rscu.rscu_values.keys()), 61, 'Dict length incorrect')
        self.assertEqual(rscu.rscu_values['ATG'], 1.0, 'Values incorrect')
        self.assertEqual(rscu.rscu_values['GCC'], 2.0, 'Values incorrect')
        self.assertEqual(rscu.rscu_values['GCT'], 2.0, 'Values incorrect')
        
    def test_extract_window(self):
        record = SeqRecord(Seq("ATGGTCAAATTA"), id='id')
        w = 3
        d = 3
        cs = 3
        window  = extract_window(record, w, d, cs)
        self.assertEqual(window, "AAA", "window indexing is wrong")

if __name__ == '__main__':
    unittest.main()