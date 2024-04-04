import unittest
import src.MutFinder as mf


class TestPosition(unittest.TestCase):

    def test_adjust_normal(self):
        sequence = 'ACTG'
        result = mf.adjust_position(sequence, 1)
        self.assertEqual(result, 0)
        
    def test_adjust_end(self):
        sequence = 'ACTG'
        result = mf.adjust_position(sequence, 4)
        self.assertEqual(result, 3)
        
    def test_adjust_gap(self):
        sequence = '-CTG'
        result = mf.adjust_position(sequence, 1)
        self.assertEqual(result, 1)
        
    def test_adjust_gap_long(self):
        sequence = '---AT--GC'
        result = mf.adjust_position(sequence, 1)
        self.assertEqual(result, 3)
        
    def test_adjust_gap_multiple(self):
        sequence = '---AT-----GC'
        result = mf.adjust_position(sequence, 4)
        self.assertEqual(result, 11)
        
    def test_adjust_gap_end(self):
        sequence = '---AT-----'
        result = mf.adjust_position(sequence, 4)
        self.assertEqual(result, None)
