import unittest
import flumut.flumut as fm


class TestPosition(unittest.TestCase):

    def test_adjust_normal(self):
        sequence = 'ACTG'
        result = fm.adjust_position(sequence, 1)
        self.assertEqual(result, 0)
        
    def test_adjust_end(self):
        sequence = 'ACTG'
        result = fm.adjust_position(sequence, 4)
        self.assertEqual(result, 3)
        
    def test_adjust_gap(self):
        sequence = '-CTG'
        result = fm.adjust_position(sequence, 1)
        self.assertEqual(result, 1)
        
    def test_adjust_gap_long(self):
        sequence = '---AT--GC'
        result = fm.adjust_position(sequence, 1)
        self.assertEqual(result, 3)
        
    def test_adjust_gap_multiple(self):
        sequence = '---AT-----GC'
        result = fm.adjust_position(sequence, 4)
        self.assertEqual(result, 11)
