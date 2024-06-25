import unittest
import flumut.flumut as fm
import re


class TestNameParser(unittest.TestCase):

    def test_groups(self):
        pattern = re.compile(r'(.+)_(.+)')
        name = 'my_sample_name_PB1'
        result = fm.parse_name(name, pattern, False)
        self.assertEqual(result, ('my_sample_name', 'PB1'))

    def test_named_groups(self):
        pattern = re.compile(r'(?P<sample>.+)_(?P<segment>.+)')
        name = 'my_sample_name_PB1'
        result = fm.parse_name(name, pattern, False)
        self.assertEqual(result, ('my_sample_name', 'PB1'))        

    def test_named_groups_inverted(self):
        pattern = re.compile(r'(?P<segment>.+?)_(?P<sample>.+)')
        name = 'PB1_my_sample_name'
        result = fm.parse_name(name, pattern, False)
        self.assertEqual(result, ('my_sample_name', 'PB1'))        

    def test_unmatching_pattern(self):
        pattern = re.compile(r'unmatching_pattern')
        name = 'my_sample_name_PB1'
        with self.assertRaises(SystemExit) as ar:
            result = fm.parse_name(name, pattern, False)
    
    def test_unmatching_pattern_force(self):
        pattern = re.compile(r'unmatching_pattern')
        name = 'my_sample_name_PB1'
        result = fm.parse_name(name, pattern, True)
        self.assertEqual(result, (None, None))
    
    def test_no_groups(self):
        pattern = re.compile(r'.+')
        name = 'my_sample_name_PB1'
        with self.assertRaises(SystemExit) as ar:
            result = fm.parse_name(name, pattern, False)
    
    def test_no_groups_force(self):
        pattern = re.compile(r'.+')
        name = 'my_sample_name_PB1'
        result = fm.parse_name(name, pattern, True)
        self.assertEqual(result, (None, None))
