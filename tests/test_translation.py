import unittest
import flumut.flumut as fm


class TestTranslate(unittest.TestCase):

    def test_translate_all(self):
        sequence = 'TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTTGTCGTAGTGTCTTCCTCATCGCCTCCCCCACCGACTACCACAACGGCTGCCGCAGCGTATTACTAATAGCATCACCAACAGAATAACAAAAAGGATGACGAAGAGTGTTGCTGATGGCGTCGCCGACGGAGTAGCAGAAGGGGTGGCGGAGGG---'
        result = fm.translate(sequence)
        self.assertEqual(result, ['F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V', 'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A', 'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E', 'C', 'C', '*', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G', '-'])

    def test_translate_degen(self):
        sequence = 'RTTYTTSTTWTTKTTMTTBTTDTTHTTVTTNTT'
        result = fm.translate(sequence)
        self.assertEqual(result, ['IV', 'FL', 'LV', 'FI', 'FV', 'IL', 'FLV', 'FIV', 'FIL', 'ILV', '?'])

    def test_translate_del_inframe(self):
        sequence = 'ATG---CGT'
        result = fm.translate(sequence)
        self.assertEqual(result, ['M', '-', 'R'])

    def test_translate_del_frameshift(self):
        sequence = 'ATG-CGTTAC--ATG'
        result = fm.translate(sequence)
        self.assertEqual(result, ['M', 'R', 'Y', '-', 'M'])

    def test_translate_delayed_start(self):
        sequence = '---ATG---CGT'
        result = fm.translate(sequence)
        self.assertEqual(result, ['-', 'M', '-', 'R'])

    def test_translate_delayed_start_trunc(self):
        sequence = '--GATG---CGT'
        result = fm.translate(sequence)
        self.assertEqual(result, ['?', 'M', '-', 'R'])
        
    def test_translate_delayed_start_trunc_2(self):
        sequence = '----G---ATGC---C-GT'
        result = fm.translate(sequence)
        self.assertEqual(result, ['-', '?', '?', 'C', '-', 'R', '?'])

    pass
