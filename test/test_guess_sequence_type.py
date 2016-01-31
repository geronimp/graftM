#!/usr/bin/env python

#=======================================================================
# Authors: Ben Woodcroft, Joel Boyd
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import os
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.unpack_sequences import UnpackRawReads

class Tests(unittest.TestCase):
    def test__guess_sequence_type(self):
        
        test_read=""">HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444
GTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGACGTACAGTTATCGTCGGTCTCGACCAGGCCCACGAGACGCTGACCAGAAGACTCGGTAAGGAAGTTACCCCTGAGACCATCAACGGCTATCTTGAGGCGTTGAACCAC"""
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(test_read)
            fasta.flush()

            urr = UnpackRawReads(fasta.name)
            
            self.assertEqual('aminoacid', urr._guess_sequence_type_from_string('P'*10))
            self.assertEqual('aminoacid', urr._guess_sequence_type_from_string('P'*10+'T'*89))
            self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('P'*10+'T'*90))
            self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('A'*300+'E'*999)) #only look at the first 300bp
            self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('a'*10+'T'*89)) #lowercase

    def test_find_slash_endings(self):
        slash_read_1 = """>FCC0WM1ACXX:2:1101:2167:2180#GTCCAGAA/1
CNTGCAGCCAAGTTGGCCGTTTCCGGCGCGATTGCAGATAAAAGCGCAGGCTGCTCCAAGGATAGGACCCCAGCTTCTGTTCAGGCCTCAGCAACTTCGC"""
        slash_read_2 = """>FCC0WM1ACXX:2:1101:2167:2180#GTCCAGAA/2
CNTGCAGCCAAGTTGGCCGTTTCCGGCGCGATTGCAGATAAAAGCGCAGGCTGCTCCAAGGATAGGACCCCAGCTTCTGTTCAGGCCTCAGCAACTTCGC"""
        false_slash_1 = """>FCC0WM1ACXX:2:1101:2167:2180#GTCCAGAA/3
CNTGCAGCCAAGTTGGCCGTTTCCGGCGCGATTGCAGATAAAAGCGCAGGCTGCTCCAAGGATAGGACCCCAGCTTCTGTTCAGGCCTCAGCAACTTCGC"""
        false_slash_2 = """>FCC0WM1ACXX:2:1101:2167:2180#GTCCAGAA/a
CNTGCAGCCAAGTTGGCCGTTTCCGGCGCGATTGCAGATAAAAGCGCAGGCTGCTCCAAGGATAGGACCCCAGCTTCTGTTCAGGCCTCAGCAACTTCGC"""
        reads = [slash_read_1, slash_read_2, false_slash_1, false_slash_2]
        truth = [True, True, False, False]
        
        for read, expected in zip(reads, truth):
            with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
                fasta.write(read)
                fasta.flush()
                urr = UnpackRawReads(fasta.name)
                self.assertEqual(expected, urr.slash_endings)
                
if __name__ == "__main__":
    unittest.main()
