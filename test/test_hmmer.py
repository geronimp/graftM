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
import tempfile
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.sequence_searcher import SequenceSearcher

class Tests(unittest.TestCase):
    def test_merg_aln(self):
        forward_reads='''>no_overlap
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_mismatch_in_reverse
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_mismatch_in_forward
-------CGTATGCAACCTTCCTT---------------------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTATGCATCCTACCTT---------------------------------------'''
        reverse_reads='''>no_overlap
---------------------------------------CGTATGCAACCTACCTT-------
>overlap_all_match
----------------CCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_reverse
----------------CCTTCCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_forward
----------------CCTACCTTCAACCTACCTT----------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTTTGCAAGCTACCTT---------------------------------------'''
        expected_aln='''>no_overlap
-------CGTATGCAACCTACCTT---------------CGTATGCAACCTACCTT-------
>overlap_all_match
-------CGTATGCAACCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_reverse
-------CGTATGCAACCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_forward
-------CGTATGCAACCTTCCTTCAACCTACCTT----------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTATGCATCCTACCTT---------------------------------------'''.split()
        with tempfile.NamedTemporaryFile(suffix='_forward.fa') as forward_file:
            with tempfile.NamedTemporaryFile(suffix='_reverse.fa') as reverse_file:
                with tempfile.NamedTemporaryFile(suffix='.fa') as output_file:
                    forward_file.write(forward_reads)
                    reverse_file.write(reverse_reads)
                    forward_file.flush()
                    reverse_file.flush()
                    SequenceSearcher(None).merge_forev_aln([forward_file.name],[reverse_file.name],[output_file.name])
                    count = 0
                    for line in open(output_file.name):
                        self.assertEqual(expected_aln[count], line.strip())
                        count += 1
                    self.assertEqual(count, len(open(output_file.name).readlines()))

if __name__ == "__main__":
    unittest.main()
