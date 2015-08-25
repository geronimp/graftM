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
from graftm.sequence_search_results import HMMSearchResult, SequenceSearchResult

class Tests(unittest.TestCase):
    def test_whacky_directions(self):
        hmmout = '''#                                                                                                       --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
2524288035           -            897 DNGNGWU00005_mingle_output_good_seqs_Tx7Fta.aln -            808   1.9e-16   49.0   0.7   1   2   2.5e-15   2.5e-15   45.4   0.0   283   406     8   148     2   156 0.78 G480DRAFT_04506 small GTP-binding protein domain [Blautia producta DSM 2950]
2524288035           -            897 DNGNGWU00005_mingle_output_good_seqs_Tx7Fta.aln -            808   1.9e-16   49.0   0.7   2   2      0.91      0.91   -2.8   0.1   142   142   706   706   601   815 0.58 G480DRAFT_04506 small GTP-binding protein domain [Blautia producta DSM 2950]
2524285235           -            530 DNGNGWU00005_mingle_output_good_seqs_Tx7Fta.aln -            808   2.1e-13   38.9   0.0   1   1   3.9e-13   3.9e-13   38.1   0.0   326   397    78   149    61   164 0.90 G480DRAFT_01694 bacterial peptide chain release factor 3 (bRF-3) [Blautia producta DSM 2950]
#
# Program:         hmmsearch
# Version:         3.1b2 (February 2015)
# Pipeline mode:   SEARCH
# Query file:      test/data/DNGNGWU00005_mingle_output_good_seqs.hmm
# Target file:     -
# Option settings: hmmsearch -o /dev/null --domtblout /tmp/single.graftm/single/single.hmmout.csv --noali -E 1e-5 --cpu 5 test/data/DNGNGWU00005_mingle_output_good_seqs.hmm - 
# Current dir:     /home/ben/git/graftM
# Date:            Sat Aug 22 22:23:06 2015
# [ok]
'''
        with tempfile.NamedTemporaryFile(prefix='graftm_test_sequence_search_result') as f:
            f.write(hmmout)
            f.flush()
            
            res = HMMSearchResult.import_from_hmmsearch_table(f.name)
            lres = list(res.each([SequenceSearchResult.QUERY_ID_FIELD,
                                    SequenceSearchResult.ALIGNMENT_DIRECTION]))
            self.assertEqual([['2524288035',True],['2524285235',True]], lres)
        
        
        
if __name__ == "__main__":
    unittest.main()
