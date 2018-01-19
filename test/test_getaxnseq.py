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
import os.path
import sys
import tempfile
from StringIO import StringIO
import string

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.getaxnseq import Getaxnseq

class Tests(unittest.TestCase):
    
    def test_hello_world(self):
        with tempfile.NamedTemporaryFile(prefix='graftm_test_getaxnseq') as tmp_seq:
            with tempfile.NamedTemporaryFile(prefix='graftm_test_getaxnseq') as tmp_tax:
                Getaxnseq().write_taxonomy_and_seqinfo_files({'seq1': ['k__me','p__you'],
                                                              'seq2': []},
                                                             tmp_tax.name,
                                                             tmp_seq.name)
                expected = "\n".join([','.join(p) for p in [['seqname','tax_id'],
                    ['seq2','Root'],
                    ['seq1','p__you']]])+"\n"
                self.assertEqual(expected, open(tmp_seq.name).read())
                expected = '\n'.join(["tax_id,parent_id,rank,tax_name,root,rank_0,rank_1",
                                      "Root,Root,root,Root,Root,,",
                                      "k__me,Root,rank_0,k__me,Root,k__me,",
                                      "p__you,k__me,rank_1,p__you,Root,k__me,p__you"])+"\n"
                self.assertEqual(expected, open(tmp_tax.name).read())

                
    def test_read_taxtastic_taxonomy_and_seqinfo(self):
        tax = StringIO('\n'.join(['tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species',
                                      'Root,Root,root,Root,Root,,,,,,,',
                                      'k__me,Root,kingdom,k__me,Root,k__me,,,,,,',
                                      'p__you,k__me,phylum,p__you,Root,k__me,p__you,,,,,'])+"\n")
        seq = StringIO("\n".join([','.join(p) for p in [['seqname','tax_id'],
                    ['seq2','Root'],
                    ['seq1','p__you']]])+"\n")
        self.assertEqual({'seq1': ['k__me','p__you'],
                          'seq2': []},
                         Getaxnseq().read_taxtastic_taxonomy_and_seqinfo(tax, seq))
        
    def test_more_than_seven_levels(self):
        with tempfile.NamedTemporaryFile(prefix='graftm_test_getaxnseq') as tmp_seq:
            with tempfile.NamedTemporaryFile(prefix='graftm_test_getaxnseq') as tmp_tax:
                Getaxnseq().write_taxonomy_and_seqinfo_files({'seq1': string.split('k__me p__you c__came over for great spaghetti extra'),
                                                              'seq1.5': string.split('k__me p__you c__came over for great spaghetti'),
                                                              'seq2': []},
                                                             tmp_tax.name,
                                                             tmp_seq.name)
                expected = "\n".join([','.join(p) for p in [['seqname','tax_id'],
                    ['seq2','Root'],
                    ['seq1','extra'],
                    ['seq1.5','spaghetti']]])+"\n"
                self.assertEqual(expected, open(tmp_seq.name).read())
                expected = '\n'.join(["tax_id,parent_id,rank,tax_name,root,rank_0,rank_1,rank_2,rank_3,rank_4,rank_5,rank_6,rank_7",
                                      "Root,Root,root,Root,Root,,,,,,,,",
                                      "k__me,Root,rank_0,k__me,Root,k__me,,,,,,,",
                                      "p__you,k__me,rank_1,p__you,Root,k__me,p__you,,,,,,",
                                      "c__came,p__you,rank_2,c__came,Root,k__me,p__you,c__came,,,,,",
                                      "over,c__came,rank_3,over,Root,k__me,p__you,c__came,over,,,,",
                                      "for,over,rank_4,for,Root,k__me,p__you,c__came,over,for,,,",
                                      "great,for,rank_5,great,Root,k__me,p__you,c__came,over,for,great,,",
                                      "spaghetti,great,rank_6,spaghetti,Root,k__me,p__you,c__came,over,for,great,spaghetti,",
                                      "extra,spaghetti,rank_7,extra,Root,k__me,p__you,c__came,over,for,great,spaghetti,extra"])+"\n"

                self.assertEqual(expected, open(tmp_tax.name).read())



if __name__ == "__main__":
    unittest.main()
