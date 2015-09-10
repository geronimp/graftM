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
import re
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
import graftm.hmmsearcher

class HmmsearcherTests(unittest.TestCase):
    path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
    
    def test_generates_correct_cmd_single_hmm(self):
        searcher = graftm.hmmsearcher.HmmSearcher(1)
        cmd = searcher._hmm_command('orfm some', [(['hmm1','out1'],1)])
        self.assertEqual('orfm some | hmmsearch  --cpu 1 -o /dev/null --noali --domtblout out1 hmm1 -', cmd)
        
    def test_generates_correct_cmd_extra_args(self):
        searcher = graftm.hmmsearcher.HmmSearcher(1, extra_args='-E 1e-99')
        cmd = searcher._hmm_command('orfm some', [(['hmm1','out1'],1)])
        self.assertEqual('orfm some | hmmsearch -E 1e-99 --cpu 1 -o /dev/null --noali --domtblout out1 hmm1 -', cmd)
        
    def test_generates_multi_hmms(self):
        searcher = graftm.hmmsearcher.HmmSearcher(1)
        cmd = searcher._hmm_command('orfm some', [(['hmm1','out1'],1), (['hmm2','out2'],2)])
        self.assertEqual('orfm some | tee >(hmmsearch  --cpu 1 -o /dev/null --noali --domtblout out1 hmm1 -) | hmmsearch  --cpu 2 -o /dev/null --noali --domtblout out2 hmm2 -', cmd)
        
    def test_munch_off_batch_single_cpu(self):
        searcher = graftm.hmmsearcher.HmmSearcher(1)
        queue = [['hmm1','out1'],['hmm2','out2']]
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm1','out1'],1]], pairs_to_run)
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm2','out2'],1]], pairs_to_run)
        
    def test_munch_off_batch_single_hmm(self):
        searcher = graftm.hmmsearcher.HmmSearcher(10)
        queue = [['hmm1','out1']]
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm1','out1'],10]], pairs_to_run)
        
    def test_munch_off_batch_multi_cpu(self):
        searcher = graftm.hmmsearcher.HmmSearcher(2)
        queue = [['hmm1','out1'],['hmm2','out2'],['hmm3','out3']]
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm1','out1'],1], [['hmm2','out2'],1]], pairs_to_run)
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm3','out3'],2]], pairs_to_run)
        
    def test_munch_off_batch_multi_cpu_heterogenous(self):
        searcher = graftm.hmmsearcher.HmmSearcher(5)
        queue = [['hmm1','out1'],['hmm2','out2']]
        pairs_to_run = searcher._munch_off_batch(queue)
        self.assertEqual([[['hmm1','out1'],3], [['hmm2','out2'],2]], pairs_to_run)
        
    def test_actually_runs(self):
        searcher = graftm.hmmsearcher.HmmSearcher(5)
        faa_file = os.path.join(self.path_to_data, 'mcrA.gpkg/mcrA_1.1.faa')
        hmm_file = os.path.join(self.path_to_data, 'mcrA.gpkg/mcrA.hmm')
        with tempfile.NamedTemporaryFile(suffix='.fa') as output:
            searcher.hmmsearch('cat %s' % faa_file, 
                               [hmm_file],
                               [output.name])
            expected = '''#                                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name         accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
# ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
example_partial_mcra8 -            162 mcrA.fasta           -            557     2e-88  286.6   5.0   1   1   2.5e-89   2.2e-88  286.4   5.0   332   487     1   162     1   162 0.99 -
#'''.split('\n')
            observed = re.sub(r'\t', ' ', open(output.name).read()).split("\n")
            self.assertEqual(expected, observed[:5])
        
    def test_nhmmer_searcher(self):
        searcher = graftm.hmmsearcher.NhmmerSearcher(20)
        cmd = searcher._hmm_command('cat some', [(['hmm1','out1'],10)])
        self.assertEqual('cat some | nhmmer  --cpu 10 -o /dev/null --noali --tblout out1 hmm1 -', cmd)        
        


if __name__ == "__main__":
    unittest.main()
