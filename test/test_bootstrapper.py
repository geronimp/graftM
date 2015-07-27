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
import subprocess
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.bootstrapper import Bootstrapper

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    def test_hello_world(self):
        boots = Bootstrapper(search_hmm_files = [os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')],
                             evalue='1e-5',
                             maximum_range=1000,
                             threads=1)
        with tempfile.NamedTemporaryFile() as tf:
            self.assertEqual(True,
                             boots.generate_hmm_from_contigs(\
                                [os.path.join(path_to_data,'bootstrapper','contigs.fna')],
                                tf.name))
            
            self.assertEqual("HMMER3/f [3.1b2 | February 2015]\n",
                             subprocess.check_output("head -n1 %s" % tf.name,
                                                     shell=True))
            self.assertEqual(subprocess.check_output("tail -n+17 %s" %
                                                     os.path.join(path_to_data,
                                                                  'bootstrapper',
                                                                  'expected.hmm'),
                                                     shell=True),
                             "".join(open(tf.name).readlines()[16:]))
            
    def test_no_hits(self):
        boots = Bootstrapper(search_hmm_files = [os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')],
                             evalue='1e-5',
                             maximum_range=1000,
                             threads=1)
        with tempfile.NamedTemporaryFile(suffix=".fasta") as contigs:
            contigs.write(">contig\n")
            contigs.write('A'*300+"\n")
            contigs.flush()
            with tempfile.NamedTemporaryFile() as tf:
                self.assertEqual(False,
                             boots.generate_hmm_from_contigs(\
                                [contigs.name],
                                tf.name))
        


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
    
