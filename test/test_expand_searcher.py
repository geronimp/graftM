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

from graftm.graftm_package import GraftMPackage
from graftm.expand_searcher import ExpandSearcher


import unittest
import os.path
import sys
import tempfile
import subprocess
import logging
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    def test_hello_world(self):
        expandsearcher = ExpandSearcher(search_hmm_files = [os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')],
                             evalue='1e-5',
                             maximum_range=1000,
                             threads=1)
        with tempfile.NamedTemporaryFile() as tf:
            self.assertEqual(True,
                             expandsearcher.generate_expand_search_database_from_contigs(\
                                [os.path.join(path_to_data,'bootstrapper','contigs.fna')],
                                tf.name,
                                "hmmsearch"))

            self.assertTrue(subprocess.check_output("head -n1 %s" % tf.name,
                                        shell=True) in
                ["HMMER3/f [3.1b2 | February 2015]\n",
                 "HMMER3/f [3.2.1 | June 2018]\n"])
            self.assertEqual('NSEQ  2\n', open(tf.name).readlines()[10])

    def test_hello_world_diamond(self):
        gpkg=os.path.join(path_to_data, "bootstrapper", "D1_gpkg_for_diamond.gpkg")
        expandsearcher = ExpandSearcher(search_hmm_files = [os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')],
                             evalue='1e-5',
                             maximum_range=1000,
                             threads=1,
                             graftm_package=GraftMPackage.acquire(gpkg))
        with tempfile.NamedTemporaryFile() as tf:
            self.assertEqual(True,
                             expandsearcher.generate_expand_search_database_from_contigs(\
                                [os.path.join(path_to_data,'bootstrapper','diamond_bootstrap_contigs.fna')],
                                tf.name,
                                "diamond"))

    def test_no_hits(self):
        expandsearcher = ExpandSearcher(search_hmm_files = [os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')],
                             evalue='1e-5',
                             maximum_range=1000,
                             threads=1)
        with tempfile.NamedTemporaryFile(suffix=".fasta") as contigs:
            contigs.write(">contig\n")
            contigs.write('A'*300+"\n")
            contigs.flush()
            with tempfile.NamedTemporaryFile() as tf:
                self.assertEqual(False,
                             expandsearcher.generate_expand_search_database_from_contigs(\
                                [contigs.name],
                                tf.name,
                                "hmmsearch"))

    def test_bootstrap_executable(self):
        with tempfile.NamedTemporaryFile() as tf:
            cmd = '%s expand_search --verbosity 5 --contigs %s --output_hmm %s --search_hmm_files %s' % (path_to_script,
                                                                                              os.path.join(path_to_data,'bootstrapper','contigs.fna'),
                                                                                              tf.name,
                                                                                              os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm'))
            extern.run(cmd)
            self.assertTrue(
                subprocess.check_output("head -n1 %s" % tf.name,
                                        shell=True) in
                ["HMMER3/f [3.1b2 | February 2015]\n",
                 "HMMER3/f [3.2.1 | June 2018]\n"])
            self.assertEqual('NSEQ  2\n', open(tf.name).readlines()[10])


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()

