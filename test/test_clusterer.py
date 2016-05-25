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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import tempdir
import unittest
import os.path
import sys
import tempfile
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    
    def test_alignment_rereplication(self):
        gpkg = os.path.join(path_to_data,'61_otus.gpkg')
        test_sequences=""">FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA_2/1
ACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGGATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGACGCCGCGTGTGGGATGAAGGCCCTCGGGT
>FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/1
ACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGGATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGACGCCGCGTGTGGGATGAAGGCCCTCGGGT
>FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/2
CGGGGTATCTAATCCCGTTCGCTCCCCTAGCTTTCGTGCCTCAGCGTCAGAAAAGACCCAGTGAGCCGCTTTCGCCCCCGGTGTTCCTTAGGATATCAAC
"""
        expected_rereplicated_alignment=""">FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/2
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTGATATCCTAAGGAACACCGGGGGCGAAAGCGGCTCACTGGGTCTTCTGACGCTGAGGCACGAAAGCTAGGGGAGCGAACGGGATTAGATACCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA_2/1
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGGATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGACGCCGCGTGTGGGATGAAGGCCCTCGGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/1
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGGATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGACGCCGCGTGTGGGATGAAGGCCCTCGGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------""".split()
        with tempfile.NamedTemporaryFile(suffix=".fa") as tf:
            tf.write(test_sequences)
            tf.flush()
            
            with tempdir.TempDir() as tmp:
                cmd = "graftM graft --forward %s --graftm_package %s --output_directory %s --force" % (tf.name,
                                                                                                       gpkg,
                                                                                                       tmp)
                extern.run(cmd)
                
                filename=os.path.splitext(os.path.basename(tf.name))[0]
                observed_rereplicated_alignment = [x.strip() for x in open(os.path.join(tmp, filename, "%s_hits.aln.fa" % filename))]
                
                self.assertEquals(expected_rereplicated_alignment, 
                                  observed_rereplicated_alignment)


if __name__ == "__main__":
    unittest.main()
