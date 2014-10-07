#!/usr/bin/env python

#=======================================================================
# Author:
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
import sys
import subprocess
import tempfile
import os.path
import tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
  def test16S(self):
    data = os.path.join(path_to_data,'1_small_16S_tree')
    with tempdir.TempDir() as tmp:
      cmd = '%s -m %s/16S.hmm -c %s/16S_gg201308_82.refpkg -f %s/two_examples.fa -t dna -g %s/16S_gg201308_82.refpkg/82_otus.fasta -o %s' % (
        path_to_script,
        data, data, data, data,
        tmp
      )
      output = subprocess.check_output(cmd, shell=True)
      otuTableFile = os.path.join(tmp, 'two_examples_otu_table.txt')
      lines = (
        "\t".join(('#OTU_ID','two_examples','ConsensusLineage')),
          "\t".join(('0','1','Root;k__Bacteria;p__Gemmatimonadetes;c__Gemm-1')),
          "\t".join(('1','1','Root;k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae')),
      )
      count = 0
      for line in open(otuTableFile):
        self.assertEqual(lines[count], line.strip())
        count += 1
      self.assertEqual(count, 3)


if __name__ == "__main__":
	unittest.main()
