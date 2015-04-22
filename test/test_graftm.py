#!/usr/bin/env python

#=======================================================================
# Author: Ben Woodcroft
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
import subprocess
import os.path
import tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    
    # Tests on searching for rRNA sequence in nucleic acid sequence
    def test_single_forward_read_run_16S(self):
        data = os.path.join(path_to_data,'16S.gpkg', 'two_examples.fa')
        package = os.path.join(path_to_data,'16S.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script, 
                                                                                               data, 
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'two_examples' , 'two_examples_count_table.txt')
            lines = ("\t".join(('#ID','two_examples','ConsensusLineage')),
                    "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)
        
    #def test_single_paired_read_run_16S(self): pass
    
    #def test_multiple_forward_read_run_16S(self): pass
    
    #def test_multiple_paired_read_run_16S(self): pass
    
    # Tests on searching for proteins in nucelic acid sequence
    def test_single_forward_read_run_McrA(self): 
        data = os.path.join(path_to_data,'mcrA.gpkg', 'eg.fa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'eg', 'eg_count_table.txt')
            lines = ("\t".join(('#ID','eg','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)
        
    #def test_single_paired_read_run_McrA(self): pass
    
    #def test_multiple_forward_read_run_McrA(self): pass
    
    #def test_multiple_paired_read_run_McrA(self): pass

    # tests on searching amino acid sequence....
    #def test_multiple_forward_read_run_McrA_aa(self): pass
    
    #def test_multiple_paired_read_run_McrA_aa(self): pass

    #def test_single_forward_read_run_McrA_aa(self): pass
    
    #def test_single_paired_read_run_McrA_aa(self): pass



if __name__ == "__main__":
    unittest.main()
