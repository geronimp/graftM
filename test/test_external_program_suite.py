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
import tempdir
import sys
import extern
import logging
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.external_program_suite import ExternalProgramSuite


def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
                return os.path.join(path)

    return None

class Tests(unittest.TestCase):

    def test_native_external_fasttree_command(self):
                
        if which("FastTree"):
        
            fasttree_original_path=which("FastTree")
            fasttree_possibilities = {'FastTreeMP':'fasttreeMP',
                                      'FastTree':'fasttree'}
        elif which("fasttree"):
            fasttree_original_path=which("fasttree")
            fasttree_possibilities = {'fasttreeMP':'FastTreeMP',
                                      'fasttree':'FastTree'}
        else:
            self.fail("No version of fasttree is installed that graftM can recognise!")
            
        if fasttree_original_path:
            for orig, possibility in fasttree_possibilities.iteritems():    
                with tempdir.TempDir() as tmp_file:
                    os.symlink(os.path.join(fasttree_original_path,
                                            orig), 
                               os.path.join(tmp_file,
                                            possibility))
                    
                    os.environ["PATH"] = os.environ["PATH"].replace(fasttree_original_path,
                                                                    tmp_file)
                    
                    commands = ExternalProgramSuite(['FastTreeMP'])
                    self.assertTrue(commands.fasttree in ['FastTreeMP',
                                                          'fasttree',
                                                          'fasttreeMP',
                                                          'FastTree'])
                    cmd = commands.fasttree
                    extern.run(cmd)
                    os.environ["PATH"] = os.environ["PATH"].replace(tmp_file,
                                                                    fasttree_original_path)
            
        
if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()