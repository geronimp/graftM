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
import subprocess
import os.path
import tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


class Tests(unittest.TestCase):

    def test_hello_world(self):
        with tempdir.TempDir() as tmp:
            with tempdir.TempDir() as tmp2:
                cmd1 = "%s create --verbosity 1 --alignment %s --taxonomy %s --rerooted_tree %s --output %s" \
                    %(path_to_script,
                      os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                      os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                      os.path.join(path_to_data,'create','homologstre.tree'),
                      tmp)
                subprocess.check_call(cmd1, shell=True)
                cmd2 = "%s graft --verbosity 1 --graftm_package %s --forward %s --output_directory %s" \
                    % (path_to_script,
                       "%s.gpkg" % tmp,
                       os.path.join(path_to_data,'create','test.faa'),
                       tmp2+"_")
                subprocess.check_call(cmd2, shell=True)

    def test_rerooted_tree_with_node_names(self):
        with tempdir.TempDir() as tmp:
            with tempdir.TempDir() as tmp2:
                cmd1 = "%s create --verbosity 1 --alignment %s --taxonomy %s --rerooted_tree %s --output %s" \
                    %(path_to_script,
                      os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                      os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                      os.path.join(path_to_data,'create','decorated.tree'),
                      tmp)
                subprocess.check_call(cmd1, shell=True)
                cmd2 = "%s graft --verbosity 1 --graftm_package %s --forward %s --output_directory %s" \
                    % (path_to_script,
                       "%s.gpkg" % tmp,
                       os.path.join(path_to_data,'create','test.faa'),
                       tmp2+"_")
                subprocess.check_call(cmd2, shell=True)


if __name__ == "__main__":
    unittest.main()
