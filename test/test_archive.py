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
from graftm.archive import Archive

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


class Tests(unittest.TestCase):

    def test_round_trip(self):
        with tempdir.in_tempdir():
            archiver = Archive()
            ingpkg = os.path.join(path_to_data, 'mcrA.gpkg')
            arc = 'compressed.gpkg.tar.gz'
            out = 'reconstituted.gpkg'
            archiver.create(ingpkg, arc)
            archiver.extract(arc, out)
            self.assertEquals(
                ['CONTENTS.json', 'mcrA.faa', 'mcrA.faa.dmnd', 'mcrA.hmm', 'mcrA.refpkg'],
                sorted(os.listdir(out)))

if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
