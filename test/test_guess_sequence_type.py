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
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.unpack_sequences import UnpackRawReads

class Tests(unittest.TestCase):
    def test__guess_sequence_type(self):
        urr = UnpackRawReads(None)
        self.assertEqual('aminoacid', urr._guess_sequence_type_from_string('P'*10))
        self.assertEqual('aminoacid', urr._guess_sequence_type_from_string('P'*10+'T'*89))
        self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('P'*10+'T'*90))
        self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('A'*300+'E'*999)) #only look at the first 300bp
        self.assertEqual('nucleotide', urr._guess_sequence_type_from_string('a'*10+'T'*89)) #lowercase

    def test_stars(self):
        urr = UnpackRawReads(None)
        self.assertEqual('aminoacid', urr._guess_sequence_type_from_string('P'*10+"*"))


if __name__ == "__main__":
    unittest.main()
