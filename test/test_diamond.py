#!/usr/bin/env python3

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
import sys
import os
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.diamond import Diamond
from graftm.unpack_sequences import UnpackRawReads

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    _protein_query = '''>seq1 637699780 METHANOSARCINA BARKERI MCRA
SKFKKDMEVKFGGDITDKTAKFLRLGPEQDPRKVEMIKAGKE
IAEKRGIAFYNPMMHSGAPLGQRAITPYTISGTDIVCEPDDLHYVNNAAMQQMWDDIRRT
CIVGLDMAHETLEKRLGKEVTPETINHYLEVLNHAMPGAAVVQEMMVETHPALVDDCYVK
VFTGDDALADEIDKQFLIDINKEFSEEQAAQIKASIGKTSWQAIHIPTIVSRTTDGAQTS
RWAAMQIGMSFISAYAMCAGEAAVADLSFAAKHAALVSMGEMLPARRARGPNEPGGLSFG
HLSDIVQTSRVSEDPAKIALEVVGAGCMLYDQIWLGSYMSGGVGFTQYATAAYTDDILD
NNTYYDVDYINDKYNKDNKVKASLEVVKDIATESTLYGIETYEKFPTALEDHFG
GSQRATVLAAAAGVACSLATGNANAGLSGWYLSMYLHKEAWGRLGFFGFDLQDQCGATNV
LSYQGDEGLPDELRGPNYPNYAMNVGHQGGYAGIAQAAHSGRGDAFTVNPLLKVCFADDL
LPFNFAEPRREFGRGAIREFVPAGERSLVIPAK
>seq2 638201361 METHANOCALDOCOCCUS JANNASCHII MRTA
DVEKKLFLKALKEKFEEDPKEKYTKFYIFGGWRQSARKREFVEFAQKLIEKRGIPF
YNPDIGVPLGQRKLMTYKISGTDAFVEGDDLHFCNNAAIQQLVDDIKRTVIVGMDTAH
AVLEKRLGVEVTPETINEYMETINHALPGGAVVQEHMVEVHPGLVWDCYAKIFTGNDELA
DEIDKRFLIDINKEFPEEQAEQIKKYIGNRTYQVSRVPTIVVRCCDGGTVSRWSAMQIGM
SFITAYKLCAGEAAIADFSYAAKHADVIQMGMILPARRARGPNEPGGVPFGIFADIIQTS
RVSDDPAQVTLEVIGAAATFYDQVWLGSYMSGGVGFTQYASATYTDDILDDFVYYGMDY
VEKKYGLCGVKPSMEVVKDIATEVTLYGLEQYDEYPALLEDHFGGSQRAGVTAAAAGCS
VAFATGNSNAGINGWYLSQILHKEYHSRLGFYGYDLQDQCGAANSLSIRSDEGLLHECRG'''

    def assertEquivalent(self, list1, list2):
        pairs = [(item1, item2) for item1, item2 in zip(list1, list2)]

        for res1, res2 in pairs:
            length = len(res1)
            self.assertEqual(
                [res1[i] for i in range(0, length) if not i in [9, 10]],
                [res2[i] for i in range(0, length) if not i in [9, 10]])

            self.assertAlmostEqual(float(res1[9]), float(res2[9]), 20)
            self.assertAlmostEqual(float(res1[10]), float(res2[10]), -2)


    def test_blastp(self):
        with tempfile.NamedTemporaryFile() as f:
            f.write(self._protein_query.encode())
            f.flush()
            d = Diamond(os.path.join(path_to_data,'diamond','mcra.faa.dmnd'))
            res = d.run(f.name, UnpackRawReads.PROTEIN_SEQUENCE_TYPE)
            # ben@ben:~/git/graftM.local/test$ diamond view -a /tmp/a
            # seq1    637699780    100.0    548    0    0    1    548    1    548    0.0e+00    1111.3
            # seq2    638201361    100.0    472    0    0    1    472    1    472    2.9e-283    963.0
            self.assertEquivalent(

                             [['seq1', '637699780', '100', '548', '0', 1, 548, 1, 548, '0.0e+00', '1103.6', True],
                              ['seq2', '638201361', '100', '472', '0', 1, 472, 1, 472, '1.1e-282', '961.1', True]],
                             list([x[:-1] for x in res.each(res.fields)])
                             )

    def test_basename(self):
        with tempfile.NamedTemporaryFile() as f:
            f.write(self._protein_query.encode())
            f.flush()
            d = Diamond(os.path.join(path_to_data,'diamond','mcra.faa.dmnd'))
            base = 'mybase'
            daa = "%s.daa" % base
            res = d.run(f.name, UnpackRawReads.PROTEIN_SEQUENCE_TYPE, daa_file_basename=base)
            # ben@ben:~/git/graftM.local/test$ diamond view -a /tmp/a
            # seq1    637699780    100.0    548    0    0    1    548    1    548    0.0e+00    1111.3
            # seq2    638201361    100.0    472    0    0    1    472    1    472    2.9e-283    963.0

            self.assertEquivalent(
                             [['seq1', '637699780', '100', '548', '0', 1, 548, 1, 548, '0.0e+00', '1103.6', True],
                              ['seq2', '638201361', '100', '472', '0', 1, 472, 1, 472, '1.1e-282', '961.1', True]],
                             list([x[:-1] for x in res.each(res.fields)])
                             )
            self.assertTrue(os.path.exists(daa))
            os.remove("%s.daa" % base)


if __name__ == "__main__":
    unittest.main()
