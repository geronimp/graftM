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
import tempfile
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.decoy_filter import DecoyFilter
from graftm.sequence_io import SequenceIO
from graftm.diamond import Diamond

class Tests(unittest.TestCase):
    eg1 = """>PROKKA_03952 Electron transfer flavoprotein-ubiquinone oxidoreductase
MSVESTTAPSPAPDLDRQTMEVDIACAGFGPAMGGFLTTLTRAWSENPADPAFESKAAPG
MPLQVLCYERADDIAAGVSGVVTRAQGIRASFPGLNPAEIPMAVEVTHERVLYLLDPIGA
SRRSLTLRAADRVLRALGPLLGARDHAFELPWTPAFLQKHGGLVLSIGQFNQWVGSQLMA
TGLVQIWPGTPVSAPIFPDAAGKDGPAEKTVLGLRLADQGVDRFGTPADGFMPGMDVRAH
LTVVGDGPVGAVSQSIDQTLGLPQGHAHREWALGMKMVIELPEDSPLQPGAVWHTFGYPE
PEIFGFLYVHPQRLASVGIFVPSWMSDPSRTAYRYLQHYIQHPALWRYLKDGILRSWGAK
SLEESGKRGEPFLAGDGYARIGESSGSTNMLTGSGVDEAWTTGTQLAEAVVELLRAGKPF
TRENLAATYETRRRASWVERGAREAQNARNGFHGGIVKGMIGMALAGLTRGHISFQAHIP
PAPKQIRPMSSRSAARLKLKNGAVLALENGRPLHDALLTARGWPEIAFDGRLLVTQQDAL
LMGGKVQAMPGFADHVVFRDRRLCIACEEKTCIAMCSGQAITEGADGVPAFEREKCVYCG
ACLWNCAQSSGGDRSNIDFLAGAGGLHSAEN
"""
    eg2 = """>PROKKA_03206_split_1 PROKKA_03206 Sulfite reductase [ferredoxin]
KREKNPWEAFDEVRAFARAGRSSVVPEWASYYFKWWGIYTQGDGVGATGGKGGEGLASDY
FMMRIAIPNGIVSAKQLRVIGSLTRKYARNLADITVRQNIQLHWLTIESLPEVVDALHAV
GLSPRGACGDVVRNVTGCPLAGVAADELIDASPFAAQIAHLLTANSAFYNLPRKFKISVT
GCPSWCAYPEINDIGLTAVKHNGEVGFSLRVGGGLSADPHLAVRLDAFILPDQALHVVRA
VTEIFRDQQGLRESRDRARLKHLFLKEGWTAGSFLAELQSRLDFALLPAAPEQPPIDVFR
DHAGIHPQRLPGLSYVGVTVLRGRLTGEQLQAIADLADRFGSGALRTTVSQNLLFIDIPS
RKAAELVRELGHIGLQVEGSPFWRGAVACTGTEFCKLAITETKGFTRWLVDELEERLPEF
DQQLKLNVTGCPNSCGQHWVADIGIEGKKIKHDGKLTDAYYFCLGGAVGLHAAIARPVGY
RCPAPLVPEAIERLLRQYLATRLRDENLRAWFSRHSNDELRAHLAGE
"""

    def test_hello_world(self):
        with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f1:
            with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f2:
                f1.write(self.eg1)
                f1.flush()
                extern.run("diamond makedb --in %s --db %s.dmnd" %\
                           (f1.name, f1.name))
                f2.write(self.eg1)
                f2.write(self.eg2)
                f2.flush()
                extern.run("diamond makedb --in %s --db %s.dmnd" %\
                           (f2.name, f2.name))
                with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f3:
                    with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f4:
                        f3.write(self.eg1)
                        f3.flush()
                        ret = DecoyFilter(
                            Diamond(f2.name+".dmnd"),
                            Diamond(f1.name+".dmnd")).filter(f1.name, f4.name)
                        self.assertEqual(True, ret)
                        seqs = SequenceIO().read_fasta_file(f4.name)
                        self.assertEqual(1, len(seqs))
                        self.assertEqual("PROKKA_03952", seqs[0].name)
                # clean up
                os.remove(f1.name+".dmnd")
                os.remove(f2.name+".dmnd")

    def test_no_decoys(self):
        with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f1:
            f1.write(self.eg1)
            f1.flush()
            extern.run("diamond makedb --in %s --db %s.dmnd" %\
                       (f1.name, f1.name))
            with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f3:
                with tempfile.NamedTemporaryFile(prefix='graftm_decoy_test') as f4:
                    f3.write(self.eg1)
                    f3.flush()
                    ret = DecoyFilter(Diamond(f1.name+".dmnd")).filter(f1.name, f4.name)
                    self.assertEqual(True, ret)
                    seqs = SequenceIO().read_fasta_file(f4.name)
                    self.assertEqual(1, len(seqs))
                self.assertEqual("PROKKA_03952", seqs[0].name)
        # clean up
        os.remove(f1.name+".dmnd")
        
if __name__ == "__main__":
    unittest.main()
