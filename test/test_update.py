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
from Bio import SeqIO
from Bio.Seq import Seq

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.create import Create
from graftm.graftm_package import GraftMPackageVersion2, GraftMPackage
from graftm.sequence_io import Sequence
from graftm.external_program_suite import ExternalProgramSuite

prerequisites = ExternalProgramSuite(['taxit', 'FastTreeMP', 'seqmagick', 'hmmalign', 'mafft'])
path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


class Tests(unittest.TestCase):

    def test_hello_world(self):
        with tempdir.in_tempdir():
            with tempfile.NamedTemporaryFile() as fasta:
                with tempfile.NamedTemporaryFile() as tax:
                    fasta.write('''>KYC55281.1 Methyl-coenzyme M reductase I subunit alpha [Arc I group archaeon ADurb1013_Bin02101]
MVYKDDKHNFMQAMKKKFEEAPDKRQTKFYVYGGYKQNKRKVEFHDAGQQIAKERGIPGYNPSVGMPQGQ
RVLMPYQLSHTDIIANMDDLHFVNNAAMQQAWDDMRRTILVGLDSPHNILEKRLGKEVTPETINHYLEVV
NHSMPGAAVIQEHMVETDPRLVKDSYVKVYSGNDELIDEIDSRFVIDINKEFPADQAAQLKKAIGKSMWQ
VVRIPTVVGRVCDGGTVSRWSAMQIGMAFIGAYNLCAGEAATADFAYAAKHGSVIQMSDMMPARRARGPN
EPGGLMYGVVSDCAQSMAKYPDDPARHSLESIALAALIYDQIYLGSYMSGGVGFTQYATAAYTDNILEDF
VYWGMEHVKDKYGSLAKQKPSVKLINDIGTDVAMYCLEQYELYPAVMETHFGGSQRATCISAAAGTSVSM
ATGNAQAGLSAWYLACNVHKEQMGRFGFYGYDLQDQIGAANTFSYRSDEGLPFELRGGNYPSYAMNVGHQ
SAYTGIVAAAHSARGDAWALSPHVKVAFADRSLPFDFANITKEFGRGAMREFVPAGERDLIIP
''')
                    fasta.flush()
                    tax.write("KYC55281.1\tRoot; mcrA; Euryarchaeota_mcrA; Methanofastidiosa\n")
                    tax.flush()
                    cmd1 = "%s update --graftm_package %s --sequences %s --taxonomy %s --output %s" %(
                        path_to_script,
                        os.path.join(path_to_data,'mcrA.gpkg'),
                        fasta.name,
                        tax.name,
                        'updated.gpkg')
                    extern.run(cmd1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
