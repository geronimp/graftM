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

from graftm.graftm_package import GraftMPackage
from graftm.external_program_suite import ExternalProgramSuite
from graftm.update import Update
from graftm.sequence_io import SequenceIO

prerequisites = ExternalProgramSuite(['taxit', 'FastTreeMP', 'hmmalign', 'mafft'])
path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


class Tests(unittest.TestCase):
    extra_mcra_fasta = '''>KYC55281.1 Methyl-coenzyme M reductase I subunit alpha [Arc I group archaeon ADurb1013_Bin02101]
MVYKDDKHNFMQAMKKKFEEAPDKRQTKFYVYGGYKQNKRKVEFHDAGQQIAKERGIPGYNPSVGMPQGQ
RVLMPYQLSHTDIIANMDDLHFVNNAAMQQAWDDMRRTILVGLDSPHNILEKRLGKEVTPETINHYLEVV
NHSMPGAAVIQEHMVETDPRLVKDSYVKVYSGNDELIDEIDSRFVIDINKEFPADQAAQLKKAIGKSMWQ
VVRIPTVVGRVCDGGTVSRWSAMQIGMAFIGAYNLCAGEAATADFAYAAKHGSVIQMSDMMPARRARGPN
EPGGLMYGVVSDCAQSMAKYPDDPARHSLESIALAALIYDQIYLGSYMSGGVGFTQYATAAYTDNILEDF
VYWGMEHVKDKYGSLAKQKPSVKLINDIGTDVAMYCLEQYELYPAVMETHFGGSQRATCISAAAGTSVSM
ATGNAQAGLSAWYLACNVHKEQMGRFGFYGYDLQDQIGAANTFSYRSDEGLPFELRGGNYPSYAMNVGHQ
SAYTGIVAAAHSARGDAWALSPHVKVAFADRSLPFDFANITKEFGRGAMREFVPAGERDLIIP
'''
    extra_mcra_taxonomy = "KYC55281.1\tmcrA; Euryarchaeota_mcrA; Methanofastidiosa\n"
    def test_hello_world(self):
        with tempdir.in_tempdir():
            with tempfile.NamedTemporaryFile() as fasta:
                with tempfile.NamedTemporaryFile() as tax:
                    fasta.write(Tests.extra_mcra_fasta)
                    fasta.flush()
                    tax.write(Tests.extra_mcra_taxonomy)
                    tax.flush()
                    prev_path = os.path.join(path_to_data,'mcrA.10seqs.gpkg')
                    cmd1 = "%s update --graftm_package %s --sequences %s --taxonomy %s --output %s" %(
                        path_to_script,
                        prev_path,
                        fasta.name,
                        tax.name,
                        'updated.gpkg')
                    extern.run(cmd1)

                    prev = GraftMPackage.acquire(prev_path)
                    up = GraftMPackage.acquire('updated.gpkg')
                    prevhash = prev.taxonomy_hash()
                    taxhash = up.taxonomy_hash()
                    self.assertEqual(len(prevhash)+1,
                                     len(taxhash))
                    self.assertEqual(['mcrA','Euryarchaeota_mcrA','Methanofastidiosa'],
                                     taxhash['KYC55281.1'])
                    self.assertEqual(prevhash['638165755'],
                                     taxhash['638165755'])
                    seqio = SequenceIO()
                    self.assertEqual(
                        len(seqio.read_fasta_file(prev.unaligned_sequence_database_path()))+1,
                        len(seqio.read_fasta_file(up.unaligned_sequence_database_path())))

    def test_autodecorate(self):
        with tempdir.in_tempdir():
            with tempfile.NamedTemporaryFile() as fasta:
                fasta.write(Tests.extra_mcra_fasta)
                fasta.flush()
                
                prev_path = os.path.join(path_to_data,'mcrA.10seqs.gpkg')
                update = Update(prerequisites)
                update.update(
                    input_sequence_path = fasta.name,
                    input_graftm_package_path = prev_path,
                    output_graftm_package_path = 'updated.gpkg')
                prev = GraftMPackage.acquire(prev_path)
                up = GraftMPackage.acquire('updated.gpkg')
                prevhash = prev.taxonomy_hash()
                taxhash = up.taxonomy_hash()
                self.assertEqual(11, len(taxhash)) #hard-code 11 because of
                                                   #https://github.com/geronimp/graftM/issues/204
                self.assertEqual(['mcrA','Euryarchaeota_mcrA', 'Methanomicrobia'],
                                 taxhash['KYC55281.1'])
                
                self.assertEqual(prevhash['638165755'],
                                 taxhash['638165755'])
                seqio = SequenceIO()
                self.assertEqual(
                    len(seqio.read_fasta_file(prev.unaligned_sequence_database_path()))+1,
                    len(seqio.read_fasta_file(up.unaligned_sequence_database_path())))
                
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
