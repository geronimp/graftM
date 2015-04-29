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
import tempfile
from pandas.util.testing import assertRaises

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    
    def test_finds_reverse_complement(self):
        reads='''>NS500333:16:H16F3BGXX:1:11101:11211:1402 1:N:0:CGAGGCTG+CTCCTTAC
GAGCGCAACCCTCGCCTTCAGTTGCCATCAGGTTTGGCTGGGCACTCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGTGACAGTG
>NS500333:16:H16F3BGXX:1:11101:25587:3521 1:N:0:CGAGGCTG+CTCCTTAC
GAGTCCGGACCGTGTCTCAGTTCCGGTGTGGCTGGTCGTCCTCTCAGACCAGCTACGGATTGTCGCCTTGGTGAGCCATTACCTCACCAACTAGCTAATCCGACTTTGGCTCATCCAATAGCGCGAGGTCTTACGATCCCCCGCTTTCT'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'16S.gpkg')
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force --search_only' % (path_to_script,
                                                                                                                 data,
                                                                                                                 package,
                                                                                                                 tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>NS500333:16:H16F3BGXX:1:11101:11211:1402
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAGCGCAACCCTCGCCTTCAGTTGCCATCATGGGCACTCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>NS500333:16:H16F3BGXX:1:11101:25587:3521
-----------------------------------------------------------------------------------------------------------------------------------------------------------GCGCTATTGGATGAGCCAAAGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACAATCCGTAGCTGGTCTGAGAGGACGACCAGCCACACCGGAACTGAGACACGGTCCGGACTC--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))
                
                
    def test_multiple_hits_on_same_contig(self):
        contig = '''>AB11.qc.1_(paired)_contig_360665
CTGCTGGCGAAACTCGGCTGAAAAAAAGAGATAAAAAGTAGGCACGAATGCCTACTCTGT
TTTAAAATGTTTAATTATTGAGTTTATTTTGCTGGGATGATGAGAGATCTCTCACCAGCA
GGCACGAACTCTCTCATGGCACCGCGACCGAACTCTCTCCTGGGTTCAGCGAAGTTGAAG
GGCATGAGATCGTCAGCGAAGCAGGCCTTGATGAGCGGGTTGACGGTGAATGCGTCGCCG
CGGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGCCTGAGCGATACCTGCGTATCCACC
CTGGTGGCCGACGTTCATTGCGTAGTTTGGATAGTTTGGACCACGGAGTTCGTCTGGGAG
ACCTTCGTCGCCCTGGTAGGACAGAACGTTTGTGGCACCGCACTGGTCCTGCAGGTCGTA
TCCGAAGAAGCCGAGACGGCCCCATGCTTCCTTGTGCAGGTACATGGAGAGGTACCAGCC
AGAGAGACCAGCATTTGCGTTTGCGGTTGCAAGAGCACATGCGACACCGGCTGCAGCTGC
GAGCACGGTTGCTCTCTGGGATCCACCGAAGTGGTCTTCAAGGGCAGTCGGGAACTTCTC
GTAGGTCTCGATACCATAGATTGTGGACTCGGTTGCGATGTCCTTTACGACTTCGAGAGT
TGCTTTTACCTTGTTG'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(contig)
            fasta.flush()

            data = fasta.name
            package = os.path.join(path_to_data,'mcrA.gpkg')

            with tempdir.TempDir() as tmp:
                cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, sample_name, '%s_count_table.txt' % sample_name)
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('0','2','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                         )
                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1
                self.assertEqual(count, 2)

                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>AB11.qc.1_(paired)_contig_360665_196_4_7
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------KVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGGSQRATVLAAAAGVACALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAMNVGHQGGYAGIAQAXX------------------------------------------------------
>AB11.qc.1_(paired)_contig_360665_87_6_1
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------XXRGDAFTVNPLIKACFADDLMPFNFAEPRREFGRGAMREFVPAGERSLIIPAK
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))
    # Tests on searching for rRNA sequence in nucleic acid sequence
    def test_single_forward_read_run_16S(self):
        data = os.path.join(path_to_data,'16S.gpkg', '16S_1.1.fa')
        package = os.path.join(path_to_data,'16S.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)



    def test_single_paired_read_run_16S(self):
        data_for = os.path.join(path_to_data,'16S.gpkg', '16S_1.1.fa')
        data_rev = os.path.join(path_to_data,'16S.gpkg', '16S_1.2.fa')
        package = os.path.join(path_to_data,'16S.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for,
                                                                                               data_rev,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1', 'reverse')))

    def test_multiple_forward_read_run_16S(self):
        data_for1 = os.path.join(path_to_data,'16S.gpkg', '16S_1.1.fa')
        data_for2 = os.path.join(path_to_data,'16S.gpkg', '16S_2.1.fa')
        package = os.path.join(path_to_data,'16S.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

            otuTableFile = os.path.join(tmp, '16S_2' , '16S_2_count_table.txt')
            lines = ("\t".join(('#ID','16S_2','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

    def test_multiple_paired_read_run_16S(self):
        data_for1 = os.path.join(path_to_data,'16S.gpkg', '16S_1.1.fa')
        data_for2 = os.path.join(path_to_data,'16S.gpkg', '16S_2.1.fa')
        data_rev1 = os.path.join(path_to_data,'16S.gpkg', '16S_1.2.fa')
        data_rev2 = os.path.join(path_to_data,'16S.gpkg', '16S_2.2.fa')
        package = os.path.join(path_to_data,'16S.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --reverse %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               data_rev1,
                                                                                               data_rev2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

            otuTableFile = os.path.join(tmp, '16S_2' , '16S_2_count_table.txt')
            lines = ("\t".join(('#ID','16S_2','ConsensusLineage')),
                     "\t".join(('0','1','k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                     "\t".join(('1','1','k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)


            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_2', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_2', 'reverse')))

    # Tests on searching for proteins in nucelic acid sequence
    def test_single_forward_read_run_McrA(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

    def test_single_paired_read_run_McrA(self):
        data_for = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        data_rev = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.2.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for,
                                                                                                           data_rev,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'reverse')))

    def test_multiple_forward_read_run_McrA(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)


    def test_multiple_paired_read_run_McrA(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.fna')
        data_rev1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.2.fna')
        data_rev2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.2.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --reverse %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2', 'reverse')))

    # tests on searching amino acid sequence....
    def test_multiple_forward_read_run_McrA_aa(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

    def test_multiple_paired_read_run_McrA_aa(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.faa')
        data_rev1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.2.faa')
        data_rev2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.2.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s,%s --reverse %s,%s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2', 'reverse')))

    def test_single_forward_read_run_McrA_aa(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

    def test_single_paired_read_run_McrA_aa(self):
        data_for = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        data_rev = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.2.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for,
                                                                                                           data_rev,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'reverse')))

    def test_search_only(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_only --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_hits.aln.fa')
            for line in open(alnFile):
                self.assertEqual(expected[count], line.strip())
                count += 1
            self.assertEqual(count, len(expected))

    def test_search_only_specifying_hmm(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_only --forward %s --search_hmm_files %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               hmm,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_hits.aln.fa')
            for line in open(alnFile):
                self.assertEqual(expected[count], line.strip())
                count += 1
            self.assertEqual(count, len(expected))

    def test_search_only_specifying_multiple_hmms_but_no_aln_file(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_only --forward %s --search_hmm_files %s %s --output_directory %s --force 2>/dev/null' % (path_to_script,
                                                                                               data,
                                                                                               hmm, hmm2,
                                                                                               tmp)
            with assertRaises(subprocess.CalledProcessError):
                subprocess.check_output(cmd, shell=True)

    def test_search_only_specifying_multiple_hmms_and_aln_file(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_only --forward %s --search_hmm_files %s %s --aln_hmm_file %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               hmm2,hmm,
                                                                                               hmm,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_hits.aln.fa')
            for line in open(alnFile):
                self.assertEqual(expected[count], line.strip())
                count += 1
            self.assertEqual(count, len(expected))


if __name__ == "__main__":
    unittest.main()
