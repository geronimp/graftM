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
from graftm.hmmer import Hmmer

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_samples = os.path.join(os.path.dirname(os.path.realpath(__file__)),'sample_runs')


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
                cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
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
                         "\t".join(('0','2','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)
            
    def test_multiple_hsps_in_same_orf_of_fastq_sequence(self):
        fq = '''@NS500333:6:H1124BGXX:2:11107:13774:3316 1:N:0:GATCAG
CGCTTCCAGGTCGTCACCGGCCAACTCGCGAACCCGTCGCGGATCAAACTCGTGCGGCGCAACATCGCCCGTGTCCGCACGCAGATCAGTAAGTTGCAGATCGACCGTGTCCGCGCTGACCTGAAGAACGAGTACCAGACGCTGATCCAGG
+
AAAAAFFFAFFFFFF<FFFFFFAAFFFFFF)FFFFAFFFFFFFFFFFFFFFFFFFFFFFF7FF7FFFFFFFF<FFFFFFFFFFFFFFFFFFFF<FFFFFFAFFAFFFFFFFFFFFFFF.AFFFAFFF<FFAFFFFF<FFFFAF<A.FFF7F'''
        hmm = os.path.join(path_to_data,'hmms','DNGNGWU00027.hmm')
        with tempdir.TempDir() as tmp_in:
            fastq = open(os.path.join(tmp_in, 'a.fq'),'w')
            fastq.write(fq)
            fastq.flush()
            fastq_gz = '%s.gz' % fastq.name
            subprocess.check_call('gzip %s' % fastq.name, shell=True)
            
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --forward %s --search_hmm_files %s --search_and_align_only --output_directory %s/out' % (path_to_script,
                                                                                                   fastq_gz,
                                                                                                   hmm,
                                                                                                   tmp)
                
                subprocess.check_call(cmd, shell=True)
                self.assertEqual(open(os.path.join(tmp,'out','a','a_hits.aln.fa')).read(),
                                 "\n".join(['>NS500333:6:H1124BGXX:2:11107:13774:3316_1_1_2',
                                  '---------------------------RFQVVTGQLANPSRIKLVRRNIARVRTQISK----',''])
                                 )



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
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
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
            cmd = '%s graft --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

            otuTableFile = os.path.join(tmp, '16S_2' , '16S_2_count_table.txt')
            lines = ("\t".join(('#ID','16S_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
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
            cmd = '%s graft --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               data_rev1,
                                                                                               data_rev2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, '16S_1' , '16S_1_count_table.txt')
            lines = ("\t".join(('#ID','16S_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 3)

            otuTableFile = os.path.join(tmp, '16S_2' , '16S_2_count_table.txt')
            lines = ("\t".join(('#ID','16S_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae')),
                     "\t".join(('1','1','Root; k__Bacteria; p__Gemmatimonadetes; c__Gemm-1')),
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
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_1', 'mcrA_1_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'mcrA_2', 'mcrA_2_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_2','ConsensusLineage')),
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
                     "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1', 'reverse')))

    def test_search_and_align_only(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_and_align_only --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
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

    def test_search_and_align_only_specifying_hmm(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_and_align_only --forward %s --search_hmm_files %s --output_directory %s --force' % (path_to_script,
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

    def test_search_and_align_only_specifying_multiple_hmms_but_no_aln_file(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_and_align_only --forward %s --search_hmm_files %s %s --output_directory %s --force 2>/dev/null' % (path_to_script,
                                                                                               data,
                                                                                               hmm, hmm2,
                                                                                               tmp)
            with self.assertRaises(subprocess.CalledProcessError):
                subprocess.check_output(cmd, shell=True)

    def test_search_and_align_only_specifying_multiple_hmms_and_aln_file(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --search_and_align_only --forward %s --search_hmm_files %s %s --aln_hmm_file %s --output_directory %s --force' % (path_to_script,
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
            
    def test_search_and_align_only_specifying_hmm_files_and_aln_file(self):
        
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')
    
        with tempfile.NamedTemporaryFile(suffix='.txt') as hmms:
            hmms.write(hmm)
            hmms.write("\n")
            hmms.write(hmm2)
            hmms.write("\n")
            hmms.flush()
            print hmm, hmm2
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --search_and_align_only --forward %s --search_hmm_list_file %s --aln_hmm_file %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   hmms.name,
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

    def test_min_orf_length(self):
        fa = '''>long_partial_mcra_488bp
GGTGGTGTCGGATTCACACAGTATGCTACAGCTGCATACACTGATGATATCCTCGACAATAACGTGTACT
ACAACATCGACTACATCAACGACAAGTACAAGGGTGCTGCAAACATCGGCACGGACAACAAGGTAAAAGC
AACTCTCGAAGTCGTAAAGGACATCGCAACCGAGTCCACAATCTATGGTATCGAGACCTACGAGAAGTTC
CCGACTGCCCTTGAAGACCACTTCGGTGGNTCCCAGAGAGCAACCGTGCTCGCAGCCGCAGCCGGTGTCN
GTAGTGCCCTCGCAACCGCAAACGCAAATGCCGGTCTCTCTGGCTGGTACCTCTCCATGTACCTGCACAA
GGAAGCATGGGGCCGTCTCGGCTTCTTCGGATACGACCTGCAGGACCAGTGCGGTGCCACAAATGTTCTG
TCCTACCAGGGCGACGAAGGTCTCCCAGACGAACTCCGTGGTCCAAACTATCCTAACTACGCAATGAA
>short_partial_mcra_235bp
GGTGGTGTCGGATTCACACAGTATGCTACAGCTGCATACACTGATGATATCCTCGACAATAACGTGTACT
ACAACATCGACTACATCAACGACAAGTACAAGGGTGCTGCAAACATCGGCACGGACAACAAGGTAAAAGC
AACTCTCGAAGTCGTAAAGGACATCGCAACCGAGTCCACAATCTATGGTATCGAGACCTACGAGAAGTTC
CCGACTGCCCTTGAAGACCACTTCG
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(fa)
            fasta.flush()

            data = fasta.name
            package = os.path.join(path_to_data,'mcrA.gpkg')

            with tempdir.TempDir() as tmp:
                cmd = '%s graft --min_orf_length 300 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, sample_name, '%s_count_table.txt' % sample_name)
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                         )
                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1
                self.assertEqual(count, 2)

                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>long_partial_mcra_488bp_1_1_3
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))

    def test_restrict_read_length(self):
        fa = '''>long_partial_mcra_488bp_with_As
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
GGTGGTGTCGGATTCACACAGTATGCTACAGCTGCATACACTGATGATATCCTCGACAATAACGTGTACT
ACAACATCGACTACATCAACGACAAGTACAAGGGTGCTGCAAACATCGGCACGGACAACAAGGTAAAAGC
AACTCTCGAAGTCGTAAAGGACATCGCAACCGAGTCCACAATCTATGGTATCGAGACCTACGAGAAGTTC
CCGACTGCCCTTGAAGACCACTTCGGTGGNTCCCAGAGAGCAACCGTGCTCGCAGCCGCAGCCGGTGTCN
GTAGTGCCCTCGCAACCGCAAACGCAAATGCCGGTCTCTCTGGCTGGTACCTCTCCATGTACCTGCACAA
GGAAGCATGGGGCCGTCTCGGCTTCTTCGGATACGACCTGCAGGACCAGTGCGGTGCCACAAATGTTCTG
TCCTACCAGGGCGACGAAGGTCTCCCAGACGAACTCCGTGGTCCAAACTATCCTAACTACGCAATGAA
>short_partial_mcra_235bp
GGTGGTGTCGGATTCACACAGTATGCTACAGCTGCATACACTGATGATATCCTCGACAATAACGTGTACT
ACAACATCGACTACATCAACGACAAGTACAAGGGTGCTGCAAACATCGGCACGGACAACAAGGTAAAAGC
AACTCTCGAAGTCGTAAAGGACATCGCAACCGAGTCCACAATCTATGGTATCGAGACCTACGAGAAGTTC
CCGACTGCCCTTGAAGACCACTTCG
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(fa)
            fasta.flush()

            data = fasta.name
            package = os.path.join(path_to_data,'mcrA.gpkg')

            with tempdir.TempDir() as tmp:
                cmd = '%s graft --restrict_read_length 102 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, sample_name, '%s_count_table.txt' % sample_name)
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('0','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                         )
                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1
                self.assertEqual(count, len(lines))

                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>short_partial_mcra_235bp_1_1_1
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYK------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))


    def test_concatenated_OTU_table(self):
        reads_1=os.path.join(path_to_samples, "sample_16S_1.1.fa")
        reads_2=os.path.join(path_to_samples, "sample_16S_2.1.fa")
        gpkg=os.path.join(path_to_data, "61_otus.gpkg")
        
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                  reads_1,
                                                                                                  reads_2,
                                                                                                  gpkg,
                                                                                                  tmp)
            subprocess.check_output(cmd, shell=True)
            comb_otu_table=os.path.join(tmp, "combined_count_table.txt")
            expected=('\t'.join(('#ID', 'sample_16S_1', 'sample_16S_2', 'ConsensusLineage')),
                      '\t'.join(('0','1','2','Root; k__Bacteria; p__Cyanobacteria; c__Chloroplast')),
                      '\t'.join(('1','35','42','Root; k__Bacteria; p__Proteobacteria')),
                      '\t'.join(('2','0','2','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__mitochondria')),
                      '\t'.join(('3','42','120','Root; k__Bacteria')),
                      '\t'.join(('4','8','6','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales')),
                      '\t'.join(('5','4','18','Root')),
                      '\t'.join(('6','6','10','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria'))
            )
            count=0
            for line in open(comb_otu_table):
                self.assertEqual(expected[count], line.strip())
                count+=1
            self.assertEqual(count, len(open(comb_otu_table).readlines()))
            
                                                                                                                          

        
    def test_fastq_gz_input(self):
        reads='''@NS500333:16:H16F3BGXX:1:11101:11211:1402 1:N:0:CGAGGCTG+CTCCTTAC
GAGCGCAACCCTCGCCTTCAGTTGCCATCAGGTTTGGCTGGGCACTCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGTGACAGTG
+
DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
@NS500333:16:H16F3BGXX:1:11101:25587:3521 1:N:0:CGAGGCTG+CTCCTTAC
GAGTCCGGACCGTGTCTCAGTTCCGGTGTGGCTGGTCGTCCTCTCAGACCAGCTACGGATTGTCGCCTTGGTGAGCCATTACCTCACCAACTAGCTAATCCGACTTTGGCTCATCCAATAGCGCGAGGTCTTACGATCCCCCGCTTTCT
+
DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
'''
        with tempfile.NamedTemporaryFile(suffix='.fq.gz') as fastq_gz:
            with tempfile.NamedTemporaryFile(suffix='.fq') as fastq:
                fastq.write(reads)
                fastq.flush()
                cmd = 'gzip -c %s > %s' % (fastq.name, fastq_gz.name)
                subprocess.check_output(cmd, shell=True)
                data = fastq_gz.name
                package = os.path.join(path_to_data,'16S.gpkg')
                with tempdir.TempDir() as tmp:
                    cmd = '%s graft --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                     data,
                                                                                                                     package,
                                                                                                                     tmp)
                    subprocess.check_output(cmd, shell=True)
                    sample_name = os.path.basename(fastq_gz.name[:-6])
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

    def test_merg_aln(self):
        forward_reads='''>no_overlap
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_mismatch_in_reverse
-------CGTATGCAACCTACCTT---------------------------------------
>overlap_mismatch_in_forward
-------CGTATGCAACCTTCCTT---------------------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTATGCATCCTACCTT---------------------------------------'''
        reverse_reads='''>no_overlap
---------------------------------------CGTATGCAACCTACCTT-------
>overlap_all_match
----------------CCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_reverse
----------------CCTTCCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_forward
----------------CCTACCTTCAACCTACCTT----------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTTTGCAAGCTACCTT---------------------------------------'''
        expected_aln='''>no_overlap
-------CGTATGCAACCTACCTT---------------CGTATGCAACCTACCTT-------
>overlap_all_match
-------CGTATGCAACCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_reverse
-------CGTATGCAACCTACCTTCAACCTACCTT----------------------------
>overlap_mismatch_in_forward
-------CGTATGCAACCTTCCTTCAACCTACCTT----------------------------
>complete_overlap_all_match
-------CGTATGCAACCTACCTT---------------------------------------
>complete_overlap_mismatch
-------CGTATGCATCCTACCTT---------------------------------------'''.split()
        with tempfile.NamedTemporaryFile(suffix='_forward.fa') as forward_file:
            with tempfile.NamedTemporaryFile(suffix='_reverse.fa') as reverse_file:
                with tempfile.NamedTemporaryFile(suffix='.fa') as output_file:
                    forward_file.write(forward_reads)
                    reverse_file.write(reverse_reads)
                    forward_file.flush()
                    reverse_file.flush()
                    Hmmer(None).merge_forev_aln([forward_file.name,reverse_file.name],[output_file.name])
                    count = 0
                    for line in open(output_file.name):
                        self.assertEqual(expected_aln[count], line.strip())
                        count += 1
                    self.assertEqual(count, len(open(output_file.name).readlines()))
                    
if __name__ == "__main__":
    unittest.main()
