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
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import subprocess
import os.path
import tempdir
import tempfile
import extern
import json

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_samples = os.path.join(os.path.dirname(os.path.realpath(__file__)),'sample_runs')


class Tests(unittest.TestCase):

    def test_euk_check(self):
        reads = """>euk1
TCAAATGTCTGCCCTATCAACTATTGATGGTAGTGTAGAGGACTACCATGGTTGCGACGGGTAACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTG
>euk2
GTTGCGACGGGTAACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAATGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTACCCAA
>euk3
AGAAATGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTACCCAATCCCGGCACGGGGAGGTAGTGACGAGAAATAACAATATGGACCTCTCTAA
>euk4
TCCCGGCACGGGGAGGTAGTGACGAGAAATAACAATATGGACCTCTCTAACGATGGTCCATAATTGGAATGAGTTGAGCATAAATCCTTTTGCAAGGATC
>bac1
CGGCTTTCATGATGTAACGGGCGGTGTGTACAAGGCCTGGGAACGTATTCACGGCGCCGTAGCTGATGCGCCATTACTAGCGATTCCAACTTCATGTCGT
>bac2
GTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGACGCAACGCGAAGAACCTTACCTGGGCTTGACATGGTAG
>bac3
CTTGTGCGGGGACCCGTCAATTCCTTTGAGTTTTAACCTTGCGGCCGTAGTCCCCAGGCGGTGAACTTATCGCGTTAGCTTCGGCACCGACGGCTTTGAA
>bac4
CCCCGCCTTCCTCCCCTTTGTCAGAGGCAGTTCTGCTAGAGTGCGCCCCGATTGCTCGGGGGTAGCAACTAACCGTAAGGGTTGCGCTCGTTGCGGGACT"""
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'61_otus.gpkg')

            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --graftm_package %s --output_directory %s --force --euk_check' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                extern.run(cmd)

                sample_name = os.path.basename(fasta.name[:-3])
                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','1','Root; k__Bacteria')),
                         "\t".join(('2','1','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales')),
                         "\t".join(('3','1','Root; k__Bacteria; p__Proteobacteria')),
                         "\t".join(('4','1','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria')),
                         )

                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1

    def test_cluster(self):
        reads = """>HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444
GTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGACGTACAGTTATCGTCGGTCTCGACCAGGCCCACGAGACGCTGACCAGAAGACTCGGTAAGGAAGTTACCCCTGAGACCATCAACGGCTATCTTGAGGCGTTGAACCAC
>HWI-ST1243:121:D1AF9ACXX:8:1101:14973:25766
CATGACGTCAGTCCCGGAGACGTTGTACGGCATAAGTTGCCGCTGACCAAGCGTCATACCGCCGAGGTGCACGTTGACGTTGTATGCGGGAATTCCTCGGTCCTTGGCGATCTTTGCGCCCCATTCCTGGAACTCGCGCTTGTCCTTGGA
>HWI-ST1243:121:D1AF9ACXX:8:1101:4414:35570
GGCCGATATGGTCCAGTCTTCGAGGAAGTACCCGAACGACCCCGCCAGGCAAGCGCTTGAGACAGTAGCACTTGGTGCGGTAATCTTCGACCAGATCTACCTCGGTTCCTGCATGTCCGGCGGTGTTGGGTTCACCCAATACGCCACCGC
>HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444_2
GTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGACGTACAGTTATCGTCGGTCTCGACCAGGCCCACGAGACGCTGACCAGAAGACTCGGTAAGGAAGTTACCCCTGAGACCATCAACGGCTATCTTGAGGCGTTGAACCAC
>HWI-ST1243:121:D1AF9ACXX:8:1101:14973:25766_2
CATGACGTCAGTCCCGGAGACGTTGTACGGCATAAGTTGCCGCTGACCAAGCGTCATACCGCCGAGGTGCACGTTGACGTTGTATGCGGGAATTCCTCGGTCCTTGGCGATCTTTGCGCCCCATTCCTGGAACTCGCGCTTGTCCTTGGA"""
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'mcrA_with_dmnd.gpkg/')
            with tempdir.TempDir() as tmp:                
                cmd = '%s graft --verbosity 2 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                                 data,
                                                                                                                 package,
                                                                                                                 tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','4','Root; p__Euryarchaeota')),
                         "\t".join(('2','1','Root; p__Euryarchaeota; c__Methanomicrobia_2; o__[Methanomassiliicoccus]; f__[Methanomassiliicoccus]; g__Methanomassiliicoccus')),
                         )
                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1
                self.assertEqual(count, 3)
                
                placements_file = os.path.join(tmp, sample_name, "placements.jplace")
                open_placements_file = json.load(open(placements_file))
                self.assertEqual(len(open_placements_file["placements"]), 3)
                placements = sorted((u'HWI-ST1243:121:D1AF9ACXX:8:1101:4414:35570_2_2_3', 
                                     u'HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444_2_1_1_2',
                                     u'HWI-ST1243:121:D1AF9ACXX:8:1101:14973:25766_1_4_3'))
                for idx, placed_read_name in enumerate(sorted(tuple([y[0] for y in [x['nm'][0] for x in open_placements_file["placements"]]]))):
                    self.assertEqual(placements[idx], placed_read_name)
                                                
    
    def test_diamond_translate(self):
        reads='''>HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444
GTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGACGTACAGTTATCGTCGGTCTCGACCAGGCCCACGAGACGCTGACCAGAAGACTCGGTAAGGAAGTTACCCCTGAGACCATCAACGGCTATCTTGAGGCGTTGAACCAC
>HWI-ST1243:121:D1AF9ACXX:8:1101:14973:25766
CATGACGTCAGTCCCGGAGACGTTGTACGGCATAAGTTGCCGCTGACCAAGCGTCATACCGCCGAGGTGCACGTTGACGTTGTATGCGGGAATTCCTCGGTCCTTGGCGATCTTTGCGCCCCATTCCTGGAACTCGCGCTTGTCCTTGGA
>HWI-ST1243:121:D1AF9ACXX:8:1101:4414:35570
GGCCGATATGGTCCAGTCTTCGAGGAAGTACCCGAACGACCCCGCCAGGCAAGCGCTTGAGACAGTAGCACTTGGTGCGGTAATCTTCGACCAGATCTACCTCGGTTCCTGCATGTCCGGCGGTGTTGGGTTCACCCAATACGCCACCGC
>HWI-ST1243:121:D1AF9ACXX:8:1101:3408:44386
CTACAAGCTCTCCGGAACTGAGTATGTAGTAGAGGGCGACGATCTGCACTTCGTCAACAATCCAGCCATTCAGCAGATGTGGGATGATGTCAGGAGGACCATTGTAGTAGGTCTGGACATGGCCCACGAGGTTCTGCAAAAGCGCCTCGG
>HWI-ST1243:121:D1AF9ACXX:8:1101:4036:46970
CAACAAGTTGTTCCCTGCAAAGCAAGCTGGGCAACTTAAGAAGGCAATCGGCAAGACGCTCTGGCAGGGTGTGCACATCCCGACGATCGTCTCACGAACGTGCGACGGCGGTAAGACGAGCAGATGGGTTGCCTGGCAGACGTGTATGTG
>HWI-ST1243:121:D1AF9ACXX:8:1101:6810:63372
CATGACGTCAGTCCCGGAGACGTTGTACGGCATAAGTTGCCGCTGACCAAGCGTCATACCGCCGAGGTGCACGTTGACGTTGTATGCGGGAATTCCTCGGTCCTTGGCGATCTTTGCGCCCAATTCCTGGAACTCGCGCTTGTCCTTGGA
>HWI-ST1243:121:D1AF9ACXX:8:1101:13819:75491
CACGAGCAGATGGTCTGCTATGCAGATCGGTATGTCGTTCATCGCTGCATACAACATGTGTGCCGGTGAGGCAGCTGTGGCAGACCTGGCGTTTGCGCCCAAGCGCGCAAGCCTCTTGGGGGTGGCGGACGTGCTTCCACACTCCCGCAC
>HWI-ST1243:121:D1AF9ACXX:8:1101:17503:91050
CGCTTGGTCAGCGGCAACTTATGCCGTACAACGTCTCCGGGACTGACGTCATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGACGTACAGTTATCGTCGGTCTCGACCAGGCCC
>HWI-ST1243:121:D1AF9ACXX:8:1101:2683:99866
TCCCGGCAAAGAGTGCAGCAGCACTGAAGGCCACAGTCGGCAAGTCGATGTACCAGGCAGTGCACATCCCCACAACCGTCAGCAGGACCTGCGACGGCGGAACCACCTCCCGGTGGTCTGCAATGCAGATCGGTATGTCCTTCATTGGAG
>HWI-ST1243:121:D1AF9ACXX:8:1102:6995:2759
GATCGTATATGTCAGGCGGTGTCGGTTTCACGCAGTACGCAACTGCAGCATACACCGATAACATCCTCGACGACTACAGCTACTACGGCAACGACTACGCCAAGAAGTACGGAGCGGACGGAGCGGCGCCAGCGACGATGGACGGCGTTC'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'mcrA_with_dmnd.gpkg/')
            with tempdir.TempDir() as tmp:                
                cmd = '%s graft --verbosity 2 --search_method diamond --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                 data,
                                                                                                                 package,
                                                                                                                 tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                orfFile = os.path.join(tmp, sample_name, '%s_orf.fa' % sample_name)
                expected_aln = '''>HWI-ST1243:121:D1AF9ACXX:8:1101:4036:46970
NKLFPAKQAGQLKKAIGKTLWQGVHIPTIVSRTCDGGKTSRWVAWQTCM
>HWI-ST1243:121:D1AF9ACXX:8:1101:6810:63372
SKDKREFQELGAKIAKDRGIPAYNVNVHLGGMTLGQRQLMPYNVSGTDVM
>HWI-ST1243:121:D1AF9ACXX:8:1101:3408:44386
YKLSGTEYVVEGDDLHFVNNPAIQQMWDDVRRTIVVGLDMAHEVLQKRL
>HWI-ST1243:121:D1AF9ACXX:8:1101:2683:99866
PAKSAAALKATVGKSMYQAVHIPTTVSRTCDGGTTSRWSAMQIGMSFIG
>HWI-ST1243:121:D1AF9ACXX:8:1101:14973:25766
SKDKREFQEWGAKIAKDRGIPAYNVNVHLGGMTLGQRQLMPYNVSGTDVM
>HWI-ST1243:121:D1AF9ACXX:8:1101:13819:75491
TSRWSAMQIGMSFIAAYNMCAGEAAVADLAFAPKRASLLGVADVLPHSR
>HWI-ST1243:121:D1AF9ACXX:8:1101:17503:91050
LGQRQLMPYNVSGTDVMCEGDDLHYVNNPAMQQMWDEIRRTVIVGLDQA
>HWI-ST1243:121:D1AF9ACXX:8:1101:4414:35570
ADMVQSSRKYPNDPARQALETVALGAVIFDQIYLGSCMSGGVGFTQYAT
>HWI-ST1243:121:D1AF9ACXX:8:1101:12684:12444
VNNPAMQQMWDEIRRTVIVGLDQAHETLTRRLGKEVTPETINGYLEALNH
>HWI-ST1243:121:D1AF9ACXX:8:1102:6995:2759
SYMSGGVGFTQYATAAYTDNILDDYSYYGNDYAKKYGADGAAPATMDGV
'''.split()
                count = 0
                for line in open(orfFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(orfFile).readlines()))
                
                self.assertTrue(os.path.exists(os.path.join(tmp, sample_name, '%s_diamond_search.daa' % sample_name)), "should keep the diamond daa file")
                



    def test_split_for_reads_hitting_same_HMM_region_but_within_merge(self):
        read='''>1000
AACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAACAAATCTATACAGTGTAACTGCGAATGGCTCATTATATCAGTTATAA
TTTATTTGATGGTCCTTGCTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCATCAAGCCCCGACTTCTGGAAGGGGTGTATTTATTAGATGGAAACCAATGCGGGGAA
ACCCGGAACCTGGTGATTCATAATAACTTTCGGATCGATTTATTCGATGCATCATTCAAGTTTCTGCCCTATCAGCTTTGGATGGTAGGGTATTGGCCTACCATGGCTTTGACGG
GTAACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTACCCAATCCTGACACAGGGAGGTAGTGACAA
TAAATAACAATGCCGGGCTTTCTAGTCTGGCACTTGGAATGAGAACAATTTAAATCCCTTATCGAGGATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXAACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAACAAATCTATACAGTGTAACTGCGAAT
GGCTCATTATATCAGTTATAATTTATTTGATGGTCCTTGCTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCATCAAGCCCCGACTTCTGGAAGGGGTGTATTTATTA
GATGGAAACCAATGCGGGGAAACCCGGAACCTGGTGATTCATAATAACTTTCGGATCGATTTATTCGATGCATCATTCAAGTTTCTGCCCTATCAGCTTTGGATGGTAGGGTATT
GGCCTACCATGGCTTTGACGGGTAACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTACCCAATCCT
GACACAGGGAGGTAGTGACAATAAATAACAATGCCGGGCTTTCTAGTCTGGCACTTGGAATGAGAACAATTTAAATCCCTTATCGAGGATCAATTGGAGGGCAAGTCTGGTGCCA
GCAGCCGCGGTAATTC'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(read)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'61_otus.gpkg')
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                 data,
                                                                                                                 package,
                                                                                                                 tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                hits_file = os.path.join(tmp, sample_name, '%s_hits.fa' % sample_name)
                expected_reads='''>1000_split_1 1000\nAGTTTCTGCCCTATCAGCTTTGGATGGTAGGGTATTGGCCTACCATGGCTTTGACGGGTA\nACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGG\nAAGGCAGCAGGCGCGTAAATTACCCAATCCTGACACAGGGAGGTAGTGACAATAAATAAC\nAATGCCGGGCTTTCTAGTCTGGCACTTGGAATGAGAACAATTTAAATCCCTTATCGAGGA\nTCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCXX\n>1000_split_2 1000\nAGTTTCTGCCCTATCAGCTTTGGATGGTAGGGTATTGGCCTACCATGGCTTTGACGGGTA\nACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGG\nAAGGCAGCAGGCGCGTAAATTACCCAATCCTGACACAGGGAGGTAGTGACAATAAATAAC\nAATGCCGGGCTTTCTAGTCTGGCACTTGGAATGAGAACAATTTAAATCCCTTATCGAGGA\nTCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAAT'''.split('\n')
                count = 0
                for line in open(hits_file):
                    self.assertEqual(expected_reads[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(hits_file).readlines()))
        
            
    def test_finds_reverse_complement(self):
        reads='''>NS500333:16:H16F3BGXX:1:11101:11211:1402 1:N:0:CGAGGCTG+CTCCTTAC
GAGCGCAACCCTCGCCTTCAGTTGCCATCAGGTTTGGCTGGGCACTCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGTGACAGTG
>NS500333:16:H16F3BGXX:1:11101:25587:3521 1:N:0:CGAGGCTG+CTCCTTAC
GAGTCCGGACCGTGTCTCAGTTCCGGTGTGGCTGGTCGTCCTCTCAGACCAGCTACGGATTGTCGCCTTGGTGAGCCATTACCTCACCAACTAGCTAATCCGACTTTGGCTCATCCAATAGCGCGAGGTCTTACGATCCCCCGCTTTCT'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'61_otus.gpkg')
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                 data,
                                                                                                                 package,
                                                                                                                 tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>NS500333:16:H16F3BGXX:1:11101:11211:1402
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAGCGCAACCCTCGCCTTCAGTTGCCATCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    >NS500333:16:H16F3BGXX:1:11101:25587:3521
    -----------------------------------------------------------------------------------------------------------------------------------------------------------CGCTATTGGATGAGCCAAAGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACAATCCGTAGCTGGTCTGAGAGGACGACCAGCCACACCGGAACTGAGACACGGTCCGGACT--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))
                
    def test_reverse_complement_protein_db(self):
        reads='''>example_partial_mcra_revcom
TTCATTGCGTAGTTAGGATAGTTTGGACCACGGAGTTCGTCTGGGAGACCTTCGTCGCCCTGGTAGGACAGAACAT
TTGTGGCACCGCACTGGTCCTGCAGGTCGTATCCGAAGAAGCCGAGACGGCCCCATGCTTCCTTGTGCAGGTACAT
GGAGAGGTACCAGCCAGAGAGACCGGCATTTGCGTTTGCGGTTGCGAGGGCACTACNGACACCGGCTGCGGCTGCG
AGCACGGTTGCTCTCTGGGANCCACCGAAGTGGTCTTCAAGGGCAGTCGGGAACTTCTCGTAGGTCTCGATACCAT
AGATTGTGGACTCGGTTGCGATGTCCTTTACGACTTCGAGAGTTGCTTTTACCTTGTTGTCCGTGCCGATGTTTGC
AGCACCCTTGTACTTGTCGTTGATGTAGTCGATGTTGTAGTACACGTTATTGTCGAGGATATCATCAGTGTATGCA
GCTGTAGCATACTGTGTGAATCCGACACCACC
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()

            data = fasta.name
            package = os.path.join(path_to_data,'mcrA.gpkg')

            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                         )
                count = 0
                for line in open(otuTableFile):
                    self.assertEqual(lines[count], line.strip())
                    count += 1
                self.assertEqual(count, 2)

                alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                expected_aln = '''>example_partial_mcra_revcom_3_6_9
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------
'''.split()
                count = 0
                for line in open(alnFile):
                    self.assertEqual(expected_aln[count], line.strip())
                    count += 1
                self.assertEqual(count, len(open(alnFile).readlines()))        
        

    def test_two_files_one_no_sequences_hit_nucleotide(self):
        reads='''>NS500333:16:H16F3BGXX:1:11101:11211:1402 1:N:0:CGAGGCTG+CTCCTTAC
GAGCGCAACCCTCGCCTTCAGTTGCCATCAGGTTTGGCTGGGCACTCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGTGACAGTG
>NS500333:16:H16F3BGXX:1:11101:25587:3521 1:N:0:CGAGGCTG+CTCCTTAC
GAGTCCGGACCGTGTCTCAGTTCCGGTGTGGCTGGTCGTCCTCTCAGACCAGCTACGGATTGTCGCCTTGGTGAGCCATTACCTCACCAACTAGCTAATCCGACTTTGGCTCATCCAATAGCGCGAGGTCTTACGATCCCCCGCTTTCT'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            package = os.path.join(path_to_data,'61_otus.gpkg')
            with tempdir.TempDir() as tmp:
                with tempfile.NamedTemporaryFile(suffix='.fa') as empty_fasta:
                    empty_fasta.write('>seq\n')
                    empty_fasta.write('A'*1000+"\n")
                    empty_fasta.flush()
                    cmd = '%s graft --verbosity 2  --forward %s %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                     data,
                                                                                                                     empty_fasta.name,
                                                                                                                     package,
                                                                                                                     tmp)
                    subprocess.check_output(cmd, shell=True)
                    sample_name = os.path.basename(fasta.name[:-3])
                    alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                    expected_aln = '''>NS500333:16:H16F3BGXX:1:11101:11211:1402
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAGCGCAACCCTCGCCTTCAGTTGCCATCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        >NS500333:16:H16F3BGXX:1:11101:25587:3521
        -----------------------------------------------------------------------------------------------------------------------------------------------------------CGCTATTGGATGAGCCAAAGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACAATCCGTAGCTGGTCTGAGAGGACGACCAGCCACACCGGAACTGAGACACGGTCCGGACT--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
                cmd = '%s graft --verbosity 2 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','2','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
        data = os.path.join(path_to_data,'16S_inputs','16S_1.1.fa')
        package = os.path.join(path_to_data,'61_otus.gpkg')
        
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','16S_1.1','ConsensusLineage')),
                     "\t".join(('1','2','Root; k__Bacteria')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, len(lines))
            
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
                cmd = '%s graft --verbosity 2  --forward %s --search_hmm_files %s --search_and_align_only --output_directory %s/out' % (path_to_script,
                                                                                                   fastq_gz,
                                                                                                   hmm,
                                                                                                   tmp)
                
                subprocess.check_call(cmd, shell=True)
                self.assertEqual(open(os.path.join(tmp,'out','a','a_hits.aln.fa')).read(),
                                 "\n".join(['>NS500333:6:H1124BGXX:2:11107:13774:3316_1_1_2',
                                  '---------------------------RFQVVTGQLANPSRIKLVRRNIARVRTQISK----',''])
                                 )



    def test_single_paired_read_run_16S(self):
        data_for = os.path.join(path_to_data,'16S_inputs','16S_1.1.fa')
        data_rev = os.path.join(path_to_data,'16S_inputs','16S_1.2.fa')
        package = os.path.join(path_to_data,'61_otus.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for,
                                                                                               data_rev,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','16S_1.1','ConsensusLineage')),
                     "\t".join(('1','2','Root; k__Bacteria')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, len(lines))

            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1.1', 'reverse')))

    def test_multiple_forward_read_run_16S(self):
        data_for1 = os.path.join(path_to_data,'16S_inputs', '16S_1.1.fa')
        data_for2 = os.path.join(path_to_data,'16S_inputs', '16S_2.1.fa')
        package = os.path.join(path_to_data,'61_otus.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','16S_1.1','16S_2.1','ConsensusLineage')),
                     "\t".join(('1','2','2','Root; k__Bacteria')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, len(lines))

    def test_multiple_paired_read_run_16S(self):
        data_for1 = os.path.join(path_to_data,'16S_inputs','16S_1.1.fa')
        data_for2 = os.path.join(path_to_data,'16S_inputs','16S_2.1.fa')
        data_rev1 = os.path.join(path_to_data,'16S_inputs','16S_1.2.fa')
        data_rev2 = os.path.join(path_to_data,'16S_inputs','16S_2.2.fa')
        package = os.path.join(path_to_data,'61_otus.gpkg')
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data_for1,
                                                                                               data_for2,
                                                                                               data_rev1,
                                                                                               data_rev2,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','16S_1.1','16S_2.1','ConsensusLineage')),
                     "\t".join(('1','2','2','Root; k__Bacteria')),
                    )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, len(lines))

            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_1.1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_2.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, '16S_2.1', 'reverse')))

    # Tests on searching for proteins in nucelic acid sequence
    def test_single_forward_read_run_McrA(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','ConsensusLineage')),
                     "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)
            
    # Tests on searching for proteins in nucelic acid sequence
    def test_two_files_one_no_sequences_hit_protein(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            with tempfile.NamedTemporaryFile(suffix='.fa') as empty_fasta:
                empty_fasta.write('>seq\n')
                empty_fasta.write('A'*1000+"\n")
                empty_fasta.flush()
                
                cmd = '%s graft --verbosity 2  --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   empty_fasta.name,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID','mcrA_1.1','ConsensusLineage')),
                         "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --verbosity 2  --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for,
                                                                                                           data_rev,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','ConsensusLineage')),
                     "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'reverse')))

    def test_multiple_forward_read_run_McrA(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','mcrA_2.1','ConsensusLineage')),
                     "\t".join(('1','1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --verbosity 2  --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','mcrA_2.1','ConsensusLineage')),
                     "\t".join(('1','1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2.1', 'reverse')))

    # tests on searching amino acid sequence....
    def test_multiple_forward_read_run_McrA_aa(self):
        data_for1 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        data_for2 = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_2.1.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','mcrA_2.1','ConsensusLineage')),
                     "\t".join(('1','1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --verbosity 2  --forward %s %s --reverse %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for1,
                                                                                                           data_for2,
                                                                                                           data_rev1,
                                                                                                           data_rev2,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','mcrA_2.1','ConsensusLineage')),
                     "\t".join(('1','1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            subprocess.check_output(cmd, shell=True)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'reverse')))
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_2.1', 'reverse')))

    def test_single_forward_read_run_McrA_aa(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.faa')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','ConsensusLineage')),
                     "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --verbosity 2  --forward %s --reverse %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                           data_for,
                                                                                                           data_rev,
                                                                                                           package,
                                                                                                           tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            lines = ("\t".join(('#ID','mcrA_1.1','ConsensusLineage')),
                     "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
                     )
            count = 0
            for line in open(otuTableFile):
                self.assertEqual(lines[count], line.strip())
                count += 1
            self.assertEqual(count, 2)

            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'forward'))) # Check forward and reverse reads exist.
            self.assertTrue(os.path.isdir(os.path.join(tmp, 'mcrA_1.1', 'reverse')))

    def test_search_and_align_only(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        package = os.path.join(path_to_data,'mcrA.gpkg')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2 --search_and_align_only --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               package,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1.1', 'mcrA_1.1_hits.aln.fa')
            for line in open(alnFile):
                self.assertEqual(expected[count], line.strip())
                count += 1
            self.assertEqual(count, len(expected))

    def test_search_and_align_only_specifying_hmm(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2 --search_and_align_only --forward %s --search_hmm_files %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               hmm,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1.1', 'mcrA_1.1_hits.aln.fa')
            for line in open(alnFile):
                self.assertEqual(expected[count], line.strip())
                count += 1
            self.assertEqual(count, len(expected))

    def test_search_and_align_only_specifying_multiple_hmms_but_no_aln_file(self):
        data = os.path.join(path_to_data,'mcrA.gpkg', 'mcrA_1.1.fna')
        hmm = os.path.join(path_to_data,'mcrA.gpkg','mcrA.hmm')
        hmm2 = os.path.join(path_to_data,'mcrA_second_half.gpkg','mcrA.300-557.aln.fasta.hmm')

        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 2 --search_and_align_only --forward %s --search_hmm_files %s %s --output_directory %s --force 2>/dev/null' % (path_to_script,
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
            cmd = '%s graft --verbosity 2  --search_and_align_only --forward %s --search_hmm_files %s %s --aln_hmm_file %s --output_directory %s --force' % (path_to_script,
                                                                                               data,
                                                                                               hmm2,hmm,
                                                                                               hmm,
                                                                                               tmp)
            subprocess.check_output(cmd, shell=True)
            otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
            # otu table should not exist
            self.assertFalse(os.path.isfile(otuTableFile))

            expected = ['>example_partial_mcra_1_1_8',
                        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
            count = 0
            alnFile = os.path.join(tmp, 'mcrA_1.1', 'mcrA_1.1_hits.aln.fa')
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
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --search_and_align_only --forward %s --search_hmm_list_file %s --aln_hmm_file %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   hmms.name,
                                                                                                   hmm,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                # otu table should not exist
                self.assertFalse(os.path.isfile(otuTableFile))
                expected = ['>example_partial_mcra_1_1_8',
                            '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGVGFTQYATAAYTDDILDNNVYYNIDYINDKYKTDNKVKATLEVVKDIATESTIYGIETYEKFPTALEDHFGXSQRATVLAAAAGVXSALATANANAGLSGWYLSMYLHKEAWGRLGFFGYDLQDQCGATNVLSYQGDEGLPDELRGPNYPNYAM----------------------------------------------------------------------']
                count = 0
                alnFile = os.path.join(tmp, 'mcrA_1.1', 'mcrA_1.1_hits.aln.fa')
                
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
                cmd = '%s graft --verbosity 2 --min_orf_length 300 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
                cmd = '%s graft --verbosity 2  --restrict_read_length 102 --forward %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                   data,
                                                                                                   package,
                                                                                                   tmp)
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])

                otuTableFile = os.path.join(tmp, 'combined_count_table.txt')
                lines = ("\t".join(('#ID',sample_name,'ConsensusLineage')),
                         "\t".join(('1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanosarcinales; Methanosarcinaceae; Methanosarcina')),
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
            cmd = '%s graft --verbosity 5  --forward %s %s --graftm_package %s --output_directory %s --force' % (path_to_script,
                                                                                                  reads_1,
                                                                                                  reads_2,
                                                                                                  gpkg,
                                                                                                  tmp)
            extern.run(cmd)
            comb_otu_table=os.path.join(tmp, "combined_count_table.txt")
            
            expected=('\t'.join(('#ID', 'sample_16S_1.1', 'sample_16S_2.1', 'ConsensusLineage')),
                      '\t'.join(('1','1','1','Root; k__Bacteria; p__Cyanobacteria; c__Chloroplast')),
                      '\t'.join(('2','9','36','Root; k__Bacteria; p__Proteobacteria')),
                      '\t'.join(('3','0','2','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__mitochondria')),
                      '\t'.join(('4','18','106','Root; k__Bacteria')),
                      '\t'.join(('5','5','6','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales')),
                      '\t'.join(('6','1','18','Root')),
                      '\t'.join(('7','2','9','Root; k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria'))
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
                package = os.path.join(path_to_data,'61_otus.gpkg')
                with tempdir.TempDir() as tmp:
                    cmd = '%s graft --verbosity 2  --forward %s --graftm_package %s --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                     data,
                                                                                                                     package,
                                                                                                                     tmp)
                    subprocess.check_output(cmd, shell=True)
                    sample_name = os.path.basename(fastq_gz.name[:-6])
                    alnFile = os.path.join(tmp, sample_name, '%s_hits.aln.fa' % sample_name)
                    expected_aln = '''>NS500333:16:H16F3BGXX:1:11101:11211:1402
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAGCGCAACCCTCGCCTTCAGTTGCCATCTGAAGGAACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCGGT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    >NS500333:16:H16F3BGXX:1:11101:25587:3521
    -----------------------------------------------------------------------------------------------------------------------------------------------------------CGCTATTGGATGAGCCAAAGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACAATCCGTAGCTGGTCTGAGAGGACGACCAGCCACACCGGAACTGAGACACGGTCCGGACT--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    '''.split()
                    count = 0
                    for line in open(alnFile):
                        self.assertEqual(expected_aln[count], line.strip())
                        count += 1
                    self.assertEqual(count, len(open(alnFile).readlines()))
                    
    def test_bootrap_contigs(self):
        # this read is picked up by bootstrap but not regular graftm at the evalue
        testing_read = '''>196339_2
AAGATGAGCGCCAGGCTCTTCCTCCCCTTGTTGTTGCAGGGGATCATGAAGTCGATGTAT
TGATGCGGCGTGTCTGTATCTACAAGCGCCACCACGGGTA
>197849_1
ATCTACAAGCGCCACCACGGGTATCCCCATCTTGGCCGCCTCGGCGACGGCTTGGGAGTC
TAGTCTCGGGTCTACTACGAATAGCAAGTCTACCTCAAGG
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            original_hmm = os.path.join(path_to_data,'bootstrapper','DNGNGWU00001.hmm')
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --forward %s --search_hmm_files %s --aln_hmm %s --evalue 1e-14 --output_directory %s --force --search_and_align_only' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 original_hmm,
                                                                                                                 original_hmm,
                                                                                                                 tmp)
                
                subprocess.check_output(cmd, shell=True)
                sample_name = os.path.basename(fasta.name[:-3])
                self.assertEqual('', open(os.path.join(tmp, sample_name, '%s_hits.fa' % sample_name)).read())
                
                cmd = "%s --bootstrap_contigs %s" % (cmd, 
                                                     os.path.join(path_to_data,'bootstrapper','contigs.fna'))

                subprocess.check_output(cmd, shell=True)
                self.assertEqual(testing_read, open(os.path.join(tmp, sample_name, '%s_hits.fa' % sample_name)).read())

                bootstrap_hmm_path = os.path.join(tmp, 'bootstrap.hmm')
                self.assertTrue(os.path.isfile(bootstrap_hmm_path))
                self.assertEqual('HMMER3/f [3.1b2 | February 2015]\n',
                                 open(bootstrap_hmm_path).readlines()[0])
                
    def test_diamond_search_with_pplacer(self):
        testing_read = '''>Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040
ATGGCTACTGAAAAAACACAAAAGATGTTCCTCGAGGCGATGAAAAAGAAGTTCGCAGAGGACCCTACTTCAAACAAGACGACCTATAAGCGCGAGGGGTGGACTCAGTCCAAGGACAAGCGCGAGTTCCAGGAATGGGGCGCAAAAATCGCCAAGGACCGTGGAATACCGGCGTACAACGTCAACGTCCACCTCGGCGGTATGACCCTCGGCCAGCGGCAACTCATGCCGTACAATGTCTCTGGGACCGACGTGATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGGCGTACGGTTATCGTAGGTCTTGACACCGCTCACGAGACGCTGACCAGGAGACTCGGCAAGGAGGTTACCCCCGAGACCATCAACGGCTATCTCGAGGCATTGAACCACACGATGCCCGGTGCGGCCATTGTCCAAGAACACATGGTGGAAACCCACCCTGCGCTCGTTGAAGACTGCTTCGTAAAAGTCTTCACCGGCGACGATGACCTCGCC'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2 --search_method diamond --forward %s --output_directory %s --force --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                subprocess.check_output(cmd, shell=True)
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                
    def test_diamond_search_with_pplacer_protein_input(self):
        testing_read = '''>2518787893 METHANOFLORENS STORDALENMIRENSIS MCRA
MATEKTQKMFLEAMKKKFAEDPTSNKTTYKREGWTQSKDKREFQEWGAKIAKDRGIPAY
NVNVHLGMTLGQRQLMPYNVSGTDVMCEGDDLHYVNNPAMQQMWDEIRRTVIVGLDTAH
ETLTRRLGKEVTPETINGYLEALNHTMPGAAIVQEHMVETHPALVEDCFVKVFTGDDDLA'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2 --search_method diamond --forward %s --output_directory %s --force --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                subprocess.check_output(cmd, shell=True)
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                
    def test_diamond_placement_method(self):
        testing_read = '''>Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040
ATGGCTACTGAAAAAACACAAAAGATGTTCCTCGAGGCGATGAAAAAGAAGTTCGCAGAGGACCCTACTTCAAACAAGACGACCTATAAGCGCGAGGGGTGGACTCAGTCCAAGGACAAGCGCGAGTTCCAGGAATGGGGCGCAAAAATCGCCAAGGACCGTGGAATACCGGCGTACAACGTCAACGTCCACCTCGGCGGTATGACCCTCGGCCAGCGGCAACTCATGCCGTACAATGTCTCTGGGACCGACGTGATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGGCGTACGGTTATCGTAGGTCTTGACACCGCTCACGAGACGCTGACCAGGAGACTCGGCAAGGAGGTTACCCCCGAGACCATCAACGGCTATCTCGAGGCATTGAACCACACGATGCCCGGTGCGGCCATTGTCCAAGAACACATGGTGGAAACCCACCCTGCGCTCGTTGAAGACTGCTTCGTAAAAGTCTTCACCGGCGACGATGACCTCGCC'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            sample_name = os.path.basename(fasta.name[:-3])
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --forward %s --output_directory %s --force --assignment_method diamond --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                subprocess.check_output(cmd, shell=True)
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales; Methanoflorentaceae; Methanoflorens']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                self.assertTrue(os.path.exists(os.path.join(tmp, sample_name, '%s.hmmout.csv' % sample_name)), "should keep the hmmer search file")
                
    def test_diamond_placement_method_protein_input(self):
        testing_read = '''>2518787893 METHANOFLORENS STORDALENMIRENSIS MCRA
MATEKTQKMFLEAMKKKFAEDPTSNKTTYKREGWTQSKDKREFQEWGAKIAKDRGIPAY
NVNVHLGMTLGQRQLMPYNVSGTDVMCEGDDLHYVNNPAMQQMWDEIRRTVIVGLDTAH
ETLTRRLGKEVTPETINGYLEALNHTMPGAAIVQEHMVETHPALVEDCFVKVFTGDDDLA'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2  --forward %s --output_directory %s --force --assignment_method diamond --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                subprocess.check_output(cmd, shell=True)
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales; Methanoflorentaceae; Methanoflorens']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                
    def test_diamond_search_and_place_method(self):
        testing_read = '''>Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040
ATGGCTACTGAAAAAACACAAAAGATGTTCCTCGAGGCGATGAAAAAGAAGTTCGCAGAGGACCCTACTTCAAACAAGACGACCTATAAGCGCGAGGGGTGGACTCAGTCCAAGGACAAGCGCGAGTTCCAGGAATGGGGCGCAAAAATCGCCAAGGACCGTGGAATACCGGCGTACAACGTCAACGTCCACCTCGGCGGTATGACCCTCGGCCAGCGGCAACTCATGCCGTACAATGTCTCTGGGACCGACGTGATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGGCGTACGGTTATCGTAGGTCTTGACACCGCTCACGAGACGCTGACCAGGAGACTCGGCAAGGAGGTTACCCCCGAGACCATCAACGGCTATCTCGAGGCATTGAACCACACGATGCCCGGTGCGGCCATTGTCCAAGAACACATGGTGGAAACCCACCCTGCGCTCGTTGAAGACTGCTTCGTAAAAGTCTTCACCGGCGACGATGACCTCGCC'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing_read)
            fasta.flush()
            sample_name = os.path.basename(fasta.name[:-3])
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --search_method diamond --forward %s --output_directory %s --force --assignment_method diamond --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                extern.run(cmd)
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales; Methanoflorentaceae; Methanoflorens']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                
                
                expected = [['#ID',os.path.basename(fasta.name)[:-3],'ConsensusLineage'],
                            ['1','1','Root; mcrA; Euryarchaeota_mcrA; Methanomicrobia; Methanocellales; Methanoflorentaceae; Methanoflorens']]
                expected = ['\t'.join(l) + '\n' for l in expected]
                self.assertEqual(expected, open(os.path.join(tmp,'combined_count_table.txt')).readlines())
                
                self.assertTrue(os.path.exists(os.path.join(tmp, sample_name, '%s_diamond_search.daa' % sample_name)), "should keep the diamond search file")
                self.assertTrue(os.path.exists(os.path.join(tmp, sample_name, '%s_diamond_assignment.daa' % sample_name)), "should keep the diamond assign file")

    def test_hit_where_sequence_evalue_is_good_but_individuals_bad(self):
        # the first one is a real hit, the second has several hits better but none better than 1e-5
        testing = '''>2509711280 Pleur7313DRAFT_05268 ribosomal protein S19, bacterial/organelle [Pleurocapsa sp. PCC 7319]
MSRSLKKGPFIADSLLKKIEKLNANNKKEVIKTWSRASTILPQMVGHTIAVHNGRQHIPVFISDQMVGHKLGEFAPTRTFRGHAKSDKKGRR
>2509712449 Pleur7313DRAFT_06440 FOG: WD40 repeat [Pleurocapsa sp. PCC 7319]
MSKSVVINLGSGDLYKGFPRVTAQLWAAGYPRPEQFIGSLPAAPALAESYRTWQSIYKALCSRLVLFSRGPDDDDDELEIDEGGITQVSEVSFDALGQQLHESLNAWLKSVEFLNIERQLRSQLDPSEEIRVIIETNDEQLRRLPWHCWDFVEDYYYSEAALSRPEYKRIKRSPSKINRKKVRILAVLGNSQGIDLEVERRFLKGLPDAETVFLVNPSRHEFNEQLWNSQGWDILFFAGHSQSEGKTGRIYLNENKTNNSLTVEQLKEALKEAIKNDLYLAIFNSCDGLGLANALGKLNLPQIIVMREPVPNRVAQEFFKHFLSAFASQRSPLYSAVRQARRQLQGLENEFPSASWLPVISQNPAVEPPTWLQLGGIPPCPYRGLSAFREEDAHLFCGREEFTANLLKEVKSKPLVAVVGASGSGKSSVVFAGLIPHLQQDTDLCWQIVSFRPGNNPVDSLAAALVPLWQQRENNQEDNLASRDNSIVELELFLALRRDELALSKLIKSLVESPKTRLLLIADQFEELYTLSNQAERQFFLTLLLNAVKFAPAFTLVLTLRADFYGYILGDRSWSDALQGAIQNLGSMSRSELQLVIEKPAALMQVKLEQGLTDKLINNVWEQPGRLPLLEFALTQLWSKQQHGLLTHQAYSEIGGVSSALANHAEAVYAQLSETDKQLVQQVFLQLVRLGEETEATRRLATRDEIKAENWDLVTRLASARLVVTNHRYSTGKETVEIVHESLIKSWGRLQQWLLLDGDFRRWQEQLRTAIGQWQWQGSSSNQDTLLRGKPLTDASNWLQKRFQELTDSERSFIQASLEQRNNHIKAEKRRQQWTILGLTLGSILAVSLMGVARWQWQKVRIGELYALTKSSEVLFADNDRLNALVKAIAAQEKLQKLGRVDTKLQNQVESILRRSIYGANESNRLSGHHGAVWGVAFSPDGQTIASTSWDNTIKLWSRNGKELKTITGHVEEVWGVAFSPDGQTIATASGDNTVKLWQHDGTLLKTLTGHGDIVYSVAFSPDGQAIATASGDNTVKLWQRDGTLVKTLTGHEASVWGVAFNPNGQTIASAGWDKKVKLWSRDGALLKTLESHQAPVNALAFSPDGQTIATASNDKTVKLWSRDGTLLKTLEGHDDDVWGISFSPDGQTFASVSGDKTVKLWSRDGTLLQTLEGHDDEVWGVSFSPDGQTLASAGDDKTIRLWQRDNKLLTTLSSHSAGVNGVAFSPDEQLIASAGWDKTVKLWNSNGSLLKTLIGHSAAINALAFEPGGNIIATASADNTVKFWQRDGTLLKTLIGHRAGVNAVVFSPDGEIVASASADTTVKLWSREGTLLKTLTGHRAGVNAVVFSPDGEMVASASADNTVKFWSRDGTELKTFKAHGDRLYALAFSFDGQMLATASADNTVKLWNRDGTELKTLNGHNGTVWGVAFSPDGQTIATASGDNTVNLWKLDGTLLTTLNAHNSAVNAIAFNSDGNRLTSASEDRTVIVWDLNRVRSLDEILAYSCDWVQHYLQNNANVNPEALLLCGKIENLNPQ
'''
        
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing)
            fasta.flush()
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 2 --search_and_align_only --forward %s --output_directory %s --force --search_hmm_files %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'s19.hmm'))
                extern.run(cmd)
                expected = '''>2509711280 Pleur7313DRAFT_05268 ribosomal protein S19, bacterial/organelle [Pleurocapsa sp. PCC 7319]
MSRSLKKGPFIADSLLKKIEKLNANNKKEVIKTWSRASTILPQMVGHTIAVHNGRQHIPV
FISDQMVGHKLGEFAPTRTFRGHAKSDKKGRR
'''
                self.assertEqual(expected, open(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "%s_hits.fa" % os.path.basename(fasta.name)[:-3])).read())
                
    def test_forward_and_reverse_slash_type_fastq(self):
        fwd = '''@FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/1
ACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGGATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGACGCCGCGTGTGGGATGAAGGCCCTCGGGT
+
^^_cccacgeecafgfghhhhheYea_c^efaff`eggfhhhhf]`_d]]X^aa\]^Y^a`[^bbaZaaaa]]a]_[]HTX`[[[_[]`_][`^^`]__a
'''
        rev = '''@FCC0WM1ACXX:2:2208:12709:74426#GTCCAGAA/2
CGGGGTATCTAATCCCGTTCGCTCCCCTAGCTTTCGTGCCTCAGCGTCAGAAAAGACCCAGTGAGCCGCTTTCGCCCCCGGTGTTCCTTAGGATATCAAC
+
\a_ccO_ceceeehdgaffZ^degfggfdefggib^ef^cecRafeefgdZecf_dd`gbcZ___b]aUZaa`aa__aaX__TT[]_bY]Y`]RG]`b_b
'''
        with tempfile.NamedTemporaryFile(suffix='.fq') as fwd_f:
            fwd_f.write(fwd)
            fwd_f.flush()
            with tempfile.NamedTemporaryFile(suffix='.fq') as rev_f:
                rev_f.write(rev)
                rev_f.flush()

                with tempdir.TempDir() as tmp:
                    cmd = '%s graft --verbosity 5 --forward %s --reverse %s --output_directory %s --force --graftm_package %s' % (path_to_script,
                                                                                                                     fwd_f.name,
                                                                                                                     rev_f.name,
                                                                                                                     tmp,
                                                                                                                     os.path.join(path_to_data,'61_otus.gpkg'))
                    extern.run(cmd)
                    raise Exception("At least it runs, but nt sure what is expected for this test. Before it didn't even run")

    @unittest.skip("known failure")
    def test_input_sequence_near_duplication(self):
        testing = '''>Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040
ATGGCTACTGAAAAAACACAAAAGATGTTCCTCGAGGCGATGAAAAAGAAGTTCGCAGAGGACCCTACTTCAAACAAGACGACCTATAAGCGCGAGGGGTGGACTCAGTCCAAGGACAAGCGCGAGTTCCAGGAATGGGGCGCAAAAATCGCCAAGGACCGTGGAATACCGGCGTACAACGTCAACGTCCACCTCGGCGGTATGACCCT
CGGCCAGCGGCAACTCATGCCGTACAATGTCTCTGGGACCGACGTGATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGGCGTACGGTTATCGTAGGTCTTGACACCGCTCACGAGACGCTGACCAGGAGACTCGGCAAGGAGGTTACCCCCGAGACCATCAACGGCTATCTCGA
GGCATTGAACCACACGATGCCCGGTGCGGCCATTGTCCAAGAACACATGGTGGAAACCCACCCTGCGCTCGTTGAAGACTGCTTCGTAAAAGTCTTCACCGGCGACGATGACCTCGCC
>another_Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040
ATGGCTACTGAAAAAACACAAAAGATGTTCCTCGAGGCGATGAAAAAGAAGTTCGCAGAGGACCCTACTTCAAACAAGACGACCTATAAGCGCGAGGGGTGGACTCAGTCCAAGGACAAGCGCGAGTTCCAGGAATGGGGCGCAAAAATCGCCAAGGACCGTGGAATACCGGCGTACAACGTCAACGTCCACCTCGGCGGTATGACCCT
CGGCCAGCGGCAACTCATGCCGTACAATGTCTCTGGGACCGACGTGATGTGTGAAGGCGATGACCTCCACTACGTCAACAACCCCGCAATGCAACAGATGTGGGATGAGATCAGGCGTACGGTTATCGTAGGTCTTGACACCGCTCACGAGACGCTGACCAGGAGACTCGGCAAGGAGGTTACCCCCGAGACCATCAACGGCTATCTCGA
GGCATTGAACCACACGATGCCCGGTGCGGCCATTGTCCAAGAACACATGGTGGAAACCCACCCTGCGCTCGTTGAAGACTGCTTCGTAAAAGTCTTCACCGGCGACGATGACCTCGCC
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing)
            fasta.flush()
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --output_directory %s --force --graftm_package %s' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'mcrA.gpkg'))
                extern.run(cmd)
                expected = '''
>Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040_1_1_7_0 Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040_1_1_7
MATEKTQKMFLEAMKKKFAEDPTSNKTTYKR-EGWTQSKDKREFQEWGAKIAKDRGIPAY
NVNVHLGMTLGQRQLMPYNVSGTDVMCEGDDLHYVNNPAMQQMWDEIRRTVIVGLDTAHE
TLTRRLGKEVTPETINGYLEALNHTMPGAAIVQEHMVETHPALVEDCFVKVFTGDDDLA-
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-----------------
>another_Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040_1_1_7_0 Methanoflorens_stordalmirensis_v4.3_scaffold3_chopped_215504-216040_1_1_7
MATEKTQKMFLEAMKKKFAEDPTSNKTTYKR-EGWTQSKDKREFQEWGAKIAKDRGIPAY
NVNVHLGMTLGQRQLMPYNVSGTDVMCEGDDLHYVNNPAMQQMWDEIRRTVIVGLDTAHE
TLTRRLGKEVTPETINGYLEALNHTMPGAAIVQEHMVETHPALVEDCFVKVFTGDDDLA-
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-----------------
'''
                self.assertEqual(expected, open(os.path.join(tmp, "combined_alignment.aln.fa")).read())
                
                jplace = json.load(open(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "placements.jplace")))
                self.assertEqual(2, len(jplace['placements']))
        
    def test_filter_minimum(self):
        testing = '''>NS500333:6:H1124BGXX:1:23310:10768:12778 1:N:0:GATCAG
CGGGAGGAACACCAGTGGCGAAGGCGGCTTCCTGGCCTGTTCTTGACGCTGAGGCGCGAA
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(testing)
            fasta.flush()
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --output_directory %s --force --graftm_package %s --filter_minimum 0' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'61_otus.gpkg'))
                extern.run(cmd)
                expected = '>NS500333:6:H1124BGXX:1:23310:10768:12778\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGAGGAACACCAGTGGCGAAGGCGGCTTCCTGGCCTGCTTGACGCTGAGGCGCGAA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'
                self.assertEqual(expected, open(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "%s_hits.aln.fa" % os.path.basename(fasta.name)[:-3])).read())
                
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --output_directory %s --force --graftm_package %s --filter_minimum 150' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'61_otus.gpkg'))
                extern.run(cmd)
                self.assertFalse(os.path.exists(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "%s_hits.aln.fa" % os.path.basename(fasta.name)[:-3])))
                
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --output_directory %s --force --graftm_package %s --filter_minimum 56' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'61_otus.gpkg'))
                extern.run(cmd)
                self.assertTrue(os.path.exists(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "%s_hits.aln.fa" % os.path.basename(fasta.name)[:-3])))
                
            with tempdir.TempDir() as tmp:
                cmd = '%s graft --verbosity 5 --forward %s --output_directory %s --force --graftm_package %s --filter_minimum 57' % (path_to_script,
                                                                                                                 fasta.name,
                                                                                                                 tmp,
                                                                                                                 os.path.join(path_to_data,'61_otus.gpkg'))
                extern.run(cmd)
                self.assertFalse(os.path.exists(os.path.join(tmp, os.path.basename(fasta.name)[:-3], "%s_hits.aln.fa" % os.path.basename(fasta.name)[:-3])))

        
    def test_max_samples_for_krona(self):
        reads_1=os.path.join(path_to_samples, "sample_16S_1.1.fa")
        reads_2=os.path.join(path_to_samples, "sample_16S_2.1.fa")
        gpkg=os.path.join(path_to_data, "61_otus.gpkg")
        
        with tempdir.TempDir() as tmp:
            cmd = '%s graft --verbosity 5  --forward %s %s --graftm_package %s --output_directory %s --force --max_samples_for_krona 2' % (path_to_script,
                                                                                                  reads_1,
                                                                                                  reads_2,
                                                                                                  gpkg,
                                                                                                  tmp)
            extern.run(cmd)
            self.assertTrue(os.path.exists(os.path.join(tmp, "combined_count_table.txt")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "krona.html")))
            
            cmd = '%s graft --verbosity 5  --forward %s %s --graftm_package %s --output_directory %s --force --max_samples_for_krona 1' % (path_to_script,
                                                                                                  reads_1,
                                                                                                  reads_2,
                                                                                                  gpkg,
                                                                                                  tmp)
            extern.run(cmd)
            self.assertTrue(os.path.exists(os.path.join(tmp, "combined_count_table.txt")))
            self.assertFalse(os.path.exists(os.path.join(tmp, "krona.html")))

        
                    
if __name__ == "__main__":
    unittest.main()
