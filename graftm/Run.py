#!/usr/bin/env python

__author__ = "Joel Boyd, Ben Woodcroft"
__copyright__ = "Copyright 2014"
__credits__ = ["Joel Boyd", "Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd, Ben Woodcroft"
__email__ = "joel.boyd near uq.net.au, b.woodcroft near uq.edu.au"
__status__ = "Development"
__version__ = "0.4.0"

from graftm.Messenger import Messenger
from graftm.GraftMFiles import GraftMFiles
from graftm.Extract_Sequences import Extract
from graftm.Hmmer import Hmmer
from graftm.HouseKeeping import HouseKeeping
from graftm.Stats_And_Summary import Stats_And_Summary
from graftm.Alignment_Manager import Alignment_Manager
from graftm.DatManip import DatManip
from graftm.Pplacer import Pplacer
from graftm.krona_from_community_profiles import KronaBuilder
from graftm.assembler import TaxoGroup

import timeit
import tempfile
import os
import IPython



class Run:
    
    def __init__(self, args):
        self.args = args
        
        self.setattributes(self.args)
        
    def setattributes(self, args):
        self.KB = KronaBuilder()
        self.AM = Alignment_Manager()
        self.DM = DatManip()
        self.HK = HouseKeeping()
        self.SAS = Stats_And_Summary()
        self.TG = TaxoGroup()
        self.ES = Extract()
        
        if args.subparser_name == 'graft':
            self.HK.set_attributes(self.args)
            self.H = Hmmer(self.args.hmm_file)
            self.sequence_pair_list, self.input_file_format = self.HK.parameter_checks(args)
            self.P = Pplacer(self.args.reference_package)
    
    def filter(self):
        
        print '''
                           FILTER
           
                   Joel Boyd, Ben Woodcroft
                                        
'''     
        hit_hash = {}
        HMM_list = os.listdir(self.args.filter_hmms)
        for hmm in HMM_list:
            hit_hash[hmm] = {}
            self.H = Hmmer(hmm)
        

        

    
    
    def protein_pipeline(self, base, summary_dict, sequence_file_list):
            
            run_stats = summary_dict[base]
                    
            # Search for reads
            
            Messenger().message('Searching %s using %s' % (os.path.basename(sequence_file_list[0]), os.path.basename(self.args.hmm_file)))
            start = timeit.default_timer()
            
            hmm_search_output = self.H.hmmsearch(self.GMF.forward_read_hmmsearch_output_path(base), 
                                                 self.GMF.reverse_read_hmmsearch_output_path(base),
                                                 sequence_file_list, 
                                                 self.input_file_format, 
                                                 self.args.input_sequence_type,
                                                 self.args.threads,
                                                 self.args.eval)
                                            
            stop = timeit.default_timer()
            run_stats['search_t'] = str(int(round((stop - start), 0)) )
            
            # Interpret the results of the hmmsearch output
            Messenger().message('Reading results')
            evals, n_total_reads, rev_true = self.DM.csv_to_titles(hmm_search_output, 
                                                                   self.args.type,
                                                                   self.GMF.readnames_output_path(base),
                                                                   base,
                                                                   self.args.input_sequence_type)
            run_stats['rev_true'] = rev_true
            run_stats['n_total_reads'] = n_total_reads
            run_stats['evals'] = evals    
            # Extract the hists from the original file
            Messenger().message('Extracting reads')
            start = timeit.default_timer() 
            
            self.DM.extract_from_raw_reads(sequence_file_list[0], 
                                           self.args.hmm_file, 
                                           self.GMF.readnames_output_path(base), 
                                           self.input_file_format, 
                                           self.args.output_directory,
                                           self.GMF.fa_output_path(base),
                                           self.GMF.orf_output_path(base),
                                           self.GMF.orf_hmmsearch_output_path(base),
                                           self.GMF.orf_titles_output_path(base),
                                           self.GMF.orf_fasta_output_path(base),
                                           False,
                                           self.args.input_sequence_type)  # And extract from the original file
            
            stop = timeit.default_timer()
            
            
            Messenger().message('Aligning reads to reference package database')
            start = timeit.default_timer()
            
            if self.args.search_only:
                Messenger().message('Stopping before alignment\n')
                
                exit(0)
            
            if self.args.input_sequence_type == 'nucleotide':
                self.H.hmmalign(base,
                                self.GMF.orf_fasta_output_path(base),
                                run_stats,
                                self.args)
            
            elif self.args.input_sequence_type == 'protein':
                self.H.hmmalign(base,
                                self.GMF.fa_output_path(base),
                                run_stats,
                                self.args)

            # Fix up the alignment file by removing the insertions
            if run_stats['rev_true']:
                self.AM.alignment_correcter([self.GMF.conv_output_for_path(base), self.GMF.conv_output_rev_path(base)], 
                                            self.GMF.aligned_fasta_output_path(base))
                
            else:
                self.AM.alignment_correcter([self.GMF.conv_output_for_path(base)], 
                                            self.GMF.aligned_fasta_output_path(base))

            stop = timeit.default_timer()
            run_stats['aln_t'] = str(int(round((stop - start), 0)) )
                
            run_stats['extract_t'] = str(int(round((stop - start), 0)) )
            run_stats['n_contamin_euks'] = 'N/A'
            run_stats['n_uniq_euks'] = 'N/A'
            run_stats['euk_check_t'] = 'N/A'
            
            summary_dict[base] = run_stats
            return self.GMF.aligned_fasta_output_path(base), summary_dict

    def dna_pipeline(self, base, summary_dict, sequence_file_list):
            run_stats = summary_dict[base]           
            
            # Search for reads
            Messenger().message('Searching %s using %s' % (os.path.basename(sequence_file_list[0]), os.path.basename(self.args.hmm_file)))
            start = timeit.default_timer()
            
            hmm_search_output = self.H.nhmmer(self.GMF.forward_read_hmmsearch_output_path(base), 
                                              self.GMF.reverse_read_hmmsearch_output_path(base), 
                                              sequence_file_list, 
                                              self.input_file_format,
                                              self.args.threads,
                                              self.args.eval)
                                        
            stop = timeit.default_timer()
            run_stats['search_t'] = str(int(round((stop - start), 0)) )
            
            
            # Interpret the results of the nhmmer output
            Messenger().message('Reading results')

            evals, n_total_reads, rev_true = self.DM.csv_to_titles(hmm_search_output, 
                                                                   self.args.type, 
                                                                   self.GMF.readnames_output_path(base),
                                                                   base,
                                                                   self.args.input_sequence_type)
            run_stats['rev_true'] = rev_true
            run_stats['n_total_reads'] = n_total_reads
            run_stats['evals'] = evals
            # Extract the hists from the original file
            Messenger().message('Extracting reads')
            start = timeit.default_timer()
            
            self.DM.extract_from_raw_reads(sequence_file_list[0], 
                                           self.args.hmm_file, 
                                           self.GMF.readnames_output_path(base), 
                                           self.input_file_format, 
                                           self.args.output_directory, 
                                           self.GMF.fa_output_path(base),
                                           self.GMF.orf_output_path(base),
                                           self.GMF.orf_hmmsearch_output_path(base),
                                           self.GMF.orf_titles_output_path(base),
                                           self.GMF.orf_fasta_output_path(base),
                                           True,
                                           self.args.input_sequence_type)
            
            stop = timeit.default_timer()
            run_stats['extract_t'] = str(int(round((stop - start), 0)) )
            
            
            # Check for Eukarytoic contamination
            Messenger().message("Checking for Eukaryotic contamination")
            start = timeit.default_timer()

            n_contamin_euks, n_uniq_euks = self.H.check_euk_contamination(self.GMF.euk_free_path(base),
                                                                          self.GMF.euk_contam_path(base), 
                                                                          self.GMF.fa_output_path(base), 
                                                                          run_stats['evals'], 
                                                                          self.SAS.check_read_length(self.GMF.fa_output_path(base), "D"), 
                                                                          self.input_file_format, 
                                                                          self.args.threads, 
                                                                          self.args.eval,
                                                                          self.args.check_total_euks,
                                                                          sequence_file_list[0])
            stop = timeit.default_timer()
            
            Messenger().message('Aligning reads to reference package database')
            start = timeit.default_timer()
            
            if self.args.search_only:
                Messenger().message('Stopping before alignment\n')
                
                exit(0)
            
            
            self.H.hmmalign(base,
                            self.GMF.euk_free_path(base),
                            run_stats,
                            self.args)
            

            # Fix up the alignment file by removing the insertions
            if run_stats['rev_true']:
                self.AM.alignment_correcter([self.GMF.conv_output_for_path(base), self.GMF.conv_output_rev_path(base)], 
                                            self.GMF.aligned_fasta_output_path(base))
                
            else:
                self.AM.alignment_correcter([self.GMF.conv_output_for_path(base)], 
                                            self.GMF.aligned_fasta_output_path(base))

            stop = timeit.default_timer()
            run_stats['aln_t'] = str(int(round((stop - start), 0)) )
            
            
            
            run_stats['euk_check_t'] = str( int(round((stop - start), 0)) )
            run_stats['n_contamin_euks'] = n_contamin_euks
            run_stats['n_uniq_euks'] = n_uniq_euks        
    
            summary_dict[base] = run_stats
            return self.GMF.aligned_fasta_output_path(base), summary_dict
                                
    def placement(self, aln_files, base_list, GM_temp, summary_dict):      

            # Check that there is a reference package to place in. If not, exit.     
            if self.args.seach_and_align_only:
                Messenger().message('Stopping before placement\n')
                for base in base_list:
                    self.GMF = GraftMFiles(base, self.args.output_directory)
                    self.HK.delete([self.GMF.for_aln_path(base), 
                                    self.GMF.rev_aln_path(base), 
                                    self.GMF.sto_for_output_path(base), 
                                    self.GMF.sto_rev_output_path(base), 
                                    self.GMF.conv_output_rev_path(base), 
                                    self.GMF.conv_output_for_path(base), 
                                    self.GMF.euk_free_path(base), 
                                    self.GMF.euk_contam_path(base), 
                                    self.GMF.readnames_output_path(base), 
                                    self.GMF.forward_read_hmmsearch_output_path(base), 
                                    self.GMF.sto_output_path(base), 
                                    self.GMF.reverse_read_hmmsearch_output_path(base),
                                    self.GMF.orf_titles_output_path(base), 
                                    self.GMF.orf_hmmsearch_output_path(base),
                                    self.GMF.orf_output_path(base)])
                exit(0)
            
            # Combine sequences (4 Xtra zpeed)
            Messenger().message("Combining alignments")
            alias_hash = self.AM.name_changer(base_list,
                                              self.GMF.comb_aln_fa(GM_temp),
                                              self.args)
            
            # Place in tree
            Messenger().message('Placing reads into reference package tree')
            start = timeit.default_timer()

            placements = self.P.pplacer(self.GMF.comb_aln_fa(GM_temp),
                                        self.args.threads,
                                        base_list,
                                        GM_temp)
            
           
            stop = timeit.default_timer()
            summary_dict['place_t'] = str( int(round((stop - start), 0)) )
            
            jplace_path_list = self.DM.jplace_split(placements, alias_hash, self.args)
            
            # Create Guppy File
            Messenger().message('Creating Guppy file')
            start = timeit.default_timer()
            n_placements = self.P.guppy_class(self.GMF.guppy_file_output_path(),
                                              jplace_path_list,
                                              GM_temp,
                                              self.args)
            
            for count in n_placements.keys():
                summary_dict[count]['reads_found'] = n_placements[count]

            # Generate Summary Table
            base_path_list = [self.args.output_directory + '/' +base+'/'+base for base in base_list]
            for idx, base in enumerate(base_list):
                self.GMF = GraftMFiles(base, self.args.output_directory)
                Messenger().message('Building summary table for %s' % base)
                self.SAS.otu_builder(base_path_list[idx] + self.GMF.guppy_file_output_path(), 
                                     base_path_list[idx] + self.GMF.summary_table_output_path(),
                                     self.args.placements_cutoff,
                                     base)
            
                
                # Generate coverage table
                Messenger().message('Building coverage table for %s' % base)
                self.SAS.coverage_of_hmm(self.args.hmm_file, 
                                         base_path_list[idx] + self.GMF.summary_table_output_path(), 
                                         base_path_list[idx] + self.GMF.coverage_table_path(), 
                                        self.SAS.check_read_length(self.GMF.fa_output_path(base), self.args.type))
                
                Messenger().message('Building krona for %s' % base)
                self.KB.otuTablePathListToKrona(base_path_list[idx] + self.GMF.summary_table_output_path(), 
                                                base_path_list[idx] + self.GMF.krona_output_path())
            
            stop = timeit.default_timer()
            summary_dict['summary_t'] = str(int(round((stop - start), 0)) )
            
            # Compile basic run statistics if they are wanted
            summary_dict['stop_all'] = timeit.default_timer()
            summary_dict['all_t'] = str(int(round((summary_dict['stop_all'] - summary_dict['start_all']), 0)) )
            
            
            self.SAS.build_basic_statistics(summary_dict, base_list, self.GMF.basic_stats_path(), self.args.type)
            
            
            # Delete unncessary files
            Messenger().message('Cleaning up')
            for base in base_list:
                self.GMF = GraftMFiles(base, self.args.output_directory)
                self.HK.delete([self.GMF.for_aln_path(base), 
                                self.GMF.rev_aln_path(base), 
                                self.GMF.sto_for_output_path(base), 
                                self.GMF.sto_rev_output_path(base), 
                                self.GMF.conv_output_rev_path(base), 
                                self.GMF.conv_output_for_path(base), 
                                self.GMF.euk_free_path(base), 
                                self.GMF.euk_contam_path(base), 
                                self.GMF.readnames_output_path(base), 
                                self.GMF.forward_read_hmmsearch_output_path(base), 
                                self.GMF.sto_output_path(base), 
                                self.GMF.reverse_read_hmmsearch_output_path(base),
                                self.GMF.orf_titles_output_path(base), 
                                self.GMF.orf_hmmsearch_output_path(base),
                                self.GMF.orf_output_path(base)])
            
            Messenger().message('Done, thanks for using graftM!\n')

    def graft(self):
        

        # Prepping - checking all necessary parameters are present, 
        # checking formats, making working directories
            
        if self.args.subparser_name == 'graft':
            print '''
                            GRAFT
           
                   Joel Boyd, Ben Woodcroft
                                                     __/__
                                              ______|
      _- - _                         ________|      |_____/
       - -            -             |        |____/_
       - _     --->  -   --->   ____|
      - _-  -         -             |      ______
         - _                        |_____|
       -                                  |______
                '''
            
            timer_hash = {}
            
            summary_table = {'euks_checked': self.args.check_total_euks,
                             'base_list': [],
                             'seqs_list': [],
                             'start_all': timeit.default_timer()}
        
            GM_temp = tempfile.mkdtemp(prefix='GM_temp_')
            
            # Check that a specific output directory has not been specified as a flag
            if not hasattr(self.args, 'output_directory'):
                self.args.output_directory = "GraftM_proc"

            self.HK.make_working_directory(self.args.output_directory, self.args.force)
                
            
            for pair in self.sequence_pair_list:  
                
                # If an input sequence type has not been specified, attempt to autocheck.
                if not hasattr(self.args, 'input_sequence_type'):
                    setattr(self.args, 'input_sequence_type', self.HK.guess_sequence_type(pair[0], self.input_file_format))
                    
                base = os.path.basename(pair[0]).split('.')[0]
                summary_table[base] = {}
                
                # Set the working path for all of the output files provided by 
                # GraftM.
                
                self.GMF = GraftMFiles(base, self.args.output_directory)
                self.HK.make_working_directory(self.args.output_directory +'/'+ base, self.args.force)
                
                Messenger().header("Working on %s" % base)
                
                # Set pipeline
                setattr(self.args, 'type', self.HK.setpipe(self.args.hmm_file))
                
                # Find input type
                
                # Protein pipeline
                if self.args.type == 'P':
                    info = self.protein_pipeline(    
                                            base, 
                                            summary_table,
                                            pair
                                                )
                
                # DNA pipeline
                elif self.args.type == 'D':
                    info = self.dna_pipeline(
                                        base, 
                                        summary_table,
                                        pair
                                            )
        
                summary_table['seqs_list'].append(info[0])
                summary_table['base_list'].append(base)
        
            Messenger().header("Grafting")
            self.placement(' '.join(summary_table['seqs_list']),
                           summary_table['base_list'],
                           GM_temp,
                           info[1])
        
        elif self.args.subparser_name == 'filter':
            contents = self.HK.read_contents(self.args.graftm_package)
            self.filter_pipe()
    
    def manage(self):
        print '''
                            MANAGE
           
                   Joel Boyd, Ben Woodcroft
                      
'''     
        
        if self.args.seq:
            self.ES.extract(self.args)
        
        
    def assemble(self):
        print '''
                           ASSEMBLE
           
                   Joel Boyd, Ben Woodcroft
                                        
                                   
          _- - _               ___            __/
           -                  /___\____      /\/   
           - _     --->   ___/       \_\     \/      
          - _-           /_/            \    /
             - _        /                \__/
                       /
'''
        self.TG.main(self.args)
        
        
        
        
    
    def main(self):

        if self.args.subparser_name == 'graft':
            self.graft()
            
        elif self.args.subparser_name == 'assemble':
            self.assemble()
            
        elif self.args.subparser_name == 'filter':
            self.filter()
        
        elif self.args.subparser_name == 'manage':
            self.manage()
            
