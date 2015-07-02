#!/usr/bin/env python

import timeit
import os
import logging

from graftm.sequence_search_results import SequenceSearchResult, HMMSearchResult
from graftm.graftm_output_paths import GraftMFiles
from graftm.extract_sequences import Extract
from graftm.hmmer import Hmmer
from graftm.housekeeping import HouseKeeping
from graftm.summarise import Stats_And_Summary
from graftm.pplacer import Pplacer
from graftm.assembler import TaxoGroup
from graftm.create import Create
from graftm.unpack_sequences import UnpackRawReads
from graftm.graftm_package import GraftMPackage
from biom.util import biom_open

PIPELINE_PROTEIN    = "P"
PIPELINE_NUCLEOTIDE = "D"

class Run:
    ### Functions that make up pipelines in GraftM

    def __init__(self, args):
        self.args = args
        self.setattributes(self.args)

    def setattributes(self, args):
        self.hk = HouseKeeping()
        self.s = Stats_And_Summary()
        self.tg = TaxoGroup()
        self.e = Extract()
        if args.subparser_name == 'graft':
            self.hk.set_attributes(self.args)
            self.hk.set_euk_hmm(self.args)
            if args.euk_check:self.args.search_hmm_files.append(self.args.euk_hmm_file)
            self.h = Hmmer(self.args.search_hmm_files, self.args.aln_hmm_file)
            self.sequence_pair_list = self.hk.parameter_checks(args)
            if hasattr(args, 'reference_package'):
                self.p = Pplacer(self.args.reference_package)
    
    def _get_sequence_directions(self, search_result):        
        directions=sum([sum(result.each([
                                         SequenceSearchResult.QUERY_ID_FIELD, 
                                         SequenceSearchResult.ALIGNMENT_DIRECTION
                                         ]
                                        ), []) \
                        for result in search_result], 
                       [])
        direction_dict={key: item for key, item in zip(directions[0::2], directions[1::2])}
        return direction_dict
        
    def protein_pipeline(self, base, sequence_file, unpack, gpkg):
        '''
        protein_pipeline - The main pipeline for GraftM searching for 
                           DNA sequence
        
        Parameters
        ----------
        base : str
            Base name of the file being searched. Used for creating file names
        sequence_file : str
            The path to the input sequences
        unpack : obj
            Object that builds the commad chunk for unpacking the raw sequences.
            First guesses file format, and unpacks appropriately. Calls 
            command_line to construct final command line string.
        gpkg : obj
            object with parameters set. graftM package object called with graftm_package.py
        
        Returns
        -------
        hit_aligned_reads : str
            Path to hit reads that are aligned to the alignemnt HMM and ready
            to be placed into the tree
        '''


        # Tell user what is being searched with what
        logging.info('Searching %s' % (os.path.basename(sequence_file)))
        
        # Search for reads using hmmsearch

        if not hit_reads:
            return False
        
        # Align the reads.
        logging.info('Aligning reads to reference package database')
        hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)
        self.h.align(
                     hit_reads,
                     hit_aligned_reads,
                     self._get_sequence_directions(search_result)
                     )

        return hit_aligned_reads

    def summarise(self, summary_dict):

        # Concatenate alignment files, place in tree, split output guppy
        # and .jplace file for the output

        # Summary steps.
        start           = timeit.default_timer()
        placements_list = []
        for base in summary_dict['base_list']:
            # First assign the hash that contains all of the trusted placements
            # to a variable to it can be passed to otu_builder, to be written
            # to a file. :)
            if summary_dict['reverse_pipe']:
                placements = summary_dict[base]['comparison_hash']['trusted_placements']
                summary_dict[base]['read_length'] = (summary_dict[base]['forward']['read_length'] + summary_dict[base]['reverse']['read_length'])/2
            elif not summary_dict['reverse_pipe']:
                placements = summary_dict[base]['trusted_placements']
            else:
                raise Exception('Programming Error: Assigning placements hash')
                
            self.s.readTax(placements, GraftMFiles(base, self.args.output_directory, False).read_tax_output_path(base))
            placements_list.append(placements)

        # Generate coverage table
        logging.info('Building coverage table for %s' % base)
        #self.s.coverage_of_hmm(self.args.aln_hmm_file,
        #                         self.gmf.summary_table_output_path(base),
        #                         self.gmf.coverage_table_path(base),
        #                         summary_dict[base]['read_length'])
        sample_names = summary_dict['base_list']
        logging.info('Writing summary table')
        with open(self.gmf.combined_summary_table_output_path(), 'w') as f:
            self.s.write_tabular_otu_table(sample_names, placements_list, f)
            
        logging.info('Writing biom file')
        with biom_open(self.gmf.combined_biom_output_path(), 'w') as f:
            biom_successful = self.s.write_biom(sample_names, placements_list, f)
        if not biom_successful:
            os.remove(self.gmf.combined_biom_output_path())
        
        logging.info('Building summary krona plot')
        self.s.write_krona_plot(sample_names, placements_list, self.gmf.krona_output_path())
        
        stop = timeit.default_timer()
        summary_dict['summary_t'] = str(int(round((stop - start), 0)) )

        # Compile basic run statistics if they are wanted
        summary_dict['stop_all'] = timeit.default_timer()
        summary_dict['all_t'] = str(int(round((summary_dict['stop_all'] - summary_dict['start_all']), 0)) )
        self.s.build_basic_statistics(summary_dict, self.gmf.basic_stats_path(), self.args.type)
        # Delete unnecessary files
        logging.info('Cleaning up')
        for base in summary_dict['base_list']:
            directions = ['forward', 'reverse']
            if summary_dict['reverse_pipe']:
                for i in range(0,2):
                    self.gmf = GraftMFiles(base, self.args.output_directory, directions[i])
                    self.hk.delete([self.gmf.for_aln_path(base),
                                    self.gmf.rev_aln_path(base),
                                    self.gmf.conv_output_rev_path(base),
                                    self.gmf.conv_output_for_path(base),
                                    self.gmf.euk_free_path(base),
                                    self.gmf.euk_contam_path(base),
                                    self.gmf.readnames_output_path(base),
                                    self.gmf.sto_output_path(base),
                                    self.gmf.orf_titles_output_path(base),
                                    self.gmf.orf_hmmsearch_output_path(base),
                                    self.gmf.hmmsearch_output_path(base),
                                    self.gmf.orf_output_path(base),
                                    self.gmf.comb_aln_fa(),
                                    self.gmf.output_for_path(base),
                                    self.gmf.output_rev_path(base)])
            elif not summary_dict['reverse_pipe']:
                self.gmf = GraftMFiles(base, self.args.output_directory, False)
                self.hk.delete([self.gmf.for_aln_path(base),
                                self.gmf.rev_aln_path(base),
                                self.gmf.conv_output_rev_path(base),
                                self.gmf.conv_output_for_path(base),
                                self.gmf.euk_free_path(base),
                                self.gmf.euk_contam_path(base),
                                self.gmf.readnames_output_path(base),
                                self.gmf.sto_output_path(base),
                                self.gmf.orf_titles_output_path(base),
                                self.gmf.hmmsearch_output_path(base),
                                self.gmf.orf_hmmsearch_output_path(base),
                                self.gmf.orf_output_path(base),
                                self.gmf.comb_aln_fa(),                                    
                                self.gmf.output_for_path(base),
                                self.gmf.output_rev_path(base)])

        logging.info('Done, thanks for using graftM!\n')

    def graft(self):
        # The Graft pipeline:
        # Searches for reads using hmmer, and places them in phylogenetic
        # trees to derive a community structure.
        if self.args.verbosity > 1:
            print '''
                                GRAFT
                                
                       Joel Boyd, Ben Woodcroft
                       
                                                         __/__
                                                  ______|
          _- - _                         ________|      |_____/
           - -            -             |        |____/_
           - _     >>>>  -   >>>>   ____|          
          - _-  -         -             |      ______
             - _                        |_____|
           -                                  |______
            ''' 
        REVERSE_PIPE   = (True if self.args.reverse else False)
        pair_direction = ['forward', 'reverse']
        gpkg           = GraftMPackage.acquire(self.args.graftm_package)
        base_list      = []
        seqs_list      = []
        search_results = []
        
        # Set the output directory if not specified and create that directory
        logging.debug('Creating working directory: %s' % self.args.output_directory)
        self.hk.make_working_directory(self.args.output_directory,
                                       self.args.force)
        
        # Set pipeline and evalue by checking HMM format
        hmm_type, hmm_tc = self.hk.setpipe(self.args.aln_hmm_file)
        logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))
        setattr(self.args, 'type', hmm_type)
        if hmm_tc:
            setattr(self.args, 'eval', '--cut_tc')
        # For each pair (or single file passed to GraftM)
        logging.debug('Working with %i file(s)' % len(self.sequence_pair_list))
        for pair in self.sequence_pair_list:
            # Set the basename, and make an entry to the summary table.
            base = os.path.basename(pair[0]).split('.')[0]
            logging.info("Working on %s" % base)
            
            # Guess the sequence file type, if not already specified to GraftM
            unpack = UnpackRawReads(pair[0])
            if hasattr(self.args, 'input_sequence_type'): 
                pass
            else:
                setattr(self.args, 'input_sequence_type',
                        self.hk.guess_sequence_type(unpack))
            logging.debug("Set sequence type of %s to %s" %(pair[0], self.args.input_sequence_type))
            
            # Make the working base subdirectory
            self.hk.make_working_directory(os.path.join(self.args.output_directory,
                                                        base),
                                           self.args.force)            

            # for each of the paired end read files
            for read_file in pair:
                unpack = UnpackRawReads(read_file)
                if not os.path.isfile(read_file): # Check file exists
                    logging.info('%s does not exist! Skipping this file..' % read_file)
                    continue
                # Set the output file_name
                if REVERSE_PIPE:
                    direction = pair_direction.pop(0)
                    logging.info("Working on %s reads" % direction)
                    self.gmf = GraftMFiles(base,
                                           self.args.output_directory,
                                           direction)
                    self.hk.make_working_directory(os.path.join(self.args.output_directory,
                                                                base,
                                                                direction),
                                                   self.args.force)
                else:
                    direction = False
                    self.gmf = GraftMFiles(base,
                                           self.args.output_directory,
                                           direction)
                
                if self.args.type == PIPELINE_PROTEIN:
                    logging.debug("Running protein pipeline")                    
                    hit_reads, search_result = self.h.p_search(
                                                               self.gmf,
                                                               self.args,
                                                               base,
                                                               unpack,
                                                               read_file,
                                                               self.args.search_method,
                                                               gpkg
                                                               )
                    
                # Or the DNA pipeline    
                elif self.args.type == PIPELINE_NUCLEOTIDE:
                    logging.debug("Running nucleotide pipeline")                   
                    hit_reads, search_result = self.h.d_search(
                                                              self.gmf,
                                                              self.args,
                                                              base,
                                                              unpack,
                                                              read_file,
                                                              self.args.euk_check,
                                                              self.args.search_method,
                                                              gpkg
                                                              )
                
                
                logging.info('Aligning reads to reference package database')
                hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)
                self.h.align(
                             hit_reads,
                             hit_aligned_reads,
                             self._get_sequence_directions(search_result)
                             )
                
                if not hit_aligned_reads:
                    continue
                else:
                    base_list.append(base)
                    seqs_list.append(hit_aligned_reads)
                    search_results.append(search_result)
                    
                    
        if self.args.merge_reads:
            merged_output=[GraftMFiles(base, self.args.output_directory, False).aligned_fasta_output_path(base) \
                           for base in base_list]
            self.h.merge_forev_aln(seqs_list, merged_output)
            seqs_list=merged_output
            REVERSE_PIPE = False
       
        # Leave the pipeline if search only was specified
        if self.args.search_and_align_only:
            logging.info('Stopping before placement\n')
            exit(0)
        elif not any(base_list):
            logging.error('No hits in any of the provided files. Cannot continue with no reads to place.\n')
            exit(0)
        self.gmf = GraftMFiles('',
                               self.args.output_directory,
                               False)
        logging.info("Placing reads into phylogenetic tree")
        placements=self.p.place(REVERSE_PIPE,
                                 seqs_list,
                                 self.args.resolve_placements,
                                 self.gmf,
                                 self.args)
        
        print 'stop here'
        exit()
        self.summarise(summary_table)


    def main(self):

        if self.args.subparser_name == 'graft':
            self.graft()

        elif self.args.subparser_name == 'assemble':
            if self.args.verbosity > 1: print '''
                           ASSEMBLE

                   Joel Boyd, Ben Woodcroft

          _- - _               ___            __/
           -                  /___\____      /\/
           - _     >>>>   ___/       \_\     \/
          - _-           /_/            \    /
             - _        /                \__/
                       /
'''
            self.tg.main(self.args)

        elif self.args.subparser_name == 'extract':
            if self.args.verbosity > 1: print '''
                           EXTRACT

                   Joel Boyd, Ben Woodcroft
                                  _
                         __/__/    |        >a
                  ______|          |->>>>   --------
         ________|      |_____/   _|        >b
        |        |____/_                    -------
    ____|                                   >c
        |      ______                       ----------
        |_____|
              |______
'''

            if self.args.seq:
                self.e.extract(self.args)

        elif self.args.subparser_name == 'create':
            if self.args.verbosity > 1: print '''
                            CREATE

                   Joel Boyd, Ben Woodcroft
                     
                                                    /                
              >a                                   /
              -------------                       /            
              >b                        |        |
              --------          >>>     |  GPKG  |
              >c                        |________|
              ----------     
'''

            if self.args.taxonomy:
                if self.args.rerooted_annotated_tree:
                    logging.error("--taxonomy is incompatible with --rerooted_annotated_tree")
                    exit(1)
                if self.args.taxtaxtic_taxonomy or self.args.taxtastic_seqinfo:
                    logging.error("--taxtastic_taxonomy is incompatible with --taxonomy")
                    exit(1)
            elif self.args.rerooted_annotated_tree:
                if self.args.taxtaxtic_taxonomy or self.args.taxtastic_seqinfo:
                    logging.error("--taxtastic_taxonomy is incompatible with --rerooted_annotated_tree")
                    exit(1)
            else:
                if not self.args.taxtastic_taxonomy and not self.args.taxtastic_seqinfo:
                    logging.error("--taxonomy, --rerooted_annotated_tree or --taxtastic_taxonomy/--taxtastic_seqinfo is required")
                    exit(1)
            if bool(self.args.taxtastic_taxonomy) ^  bool(self.args.taxtastic_seqinfo):
                logging.error("Both or neither of --taxtastic_taxonomy and --taxtastic_seqinfo must be defined")
                exit(1)
            self.hk.checkCreatePrerequisites()

            Create().main(sequences = self.args.sequences, 
                          alignment=self.args.alignment, 
                          taxonomy=self.args.taxonomy,
                          rerooted_tree=self.args.rerooted_tree,
                          tree_log=self.args.tree_log,
                          prefix=self.args.output,
                          rerooted_annotated_tree=self.args.rerooted_annotated_tree, 
                          min_aligned_percent=float(self.args.min_aligned_percent)/100,
                          taxtastic_taxonomy = self.args.taxtastic_taxonomy,
                          taxtastic_seqinfo = self.args.taxtastic_seqinfo
                          )

    
