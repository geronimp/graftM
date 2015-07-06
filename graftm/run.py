#!/usr/bin/env python

import timeit
import os
import logging

from graftm.graftm_output_paths import GraftMFiles
from graftm.extract_sequences import Extract
from graftm.hmmer import Hmmer
from graftm.housekeeping import HouseKeeping
from graftm.summarise import Stats_And_Summary
from graftm.pplacer import Pplacer
from graftm.assembler import TaxoGroup
from graftm.create import Create
from graftm.unpack_sequences import UnpackRawReads

from biom.util import biom_open
from _struct import unpack

class Run:
    _MIN_VERBOSITY_FOR_ART = 3 # with 2 then, only errors are printed

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

    def protein_pipeline(self, base, summary_dict, sequence_file, direction, unpack):
        'The main pipeline for GraftM finding protein hits in the input sequence'
        # Set a variable to store the run statistics, to be added later to
        # the summary_dict
        if direction:
            run_stats = summary_dict[base][direction]
        elif not direction:
            run_stats = summary_dict[base]
        else:
            raise Exception('Programming Error: Assigning run_stats hash')
        # Tell user what is being searched with what
        logging.info('Searching %s' % (os.path.basename(sequence_file)))
        # Search for reads using hmmsearch
        hit_reads, run_stats = self.h.p_search(self.gmf,
                                               self.args,
                                               run_stats,
                                               base,
                                               unpack,
                                               sequence_file)
        if not hit_reads:
            return summary_dict, False
        # Align the reads.
        logging.info('Aligning reads to reference package database')

        hit_aligned_reads, run_stats = self.h.align(self.gmf,
                                                    self.args,
                                                    run_stats,
                                                    base,
                                                    hit_reads)
        # Set these paramaters as N/A 'cos they don't apply to the protein pipeline
        run_stats['n_contamin_euks'] = 'N/A'
        run_stats['n_uniq_euks'] = 'N/A'
        run_stats['euk_check_t'] = 'N/A'
        if direction:
            summary_dict[base][direction] = run_stats
        elif not direction:
            summary_dict[base] = run_stats
        else:
            raise Exception('Programming Error: Logging %s hash' % direction)

        return summary_dict, hit_aligned_reads

    def dna_pipeline(self, base, summary_dict, sequence_file, direction, unpack):
        'The main pipeline for GraftM searching for DNA sequence'
        # Set a variable to store the run statistics, to be added later to
        # the summary_dict
        if direction:
            run_stats = summary_dict[base][direction]
        elif not direction:
            run_stats = summary_dict[base]
        else:
            raise Exception('Programming Error: Assigning run_stats hash')
        # Search for reads using nhmmer
        logging.info('Searching %s' % os.path.basename(sequence_file))
        hit_reads, run_stats = self.h.d_search(self.gmf,
                                               self.args,
                                               run_stats,
                                               base,
                                               unpack,
                                               sequence_file,
                                               summary_dict['euks_checked'])
        
        if not hit_reads:
            return summary_dict, False
        
        # Otherwise, run through the alignment
        logging.info('Aligning reads to reference package database')
        hit_aligned_reads, run_stats = self.h.align(self.gmf,
                                                    self.args,
                                                    run_stats,
                                                    base,
                                                    hit_reads)
        if direction:
            summary_dict[base][direction] = run_stats
        elif not direction:
            summary_dict[base] = run_stats
        else:
            raise Exception('Programming Error: Logging %s hash' % direction)
        return summary_dict, hit_aligned_reads

    def placement(self, summary_dict):
        ## This is the placement pipeline in GraftM, in aligned reads are
        ## placed into phylogenetic trees, and the results interpreted.
        ## If reverse reads are used, this is where the comparisons are made
        ## between placements, for the summary tables to be build in the
        ## next stage.
        # Concatenate alignment files, place in tree, split output guppy
        # and .jplace file for the output
        
        summary_dict = self.p.place(summary_dict,
                                    self.gmf,
                                    self.args)
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
        if self.args.verbosity >= self._MIN_VERBOSITY_FOR_ART:
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
        
        if self.args.merge_reads and not hasattr(self.args, 'reverse'):
            logging.error("--merge requires --reverse to be specified")
            exit(1)
        
        readstoplace=False # An extra check to make sure there are reads to place with
        # Set up a dictionary that will record stats as the pipeline is running
        summary_table = {'euks_checked'      : self.args.euk_check,
                         'resolve_placements': self.args.resolve_placements,
                         'base_list'         : [],
                         'seqs_list'         : [],
                         'start_all'         : timeit.default_timer(),
                         'reverse_pipe'      : False,
                         'merge_reads'       : self.args.merge_reads}
        
        # Set the output directory if not specified and create that directory
        logging.debug('Creating working directory: %s' % self.args.output_directory)
        self.hk.make_working_directory(self.args.output_directory,
                                       self.args.force)

        # For each pair (or single file passed to GraftM)
        logging.debug('Working with %i file(s)' % len(self.sequence_pair_list))
        for pair in self.sequence_pair_list:
            # Set the basename, and make an entry to the summary table.
            base = os.path.basename(pair[0]).split('.')[0]
            # Set reverse pipe if more than one pair
            if hasattr(self.args, 'reverse'):
                summary_table['reverse_pipe'] = True
                summary_table[base]           = {'reverse':{}, 
                                                 'forward':{}}
                pair_direction                = ['forward', 'reverse']
            else:
                summary_table[base]           = {}
            
            # Set pipeline and evalue by checking HMM format
            hmm_type, hmm_tc = self.hk.setpipe(self.args.aln_hmm_file)
            logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))
            setattr(self.args, 'type', hmm_type)
            if hmm_tc:
                setattr(self.args, 'eval', '--cut_tc')
                
            # Guess the sequence file type, if not already specified to GraftM
            unpack = UnpackRawReads(pair[0])
            if hasattr(self.args, 'input_sequence_type'): 
                pass
            else:
                setattr(self.args, 'input_sequence_type',
                        self.hk.guess_sequence_type(unpack))
            logging.debug("Set sequence type of %s to %s" %(pair[0], self.args.input_sequence_type))
            # Make the working base directory
            self.hk.make_working_directory(os.path.join(self.args.output_directory,
                                                        base),
                                           self.args.force)
            # tell the user which file/s is being processed
            logging.info("Working on %s" % base)

            # for each of the paired end read files
            for read_file in pair:
                if not os.path.isfile(read_file): # Check file exists
                    logging.info('%s does not exist! Skipping this file..' % read_file)
                    continue
                # Set the output file_name
                if summary_table['reverse_pipe']:
                    direction = pair_direction.pop(0)
                    logging.info("Working on %s reads" % direction)
                    self.gmf = GraftMFiles(base,
                                           self.args.output_directory,
                                           direction)
                    self.hk.make_working_directory(os.path.join(self.args.output_directory,
                                                                base,
                                                                direction),
                                                   self.args.force)
                elif not summary_table['reverse_pipe']:
                    direction = False
                    self.gmf = GraftMFiles(base,
                                           self.args.output_directory,
                                           direction)
                else:
                    raise Exception('Programming Error')
                
                if self.args.type == 'P':
                    logging.debug("Running protein pipeline")
                    summary_table, hit_aligned_reads = self.protein_pipeline(base,
                                                                             summary_table,
                                                                             read_file,
                                                                             direction,
                                                                             unpack)
                # Or the DNA pipeline
                elif self.args.type == 'D':
                    logging.debug("Running nucleotide pipeline")
                    summary_table, hit_aligned_reads = self.dna_pipeline(base,
                                                                         summary_table,
                                                                         read_file,
                                                                         direction,
                                                                         unpack)
                if not hit_aligned_reads:
                    continue
                else:
                    readstoplace=True
                    
                # Add the run stats and the completed run to the summary table
                summary_table['seqs_list'].append(hit_aligned_reads)
                if base not in summary_table['base_list']:
                    summary_table['base_list'].append(base)
        if summary_table['merge_reads']:
            merged_output=[GraftMFiles(base, self.args.output_directory, False).aligned_fasta_output_path(base) \
                           for base in summary_table['base_list']]
            self.h.merge_forev_aln(summary_table['seqs_list'], merged_output)
            summary_table['seqs_list']=[GraftMFiles(base, self.args.output_directory, False).aligned_fasta_output_path(base) \
                                       for base in summary_table['base_list']]
            summary_table['reverse_pipe']=False
        # Leave the pipeline if search only was specified
        if self.args.search_and_align_only:
            logging.info('Stopping before placement\n')
            exit(0)
        elif not readstoplace:
            logging.error('No hits in any of the provided files. Cannot continue with no reads to place.\n')
            exit(0)
        # Tell the user we're on to placing the sequences into the tree.
        self.gmf = GraftMFiles('',
                               self.args.output_directory,
                               False)
        logging.info("Placing reads into phylogenetic tree")
        self.placement(summary_table)


    def main(self):

        if self.args.subparser_name == 'graft':
            self.graft()

        elif self.args.subparser_name == 'assemble':
            if self.args.verbosity >= self._MIN_VERBOSITY_FOR_ART: print '''
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
            if self.args.verbosity >= self._MIN_VERBOSITY_FOR_ART: print '''
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
            if self.args.verbosity >= self._MIN_VERBOSITY_FOR_ART: print '''
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
                if self.args.taxtastic_taxonomy or self.args.taxtastic_seqinfo:
                    logging.error("--taxtastic_taxonomy and --taxtastic_seqinfo are incompatible with --taxonomy")
                    exit(1)
            elif self.args.rerooted_annotated_tree:
                if self.args.taxtastic_taxonomy or self.args.taxtastic_seqinfo:
                    logging.error("--taxtastic_taxonomy and --taxtastic_seqinfo are incompatible with --rerooted_annotated_tree")
                    exit(1)
            else:
                if not self.args.taxtastic_taxonomy or not self.args.taxtastic_seqinfo:
                    logging.error("--taxonomy, --rerooted_annotated_tree or --taxtastic_taxonomy/--taxtastic_seqinfo is required")
                    exit(1)
            if bool(self.args.taxtastic_taxonomy) ^  bool(self.args.taxtastic_seqinfo):
                logging.error("Both or neither of --taxtastic_taxonomy and --taxtastic_seqinfo must be defined")
                exit(1)
            self.hk.checkCreatePrerequisites()

            Create().main(alignment=self.args.alignment, 
                          taxonomy=self.args.taxonomy,
                          rerooted_tree=self.args.rerooted_tree,
                          tree_log=self.args.tree_log,
                          prefix=self.args.output,
                          rerooted_annotated_tree=self.args.rerooted_annotated_tree, 
                          min_aligned_percent=float(self.args.min_aligned_percent)/100,
                          taxtastic_taxonomy = self.args.taxtastic_taxonomy,
                          taxtastic_seqinfo = self.args.taxtastic_seqinfo
                          )

    
