#!/usr/bin/env python

import os
import logging

from graftm.sequence_search_results import SequenceSearchResult
from graftm.graftm_output_paths import GraftMFiles
from graftm.hmmer import Hmmer
from graftm.housekeeping import HouseKeeping
from graftm.summarise import Stats_And_Summary
from graftm.pplacer import Pplacer
from graftm.create import Create
from graftm.unpack_sequences import UnpackRawReads
from graftm.graftm_package import GraftMPackage

from biom.util import biom_open
from graftm.bootstrapper import Bootstrapper

PIPELINE_AA = "P"
PIPELINE_NT = "D"

class Run:
    _MIN_VERBOSITY_FOR_ART = 3 # with 2 then, only errors are printed

    def __init__(self, args):
        self.args = args
        self.setattributes(self.args)

    def setattributes(self, args):
        self.hk = HouseKeeping()
        self.s = Stats_And_Summary()
        if args.subparser_name == 'graft':
            self.hk.set_attributes(self.args)
            self.hk.set_euk_hmm(self.args)
            if args.euk_check:self.args.search_hmm_files.append(self.args.euk_hmm_file)
            self.h = Hmmer(self.args.search_hmm_files, self.args.aln_hmm_file)
            self.sequence_pair_list = self.hk.parameter_checks(args)
            if hasattr(args, 'reference_package'):
                self.p = Pplacer(self.args.reference_package)

    def _get_sequence_directions(self, search_result):
        complement_strand = {}
        for result in search_result:
            result_directions = dict(
                                     result.each([SequenceSearchResult.QUERY_ID_FIELD,
                                                  SequenceSearchResult.ALIGNMENT_DIRECTION])
                                       )
            complement_strand.update(result_directions)

        return complement_strand

    def summarise(self, base_list, trusted_placements, reverse_pipe, times, hit_read_count_list):
        '''
        summarise - write summary information to file, including otu table, biom
                    file, krona plot, and timing information

        Parameters
        ----------
        base_list : array
            list of each of the files processed by graftm, with the path and
            and suffixed removed
        trusted_placements : dict
            dictionary of placements with entry as the key, a taxonomy string
            as the value
        reverse_pipe : bool
            True = run reverse pipe, False = run normal pipeline
        times : array
            list of the recorded times for each step in the pipeline in the
            format: [search_step_time, alignment_step_time, placement_step_time]
        hit_read_count_list : array
            list containing sublists, one for each file run through the GraftM
            pipeline, each two entries, the first being the number of putative
            eukaryotic reads (when searching 16S), the second being the number
            of hits aligned and placed in the tree.
        Returns
        -------
        '''

        # Summary steps.
        placements_list = []
        for base in base_list:
            # First assign the hash that contains all of the trusted placements
            # to a variable to it can be passed to otu_builder, to be written
            # to a file. :)
            placements = trusted_placements[base]
            self.s.readTax(placements, GraftMFiles(base, self.args.output_directory, False).read_tax_output_path(base))
            placements_list.append(placements)

        #Generate coverage table
        #logging.info('Building coverage table for %s' % base)
        #self.s.coverage_of_hmm(self.args.aln_hmm_file,
        #                         self.gmf.summary_table_output_path(base),
        #                         self.gmf.coverage_table_path(base),
        #                         summary_dict[base]['read_length'])

        logging.info('Writing summary table')
        with open(self.gmf.combined_summary_table_output_path(), 'w') as f:
            self.s.write_tabular_otu_table(base_list, placements_list, f)

        logging.info('Writing biom file')
        with biom_open(self.gmf.combined_biom_output_path(), 'w') as f:
            biom_successful = self.s.write_biom(base_list, placements_list, f)
        if not biom_successful:
            os.remove(self.gmf.combined_biom_output_path())

        logging.info('Building summary krona plot')
        self.s.write_krona_plot(base_list, placements_list, self.gmf.krona_output_path())

        # Basic statistics
        placed_reads=[len(trusted_placements[base]) for base in base_list]
        self.s.build_basic_statistics(times, hit_read_count_list, placed_reads, \
                                      base_list, self.gmf.basic_stats_path())

        # Delete unnecessary files
        logging.info('Cleaning up')
        for base in base_list:
            directions = ['forward', 'reverse']
            if reverse_pipe:
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
            else:
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

        if self.args.graftm_package:
            gpkg = GraftMPackage.acquire(self.args.graftm_package)
        else:
            gpkg = None

        REVERSE_PIPE        = (True if self.args.reverse else False)
        base_list           = []
        seqs_list           = []
        search_results      = []
        clusters_list       = []
        hit_read_count_list = []
        
        # Get the maximum range, if none exists, make one from the HMM profile
        if gpkg:
            maximum_range = gpkg.maximum_range()
            diamond_db    = gpkg.diamond_database_path()
            if self.args.search_method == 'diamond':
                if not diamond_db:
                    if self.args.search_diamond_file:
                        diamond_db=self.args.search_diamond_file
                    else:
                        logging.error("%s search method selected, but no diamond database specified. \
                        Please either provide a gpkg to the --graftm_package flag, or a diamond \
                        database to the --search_diamond_file flag." % self.args.search_method)
        else:
            if self.args.maximum_range:
                maximum_range = self.args.maximum_range
            else:
                if self.args.search_method=='hmmsearch':
                    maximum_range = self.hk.get_maximum_range(self.args.aln_hmm_file)
                else:
                    logging.warning('Cannot determine maximum range when using %s pipeline and with no GraftM package specified' % self.args.search_method)
                    logging.warning('Setting maximum_range to None (linked hits will not be detected)')
                    maximum_range = None 
            if self.args.search_diamond_file:
                diamond_db = self.args.search_diamond_file
            else:
                if self.args.search_method == 'hmmsearch':
                    diamond_db = None
                else:
                    logging.error("%s search method selected, but no gpkg or diamond database selected" % self.args.search_method)

            
        # If merge reads is specified, check that there are reverse reads to merge with
        if self.args.merge_reads and not hasattr(self.args, 'reverse'):
            logging.error("--merge requires --reverse to be specified")
            exit(1)

        # Set the output directory if not specified and create that directory
        logging.debug('Creating working directory: %s' % self.args.output_directory)
        self.hk.make_working_directory(self.args.output_directory,
                                       self.args.force)

        # Set pipeline and evalue by checking HMM format
        hmm_type, hmm_tc = self.hk.setpipe(self.args.aln_hmm_file)
        logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))
        setattr(self.args, 'type', hmm_type)
        if hmm_tc:
            setattr(self.args, 'evalue', '--cut_tc')
        
        # For each pair (or single file passed to GraftM)
        logging.debug('Working with %i file(s)' % len(self.sequence_pair_list))
        for pair in self.sequence_pair_list:
            # Set the basename, and make an entry to the summary table.
            base = os.path.basename(pair[0]).split('.')[0]
            pair_direction = ['forward', 'reverse']
            logging.info("Working on %s" % base)

            # Guess the sequence file type, if not already specified to GraftM
            unpack = UnpackRawReads(pair[0])

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

                if self.args.type == PIPELINE_AA:
                    logging.debug("Running protein pipeline")
                    if self.args.bootstrap_contigs:
                        new_hmm = self.gmf.bootstrap_hmm_path()
                        if self.args.graftm_package:
                            pkg = GraftMPackage.acquire(self.args.graftm_package)
                        else:
                            pkg = None
                        boots = Bootstrapper(
                            search_hmm_files = self.args.search_hmm_files,
                            maximum_range = self.args.maximum_range,
                            threads = self.args.threads,
                            evalue = self.args.evalue,
                            min_orf_length = self.args.min_orf_length,
                            graftm_package = pkg)
                        if boots.generate_hmm_from_contigs(
                                                 self.args.bootstrap_contigs,
                                                 new_hmm):
                            self.h.search_hmm.append(new_hmm)

                    search_time, result = self.h.aa_db_search(
                                                              self.gmf,
                                                              base,
                                                              unpack,
                                                              self.args.search_method,
                                                              maximum_range,
                                                              self.args.threads,
                                                              self.args.evalue,
                                                              self.args.min_orf_length,
                                                              self.args.restrict_read_length,
                                                              diamond_db
                                                              )

                # Or the DNA pipeline
                elif self.args.type == PIPELINE_NT:
                    logging.debug("Running nucleotide pipeline")
                    search_time, result = self.h.nt_db_search(
                                                              self.gmf,
                                                              base,
                                                              unpack,
                                                              self.args.euk_check,
                                                              self.args.search_method,
                                                              maximum_range,
                                                              self.args.threads,
                                                              self.args.evalue
                                                              )

                if not result.hit_fasta():
                    logging.info('No reads found to align in %s' % base)
                    continue
                logging.info('Aligning reads to reference package database')
                hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)
                
                if self.args.cluster:
                    result.cluster()
                    clusters_list.append(result.clusters)
                    
                aln_time = self.h.align(
                                        result.hit_fasta(),
                                        hit_aligned_reads,
                                        self._get_sequence_directions(result.search_result)
                                        )

                if not hit_aligned_reads:
                    continue
                else:
                    base_list.append(base)
                    seqs_list.append(hit_aligned_reads)
                    search_results.append(result.search_result)
                    hit_read_count_list.append(result.hit_count)

        if self.args.merge_reads:
            base_list=base_list[0::2]
            merged_output=[GraftMFiles(base, self.args.output_directory, False).aligned_fasta_output_path(base) \
                           for base in base_list]
            self.h.merge_forev_aln(seqs_list[0::2], seqs_list[1::2], merged_output)
            seqs_list=merged_output
            REVERSE_PIPE = False
        
        elif REVERSE_PIPE:
            base_list=base_list[0::2]

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
        place_time, placements=self.p.place(
                                            REVERSE_PIPE,
                                            seqs_list,
                                            self.args.resolve_placements,
                                            self.gmf,
                                            self.args,
                                            result.slash_endings,
                                            gpkg.taxonomy_info_path()
                                            )

        if self.args.cluster:
            for idx, base in enumerate(base_list):
                clusters_taxonomy={}
                place    = placements[base]
                clusters = clusters_list[idx]
                for read_name in place.keys():
                    for record in clusters[read_name]:
                        clusters_taxonomy[record.name]=place[read_name]            
                placements[base].update(clusters_taxonomy)

        self.summarise(base_list, placements, REVERSE_PIPE, 
                       [search_time,aln_time[0],place_time], 
                       hit_read_count_list)


    def main(self):

        if self.args.subparser_name == 'graft':
            if self.args.verbosity >= self._MIN_VERBOSITY_FOR_ART: print '''
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
            self.graft()

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

            Create().main(
                          sequences = self.args.sequences,
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
        elif self.args.subparser_name == 'bootstrap':
            args = self.args
            if not args.graftm_package and not args.search_hmm_files:
                logging.error("bootstrap mode requires either --graftm_package or --search_hmm_files")
                exit(1)
                
            if args.graftm_package:
                pkg = GraftMPackage.acquire(args.graftm_package)
            else:
                pkg = None
                
            strapper = Bootstrapper(search_hmm_files = args.search_hmm_files,
                maximum_range = args.maximum_range,
                threads = args.threads,
                evalue = args.evalue,
                min_orf_length = args.min_orf_length,
                graftm_package = pkg)
            strapper.generate_hmm_from_contigs([args.contigs], args.output_hmm)

