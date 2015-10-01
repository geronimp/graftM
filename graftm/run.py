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
from graftm.bootstrapper import Bootstrapper
from graftm.diamond import Diamond
from graftm.getaxnseq import Getaxnseq
from graftm.sequence_io import SequenceIO
from graftm.timeit import Timer
from graftm.clusterer import Clusterer
from graftm.decorator import Decorator

from biom.util import biom_open

T=Timer()

PIPELINE_AA = "P"
PIPELINE_NT = "D"

class Run:
    _MIN_VERBOSITY_FOR_ART = 3 # with 2 then, only errors are printed
    PPLACER_TAXONOMIC_ASSIGNMENT = 'pplacer'
    DIAMOND_TAXONOMIC_ASSIGNMENT = 'diamond'

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
        db_search_results   = []

        if gpkg:
            maximum_range = gpkg.maximum_range()
            diamond_db    = gpkg.diamond_database_path()
            if self.args.search_method == 'diamond':                    
                if self.args.search_diamond_file:
                    diamond_db=self.args.search_diamond_file[0]
                if not diamond_db:
                    logging.error("%s search method selected, but no diamond database specified. \
                    Please either provide a gpkg to the --graftm_package flag, or a diamond \
                    database to the --search_diamond_file flag." % self.args.search_method)
                    raise Exception()
        else:
            # Get the maximum range, if none exists, make one from the HMM profile
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

        if self.args.assignment_method == Run.DIAMOND_TAXONOMIC_ASSIGNMENT:
            if self.args.reverse:
                logging.warn("--reverse reads specified with --assignment_method diamond. Reverse reads will be ignored.")
                self.args.reverse = None

       
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
            
        # Generate bootstrap database if required
        if self.args.bootstrap_contigs:
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
            
            #this is a hack, it should really use GraftMFiles but that class isn't currently flexible enough
            new_database = (os.path.join(self.args.output_directory, "bootstrap.hmm") \
                            if self.args.search_method == "hmmsearch" \
                            else os.path.join(self.args.output_directory, "bootstrap")
                            )
            
            if boots.generate_bootstrap_database_from_contigs(
                                     self.args.bootstrap_contigs,
                                     new_database,
                                     self.args.search_method):
                if self.args.search_method == "hmmsearch":
                    self.h.search_hmm.append(new_database)
                else:
                    diamond_db = new_database
            

        # For each pair (or single file passed to GraftM)
        logging.debug('Working with %i file(s)' % len(self.sequence_pair_list))
        for pair in self.sequence_pair_list:
            # Guess the sequence file type, if not already specified to GraftM
            unpack = UnpackRawReads(pair[0])

            # Set the basename, and make an entry to the summary table.
            base = unpack.basename()
            pair_direction = ['forward', 'reverse']
            logging.info("Working on %s" % base)

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

                
                    search_time, (result, complement_information) = self.h.aa_db_search(
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
                    search_time, (result, complement_information)  = self.h.nt_db_search(
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
                    logging.info('No reads found in %s' % base)
                    continue

                if self.args.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
                    logging.info('Aligning reads to reference package database')
                    hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)

                    aln_time, aln_result = self.h.align(
                                                        result.hit_fasta(),
                                                        hit_aligned_reads,
                                                        complement_information,
                                                        self.args.type
                                                        )
                    if aln_result:
                        seqs_list.append(hit_aligned_reads)
                    else:
                        logging.info("No more aligned sequences to place!")
                        continue
                
                db_search_results.append(result)
                base_list.append(base)
                search_results.append(result.search_result)
                hit_read_count_list.append(result.hit_count)

        if self.args.merge_reads: # not run when diamond is the assignment mode- enforced by argparse grokking
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
            logging.info('Stopping before taxonomic assignment phase\n')
            exit(0)
        elif not any(base_list):
            logging.error('No hits in any of the provided files. Cannot continue with no reads to assign taxonomy to.\n')
            exit(0)
        self.gmf = GraftMFiles('',
                               self.args.output_directory,
                               False)

        if self.args.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
            # Classification steps        
            if not self.args.no_clustering:
                C=Clusterer()
                seqs_list=C.cluster(seqs_list, REVERSE_PIPE)
            logging.info("Placing reads into phylogenetic tree")
            taxonomic_assignment_time, assignments=self.p.place(REVERSE_PIPE,
                                                                seqs_list,
                                                                self.args.resolve_placements,
                                                                self.gmf,
                                                                self.args,
                                                                result.slash_endings,
                                                                gpkg.taxtastic_taxonomy_path())
            if not self.args.no_clustering:
                assignments = C.uncluster_annotations(assignments, REVERSE_PIPE)
        
        elif self.args.assignment_method == Run.DIAMOND_TAXONOMIC_ASSIGNMENT:
            logging.info("Assigning taxonomy with diamond")
            taxonomic_assignment_time, assignments = self._assign_taxonomy_with_diamond(\
                        base_list,
                        db_search_results,
                        gpkg,
                        self.gmf)
            aln_time = 'n/a'
        else: raise Exception("Unexpected assignment method encountered: %s" % self.args.placement_method)
        self.summarise(base_list, assignments, REVERSE_PIPE,
                       [search_time,aln_time,taxonomic_assignment_time],
                       hit_read_count_list)

    @T.timeit
    def _assign_taxonomy_with_diamond(self, base_list, db_search_results,
                                      graftm_package, graftm_files):
        '''Run diamond to assign taxonomy

        Parameters
        ----------
        base_list: list of str
            list of sequence block names
        db_search_results: list of DBSearchResult
            the result of running hmmsearches
        graftm_package: GraftMPackage object
            Diamond is run against this database
        graftm_files: GraftMFiles object
            Result files are written here

        Returns
        -------
        list of
        1. time taken for assignment
        2. assignments i.e. dict of base_list entry to dict of read names to
            to taxonomies, or None if there was no hit detected.
        '''
        runner = Diamond(graftm_package.diamond_database_path(),
                         self.args.threads,
                         self.args.evalue)
        taxonomy_definition = Getaxnseq().read_taxtastic_taxonomy_and_seqinfo\
                (open(graftm_package.taxtastic_taxonomy_path()),
                 open(graftm_package.taxtastic_seqinfo_path()))
        results = {}

        # For each of the search results,
        for i, search_result in enumerate(db_search_results):
            sequence_id_to_hit = {}
            # Run diamond
            logging.debug("Running diamond on %s" % search_result.hit_fasta())
            diamond_result = runner.run(search_result.hit_fasta(), UnpackRawReads.PROTEIN_SEQUENCE_TYPE)
            for res in diamond_result.each([SequenceSearchResult.QUERY_ID_FIELD,
                                            SequenceSearchResult.HIT_ID_FIELD]):
                if res[0] in sequence_id_to_hit:
                    # do not accept duplicates
                    if sequence_id_to_hit[res[0]] != res[1]:
                        raise Exception("Diamond unexpectedly gave two hits for a single query sequence for %s" % res[0])
                else:
                    sequence_id_to_hit[res[0]] = res[1]

            # Extract taxonomy of the best hit, and add in the no hits
            sequence_id_to_taxonomy = {}
            for seqio in SequenceIO().read_fasta_file(search_result.hit_fasta()):
                name = seqio.name
                if name in sequence_id_to_hit:
                    # Add Root; to be in line with pplacer assignment method
                    sequence_id_to_taxonomy[name] = ['Root']+taxonomy_definition[sequence_id_to_hit[name]]
                else:
                    # picked up in the initial search (by hmmsearch, say), but diamond misses it
                    sequence_id_to_taxonomy[name] = ['Root']

            results[base_list[i]] = sequence_id_to_taxonomy
        return results

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
            if self.args.alignment and self.args.hmm:
                logging.error("--alignment and --hmm cannot both be set")
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
                          taxtastic_seqinfo = self.args.taxtastic_seqinfo,
                          hmm = self.args.hmm,
                          search_hmm_files = self.args.search_hmm_files,
                          force = self.args.force
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
            strapper.generate_bootstrap_database_from_contigs(args.contigs, args.output_hmm)

        
        
        
        elif self.args.subparser_name == 'decorate':
            
            if self.args.rooted_tree:
                if self.args.unrooted_tree:
                    logging.error("Both a rooted tree and an un-rooted tree were provided, so it's unclear what you are asking graftM to do. \
If you're unsure how to use GraftM decorate use graftM decorate -h")
                    exit(1)
                elif self.args.reference_tree:
                    logging.error("Both a rooted tree and reference tree were provided, so it's unclear what you are asking graftM to do. \
If you're unsure how to use GraftM decorate use graftM decorate -h")
                    exit(1)
                
                dec = Decorator(tree_path = self.args.rooted_tree)
            
            elif self.args.unrooted_tree and self.args.reference_tree:
                logging.debug("Using provided reference tree %s to reroot %s and decorating with %s" % (self.args.reference_tree, 
                                                                                                        self.args.unrooted_tree, 
                                                                                                        self.args.input_taxonomy))
                dec = Decorator(reference_tree_path = self.args.reference_tree,
                                tree_path = self.args.unrooted_tree)
            
            if self.args.input_greengenes_taxonomy:
                if self.args.input_taxtastic_seqinfo or self.args.input_taxtastic_taxonomy:
                    logging.error("Both taxtastic and greengenes taxonomy were provided, so its unclear what taxonomy you want graftM to decorate with")
                    exit(1)
                dec.main(self.args.input_greengenes_taxonomy, self.args.output_tree, self.args.output_taxonomy) 
            elif self.args.input_taxtastic_seqinfo and self.args.input_taxtastic_taxonomy:
                dec.main(self.args.input_taxtastic_taxonomy, self.args.output_tree, self.args.output_taxonomy, self.args.input_taxtastic_seqinfo) 
            else:
                logging.error("Only a taxtastic file or a seqinfo file were provided. graftM cannot continue without both.")
                exit(1)        
            
            
            
            
