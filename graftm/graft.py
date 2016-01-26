#!/usr/bin/env python
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
###
import os
import logging
import shutil

from graftm.timeit import Timer
from graftm.sequence_search_results import SequenceSearchResult
from graftm.graftm_output_paths import GraftMFiles
from graftm.search_table import SearchTableWriter
from graftm.unpack_sequences import UnpackRawReads
from graftm.graftm_package import GraftMPackage
from graftm.bootstrapper import Bootstrapper
from graftm.diamond import Diamond
from graftm.getaxnseq import Getaxnseq
from graftm.sequence_io import SequenceIO
from graftm.clusterer import Clusterer
from graftm.hmmer import Hmmer

from biom.util import biom_open

T = Timer()

class Graft:
    
    PPLACER_TAXONOMIC_ASSIGNMENT = 'pplacer'
    DIAMOND_TAXONOMIC_ASSIGNMENT = 'diamond'
    
    PIPELINE_AA = "P"
    PIPELINE_NT = "D"
    
    def __init__(self, **kwargs):
        self.assignment_method = kwargs.pop('assignment_method', None)
        self.bootstrap_contigs = kwargs.pop('bootstrap_contigs', None)
        self.euk_check = kwargs.pop('euk_check', None)
        self.evalue = kwargs.pop('evalue', None)
        self.filter_minimum = kwargs.pop('filter_minimum', None)
        self.force = kwargs.pop('force', None)
        self.graftm_packages = kwargs.pop('graftm_packages', None)
        self.maximum_range = kwargs.pop('maximum_range', None)
        self.merge_reads = kwargs.pop('merge_reads', None)
        self.min_orf_length = kwargs.pop('min_orf_length', None)
        self.no_clustering = kwargs.pop('no_clustering', None)
        self.output_directory = kwargs.pop('output_directory', None)
        self.placement_method = kwargs.pop('placement_method', None)
        self.resolve_placements = kwargs.pop('resolve_placements', None)
        self.restrict_read_length = kwargs.pop('restrict_read_length', None)
        self.reverse = kwargs.pop('reverse', None)
        self.search_and_align_only = kwargs.pop('search_and_align_only', None)
        self.search_diamond_file = kwargs.pop('search_diamond_file', None)
        self.search_hmm_files = kwargs.pop('search_hmm_files', None)
        self.search_method = kwargs.pop('search_method', None)
        self.search_only = kwargs.pop('search_only', None)
        self.threads = kwargs.pop('threads', None)
        self.sequence_pair_list = kwargs.pop('sequence_pair_list', None)
        self.type = kwargs.pop("type", None)
        
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        # Set the output directory if not specified and create that directory
        logging.debug('Creating working directory: %s' % self.output_directory)
        self.make_working_directory(self.output_directory,
                                    self.force)
        
        
        # Generate bootstrap database if required
        if self.bootstrap_contigs:
            boots_dict = {}
            if self.graftm_packages:
                for gpkg in self.graftm_packages.values():
                    boots_dict[gpkg.name] = Bootstrapper(
                                    search_hmm_files = gpkg.search_hmm_paths(),
                                    maximum_range = self.maximum_range,
                                    threads = self.threads,
                                    evalue = self.evalue,
                                    min_orf_length = self.min_orf_length,
                                    graftm_package = gpkg)
            else:
                boots_dict[None] = Bootstrapper(
                                                search_hmm_files = self.search_hmm_files,
                                                maximum_range = self.maximum_range,
                                                threads = self.threads,
                                                evalue = self.evalue,
                                                min_orf_length = self.min_orf_length,
                                                graftm_package = None)          
            
            for gene, boots in boots_dict.iteritems():
                
                # this is a hack, it should really use GraftMFiles but that class isn't currently flexible enough
                if gene:
                    name = "%s_bootstrap" % (gene)
                else:
                    name = "bootstrap"
                    
                new_database = (os.path.join(self.output_directory, "%s.hmm" % name) \
                                if self.search_method == Hmmer.HMM_SEARCH_METHOD \
                                else os.path.join(self.output_directory, name)
                                )
                
                if boots.generate_bootstrap_database_from_contigs(
                                         self.bootstrap_contigs,
                                         new_database,
                                         self.search_method):
                    
                    if self.search_method == Hmmer.HMM_SEARCH_METHOD:
                        self.h.search_hmm.append(new_database)
                        
                    else:
                        diamond_db = new_database
    
        # Set up hmmer
        self.h = Hmmer(self.search_hmm_files)
        
        import IPython ; IPython.embed()
        
    def graft(self):
     
        REVERSE_PIPE        = (True if self.reverse else False)
        base_list           = []
        seqs_list           = []
        search_results      = []
        clusters_list       = []
        hit_read_count_list = []
        db_search_results   = []
        
        # Set pipeline and evalue by checking HMM format       
        
#        if self.graftm_package:
#            gpkg = GraftMPackage.acquire(self.graftm_package)
#        else:
#            gpkg = None


#        if gpkg:
#            maximum_range = gpkg.maximum_range()
#            diamond_db    = gpkg.diamond_database_path()
#            if self.search_method == 'diamond':                    
#                if self.search_diamond_file:
#                    diamond_db=self.search_diamond_file[0]
#                if not diamond_db:
#                    raise Exception("%s search method selected, but no diamond database specified. \
#                    Please either provide a gpkg to the --graftm_package flag, or a diamond \
#                    database to the --search_diamond_file flag." % self.search_method)
#        else:
#            # Get the maximum range, if none exists, make one from the HMM profile
#            if self.maximum_range:
#                maximum_range = self.maximum_range
#            else:
#                if self.search_method=='hmmsearch':
#                    if not self.search_only:
#                        maximum_range = self.hk.get_maximum_range(self.aln_hmm_file)
#                    else:
#                        logging.debug("Running search only pipeline. maximum_range not configured.")
#                        maximum_range = None
#                else:
#                    logging.warning('Cannot determine maximum range when using %s pipeline and with no GraftM package specified' % self.search_method)
#                    logging.warning('Setting maximum_range to None (linked hits will not be detected)')
#                    maximum_range = None
#            if self.search_diamond_file:
#                diamond_db = self.search_diamond_file
#            else:
#                if self.search_method == 'hmmsearch':
#                    diamond_db = None
#                else:
#                    logging.error("%s search method selected, but no gpkg or diamond database selected" % self.search_method)



        



            

            

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
            self.make_working_directory(os.path.join(self.output_directory,
                                                        base),
                                           self.force)

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
                                           self.output_directory,
                                           direction)
                    self.make_working_directory(os.path.join(self.output_directory,
                                                                base,
                                                                direction),
                                                   self.force)
                else:
                    direction = False
                    self.gmf = GraftMFiles(base,
                                           self.output_directory,
                                           direction)

                if self.type == self.PIPELINE_AA:
                    logging.debug("Running protein pipeline")
                    maximum_range = None
                    search_time, (result, complement_information) = self.h.aa_db_search(
                                                              self.gmf,
                                                              base,
                                                              unpack,
                                                              self.search_method,
                                                              maximum_range,
                                                              self.threads,
                                                              self.evalue,
                                                              self.min_orf_length,
                                                              self.restrict_read_length,
                                                              diamond_db
                                                              )


                # Or the DNA pipeline
                elif self.type == self.PIPELINE_NT:
                    logging.debug("Running nucleotide pipeline")
                    search_time, (result, complement_information)  = self.h.nt_db_search(
                                                              self.gmf,
                                                              base,
                                                              unpack,
                                                              self.euk_check,
                                                              self.search_method,
                                                              maximum_range,
                                                              self.threads,
                                                              self.evalue
                                                              )

                if not result.hit_fasta():
                    logging.info('No reads found in %s' % base)
                    continue
                elif self.search_only:
                    db_search_results.append(result)
                    base_list.append(base)

                    continue
                
                if self.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
                    logging.info('aligning reads to reference package database')
                    hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)

                    aln_time, aln_result = self.h.align(
                                                        result.hit_fasta(),
                                                        hit_aligned_reads,
                                                        complement_information,
                                                        self.type,
                                                        filter_minimum
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
        
        # Write summary table
        
        srchtw = SearchTableWriter()
        srchtw.build_search_otu_table([x.search_objects for x in db_search_results],
                                      base_list,
                                      self.gmf.search_otu_table())
        
        if self.search_only:
            logging.info('Stopping before alignment and taxonomic assignment phase\n')
            exit(0)
        
        
        if self.merge_reads: # not run when diamond is the assignment mode- enforced by argparse grokking
            base_list=base_list[0::2]
            merged_output=[GraftMFiles(base, self.output_directory, False).aligned_fasta_output_path(base) \
                           for base in base_list]
            self.h.merge_forev_aln(unpack.slash_endings,
                                   seqs_list[0::2], seqs_list[1::2], 
                                   merged_output)
            seqs_list=merged_output
            REVERSE_PIPE = False
        elif REVERSE_PIPE:
            base_list=base_list[0::2]

        # Leave the pipeline if search only was specified
        if self.search_and_align_only:
            logging.info('Stopping before taxonomic assignment phase\n')
            exit(0)
        elif not any(base_list):
            logging.error('No hits in any of the provided files. Cannot continue with no reads to assign taxonomy to.\n')
            exit(0)
        self.gmf = GraftMFiles('',
                               self.output_directory,
                               False)

        if self.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
            # Classification steps        
            if not self.no_clustering:
                C=Clusterer(unpack.slash_endings)
                cluster_dictionary = C.seq_library

                seqs_list=C.cluster(seqs_list, REVERSE_PIPE)
            logging.info("Placing reads into phylogenetic tree")
            taxonomic_assignment_time, assignments=self.p.place(REVERSE_PIPE,
                                                                seqs_list,
                                                                self.resolve_placements,
                                                                self.gmf,
                                                                self,
                                                                result.slash_endings,
                                                                gpkg.taxtastic_taxonomy_path(),
                                                                cluster_dictionary)
            if not self.no_clustering:
                assignments = C.uncluster_annotations(assignments, REVERSE_PIPE,
                                                      self.merge_reads)
        
        elif self.assignment_method == Run.DIAMOND_TAXONOMIC_ASSIGNMENT:
            logging.info("Assigning taxonomy with diamond")
            taxonomic_assignment_time, assignments = self._assign_taxonomy_with_diamond(\
                        base_list,
                        db_search_results,
                        gpkg,
                        self.gmf)
            aln_time = 'n/a'
        else: raise Exception("Unexpected assignment method encountered: %s" % self.assignment_method)
        self._summarise(base_list, assignments, REVERSE_PIPE,
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
                         self.threads,
                         self.evalue)
        taxonomy_definition = Getaxnseq().read_taxtastic_taxonomy_and_seqinfo\
                (open(graftm_package.taxtastic_taxonomy_path()),
                 open(graftm_package.taxtastic_seqinfo_path()))
        results = {}

        # For each of the search results,
        for i, search_result in enumerate(db_search_results):
            sequence_id_to_hit = {}
            # Run diamond
            logging.debug("Running diamond on %s" % search_result.hit_fasta())
            diamond_result = runner.run(search_result.hit_fasta(),
                                        UnpackRawReads.PROTEIN_SEQUENCE_TYPE,
                                        daa_file_basename=graftm_files.diamond_assignment_output_basename(base_list[i]))
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



    def _summarise(self, base_list, trusted_placements, reverse_pipe, times, 
                  hit_read_count_list):
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
            self.s.readTax(placements, GraftMFiles(base, self.output_directory, False).read_tax_output_path(base))
            placements_list.append(placements)

        #Generate coverage table
        #logging.info('Building coverage table for %s' % base)
        #self.s.coverage_of_hmm(self.aln_hmm_file,
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
                    self.gmf = GraftMFiles(base, self.output_directory, directions[i])
                    self.hk.delete([self.gmf.for_aln_path(base),
                                    self.gmf.rev_aln_path(base),
                                    self.gmf.conv_output_rev_path(base),
                                    self.gmf.conv_output_for_path(base),
                                    self.gmf.euk_free_path(base),
                                    self.gmf.euk_contam_path(base),
                                    self.gmf.readnames_output_path(base),
                                    self.gmf.sto_output_path(base),
                                    self.gmf.orf_titles_output_path(base),
                                    self.gmf.orf_output_path(base),
                                    self.gmf.output_for_path(base),
                                    self.gmf.output_rev_path(base),
                                    self.gmf.comb_aln_fa()])
            else:
                self.gmf = GraftMFiles(base, self.output_directory, False)
                self.hk.delete([self.gmf.for_aln_path(base),
                                self.gmf.rev_aln_path(base),
                                self.gmf.conv_output_rev_path(base),
                                self.gmf.conv_output_for_path(base),
                                self.gmf.euk_free_path(base),
                                self.gmf.euk_contam_path(base),
                                self.gmf.readnames_output_path(base),
                                self.gmf.sto_output_path(base),
                                self.gmf.orf_titles_output_path(base),
                                self.gmf.orf_output_path(base),
                                self.gmf.output_for_path(base),
                                self.gmf.output_rev_path(base),
                                self.gmf.comb_aln_fa()])

        logging.info('Done, thanks for using graftM!\n')
        
        
    def make_working_directory(self, directory_path, force):
        if force:
            shutil.rmtree(directory_path, ignore_errors=True)
            os.mkdir(directory_path)
        else:
            try:
                os.mkdir(directory_path)
            except:
                logging.error('Directory %s already exists. Exiting to prevent over-writing' % directory_path)
                raise Exception('Directory %s already exists. Exiting to prevent over-writing'% directory_path)

    def delete(self, delete_list):
        for item in delete_list:
            try:
                os.remove(item)
            except:
                pass

    def get_maximum_range(self, hmm):
        '''
        If no maximum range has been specified, and if using a hmm search, a
        maximum range can be determined by using the length of the HMM

        Parameters
        ----------
        hmm : str
            path to hmm profile

        Returns
        -------
        Length to search to when linking hits on a single contig
        '''
        length=int([x for x in open(hmm) if x.startswith("LENG")][0].split()[1])
        max_length=round(length*1.5, 0)

        return max_length


