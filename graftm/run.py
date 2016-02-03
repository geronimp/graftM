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

import os
import logging
import inspect

from graftm.create import Create
from graftm.graft import Graft

from graftm.pplacer import Pplacer
from graftm.hmmer import Hmmer
from graftm.graftm_package import GraftMPackage
from graftm.bootstrapper import Bootstrapper
from graftm.decorator import Decorator
from graftm.external_program_suite import ExternalProgramSuite

class UnrecognisedSuffixError(Exception): 
    pass

class Run:
    
    _MIN_VERBOSITY_FOR_ART = 3 # with 2 then, only errors are printed

    MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES = 95
    MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES = 30
    
    _DNA = "DNA"
    _RNA = "RNA"
    _AMINO = "amino"

    def __init__(self, args):
        self.args = args

    def _check_file_existence(self, files):
        '''Iterate through files and exit(1) if any do not pass
        os.path.isfile'''
        for f in files:
            if not os.path.isfile(f):
                logging.error("The file '%s' does not appear to exist, stopping" % f)
                exit(1)
            
    def _setpipe(self, hmm_list):

        hmm_tc = False
        
        type_list = []
        tc_list = []
        
        for hmm in hmm_list:
            tc = False
            for line in open(hmm):
                if line.startswith('ALPH'):
                    code = line.split()[1]
                    if code == self._DNA or code == self._RNA:
                        type = Graft.PIPELINE_NT
                    elif code == self._AMINO:
                        type = Graft.PIPELINE_AA

                elif line.startswith('TC'):
                    tc = True
                    
            type_list.append(type)
            if tc:
                tc_list.append(tc)
        
        nr_type_list = set(type_list)
        if len(nr_type_list) > 1:
            logging.error("Search HMMs are not all %s type" % nr_type_list[0])
        else:
            type = list(nr_type_list)[0]
        
        if tc_list:
            if len(tc_list) == len(hmm_list):
                tc = True
            else:
                tc = False
                logging.warning("Not all HMMs have trusted cutoffs (tc). \
choosing the default evalue as a search cutoff.")
        
        return type, tc

                
    def main(self):
################################################################################
#                                  GRAFT                                       #
################################################################################
        if self.args.subparser_name == 'graft':
            
            commands = ExternalProgramSuite(['orfm', 'nhmmer', 'hmmsearch', 
                                                 'fxtract', 'pplacer', 
                                                 'seqmagick', 'ktImportText', 
                                                 'diamond'])
                        
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
            
            if self.args.graftm_package:
                if not hasattr(self.args, 'search_hmm_files'):
                    setattr(self.args, 'search_hmm_files', [])
                    
                    for gpkg in self.args.graftm_package:
                        if not os.path.isdir(gpkg):
                            raise Exception("%s does not exist. Are you sure you provided the correct path?" % gpkg)
                        else:
                            gpkg = GraftMPackage.acquire(gpkg)                            
                            for hmm in gpkg.search_hmm_paths():
                                self.args.search_hmm_files.append(hmm)
    
            elif hasattr(self.args, 'search_diamond_files'):
                if self.args.search_method == Hmmer.DIAMOND_SEARCH_METHOD:
                    if hasattr(self.args, 'aln_hmm_file'):
                        pass
                    else:
                        raise Exception("aln_hmm_file not specified")
                else:
                    raise Exception("Specified HMM databases when not using the diamond search pipeline. Using: %s" % (self.args.search_method))
    
            elif hasattr(self.args, 'search_hmm_files'):
                if self.args.search_method == Hmmer.HMM_SEARCH_METHOD:
                    if not hasattr(self.args, 'aln_hmm_file'):
                        if len(self.args.search_hmm_files) == 1:
                            if not self.args.search_only:
                                setattr(self.args, 'aln_hmm_file', self.args.search_hmm_files[0])
                        else:
                            raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")
    
                else:
                    raise Exception("Specified HMM search_hmm_files when not using the hmmsearch pipeline. Using: %s" % (self.args.search_method))
    
            elif hasattr(self.args, 'search_hmm_list_file'):
                if self.args.search_method == Hmmer.HMM_SEARCH_METHOD:
                    setattr(self.args, 'search_hmm_files', [x.rstrip() for x in open(self.args.search_hmm_list_file).readlines()])
                    if not hasattr(self.args, 'aln_hmm_file'):
                        if not self.args.search_only:
                            raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")
                else:
                    raise Exception("Specified HMM search_hmm_files when not using the hmmsearch pipeline. Using: %s" % (self.args.search_method))
    
            else:
                raise Exception('No graftM package or HMM specified: Do not know what to search with.')
            
            if hasattr(self.args, 'euk_hmm_file'):
                pass
            elif not hasattr(self.args, 'euk_hmm_file'):
                setattr(self.args, 'euk_hmm_file', os.path.join(os.path.dirname(inspect.stack()[-1][1]),'..','share', '18S.hmm'))
            else:
                raise Exception('Programming Error: setting the euk HMM')
            if self.args.euk_check:
                self.args.search_hmm_files.append(self.args.euk_hmm_file)

            
            # Check that the placement cutoff is between 0.5 and 1
            if float(self.args.placements_cutoff) < float(0.5) or float(self.args.placements_cutoff) > float(1.0):
                logging.info('Please specify a confidence level (-d) between 0.5 and 1.0! Found: %s' % self.args.placements_cutoff)
                exit(1)

            self._check_file_existence(self.args.forward)
            if self.args.reverse:
                self._check_file_existence(self.args.reverse)
            
            sequence_pair_list = []
            if self.args.reverse:
                if len(args.forward) != len(self.args.reverse):
                    logging.error('Confusing input. There appears to be different numbers of forward and reverse files specified')
                for i, forward_file in enumerate(self.args.forward):
                    sequence_pair_list.append([forward_file, args.reverse[i]])
            else:
                sequence_pair_list = [[f] for f in self.args.forward]

            if hasattr(self.args, 'reference_package'):
                self.p = Pplacer(self.args.reference_package)   
            
            # If merge reads is specified, check that there are reverse reads to merge with
            if self.args.merge_reads and not hasattr(self.args, 'reverse'):
                logging.error("--merge_reads requires --reverse to be specified")
                exit(1)

            if self.args.assignment_method == Graft.DIAMOND_TAXONOMIC_ASSIGNMENT:
                if self.args.reverse:
                    logging.warn("--reverse reads specified with --assignment_method diamond. Reverse reads will be ignored.")
                    self.reverse = None
                
            hmm_type, hmm_tc = self._setpipe(self.args.search_hmm_files)
            logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))
            
            setattr(self.args, 'type', hmm_type)
            if hmm_tc:
                setattr(self.args, 'evalue', '--cut_tc')
    
            if self.args.filter_minimum is not None:
                filter_minimum = self.filter_minimum
            else:
                if self.args.type == Graft.PIPELINE_NT:
                    self.args.filter_minimum = self.MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES
                elif self.args.type == Graft.PIPELINE_AA:
                    self.args.filter_minimum = self.MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES
            
            graftm_packages = {}
            search_dict = {}
            if self.args.graftm_package:
                
                for gpkg in self.args.graftm_package:
                    graftm_packages[gpkg] = GraftMPackage.acquire(gpkg)
            elif self.args.search_method:
                pass 
            else:
                raise Exception("No GraftM package, search hmms or diamond \
databases provided. GraftM can only run with at least one of those.")    
                    
            Graft(assignment_method = self.args.assignment_method,
                  bootstrap_contigs = self.args.bootstrap_contigs, 
                  euk_check = self.args.euk_check,
                  evalue = self.args.evalue, 
                  filter_minimum = self.args.filter_minimum,
                  force = self.args.force, 
                  graftm_packages = graftm_packages,
                  maximum_range = self.args.maximum_range,
                  merge_reads = self.args.merge_reads,
                  min_orf_length = self.args.min_orf_length, 
                  no_clustering = self.args.no_clustering,
                  output_directory = self.args.output_directory, 
                  resolve_placements = self.args.resolve_placements, 
                  restrict_read_length = self.args.restrict_read_length, 
                  reverse = self.args.reverse,
                  search_and_align_only = self.args.search_and_align_only, 
                  search_diamond_file = self.args.search_diamond_file, 
                  search_hmm_files = self.args.search_hmm_files, 
                  search_method = self.args.search_method,
                  search_only = self.args.search_only, 
                  sequence_pair_list = sequence_pair_list,
                  threads = self.args.threads,
                  type = self.args.type).graft()

################################################################################
#                                 CREATE                                       #
################################################################################
        elif self.args.subparser_name == 'create':

            commands = ExternalProgramSuite(['taxit', 'FastTreeMP', 
                                                 'seqmagick', 'hmmalign', 
                                                 'mafft'])
            self.create = Create(commands)
            
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
            if self.args.dereplication_level not in range(0,8):
                logging.error("Invalid dereplication level selected! please enter an integer from 0 - 7")
                exit(1)
            
            if self.args.graftm_package:
                if self.args.taxonomy and self.args.sequences:
                    logging.info("GraftM package %s specified to update with sequences in %s" % (self.args.graftm_package, self.args.sequences))
                    if not self.args.output:
                        if self.args.graftm_package.endswith(".gpkg"):
                            self.args.output = self.args.graftm_package.replace(".gpkg", "-updated.gpkg")
                        else:
                            raise UnrecognisedSuffixError("Unrecognised suffix on GraftM package %s. Please provide a graftM package with the correct suffix (.gpkg)" % (self.args.graftm_package))
                    self.create.update(self.args.sequences, self.args.taxonomy, self.args.graftm_package,
                                       self.args.output)

                else:
                    if self.args.taxonomy:
                        logging.error("--sequences must be specified to update a GraftM package")
                        exit(1)
                    elif self.args.sequences:
                        logging.error("--taxonomy must be specified to update a GraftM package")
                        exit(1)
            else:
                if not self.args.sequences:
                    if not self.args.alignment and not self.args.rerooted_annotated_tree \
                                               and not self.args.rerooted_tree:
                        logging.error("Some sort of sequence data must be provided to run graftM create")
                        exit(1)
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
                
    
                self.create.main(
                              dereplication_level = self.args.dereplication_level,
                              sequences = self.args.sequences,
                              alignment = self.args.alignment,
                              taxonomy = self.args.taxonomy,
                              rerooted_tree = self.args.rerooted_tree,
                              tree_log = self.args.tree_log,
                              prefix = self.args.output,
                              rerooted_annotated_tree = self.args.rerooted_annotated_tree,
                              min_aligned_percent = float(self.args.min_aligned_percent)/100,
                              taxtastic_taxonomy = self.args.taxtastic_taxonomy,
                              taxtastic_seqinfo = self.args.taxtastic_seqinfo,
                              hmm = self.args.hmm,
                              search_hmm_files = self.args.search_hmm_files,
                              force = self.args.force,
                              graftm_package = self.args.graftm_package,
                              threads = self.args.threads                              
                              )

################################################################################
#                              BOOTSTRAP                                       #
################################################################################
        elif self.args.subparser_name == 'bootstrap':
            if not self.args.graftm_package and not self.args.search_hmm_files:
                logging.error("bootstrap mode requires either --graftm_package or --search_hmm_files")
                exit(1)

            if self.args.graftm_package:
                pkg = GraftMPackage.acquire(self.args.graftm_package)
            else:
                pkg = None

            strapper = Bootstrapper(search_hmm_files = self.args.search_hmm_files,
                maximum_range = self.args.maximum_range,
                threads = self.args.threads,
                evalue = self.args.evalue,
                min_orf_length = self.args.min_orf_length,
                graftm_package = pkg)
            strapper.generate_bootstrap_database_from_contigs(self.args.contigs,
                                                              self.args.output_hmm,
                                                              search_method=Hmmer.HMM_SEARCH_METHOD)

################################################################################
#                                   TREE                                       #
################################################################################
        elif self.args.subparser_name == 'tree':
            if self.args.graftm_package:
                # shim in the paths from the graftm package, not overwriting
                # any of the provided paths.
                gpkg = GraftMPackage.acquire(self.args.graftm_package)
                if not self.args.rooted_tree: self.args.rooted_tree = gpkg.reference_package_tree_path()
                if not self.args.input_greengenes_taxonomy:
                    if not self.args.input_taxtastic_seqinfo:
                        self.args.input_taxtastic_seqinfo = gpkg.taxtastic_seqinfo_path()
                    if not self.args.input_taxtastic_taxonomy:
                        self.args.input_taxtastic_taxonomy = gpkg.taxtastic_taxonomy_path()
                        
            if self.args.rooted_tree:
                if self.args.unrooted_tree:
                    logging.error("Both a rooted tree and an un-rooted tree were provided, so it's unclear what you are asking GraftM to do. \
If you're unsure see graftM tree -h")
                    exit(1)
                elif self.args.reference_tree:
                    logging.error("Both a rooted tree and reference tree were provided, so it's unclear what you are asking GraftM to do. \
If you're unsure see graftM tree -h")
                    exit(1)
                    
                if not self.args.decorate:
                    logging.error("It seems a rooted tree has been provided, but --decorate has not been specified so it is unclear what you are asking graftM to do.")
                    exit(1)
                
                dec = Decorator(tree_path = self.args.rooted_tree)
            
            elif self.args.unrooted_tree and self.args.reference_tree:
                logging.debug("Using provided reference tree %s to reroot %s" % (self.args.reference_tree, 
                                                                                 self.args.unrooted_tree))
                dec = Decorator(reference_tree_path = self.args.reference_tree,
                                tree_path = self.args.unrooted_tree)
            else:
                logging.error("Some tree(s) must be provided, either a rooted tree or both an unrooted tree and a reference tree")
                exit(1)
                
            if self.args.output_taxonomy is None and self.args.output_tree is None:
                logging.error("Either an output tree or taxonomy must be provided")
                exit(1)
            if self.args.input_greengenes_taxonomy:
                if self.args.input_taxtastic_seqinfo or self.args.input_taxtastic_taxonomy:
                    logging.error("Both taxtastic and greengenes taxonomy were provided, so its unclear what taxonomy you want graftM to decorate with")
                    exit(1)
                logging.debug("Using input GreenGenes style taxonomy file")
                dec.main(self.args.input_greengenes_taxonomy, 
                         self.args.output_tree, self.args.output_taxonomy, 
                         self.args.no_unique_tax, self.args.decorate, None) 
            elif self.args.input_taxtastic_seqinfo and self.args.input_taxtastic_taxonomy:
                logging.debug("Using input taxtastic style taxonomy/seqinfo")
                dec.main(self.args.input_taxtastic_taxonomy, self.args.output_tree, 
                         self.args.output_taxonomy, self.args.no_unique_tax, 
                         self.args.decorate, self.args.input_taxtastic_seqinfo)
            else:
                logging.error("Either a taxtastic taxonomy or seqinfo file was provided. GraftM cannot continue without both.")
                exit(1)
        else:
            raise Exception("Unexpected subparser name %s" % self.args.subparser_name)
            
