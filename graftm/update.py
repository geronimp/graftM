
import logging
import extern
import tempfile
from dendropy import Tree

from graftm.create import Create
from graftm.graftm_package import GraftMPackageVersion3, GraftMPackage
from graftm.greengenes_taxonomy import GreenGenesTaxonomy
from graftm.decorator import Decorator
from graftm.getaxnseq import Getaxnseq
from graftm.rerooter import Rerooter
from graftm.tree_decorator import TreeDecorator

class UpdateDefaultOptions:
    threads=5

class UpdatedGraftMPackage:
    '''Placeholder class for a package being built. For use internal to this file
    only.'''
    pass

class Update(Create):
    def _concatenate_file(self, file_list, output):
        '''
        Call unix "cat" to concatenate a list of files

        Parameters
        ----------
        file_list: list
            List of strings, each leading to a file. These files are the ones to
            be concatenate together. E.g.:
                ["/path/to/file1", "/path/to/file2"]
        output: str
            Path to file to which to the files in file_list will be concatenated
            into.

        '''
        to_cat = ' '.join(file_list)
        logging.debug("Concatenating files: %s" % (to_cat))
        cmd = "cat %s > %s" % (to_cat, output)
        extern.run(cmd)

    def update(self, **kwargs):
        '''
        Update an existing GraftM package with new sequences and taxonomy. If no
        taxonomy is provided, attempt to decorate the new sequences with
        pre-existing taxonomy.

        Parameters
        ----------
        input_sequence_path: str
            Path to FASTA file containing sequences to add to the update GraftM
            package
        input_taxonomy_path: str
            Taxonomy corresponding to the sequences in input_sequence_path. If None,
            then attempt to assign taxonomy by decorating the tree made out of all
            sequences.
        input_graftm_package_path: str
            Path to the directory of the GraftM package that is to be updated
        output_graftm_package_path: str
            Path to the directory to which the new GraftM package will be
            written to
        '''
        input_sequence_path = kwargs.pop('input_sequence_path')
        input_taxonomy_path = kwargs.pop('input_taxonomy_path', None)
        input_graftm_package_path = kwargs.pop('input_graftm_package_path')
        output_graftm_package_path = kwargs.pop('output_graftm_package_path')
        threads = kwargs.pop('threads', UpdateDefaultOptions.threads) #TODO: add to user options
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Reading previous GraftM package")
        old_gpkg = GraftMPackage.acquire(input_graftm_package_path)
        min_input_version = 3
        if old_gpkg.version < min_input_version:
            raise InsufficientGraftMPackageVersion(
                "GraftM below version %s cannot be updated using the update function." % min_input_version +
                " Unaligned sequences are not included in these packages, therefore no new"
                " alignment/HMM/Tree can be created")

        new_gpkg = UpdatedGraftMPackage()
        new_gpkg.output = output_graftm_package_path
        new_gpkg.name = output_graftm_package_path.replace(".gpkg", "")

        #######################################
        ### Collect all unaligned sequences ###
        logging.info("Concatenating unaligned sequence files")
        new_gpkg.unaligned_sequences = "%s_sequences.fa" % (new_gpkg.name) #TODO: replace hard-coded paths like this with tempfiles
        self._concatenate_file([old_gpkg.unaligned_sequence_database_path(),
                                input_sequence_path],
                               new_gpkg.unaligned_sequences)

        #########################################################
        ### Parse taxonomy info up front so errors come early ###
        if input_taxonomy_path:
            logging.info("Reading new taxonomy information")
            input_taxonomy = GreenGenesTaxonomy.read_file(input_taxonomy_path)
            original_taxonomy_hash = old_gpkg.taxonomy_hash()
            total_taxonomy_hash = original_taxonomy_hash.copy()
            total_taxonomy_hash.update(input_taxonomy.taxonomy)
            num_duplicate_taxonomies = len(total_taxonomy_hash) - \
                                       len(input_taxonomy.taxonomy) - \
                                       len(original_taxonomy_hash)
            logging.debug("Found %i taxonomic definitions in common between the previous and updated taxonomies" % num_duplicate_taxonomies)
            if num_duplicate_taxonomies > 0:
                logging.warn("Found %i taxonomic definitions in common between the previous and updated taxonomies. Using the updated taxonomy in each case." % num_duplicate_taxonomies)

        ###############################
        ### Re-construct alignments ###
        logging.info("Multiple sequence aligning all sequences")
        new_gpkg.aligned_sequences = "%s_mafft_alignment.fa" % (new_gpkg.name)
        self._align_sequences(new_gpkg.unaligned_sequences, new_gpkg.aligned_sequences, threads)

        ########################
        ### Re-construct HMM ###
        logging.info("Creating HMM from alignment")
        new_gpkg.hmm = "%s.hmm" % (new_gpkg.name)
        new_gpkg.hmm_alignment = "%s_hmm_alignment.fa" % (new_gpkg.name)
        self._get_hmm_from_alignment(new_gpkg.aligned_sequences, new_gpkg.hmm, new_gpkg.hmm_alignment)

        #########################
        ### Re-construct tree ###
        logging.info("Generating phylogenetic tree")
        new_gpkg.unrooted_tree = "%s.tre" % (new_gpkg.name)
        new_gpkg.unrooted_tree_log = "%s.tre.log" % (new_gpkg.name)
        new_gpkg.package_type, new_gpkg.hmm_length = self._pipe_type(old_gpkg.alignment_hmm_path())
        new_gpkg.unrooted_gpkg_tree_log, new_gpkg.unrooted_gpkg_tree = \
            self._build_tree(new_gpkg.hmm_alignment, new_gpkg.name,
                             new_gpkg.package_type, self.fasttree)

        ##############################################
        ### Re-root and decorate tree if necessary ###
        if input_taxonomy_path:
            new_gpkg.gpkg_tree_log = new_gpkg.unrooted_tree_log
            new_gpkg.gpkg_tree = new_gpkg.unrooted_gpkg_tree
        else:
            logging.info("Finding taxonomy for new sequences")
            rerooter = Rerooter()
            
            old_tree = Tree.get(path=old_gpkg.reference_package_tree_path(),
                                schema='newick')
            new_tree = Tree.get(path=new_gpkg.unrooted_gpkg_tree,
                                schema='newick')
            old_tree = rerooter.reroot(old_tree)
            new_tree = rerooter.reroot(new_tree)
            # TODO: Shouldn't call an underscore method, eventually use
            # Rerooter instead.
            rerooted_tree = rerooter.reroot_by_tree(old_tree, new_tree)
            new_gpkg.gpkg_tree = "%s_gpkg.tree" % new_gpkg.name
            td = TreeDecorator(
                rerooted_tree,
                old_gpkg.taxtastic_taxonomy_path(),
                old_gpkg.taxtastic_seqinfo_path())
            
            with tempfile.NamedTemporaryFile(suffix='tsv') as taxonomy:
                td.decorate(new_gpkg.gpkg_tree, taxonomy.name, True) 
                total_taxonomy_hash = GreenGenesTaxonomy.read_file(taxonomy.name).taxonomy

            ################################
            ### Generating tree log file ###
            logging.info("Generating phylogenetic tree log file")
            new_gpkg.gpkg_tree = "%s_gpkg.tree" % new_gpkg.name
            new_gpkg.gpkg_tree_log = "%s_gpkg.tree.log" % new_gpkg.name
            self._generate_tree_log_file(new_gpkg.unrooted_tree,
                                         new_gpkg.hmm_alignment,
                                         new_gpkg.gpkg_tree,
                                         new_gpkg.gpkg_tree_log,
                                         new_gpkg.package_type,
                                         self.fasttree)

        ################################
        ### Creating taxtastic files ###
        logging.info("Writing new taxonomy files")
        new_gpkg.tt_seqinfo = "%s_seqinfo.csv" % new_gpkg.name
        new_gpkg.tt_taxonomy = "%s_taxonomy.csv" % new_gpkg.name
        gtns = Getaxnseq()

        gtns.write_taxonomy_and_seqinfo_files(
            total_taxonomy_hash,
            new_gpkg.tt_taxonomy,
            new_gpkg.tt_seqinfo)
        
        ######################
        ### Compile refpkg ###
        logging.info("Compiling pplacer refpkg")
        new_gpkg.refpkg = "%s.refpkg" % (new_gpkg.name)
        refpkg = self._taxit_create(new_gpkg.name,
                                    new_gpkg.hmm_alignment,
                                    new_gpkg.gpkg_tree,
                                    new_gpkg.gpkg_tree_log,
                                    new_gpkg.tt_taxonomy,
                                    new_gpkg.tt_seqinfo,
                                    new_gpkg.refpkg,
                                    True)

        #####################################
        ### Re-construct diamond database ###
        logging.info("Recreating DIAMOND DB")
        new_gpkg.diamond_database = "%s.dmnd" % (new_gpkg.name)
        self._create_dmnd_database(new_gpkg.unaligned_sequences, new_gpkg.name)

        ####################
        ### Compile gpkg ###
        logging.info("Compiling GraftM package")
        new_gpkg.name = "%s.gpkg" % new_gpkg.name
        GraftMPackageVersion3.compile(new_gpkg.name, new_gpkg.refpkg,
                                      new_gpkg.hmm, new_gpkg.diamond_database,
                                      self._define_range(new_gpkg.unaligned_sequences),
                                      new_gpkg.unaligned_sequences,
                                      search_hmm_files=old_gpkg.search_hmm_paths())

        ###################
        ### Test it out ###
        logging.info("Testing newly updated GraftM package works")
        self._test_package(new_gpkg.name)

        logging.info("Finished")
    
