from graftm.create import Create

class UpdatedGraftMPackage:
    pass

class Update(Create):
    def update(self, input_sequence_path, input_taxonomy_path,
               input_graftm_package_path, output_graftm_package_path):
        '''
        Update an existing GraftM pacakge with new sequences and taxonomy. If no
        taxonomy is provided, attempt to decorate the new sequences with
        pre-existing taxonomy.

        Parameters
        ----------
        input_sequence_path: str
            Path to FASTA file containing sequences to add to the update GraftM
            package
        input_taxonomy_path: str
            Taxonomy corresponding to the sequences within input_sequence_path
        input_graftm_package_path: str
            Path to the directory of the GraftM package that is to be updated
        output_graftm_package_path: str
            Path to the directory to which the new GraftM package will be
            written to
        '''
        gtns = Getaxnseq()
        old_gpkg = GraftMPackage.acquire(input_graftm_package_path)
        new_gpkg = UpdatedGraftMPackage

        if old_gpkg.version < 3:
            raise InsufficientGraftMPackageVersion("""
GraftM below version 3 cannot be updated using the update function. Unaligned
sequences are not included in these packages, therefore no new
alignment/HMM/Tree can be created""")

        new_gpkg.output = output_graftm_package_path
        new_gpkg.name = output_graftm_package_path.replace(".gpkg", "")


        #####################################
        ### Re-construct diamond database ###
        new_gpkg.unaligned_sequences = "%s_sequences.fa" % (new_gpkg.name)
        new_gpkg.diamond_database = "%s.dmnd" % (new_gpkg.name)

        self._concatenate_file([old_gpkg.unaligned_sequence_database_path(),
                                input_sequence_path],
                               new_gpkg.unaligned_sequences)
        self._create_dmnd_database(new_gpkg.unaligned_sequences, new_gpkg.name)


        ###############################
        ### Re-construct alignments ###
        new_gpkg.aligned_sequences = "%s_mafft_alignment.fa" % (new_gpkg.name)
        self._align_sequences(new_gpkg.unaligned_sequences, new_gpkg.aligned_sequences)


        ########################
        ### Re-construct HMM ###
        new_gpkg.hmm = "%s.hmm" % (new_gpkg.name)
        new_gpkg.hmm_alignment = "%s_hmm_alignment.fa" % (new_gpkg.name)
        self._get_hmm_from_alignment(new_gpkg.aligned_sequences, new_gpkg.hmm, new_gpkg.hmm_alignment)

        ########################
        ### Parse taxonomy info before tree making so errors come faster ###
        if input_taxonomy_path:
            old_tax = Getaxnseq().read_taxtastic_taxonomy_and_seqinfo(open(old_gpkg.taxtastic_taxonomy_path()),
                                                                      open(old_gpkg.taxtastic_seqinfo_path()))

        #########################
        ### Re-construct tree ###
        new_gpkg.unrooted_tree = "%s.tre" % (new_gpkg.name)
        new_gpkg.unrooted_tree_log = "%s.tre.log" % (new_gpkg.name)
        new_gpkg.ptype, new_gpkg.hmm_length = self._pipe_type(old_gpkg.alignment_hmm_path())
        self._build_tree(new_gpkg.hmm_alignment, new_gpkg.name, new_gpkg.ptype, self.fasttree)


        #################################
        ### Re-root and decorate tree ###
        new_gpkg.rooted_tree = "%s_rooted.tre" % (new_gpkg.name)
        new_gpkg.decorate_tax = "%s_decorate_tax.tsv" % (new_gpkg.name)
        new_gpkg.gg_taxonomy = "%s_greengenes_taxonomy.tsv" % (new_gpkg.name)
        reference_tree = old_gpkg.reference_package_tree_path()

        decorator = Decorator(reference_tree_path = reference_tree,
                      tree_path = new_gpkg.unrooted_tree)
        if input_taxonomy_path:
            with tempfile.NamedTemporaryFile(suffix='.tsv') as old_taxonomy_path:
                with open(new_gpkg.gg_taxonomy, "w") as decoration_taxonomy:
                    for line in ["%s\t%s\n" % \
                                 (id,
                                  '; '.join(tax + self._PREFIX_LIST[len(tax):])) \
                                 for id, tax in old_tax.iteritems()]:
                        old_taxonomy_path.write(line)
                    old_taxonomy_path.flush()

                    self._concatenate_file( [old_taxonomy_path.name,
                                             input_taxonomy_path],
                                             decoration_taxonomy.name )

                    decorator.main(decoration_taxonomy.name,
                                     new_gpkg.rooted_tree,
                                     new_gpkg.decorate_tax,
                                     False,
                                     False)
        else:
            decorator.main(old_gpkg.taxtastic_taxonomy_path,
                             new_gpkg.rooted_tree,
                             new_gpkg.decorate_tax,
                             False,
                             False,
                             old_gpkg.taxtastic_seqinfo_path)

        ################################
        ### Generating tree log file ###
        new_gpkg.gpkg_tree = "%s_gpkg.tree" % new_gpkg.name
        new_gpkg.gpkg_tree_log = "%s_gpkg.tree.log" % new_gpkg.name

        self._generate_tree_log_file(new_gpkg.rooted_tree,
                                     new_gpkg.hmm_alignment,
                                     new_gpkg.gpkg_tree,
                                     new_gpkg.gpkg_tree_log,
                                     new_gpkg.ptype,
                                     self.fasttree)

        ################################
        ### Creating taxtastic files ###
        new_gpkg.tt_seqinfo = "%s_seqinfo.csv" % new_gpkg.name
        new_gpkg.tt_taxonomy = "%s_taxonomy.csv" % new_gpkg.name
        gtns.write_taxonomy_and_seqinfo_files(GreenGenesTaxonomy.read_file(new_gpkg.gg_taxonomy).taxonomy,
                                              new_gpkg.tt_taxonomy,
                                              new_gpkg.tt_seqinfo)
        ######################
        ### Compile refpkg ###
        new_gpkg.refpkg = "%s.refpkg" % (new_gpkg.name)
        refpkg = self._taxit_create(new_gpkg.name,
                                    new_gpkg.hmm_alignment,
                                    new_gpkg.gpkg_tree,
                                    new_gpkg.gpkg_tree_log,
                                    new_gpkg.tt_taxonomy,
                                    new_gpkg.tt_seqinfo,
                                    new_gpkg.refpkg,
                                    True)

        ####################
        ### Compile gpkg ###
        new_gpkg.name = "%s.gpkg" % new_gpkg.name
        GraftMPackageVersion3.compile(new_gpkg.name, new_gpkg.refpkg,
                                      new_gpkg.hmm, new_gpkg.diamond_database,
                                      self._define_range(new_gpkg.unaligned_sequences),
                                      new_gpkg.unaligned_sequences,
                                      search_hmm_files=old_gpkg.search_hmm_paths())
    
