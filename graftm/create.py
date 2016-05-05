#!/usr/bin/env python

import os
import shutil
import tempfile
import logging
import subprocess32
import extern
import tempdir
import json
import signal
import re

from Bio import SeqIO

from graftm.hmmer import Hmmer
from graftm.tree_cleaner import TreeCleaner
from graftm.taxonomy_extractor import TaxonomyExtractor
from graftm.getaxnseq import Getaxnseq
from graftm.deduplicator import Deduplicator
from skbio.tree import TreeNode
from graftm.sequence_io import SequenceIO
from graftm.graftm_package import GraftMPackageVersion3, GraftMPackage
from graftm.decorator import Decorator
from graftm.greengenes_taxonomy import GreenGenesTaxonomy

class InsufficientGraftMPackageVersion(Exception):
    pass

class UpdatedGraftMPackage:
    pass

class Create:
    _PROTEIN_PACKAGE_TYPE = 'protein_package_type'
    _NUCLEOTIDE_PACKAGE_TYPE = 'nucleotide_package_type'
    _PREFIX_LIST = ["d__",
                    "p__",
                    "c__",
                    "o__",
                    "f__",
                    "g__",
                    "s__"]


    def __init__(self, commands):
        '''
        Class for managing the creation of GraftM packages.
        
        Parameters
        ----------
        commands: ExternalProgramSuite object
            Contains each of the command line commands as attributes. e.g.
            fasttree:
                commands.fasttree = "FastTreeMP" or "fasttree" depending on
                what is installed.
        '''
        self.h=Hmmer(None, None)
        self.the_trash=[]
        self.fasttree = commands.fasttree

    def _parse_contents(self, contents_file_path):
        '''
        Parse the contents .json file and return the dictionary

        Parameters
        ----------
        contents_file_path: str
            Path to file containing .json file to parse.

        Returns
        -------
        contents_dict: dict
            parsed .json file
        '''
        logging.debug("Parsing %s" % (contents_file_path))
        contents_dict = json.load(open(contents_file_path))
        return contents_dict

    def _check_reads_hit(self, alignment_io, min_aligned_fraction):
        '''Given an alignment return a list of sequence names that are less
        than the min_aligned_fraction'''
        to_return = []
        alignment_length = None
        for s in SeqIO.parse(alignment_io, "fasta"):
            if not alignment_length:
                alignment_length = len(s.seq)
                min_length = int(min_aligned_fraction * alignment_length)
                logging.debug("Determined min number of aligned bases to be %s" % min_length)
            elif len(s.seq) != alignment_length:
                raise Exception("Alignment file appears to not be of uniform length")

            num_unaligned = s.seq.count('-')
            num_aligned = alignment_length-num_unaligned
            logging.debug("Sequence %s has %d aligned positions" % (s.name, alignment_length-num_unaligned))
            if num_aligned <= min_length:
                to_return.append(s.name)
        return to_return

    def _align_sequences(self, input_sequences_path, output_alignment_path,
                         threads):
        '''Align sequences into alignment_file

        Parameters
        ----------
        input_sequences_path: str
            path to input sequences in fasta format
        output_alignment_path: str
            path to output alignment path
        threads: str
            number of threads to use
        Returns
        -------
        Nothing
        '''
        logging.debug("Aligning sequences using mafft")
        cmd = "mafft --anysymbol --thread %s --auto '%s' > %s" % (threads,
                                                                  input_sequences_path,
                                                                  output_alignment_path)
        extern.run(cmd)

    def _get_hmm_from_alignment(self, alignment, hmm_filename, output_alignment_filename):
        '''Return a HMM file and alignment of sequences to that HMM

        Parameters
        ----------
        alignment: str
            path to aligned proteins
        hmm_filename: str
            write the hmm to this file path
        output_alignment_filename: str
            write the output alignment to this file path

        Returns
        -------
        Nothing
        '''
        logging.debug("Building HMM from alignment")

        sto = tempfile.NamedTemporaryFile(suffix='.sto',prefix='graftm')
        tempaln = tempfile.NamedTemporaryFile(suffix='.fasta',prefix='graftm')
        cmd = "hmmbuild -O %s '%s' '%s'" % (sto.name,
                                              hmm_filename,
                                              alignment)
        out = extern.run(cmd)
        logging.debug("Got STDOUT from hmmbuild: %s" % out)

        cmd = "seqmagick convert --input-format stockholm %s %s" % (sto.name,
                                                              tempaln.name)
        out = extern.run(cmd)
        logging.debug("Got STDOUT from seqmagick: %s" % out)
        ptype, _ = self._pipe_type(hmm_filename)
        Hmmer(hmm_filename).alignment_correcter([tempaln.name],
                                                output_alignment_filename)
        return ptype

    def _align_sequences_to_hmm(self, hmm_file, sequences_file, output_alignment_file):
        '''Align sequences to an HMM, and return a path to an alignment of
        these proteins after cleanup so that they can be used for tree-making

        Parameters
        ----------
        sequences_file: str
            path to file of unaligned protein sequences
        hmm_file: str
            path to hmm file
        output_alignment_file: str
            write alignment to this file

        Returns
        -------
        nothing
        '''
        hmmer = Hmmer(hmm_file)
        tempalign = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.fasta')
        hmmer.hmmalign_sequences(hmm_file, sequences_file, tempalign.name)
        hmmer.alignment_correcter([tempalign.name], output_alignment_file)

    def _pipe_type(self, hmm):
        logging.debug("Setting pipeline type.")
        hmm_type=[x.split() for x in open(hmm).readlines() if x.startswith('ALPH') or x.startswith('LENG')]
        for item in hmm_type:
            if item[0]=='ALPH':
                if item[1]=='amino':
                    ptype=Create._PROTEIN_PACKAGE_TYPE
                elif item[1]=='DNA' or 'RNA':
                    ptype=Create._NUCLEOTIDE_PACKAGE_TYPE
                else:
                    raise Exception("Unfamiliar HMM type: %s" % (item[1]))
            elif item[0]=='LENG':
                leng=item[1]
            else:
                raise Exception("Programming Error: Misread HMM file")
        logging.debug("Set pipeline type as: %s " % ptype)
        logging.debug("Found alignment length as: %s" % leng)
        return ptype, leng

    def _check_aln_length(self, alignment):
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)

    def _build_tree(self, alignment, base, ptype, fasttree):
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == Create._NUCLEOTIDE_PACKAGE_TYPE: # If it's a nucleotide sequence
            cmd = "%s -quiet -gtr -nt -log %s -out %s %s" % (fasttree,
                                                             log_file, 
                                                             tre_file, 
                                                             alignment)
            extern.run(cmd)
        else: # Or if its an amino acid sequence
            cmd = "%s -quiet -log %s -out %s %s" % (fasttree,
                                                    log_file, 
                                                    tre_file, 
                                                    alignment)
            extern.run(cmd)

        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def _taxit_create(self, base, aln_file, tre, tre_log, tax, seq,
                      refpkg, no_reroot):
        cmd = "taxit create -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s"\
            % (aln_file, refpkg, tre, tre_log, base, tax, seq)

        if no_reroot:
            cmd += ' --no-reroot'
            logging.debug("Calling command assuming pre-rerooting: %s" % cmd)
            extern.run(cmd)
        else:
            logging.debug("Calling command: %s" % cmd)
            logging.info("Attempting to run taxit create with rerooting capabilities")
            def exit_gracefully(base):
                local_base = os.path.basename(base)
                nontemp_tax_file = 'graftm_create_taxonomy.%s.csv' % local_base
                shutil.copy(tax, nontemp_tax_file)
                nontemp_seqinfo_file = 'graftm_create_seqinfo.%s.csv' % local_base
                shutil.copy(seq, nontemp_seqinfo_file)
                nontemp_aln_file = 'graftm_create_alignment.%s.faa' % local_base
                shutil.copy(aln_file, nontemp_aln_file)
                nontemp_tre_file = 'graftm_create_tree.%s.tree' % local_base
                shutil.copy(tre, nontemp_tre_file)
                logging.error('''
taxit create failed to run in a small amount of time suggesting that
rerooting was unsuccessful. Unfortunately this tree will need to be rerooted
manually yourself using a tree editor such as ARB or FigTree.
Once you have a rerooted newick format tree, rerun graftm create
specifying the new tree with --rerooted_tree. The tree file to be rerooted is \'%s\'''' % nontemp_tre_file)

                logging.error('''
When rerunning, please use the following flags for the command line to account
for the fact that some sequences may have been removed during the deduplication
process.

graftM create --taxtastic_taxonomy %s --taxtastic_seqinfo %s --alignment %s  --rerooted_tree <REROOTED_TREE>

(plus other relevant arguments).
''' % (
       #nontemp_align_hmm_file,
       #nontemp_search_hmm_file,
       nontemp_tax_file,
       nontemp_seqinfo_file,
       nontemp_aln_file,
       #nontemp_seq_file
       )
                              )
                exit(2)
                
                
            try:
                process = subprocess32.Popen(['bash','-c',cmd],
                                           stdout=subprocess32.PIPE,
                                           stderr=subprocess32.PIPE,
                                           preexec_fn=os.setsid)
                stdout, stderr = process.communicate(timeout=60)
                logging.debug("taxit create STDOUT was: %s" % stdout)
                logging.debug("taxit create STDERR was: %s" % stderr)

                if process.returncode != 0 or re.match(r'^Error: ',stdout):
                    raise extern.ExternCalledProcessError(cmd,
                                                          process.returncode,
                                                          stderr,
                                                          stdout)
            except subprocess32.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                exit_gracefully(base)
            except extern.ExternCalledProcessError:
                exit_gracefully(base)

        return refpkg

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

    def _cleanup(self, the_trashcan):
        for every_piece_of_junk in the_trashcan:
            if os.path.isdir(every_piece_of_junk):
                shutil.rmtree(every_piece_of_junk)
            else:
                os.remove(every_piece_of_junk)

    def _define_range(self, sequences):
        '''
        define_range - define the maximum range within which two hits in a db
        search can be linked. This is defined as 1.5X the average length of all
        reads in the database.

        Parameters
        ----------
        sequences : str
            A path to the sequences in FASTA format. This fasta file is assumed
            to be in the correct format. i.e. headers start with '>'

        Returns
        -------
        max_range : int
            As described above, 1.5X the size of the average length of genes
            within the database
        '''
        sequence_count = 0
        total_sequence = 0

        for record in SeqIO.parse(open(sequences), 'fasta'): # For every sequence in the output
            total_sequence+=1 # increment the total sequence count
            sequence_count+=len(record.seq) # add the length of that sequence to the total string count

        # Get the average, multiply by 1.5, and return
        max_range = (sequence_count/total_sequence)*1.5
        return max_range

    def _generate_tree_log_file(self, tree, alignment, output_tree_file_path,
                               output_log_file_path, residue_type, fasttree):
        '''Generate the FastTree log file given a tree and the alignment that
        made that tree

        Returns
        -------
        Nothing. The log file as parameter is written as the log file.
        '''
        if residue_type==Create._NUCLEOTIDE_PACKAGE_TYPE:
            cmd = "%s -quiet -gtr -nt -nome -mllen -intree '%s' -log %s -out %s %s" %\
                                       (fasttree, tree, output_log_file_path,
                                        output_tree_file_path, alignment)
        elif residue_type==Create._PROTEIN_PACKAGE_TYPE:
            cmd = "%s -quiet -nome -mllen -intree '%s' -log %s -out %s %s" %\
                                       (fasttree, tree, output_log_file_path,
                                        output_tree_file_path, alignment)
        extern.run(cmd)

    def _remove_sequences_from_alignment(self, sequence_names, input_alignment_file, output_alignment_file):
        '''Remove sequences from the alignment file that have names in
        sequence_names

        Parameters
        ----------
        sequence_names: list of str
            names of sequences to remove
        input_alignment_file: str
            path to alignment file to remove from
        output_alignment_file: str
            path to alignment file to write to

        Returns
        -------
        int: number of sequences written to file'''
        nameset = set(sequence_names)
        num_written = 0
        with open(output_alignment_file, 'w') as f:
            for s in SeqIO.parse(open(input_alignment_file), "fasta"):
                if s.name not in nameset:
                    SeqIO.write(s, f, "fasta")
                    num_written += 1
        return num_written

    def _create_search_hmm(self, sequences, taxonomy_definition,
                           search_hmm, dereplication_level, threads):
        logging.debug("Constructing search HMM")
        with tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.faa') as temporary_alignment_fh:
            temporary_alignment = temporary_alignment_fh.name
            with tempfile.NamedTemporaryFile(prefix='graftm', suffix='.faa') as temporary_dereplicated_fh:
                temporary_dereplicated = temporary_dereplicated_fh.name
                if dereplication_level > 0:
                    dereplication_index = dereplication_level - 1

                    #####################
                    ### Dereplication ###
                    logging.debug("Dereplicating sequences at rank %i in the taxonomy file provided" % (dereplication_level))
                    dereplicated_sequence_ids = []
                    seen_taxons = set()
                    prefix_re = re.compile(r'^[a-zA-Z]__$')
                    for sequence_id, taxonomy in taxonomy_definition.iteritems():
                        try:
                            taxon = taxonomy[dereplication_index]
                        except IndexError:
                            taxon = None
                        if taxon is None or \
                           taxon == "" or \
                           prefix_re.match(taxon):
                            pass
                        
                        elif taxon in seen_taxons:
                            logging.debug("Sequence %s redundant at %i rank in the taxonomy file level: %s" % (sequence_id, dereplication_level, taxon) )
                            dereplicated_sequence_ids.append(sequence_id)
                        else:
                            seen_taxons.add(taxon)
                    if len(taxonomy_definition) - len(dereplicated_sequence_ids) < 2:
                        raise Exception("Insufficient sequences available after dereplication to create a HMM. One solution might be to modify the dereplication level to do less dereplication.")
                    logging.info("Removing %i sequences from the search HMM that are redundant at the %i rank in the taxonomy file" \
                                            % (len(dereplicated_sequence_ids),
                                               dereplication_level))
                    self._remove_sequences_from_alignment(dereplicated_sequence_ids, sequences, temporary_dereplicated)
                    self._align_sequences(temporary_dereplicated, temporary_alignment, threads)
                else:
                    logging.debug("Skipping dereplication step and using all sequences to build search HMM")
                    self._align_sequences(sequences, temporary_alignment, threads)

                self._get_hmm_from_alignment(temporary_alignment, search_hmm, temporary_alignment)


    def _create_dmnd_database(self, unaligned_sequences_path, daa_output):
        '''
        Build a diamond database using diamond makedb

        Parameters
        ----------
        unaligned_sequences_path: str
            path to a FASTA file containing unaligned sequences
        daa_output: str
            Name of output database.
        '''
        logging.debug("Building diamond database")

        cmd = "diamond makedb --in '%s' -d '%s'" % (unaligned_sequences_path, daa_output)
        extern.run(cmd)

    def _align_and_create_hmm(self, sequences, alignment, user_hmm,
                              output_align_hmm, output_alignment, threads):

        # align sequences to HMM (and potentially build hmm from alignment)
        if user_hmm:
            output_hmm = user_hmm
            ptype, _ = self._pipe_type(output_hmm)
            self._align_sequences_to_hmm(output_hmm, sequences, output_alignment)
        else:
            if not alignment:
                self._align_sequences(sequences, output_alignment, threads)
            else:
                output_alignment = alignment
            ptype = self._get_hmm_from_alignment(output_alignment,
                                                 output_align_hmm,
                                                 output_alignment)
        return ptype, output_alignment
    
    def _mask_strange_sequence_letters(self, sequences, package_type):
        '''Replace strange characters like selenocysteine (U) with X or N for
        protein and nucleotide sequences, respectively. Sequences are 
        modified in place.
        
        Parameters
        ----------
        aligned_sequences: array of Sequence objects
            sequences to mask
        package_type: _PROTEIN_PACKAGE_TYPE or _NUCLEOTIDE_PACKAGE_TYPE
            type of sequences these are
        
        Returns
        -------
        None
        '''
        if package_type == Create._PROTEIN_PACKAGE_TYPE:
            search_re = re.compile(r'([^ACDEFGHIKLMNPQRSTVWY\-X])',re.IGNORECASE)
            replace_char = 'X'
        elif package_type == Create._NUCLEOTIDE_PACKAGE_TYPE:
            search_re = re.compile(r'([^ATGC\-N])',re.IGNORECASE)
            replace_char = 'N'
        
        for s in sequences:
            if package_type == Create._NUCLEOTIDE_PACKAGE_TYPE:
                s.seq = s.seq.replace('U','T')
                s.seq = s.seq.replace('u','t')
            newseq = s.seq
            m = search_re.search(newseq)
            while m:
                logging.warn("Found a non-standard character in the sequence of %s: e.g. '%s'"\
                             % (s.name, m.group()))
                newseq = newseq[:m.start()] + replace_char + newseq[(m.start()+1):]
                m = search_re.search(newseq)
            s.seq = newseq

    def _check_for_duplicate_sequence_names(self, fasta_file_path):
        """Test if the given fasta file contains sequences with duplicate
        sequence names.

        Parameters
        ----------
        fasta_file_path: string
            path to file that is to be checked

        Returns
        -------
        The name of the first duplicate sequence found, else False.

        """
        found_sequence_names = set()
        for record in SeqIO.parse(fasta_file_path, 'fasta'):
            name = record.name
            if name in found_sequence_names:
                return name
            found_sequence_names.add(name)
        return False



    def main(self, **kwargs):
        alignment = kwargs.pop('alignment',None)
        sequences = kwargs.pop('sequences',None)
        taxonomy = kwargs.pop('taxonomy',None)
        rerooted_tree = kwargs.pop('rerooted_tree',None)
        tree_log = kwargs.pop('tree_log', None)
        prefix = kwargs.pop('prefix', None)
        rerooted_annotated_tree = kwargs.pop('rerooted_annotated_tree', None)
        user_hmm = kwargs.pop('hmm', None)
        search_hmm_files = kwargs.pop('search_hmm_files',None)
        min_aligned_percent = kwargs.pop('min_aligned_percent',0.01)
        taxtastic_taxonomy = kwargs.pop('taxtastic_taxonomy', None)
        taxtastic_seqinfo = kwargs.pop('taxtastic_seqinfo', None)
        force_overwrite = kwargs.pop('force',False)
        graftm_package = kwargs.pop('graftm_package',False)
        dereplication_level = kwargs.pop('dereplication_level',False)
        threads = kwargs.pop('threads',False)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        seqio = SequenceIO()
        locus_name = (os.path.basename(sequences).split('.')[0]
                      if sequences
                      else os.path.basename(alignment).split('.')[0])
        tmp = tempdir.TempDir()
        base = os.path.join(tmp.name, locus_name)
        insufficiently_aligned_sequences = [None]
        removed_sequence_names = []

        if prefix:
            output_gpkg_path = prefix
        else:
            output_gpkg_path = "%s.gpkg" % locus_name

        if os.path.exists(output_gpkg_path):
            if force_overwrite:
                logging.warn("Deleting previous directory %s" % output_gpkg_path)
                shutil.rmtree(output_gpkg_path)
            else:
                raise Exception("Cowardly refusing to overwrite gpkg to already existing %s" % output_gpkg_path)
        logging.info("Building gpkg for %s" % output_gpkg_path)

        # Read in taxonomy somehow
        gtns = Getaxnseq()
        if rerooted_annotated_tree:
            logging.info("Building seqinfo and taxonomy file from input annotated tree")
            taxonomy_definition = TaxonomyExtractor().taxonomy_from_annotated_tree(\
                    TreeNode.read(rerooted_annotated_tree))
        elif taxonomy:
            logging.info("Building seqinfo and taxonomy file from input taxonomy")
            taxonomy_definition = GreenGenesTaxonomy.read_file(taxonomy).taxonomy
        elif taxtastic_seqinfo and taxtastic_taxonomy:
            logging.info("Reading taxonomy from taxtastic taxonomy and seqinfo files")
            taxonomy_definition = gtns.read_taxtastic_taxonomy_and_seqinfo\
                (open(taxtastic_taxonomy),
                 open(taxtastic_seqinfo))
        else:
            raise Exception("Taxonomy is required somehow e.g. by --taxonomy or --rerooted_annotated_tree")

        # Check for duplicates
        if sequences:
            dup = self._check_for_duplicate_sequence_names(sequences)
            if dup:
                raise Exception("Found duplicate sequence name '%s' in sequences input file" % dup)
        if alignment:
            dup = self._check_for_duplicate_sequence_names(alignment)
            if dup:
                raise Exception("Found duplicate sequence name '%s' in alignment input file" % dup)

        output_alignment = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.faa').name
        align_hmm = (user_hmm if user_hmm else tempfile.NamedTemporaryFile(prefix='graftm', suffix='_align.hmm').name)
        if sequences:
            if alignment:
                ptype = self._get_hmm_from_alignment(alignment,
                                                     align_hmm,
                                                     output_alignment)
            else:
                ptype, output_alignment = self._align_and_create_hmm(sequences, alignment, user_hmm,
                                                   align_hmm, output_alignment, threads)

            logging.info("Checking for incorrect or fragmented reads")
            insufficiently_aligned_sequences = self._check_reads_hit(open(output_alignment),
                                                                     min_aligned_percent)
            while len(insufficiently_aligned_sequences) > 0:
                logging.warn("One or more alignments do not span > %.2f %% of HMM" % (min_aligned_percent*100))
                for s in insufficiently_aligned_sequences:
                    logging.warn("Insufficient alignment of %s, not including this sequence" % s)

                _, sequences2 = tempfile.mkstemp(prefix='graftm', suffix='.faa')
                num_sequences = self._remove_sequences_from_alignment(insufficiently_aligned_sequences,
                                                                      sequences,
                                                                      sequences2)
                sequences = sequences2
                if alignment:

                    _, alignment2 = tempfile.mkstemp(prefix='graftm', suffix='.aln.faa')

                    num_sequences = self._remove_sequences_from_alignment(insufficiently_aligned_sequences,
                                                                          alignment,
                                                                          alignment2)
                    alignment = alignment2
                    for name in insufficiently_aligned_sequences:
                        if rerooted_tree or rerooted_annotated_tree:
                            logging.warning('''Sequence %s in provided alignment does not meet the --min_aligned_percent cutoff. This sequence will be removed from the tree
in the final GraftM package. If you are sure these sequences are correct, turn off the --min_aligned_percent cutoff, provide it with a 0 (e.g. --min_aligned_percent 0) ''' % name)
                        removed_sequence_names.append(name)


                logging.info("After removing %i insufficiently aligned sequences, left with %i sequences" % (len(insufficiently_aligned_sequences), num_sequences))
                if num_sequences < 4:
                    raise Exception("Too few sequences remaining in alignment after removing insufficiently aligned sequences: %i" % num_sequences)
                else:
                    logging.info("Reconstructing the alignment and HMM from remaining sequences")
                    output_alignment = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.faa').name
                    if not user_hmm:
                        align_hmm = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.hmm').name
                    ptype, output_alignment= self._align_and_create_hmm(sequences, alignment, user_hmm,
                                                       align_hmm, output_alignment, threads)
                    logging.info("Checking for incorrect or fragmented reads")
                    insufficiently_aligned_sequences = self._check_reads_hit(open(output_alignment),
                                                                             min_aligned_percent)
            if not search_hmm_files:
                search_hmm = tempfile.NamedTemporaryFile(prefix='graftm', suffix='_search.hmm').name
                self._create_search_hmm(sequences, taxonomy_definition, search_hmm, dereplication_level, threads)
                search_hmm_files = [search_hmm]

            # Make sure each sequence has been assigned a taxonomy:
            aligned_sequence_objects = seqio.read_fasta_file(output_alignment)
            unannotated = []
            for s in aligned_sequence_objects:
                if s.name not in taxonomy_definition:
                    unannotated.append(s.name)
            if len(unannotated) > 0:
                for s in unannotated:
                    logging.error("Unable to find sequence '%s' in the taxonomy definition" % s)
                raise Exception("All sequences must be assigned a taxonomy, cannot continue")
        else:
            ptype = self._get_hmm_from_alignment(alignment,
                                                 align_hmm,
                                                 output_alignment)

            aligned_sequence_objects = seqio.read_fasta_file(output_alignment)
            unannotated = []
            for s in aligned_sequence_objects:
                if s.name not in taxonomy_definition:
                    unannotated.append(s.name)
            if len(unannotated) > 0:
                for s in unannotated:
                    logging.error("Unable to find sequence '%s' in the taxonomy definition" % s)
                raise Exception("All sequences must be assigned a taxonomy, cannot continue")
            
        logging.debug("Looking for non-standard characters in aligned sequences")
        self._mask_strange_sequence_letters(aligned_sequence_objects, ptype)

        # Deduplicate sequences - pplacer cannot handle these
        logging.info("Deduplicating sequences")
        dedup = Deduplicator()
        deduplicated_arrays = dedup.deduplicate(aligned_sequence_objects)
        deduplicated_taxonomy = dedup.lca_taxonomy(deduplicated_arrays, taxonomy_definition)
        deduplicated_taxonomy_hash = {}
        for i, tax in enumerate(deduplicated_taxonomy):
            deduplicated_taxonomy_hash[deduplicated_arrays[i][0].name] = tax
        deduplicated_alignment_file = base+"_deduplicated_aligned.fasta"
        seqio.write_fasta_file([seqs[0] for seqs in deduplicated_arrays],
                               deduplicated_alignment_file)

        logging.info("Removed %i sequences as duplicates, leaving %i non-identical sequences"\
                     % ((len(aligned_sequence_objects)-len(deduplicated_arrays)),
                        len(deduplicated_arrays)))

        # Get corresponding unaligned sequences
        filtered_names=[]
        for list in [x for x in [x[1:] for x in deduplicated_arrays] if x]:
            for seq in list:
                filtered_names.append(seq.name)
        _, sequences2 = tempfile.mkstemp(prefix='graftm', suffix='.faa')


        # Create tree unless one was provided
        if not rerooted_tree and not rerooted_annotated_tree:
            logging.debug("No tree provided")
            logging.info("Building tree")
            log_file, tre_file = self._build_tree(deduplicated_alignment_file, 
                                                  base, ptype, 
                                                  self.fasttree)
            no_reroot = False
        else:
            if rerooted_tree:
                logging.debug("Found unannotated pre-rerooted tree file %s" % rerooted_tree)
                tre_file=rerooted_tree
            elif rerooted_annotated_tree:
                logging.debug("Found annotated pre-rerooted tree file %s" % rerooted_tree)
                tre_file=rerooted_annotated_tree
            else:
                raise
            no_reroot = True


            # Remove any sequences from the tree that are duplicates
            cleaner = TreeCleaner()
            tree = TreeNode.read(tre_file)
            for group in deduplicated_arrays:
                [removed_sequence_names.append(s.name) for s in group[1:]]
            cleaner.remove_sequences(tree, removed_sequence_names)

            # Ensure there is nothing amiss now as a user-interface thing
            cleaner.match_alignment_and_tree_sequence_ids(\
                [g[0].name for g in deduplicated_arrays], tree)

            if tree_log:
                # User specified a log file, go with that
                logging.debug("Using user-specified log file %s" % tree_log)
                log_file = tree_log
            else:
                logging.info("Generating log file")
                log_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree_log', prefix='graftm')
                log_file = log_file_tempfile.name
                # Make the newick file simple (ie. un-arb it) for fasttree
                cleaner.clean_newick_file_for_fasttree_input(tree)
                tre_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file = tre_file_tempfile.name
                with tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm') as f:
                    tree.write(f.name)
                    self._generate_tree_log_file(f.name, deduplicated_alignment_file,
                                            tre_file, log_file, ptype, self.fasttree)

        # Create tax and seqinfo .csv files
        taxonomy_to_keep=[
                          seq.name for seq in
                                [x for x in [x[0] for x in deduplicated_arrays]
                                 if x]
                          ]
        refpkg = "%s.refpkg" % output_gpkg_path
        self.the_trash.append(refpkg)
        if taxtastic_taxonomy and taxtastic_seqinfo:
            logging.info("Creating reference package")
            refpkg = self._taxit_create(base, deduplicated_alignment_file,
                                        tre_file, log_file, taxtastic_taxonomy,
                                        taxtastic_seqinfo, refpkg, no_reroot)
        else:
            gtns = Getaxnseq()
            seq = base+"_seqinfo.csv"
            tax = base+"_taxonomy.csv"
            self.the_trash += [seq, tax]
            if rerooted_annotated_tree:
                logging.info("Building seqinfo and taxonomy file from input annotated tree")
                taxonomy_definition = TaxonomyExtractor().taxonomy_from_annotated_tree(\
                        TreeNode.read(rerooted_annotated_tree))
            elif taxonomy:
                logging.info("Building seqinfo and taxonomy file from input taxonomy")
                taxonomy_definition = GreenGenesTaxonomy.read_file(taxonomy).taxonomy
            else:
                raise Exception("Programming error: Taxonomy is required somehow e.g. by --taxonomy or --rerooted_annotated_tree")

            taxonomy_definition = {x:taxonomy_definition[x]
                                   for x in taxonomy_definition
                                   if x in taxonomy_to_keep}

            gtns.write_taxonomy_and_seqinfo_files(taxonomy_definition,
                                                  tax,
                                                  seq)

            # Create the reference package
            logging.info("Creating reference package")
            refpkg = self._taxit_create(base, deduplicated_alignment_file,
                                        tre_file, log_file, tax, seq, refpkg,
                                        no_reroot)
        if sequences:
            # Run diamond makedb
            logging.info("Creating diamond database")
            if ptype == Create._PROTEIN_PACKAGE_TYPE:
                cmd = "diamond makedb --in '%s' -d '%s'" % (sequences, base)
                extern.run(cmd)
                diamondb = '%s.dmnd' % base
            elif ptype == Create._NUCLEOTIDE_PACKAGE_TYPE:
                diamondb = None
            else: raise Exception("Programming error")
        else:
            diamondb = None

        if sequences:
            # Get range
            max_range = self._define_range(sequences)
        else:
            max_range = self._define_range(alignment)

        # Compile the gpkg
        logging.info("Compiling gpkg")
        
        GraftMPackageVersion3.compile(output_gpkg_path, refpkg, align_hmm, diamondb,
                                      max_range, sequences, search_hmm_files=search_hmm_files)
        
        logging.info("Cleaning up")
        self._cleanup(self.the_trash)

        # Test out the gpkg just to be sure.
        #
        # TODO: Use graftM through internal means rather than via extern. This
        # requires some refactoring so that graft() can be called easily with
        # sane defaults.
        logging.info("Testing gpkg package works")
        temp_output = tempdir.TempDir()
        graftM_graft_test_dir_name = temp_output.name
        temp_output.dissolve()
        cmd = "graftM graft --forward %s --graftm_package %s --output_directory %s" %(
            sequences, output_gpkg_path, graftM_graft_test_dir_name)
        extern.run(cmd)
        

        logging.info("Finished\n")

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






