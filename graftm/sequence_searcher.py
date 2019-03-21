
import extern
import os
import itertools
import logging
import tempfile
import subprocess

from Bio import SeqIO
from collections import OrderedDict
from StringIO import StringIO

from graftm.timeit import Timer
from graftm.hmmsearcher import HmmSearcher, NhmmerSearcher
from graftm.orfm import OrfM, ZcatOrfM
from graftm.diamond import Diamond
from graftm.sequence_search_results import SequenceSearchResult, HMMSearchResult
from graftm.readHmmTable import HMMreader
from graftm.db_search_results import DBSearchResult

FORMAT_FASTA = "FORMAT_FASTA"
FORMAT_FASTQ = "FORMAT_FASTQ"
FORMAT_FASTQ_GZ = "FORMAT_FASTQ_GZ"
FORMAT_FASTA_GZ = "FORMAT_FASTA_GZ"
PIPELINE_AA = "P"
PIPELINE_NT = "D"
PREVIOUS_SPAN_CUTOFF = 0.25
T = Timer()

class InterleavedFileError(Exception):
    pass

class SequenceSearcher:

    def __init__(self, search_hmm, aln_hmm=None):
        self.search_hmm = search_hmm
        self.aln_hmm = aln_hmm

    def _get_sequence_directions(self, search_result):
        sequence_directions = {}
        for result in search_result:
            for hit in result.each([SequenceSearchResult.QUERY_ID_FIELD, SequenceSearchResult.ALIGNMENT_DIRECTION]):
                sequence_directions[hit[0]]  ={"strand":[hit[1]], "entry": []}
        return sequence_directions


    def _hmmalign(self, input_path, directions, pipeline,
                  forward_reads_output_path, reverse_reads_output_path):
        '''
        Align reads to the aln_hmm. Receives unaligned sequences and
        aligns them.

        Parameters
        ----------
        input_path : str
            Filename of unaligned hits to be aligned
        directions : dict
            dictionary containing read names as keys, and complement
            as the entry (True=Forward, False=Reverse)
        pipeline: str
            either PIPELINE_AA = "P" or PIPELINE_NT = "D"
        forward_reads_output_fh: str
            Where to write aligned forward reads
        reverse_reads_output_fh: str
            Where to write aligned reverse reads
        Returns
        -------
        Nothing.
        '''

        if pipeline == PIPELINE_AA:
            reverse_direction_reads_present=False
        else:
            reverse_direction_reads_present=False in directions.values()

        with tempfile.NamedTemporaryFile(prefix='for_file', suffix='.fa') as for_file_fh:
            for_file = for_file_fh.name
            with tempfile.NamedTemporaryFile(prefix='rev_file', suffix='.fa') as rev_file_fh:
                rev_file = rev_file_fh.name
                # Align input reads to a specified hmm.
                if reverse_direction_reads_present:  # Any that are in the reverse direction would be True
                    reverse = []
                    forward = []
                    records = list(SeqIO.parse(open(input_path), 'fasta'))

                    # Split the reads into reverse and forward lists
                    for record in records:
                        read_id = record.id

                        if directions[read_id] == True:
                            forward.append(record)
                        elif directions[read_id] == False:
                            reverse.append(record)
                        else:
                            raise Exception(logging.error('Programming error: hmmalign'))
                            exit(1)

                    logging.debug("Found %i forward direction reads" % len(forward))
                    logging.debug("Found %i reverse direction reads" % len(reverse))

                    # Write reverse complement and forward reads to files
                    with open(for_file, 'w') as for_aln:
                        logging.debug("Writing forward direction reads to %s" % for_file)
                        for record in forward:
                            for_aln.write('>' + record.id + '\n')
                            for_aln.write(str(record.seq) + '\n')
                            # HMMalign and convert to fasta format
                    if any(forward):
                        self.hmmalign_sequences(self.aln_hmm, for_file, forward_reads_output_path)
                    else:
                        cmd = 'touch %s' % (forward_reads_output_path)
                        extern.run(cmd)
                    with open(rev_file, 'w') as rev_aln:
                        logging.debug("Writing reverse direction reads to %s" % rev_file)
                        for record in reverse:
                            if record.id and record.seq:
                                rev_aln.write('>' + record.id + '\n')
                                rev_aln.write(str(record.seq.reverse_complement()) + '\n')
                    self.hmmalign_sequences(self.aln_hmm, rev_file, reverse_reads_output_path)
                    conv_files = [forward_reads_output_path, reverse_reads_output_path]
                    return conv_files

                else:
                    # If there are only forward reads, just hmmalign and be done with it.
                    self.hmmalign_sequences(self.aln_hmm, input_path, forward_reads_output_path)
                    conv_files = [forward_reads_output_path]

                    return conv_files

    def hmmalign_sequences(self, hmm, sequences, output_file):
        '''Run hmmalign and convert output to aligned fasta format

        Parameters
        ----------
        hmm: str
            path to hmm file
        sequences: str
            path to file of sequences to be aligned
        output_file: str
            write sequences to this file

        Returns
        -------
        nothing
        '''

        cmd = 'hmmalign --trim %s %s' % (hmm, sequences)
        output = extern.run(cmd)
        with open(output_file, 'w') as f:
            SeqIO.write(SeqIO.parse(StringIO(output), 'stockholm'), f, 'fasta')

    def makeSequenceBinary(self, sequences, fm):
        cmd = 'makehmmerdb %s %s' % (sequences, fm)
        extern.run(cmd)

    def hmmsearch(self, output_path, input_path, unpack, seq_type, threads, cutoff, orfm):
        '''
        hmmsearch - Search raw reads for hits using search_hmm list

        Parameters
        ----------
        output_path : str
            path to output domtblout table
        input_path : str
            path to input sequences to search
        unpack : UnpackRawReads
            UnpackRawReads object, returns string command that will output
            sequences to stdout when called on command line
            (use: unpack.command_line())
        seq_type : str
            variable containing a string, either 'nucleotide' or 'aminoacid'.
            Tells the pipeline whether or not to call ORFs on the sequence.
            If sequence is 'nucleotide', ORFs are called. If not, no ORFs.
        threads : Integer
            Number of threads to use. Passed to HMMsearch command.
        cutoff : str
            cutoff for HMMsearch to use, either an evalue or --cut_tc, meaning
            use the TC score cutoff specified within the HMM. Passed to 
            HMMsearch command.
        orfm : OrfM
            Object that builds the command chunk for calling ORFs on sequences
            coming through as stdin. Outputs to stdout. Calls command_line
            to construct final command line string.

        Returns
        -------
        output_table_list : array of HMMSearchResult
            Includes the name of the output domtblout table given by hmmer

        Raises
        ------
        hmmsearcher.NoInputSequencesException
            Raised if there are no sequences fed into the HMM.
        '''

        # Define the base hmmsearch command.
        logging.debug("Using %i HMMs to search" % (len(self.search_hmm)))
        output_table_list = []
        if len(self.search_hmm) > 1:
            for hmm in self.search_hmm:
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] + '_' + os.path.split(output_path)[1])
                output_table_list.append(out)
        elif len(self.search_hmm) == 1:
            output_table_list.append(output_path)
        else:
            raise Exception("Programming error: expected 1 or more HMMs")

        # Choose an input to this base command based off the file format found.
        if seq_type == 'nucleotide':  # If the input is nucleotide sequence
            input_cmd = orfm.command_line(input_path)
        elif seq_type == 'aminoacid':  # If the input is amino acid sequence
            input_cmd = unpack.command_line()
        else:
            raise Exception('Programming Error: error guessing input sequence type')
        
        # Run the HMMsearches
        if cutoff == "--cut_tc":
            searcher = HmmSearcher(threads, cutoff) 
        else:
            searcher = HmmSearcher(threads, '--domE %s' % cutoff)
        searcher.hmmsearch(input_cmd, self.search_hmm, output_table_list)

        hmmtables = [HMMSearchResult.import_from_hmmsearch_table(x) for x in output_table_list]
        return hmmtables

    def merge_forev_aln(self, forward_aln_list, reverse_aln_list, outputs):
        '''
        merge_forev_aln - Merges forward and reverse alignments for a given run

        Parameters
        ----------
        aln_list : array
            List of the forward and reverse alignments for each of the runs
            given to graftM. **MUST** be the following pattern:
            [forward_run1, reverse_run1, forward_run2, reverse_run2 ...]
        outputs : array
            List of paths to output file to which the merged alignments from the
            aln_list will go into. Must be exactly half the size of the aln_list
            (i.e. one output file for every forward and reverse file provided)

        Returns
        -------
        Nothing - output files are known.
        '''

        orfm_regex = OrfM.regular_expression()
        def remove_orfm_end(records):

            new_dict = {}

            for key, record in records.iteritems():
                orfmregex = orfm_regex.match(key)
                if orfmregex:
                    new_dict[orfmregex.groups(0)[0]] = record
                else:
                    new_dict[key] = record
            return new_dict

        for idx, (forward_path, reverse_path) in enumerate(zip(forward_aln_list, reverse_aln_list)):
            output_path = outputs[idx]
            logging.info('Merging pair %s, %s' % (os.path.basename(forward_path), os.path.basename(reverse_path)))
            forward_reads = SeqIO.parse(forward_path, 'fasta')
            reverse_reads = remove_orfm_end(SeqIO.to_dict(SeqIO.parse(reverse_path, 'fasta')))

            with open(output_path, 'w') as out:
                for forward_record in forward_reads:
                    regex_match = orfm_regex.match(forward_record.id)
                    if regex_match:
                        id = regex_match.groups(0)[0]
                    else:
                        id = forward_record.id
                    forward_sequence = str(forward_record.seq)
                    try:
                        reverse_sequence = str(reverse_reads[id].seq)
                        new_seq = ''
                        if len(forward_sequence) == len(reverse_sequence):
                            for f, r in zip(forward_sequence, reverse_sequence):
                                if f == r:
                                    new_seq += f
                                elif f == '-' and r != '-':
                                    new_seq += r
                                elif r == '-' and f != '-':
                                    new_seq += f
                                elif f != '-' and r != '-':
                                    if f != r:
                                        new_seq += f
                                else:
                                    new_seq += '-'
                        else:
                            logging.error('Alignments do not match')
                            raise Exception('Merging alignments failed: Alignments do not match')
                        out.write('>%s\n' % forward_record.id)
                        out.write('%s\n' % (new_seq))
                        del reverse_reads[id]
                    except:

                        out.write('>%s\n' % forward_record.id)
                        out.write('%s\n' % (forward_sequence))
                for record_id, record in reverse_reads.iteritems():
                    out.write('>%s\n' % record.id)
                    out.write('%s\n' % (str(record.seq)))

    def nhmmer(self, output_path, unpack, threads, evalue):
        '''
        nhmmer - Search input path using nhmmer

        Parameters
        ----------
        output_path : str
            A string containing the path to the input sequences
        unpack : obj
            UnpackRawReads object, returns string command that will output
            sequences to stdout when called on command line
            (use: unpack.command_line())
        threads : str
            Number of threads to run. For compiling command line.
        evalue : str
            evalue to use. For compiling commmand line.

        Returns
        -------
        output_table_list : array
            Includes the name of the output domtblout table given by hmmer
        '''
        logging.debug("Using %i HMMs to search" % (len(self.search_hmm)))
        output_table_list = []
        if len(self.search_hmm) > 1:
            for hmm in self.search_hmm:
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] + '_' + os.path.split(output_path)[1])
                output_table_list.append(out)
        elif len(self.search_hmm) == 1:
            output_table_list.append(output_path)
        else:
            raise Exception("Programming error: Expected 1 or more HMMs")
        input_pipe = unpack.command_line()

        searcher = NhmmerSearcher(threads, extra_args='--incE %s -E %s' % (evalue, evalue))
        searcher.hmmsearch(input_pipe, self.search_hmm, output_table_list)

        hmmtables = [HMMSearchResult.import_from_nhmmer_table(x) for x in output_table_list]

        return hmmtables, output_table_list

    def _check_euk_contamination(self, hmm_hit_tables):
        '''
        check_euk_contamination - Check output HMM tables hits reads that hit
                                  the 18S HMM with a higher bit score.

        Parameters
        ----------
        hmm_hit_tables : array
            Array of paths to the output files produced by hmmsearch or
            nhmmer.
        run_stats : dict
            A dictionary to updatewith the number of unique 18S reads and reads
            detected by both 18S and non-18S HMMs

        Returns
        -------
        euk_reads : set
            Non-redundant set of all read names deemed to be eukaryotic
        '''
        
        euk_hit_table = HMMreader(hmm_hit_tables.pop(-1))
        other_hit_tables = [HMMreader(x) for x in hmm_hit_tables]
        reads_unique_to_eukaryotes = []
        reads_with_better_euk_hit = []

        for hit in euk_hit_table.names():
            bits = []
            for hit_table in other_hit_tables:
                if hit in hit_table.names():
                    bits.append(hit_table.bit(hit))
                else:
                    reads_unique_to_eukaryotes.append(hit)
            if bits:
                if any([x for x in bits if x > euk_hit_table.bit(hit)]):
                    continue
                else:
                    reads_with_better_euk_hit.append(hit)
            else:
                continue

        if len(reads_with_better_euk_hit) == 0:
            logging.info("No contaminating eukaryotic reads detected")
        else:
            logging.info("Found %s read(s) that may be eukaryotic" % len(reads_with_better_euk_hit + reads_unique_to_eukaryotes))

        euk_reads = set(reads_with_better_euk_hit + reads_unique_to_eukaryotes)

        return euk_reads


    def _extract_multiple_hits(self, hits, reads_path, output_path):
        '''
        splits out regions of a read that hit the HMM. For example when two of
        same gene are identified within the same contig, The regions mapping to
        the HMM will be split out and written out to a new file as a new record.

        Parameters
        ----------
        hits : dict
            A dictionary where the keys are the read names, the entry for each
            is a list of lists, each containing the range within the contig
            (or read) that mapped to the HMM. e.g.:

            {'read1': [[3, 126], [305, 413]],
             'read2': [[1, 125]],
             ...

        reads_path : str
            path to reads file containing each read or contig in FASTA format to
            be opened, and split.
        output_path : str
            path to file to which split reads will be written to in FASTA
            format.

        Returns
        -------
        Nothing, output path is known.
        '''

        complement_information = {}
        try:
            reads = SeqIO.to_dict(SeqIO.parse(reads_path, "fasta"))  # open up reads as dictionary
        except:
            logging.error("Multiple sequences found with the same ID. The input sequences are either ill formated or are interleaved. \
If you provided GraftM with an interleaved sequence file, please split them into forward and reverse reads, and provide to the the appropriate \
flags (--forward, --reverse). Otherwise, it appears that you have provided sequences with redundant IDs. GraftM doesn't know how to \
deal with these, so please remove/rename sequences with duplicate keys.")
            raise InterleavedFileError()
        with open(output_path, 'w') as out:
            for read_name, entry in hits.iteritems():  # For each contig
                ranges = entry["entry"]
                complements = entry["strand"]
                index = 1
                if len(ranges) > 1:  # if there are multiple hits in that contig
                    for r, c in zip(ranges, complements):  # for each of those hits

                        new_record = reads[read_name][r[0] - 1:r[1]]  # subset the record by the span of that hit
                        new_record.id = new_record.id + '_split_%i' % index  # give that subset record a new header
                        SeqIO.write(new_record, out, "fasta")  # and write it to output
                        index += 1  # increment the split counter
                        complement_information[new_record.id]=c

                else:  # Otherwise, just write the read back to the file
                    complement_information[read_name] = entry["strand"][0] 
                    SeqIO.write(reads[read_name], out, "fasta")
        
        return complement_information

    def _extract_from_raw_reads(self, output_path, input_reads, raw_sequences_path, input_file_format, hits):
        '''
        _extract_from_raw_reads - call fxtract to extract the hit sequences
        of the hmm/diamond search from the raw sequences file. Output into
        specified file

        Parameters
        ----------
        output_path : str
            Path of the desired output file
        input_reads : list
            list of read IDS to extract.
        raw_sequences_path : str
            Path to the raw sequences
        input_file_format : var
            Variable, either FORMAT_FASTA_GZ, FORMAT_FASTQ_GZ, FORMAT_FASTQ or
            FORMAT_FASTA, denoting the format of the input sequence
        hits : dict
            A hash with the readnames as the keys and the spans as the values

        Returns
        -------
        output_path: str
            Path to file containing extracted reads.
        
        direction_information: dict
            
            {read_1: False
             ...
             read n: True}
            
            where True = Forward direction
            and False = Reverse direction
        '''

        with tempfile.NamedTemporaryFile(prefix='_raw_extracted_reads.fa') as tmp:
            # Run fxtract to obtain reads form original sequence file
            fxtract_cmd = "fxtract -H -X -f /dev/stdin " 
            if input_file_format == FORMAT_FASTA:
                cmd = "%s %s > %s" % (fxtract_cmd, raw_sequences_path, tmp.name)
            elif input_file_format == FORMAT_FASTQ_GZ:
                cmd = "%s -z %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, tmp.name)
            elif input_file_format == FORMAT_FASTA_GZ:
                cmd = "%s -z %s > %s" % (fxtract_cmd, raw_sequences_path, tmp.name)
            elif input_file_format == FORMAT_FASTQ:
                cmd = "%s %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, tmp.name)
            else:
                raise Exception("Programming error")

            process = subprocess.Popen(["bash", "-c", cmd], 
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE)
            process.communicate('\n'.join(input_reads))
            complement_info = self._extract_multiple_hits(hits, tmp.name, output_path)  # split them into multiple reads
        
        return output_path, complement_info
        

    def alignment_correcter(self, alignment_file_list, output_file_name,
                            filter_minimum=None):
        '''
        Remove lower case insertions in alignment outputs from HMM align. Give
        a list of alignments, and an output file name, and each alignment will
        be corrected, and written to a single file, ready to be placed together
        using pplacer.

        Parameters
        ----------
        alignment_file_list : array
            List of strings, each the path to different alignments from the
            inputs provided to GraftM
        output_file_name : str
            The path and filename of the output file desired.
        filter_minimum : int
            minimum number of positions that must be aligned for each sequence
        Returns
        -------
        True or False, depending if reads were written to file

        '''

        corrected_sequences = {}
        for alignment_file in alignment_file_list:
            insert_list = []  # Define list containing inserted positions to be removed (lower case characters)
            sequence_list = list(SeqIO.parse(open(alignment_file, 'r'), 'fasta'))
            for sequence in sequence_list:  # For each sequence in the alignment
                for idx, nt in enumerate(list(sequence.seq)):  # For each nucleotide in the sequence
                    if nt.islower():  # Check for lower case character
                        insert_list.append(idx)  # Add to the insert list if it is
            insert_list = list(OrderedDict.fromkeys(sorted(insert_list, reverse=True)))  # Reverse the list and remove duplicate positions
            for sequence in sequence_list:  # For each sequence in the alignment
                new_seq = list(sequence.seq)  # Define a list of sequences to be iterable list for writing
                for position in insert_list:  # For each position in the removal list
                    del new_seq[position]  # Delete that inserted position in every sequence
                corrected_sequences['>' + sequence.id + '\n'] = (''.join(new_seq) + '\n').replace('~', '-')
        
        pre_filter_count=len(corrected_sequences)
        
        if filter_minimum:
            # Use '>' not '>=' here because the sequence is on a single line, 
            # but also includes a newline character at the end of the sequence
            corrected_sequences={key:item for key, item in corrected_sequences.iteritems() if len(item.replace('-', '')) > filter_minimum}
        
        post_filter_count=len(corrected_sequences)
        logging.info("Filtered %i short sequences from the alignment" % \
                        (pre_filter_count-post_filter_count)
                    )
        logging.info("%i sequences remaining" % post_filter_count)
        
        if len(corrected_sequences) >= 1:
            with open(output_file_name, 'w') as output_file:  # Create an open file to write the new sequences to
                for fasta_id, fasta_seq in corrected_sequences.iteritems():
                    output_file.write(fasta_id)
                    output_file.write(fasta_seq)
            return True
        else:
            return False


    def _extract_orfs(self, input_path, orfm, hit_readnames, output_path, search_method, sequence_frame_info_list=None):
        '''
        Call ORFs on a file with nucleotide sequences and extract the proteins
        whose name is in `hit_readnames`.

        Parameters
        ----------
        input_path : str
            Path to input nucleotide sequences in FASTA format.
        orfm: obj
            graftm.OrfM object with parameters already set
        hit_readnames : str
            path to a file containin the readnames of hits to the HMM, one per
            line.
        output_path : str
            Path to output orfs into, in FASTA format.
        search_method : str
            The method for searching, either 'hmmsearch' or 'diamond'
        sequence_frame_info : list
            A dataframe (list of lists) containing readname, alignment direction
            and alignment start point information
        '''
        if search_method == "hmmsearch":
            # Build and run command to extract ORF sequences:
            orfm_cmd = orfm.command_line()
            cmd = 'fxtract -H -X -f /dev/stdin <(%s %s) > %s' % (orfm_cmd, input_path, output_path)
            process = subprocess.Popen(["bash", "-c", cmd],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE)
            process.communicate('\n'.join(hit_readnames))

        elif search_method == "diamond":
            sequence_frame_info_dict = {x[0]:[x[1], x[2], x[3]] for x in sequence_frame_info_list}
            records = SeqIO.parse(input_path, "fasta")
            with open(output_path, 'w') as open_output_path:
                for record in records:
                    entry=sequence_frame_info_dict[record.id]
                    indfrom=(min(entry[2], entry[1])-1)
                    indto=max(entry[2], entry[1])
                    if entry[0] == False:
                        record.seq = max([x for x in record.seq[indfrom:indto].reverse_complement().translate().split("*")], key=len)
                    else:
                        record.seq = max([x for x in record.seq[indfrom:indto].translate().split("*")], key=len)
                    SeqIO.write(record, open_output_path, "fasta")
                open_output_path.flush()


    def _get_read_names(self, search_result, max_range):
        '''
        _get_read_names - loops through hmm hits and their alignment spans to
        determine if they are potentially linked (for example, if one gene in a 
        contig hits a hmm more than once, in two different conserved regions
        of that gene) and combines them into one 'hit' because they are
        technically the same. The total span of the hits deemed to be linked is
        returned.

        Parameters
        ----------
        search_result : obj
            SequenceSearchResult object with all paramaters defined. Used here
            to create rows containing information on alignment direction and
            alignment span.
        max_range : int
            Maximum range that a gene can extend within a contig. Any hits
            that extend beyond this length cannot be linked. max_range is
            set as 1.5 X the average length of all full length genes used
            in the search database. This is defined in the CONTENTS.json file
            within a gpkg.
        Returns
        -------
            Dictionary where keys are the contig/read name. The value for each
            entry is an array lists, one per hit in each contig, each with the
            span (min and max) of the alignment.
        '''

        splits = {}  # Define an output dictionary to be filled
        spans = []
        for result in search_result:  # Create a table (list of rows contain span, and complement information
            spans += list(
                         result.each(
                                    [SequenceSearchResult.QUERY_ID_FIELD,
                                        SequenceSearchResult.ALIGNMENT_DIRECTION,
                                        SequenceSearchResult.HIT_FROM_FIELD,
                                        SequenceSearchResult.HIT_TO_FIELD,
                                        SequenceSearchResult.QUERY_FROM_FIELD,
                                        SequenceSearchResult.QUERY_TO_FIELD]
                                       )
                         )

        for hit in spans:  # For each of these rows (i.e. hits)
            i = hit[0]  # set id to i
            c = hit[1]  # set complement to c
            ft = [min(hit[2:4]), max(hit[2:4])]  # set span as ft (i.e. from - to) - This is the amount covering the query
            qs = [min(hit[4:6]), max(hit[4:6])]  # seq the query span to qs - This is the amount covering the HMM

            if ft[0] == ft[1]: continue  # if the span covers none of the query, skip that entry (seen this before)

            if i not in splits:  # If the hit hasnt been seen yet
                splits[i] = {'span'       : [ft],
                             'strand'     : [c],
                             'query_span' : [qs]}  # add the span and complement as new entry

            else:  # otherwise (if it has been seen)
                for idx, entry in enumerate(splits[i]['span']):  # for each previously added entry
                    if splits[i]['strand'][idx] != c:  # If the hit is on the same complement strand
                        splits[i]['span'].append(ft)  # Add the new range to be split out in the future
                        splits[i]['strand'].append(c)  # Add the complement strand as well
                        splits[i]['query_span'].append(qs)
                        break

                    previous_qs      = splits[i]['query_span'][idx] # Get the query span of the previous hit

                    previous_q_range = set(range(previous_qs[0], previous_qs[1])) # Get the range of each
                    current_q_range  = set(range(qs[0], qs[1]))
                    query_overlap    = set(previous_q_range).intersection(current_q_range) # Find the intersection between the two ranges

                    previous_ft_span = set(range(entry[0], entry[1]))
                    current_ft_span  = set(range(ft[0], ft[1]))
                    if any(query_overlap): # If there is an overlap
                        ####################################################
                        # if the span over the actual read that hit the HMM
                        # for each hit overlap by > 25%, they are considered
                        # the same hit, and ignored
                        ####################################################
                        intersection_fraction = float(len(previous_ft_span.intersection(current_ft_span)))

                        if intersection_fraction / float(len(previous_ft_span)) >= PREVIOUS_SPAN_CUTOFF:
                            break
                        elif intersection_fraction / float(len(current_ft_span)) >= PREVIOUS_SPAN_CUTOFF:
                            break
                        else: # else (i.e. if the hit covers less that 25% of the sequence of the previous hit)
                            ####################################################
                            # If they made it this far, it means that the hits do not overlap.
                            # But one last check must be made to ensure they do not cover the same
                            # region in the HMM.
                            ####################################################
                            if len(query_overlap) > (len(current_q_range)*PREVIOUS_SPAN_CUTOFF): # if the overlap on the query HMM does not span over 25%
                                if (idx+1) == len(splits[i]['span']):
                                    splits[i]['span'].append(ft) # Add from-to as another entry, this is another hit.
                                    splits[i]['strand'].append(c) # Add strand info as well
                                    splits[i]['query_span'].append(qs)
                                    break # And break

                    if min(entry) < min(ft):  # if/else to determine which entry comes first (e.g. 1-5, 6-10 not 6-10, 1-5)
                        if max(ft) - min(entry) < max_range:  # Check if they lie within range of eachother
                            entry[1] = max(ft)  # ammend the entry if they are
                            break  # And break the loop
                    else:
                        if max(entry) - min(ft) < max_range:  # Check if they lie within range of eachother
                            entry[0] = min(ft)  # ammend the entry if they are
                            break  # And break the loop

                else:  # if no break occured (no overlap)
                    splits[i]['span'].append(ft)  # Add the new range to be split out in the future
                    splits[i]['strand'].append(c)  # Add the complement strand as well
                    splits[i]['query_span'].append(qs)
        return {key: {"entry":entry['span'], 'strand': entry['strand']} for key, entry in splits.iteritems()}  # return the dict, without strand information which isn't required.

    def _check_for_slash_endings(self, readnames):
        '''
        Provide a list of read names to be checked for the /1 or /2 endings
        which cause troubles downstream when working with paired reads. This
        function checks the list for ANY read ending in /1 or /2. The reads are
        first stripped of any comment lines, using the 'split' function

        Parameters
        ----------
        readnames : list
            list of strings, each a read name to be checked for the /1 or /2
            ending

        Returns
        -------
        Boolean. True if /1 or /2 endings are found, False if not.

        '''

        # Strip comment lines, if any
        readnames = [x.split()[0] for x in readnames]

        # Check for /1 and /2 endings:
        if any(x for x in readnames if x.endswith('/1') or x.endswith('/2')):
            result = True
        else:
            result = False

        return result

    @T.timeit
    def aa_db_search(self, files, base, unpack, search_method,
                     maximum_range, threads, evalue, min_orf_length,
                     restrict_read_length, diamond_database):
        '''
        Amino acid database search pipeline - pipeline where reads are searched
        as amino acids, and hits are identified using hmmsearch or diamond
        searches

        Parameters
        ----------
        files : obj
            graftm_output_paths object.
        base : str
            The name of the input file, stripped of all suffixes, and paths.
            Used for creating file names with 'files' object.
        unpack : obj
            UnpackRawReads object, returns string command that will output
            sequences to stdout when called on command line
            (use: unpack.command_line())
        search_method : str
            The method for searching, either 'hmmsearch' or 'diamond'
        maximum_range : int
            Maximum range that a gene can extend within a contig. Any hits
            that extend beyond this length cannot be linked. max_range is defined
            as 1.5 X the average length of all full length genes used in the
            search database. This is defined in the CONTENTS.json file within a
            gpkg.
        threads : int
            Number of threads for hmmer to use
        evalue : str
            evalue cutoff for hmmer to use
        min_orf_length : int
            minimum orf length for orfm to use
        restrict_read_length : int
            orf length to retrict orfm to.
        diamond_database : str
            Path to diamond database to use when searching. Set to 'None' if not
            using diamond pipeline
        Returns
        -------
        String path to amino acid fasta file of reads that hit
        '''
        # Define outputs
        if search_method == 'hmmsearch':
            output_search_file = files.hmmsearch_output_path(base)
        elif search_method == 'diamond':
            output_search_file = files.diamond_search_output_basename(base)
        hit_reads_fasta = files.fa_output_path(base)
        hit_reads_orfs_fasta = files.orf_fasta_output_path(base)

        return self.search_and_extract_orfs_matching_protein_database(\
                                                    unpack,
                                                    search_method,
                                                    maximum_range,
                                                    threads,
                                                    evalue,
                                                    min_orf_length,
                                                    restrict_read_length,
                                                    diamond_database,
                                                    output_search_file,
                                                    hit_reads_fasta,
                                                    hit_reads_orfs_fasta)

    def search_and_extract_orfs_matching_protein_database(self,
                                                      unpack,
                                                      search_method,
                                                      maximum_range,
                                                      threads,
                                                      evalue,
                                                      min_orf_length,
                                                      restrict_read_length,
                                                      diamond_database,
                                                      output_search_file,
                                                      hit_reads_fasta,
                                                      hit_reads_orfs_fasta):
        '''As per aa_db_search() except slightly lower level. Search an
        input read set (unpack) and then extract the proteins that hit together
        with their containing nucleotide sequences.

        Parameters
        ----------
        output_search_file: str
            path to hmmsearch output table or diamond basename
        hit_reads_fasta: str
            path to nucleotide sequences containing hit proteins
        hit_reads_orfs_fasta: str
            path to hit proteins, unaligned
            
        Returns
        -------
        direction_information: dict
                    
            {read_1: False
             ...
             read n: True}
            
            where True = Forward direction
            and False = Reverse direction
        
        result: DBSearchResult object containing file locations and hit 
        information
        
        '''

        # Define method of opening sequence files to stdout
        if unpack.is_zcattable():
            clazz = ZcatOrfM
        else:
            clazz = OrfM
        orfm = clazz(min_orf_length=min_orf_length,
                     restrict_read_length=restrict_read_length)
        extracting_orfm = OrfM(min_orf_length=min_orf_length,
                      restrict_read_length=restrict_read_length)

        if search_method == 'hmmsearch':
            # run hmmsearch
            search_result = self.hmmsearch(
                                           output_search_file,
                                           unpack.read_file,
                                           unpack,
                                           unpack.sequence_type(),
                                           threads,
                                           evalue,
                                           orfm
                                           )

        elif search_method == 'diamond':
            # run diamond
            search_result = Diamond(
                                     database=diamond_database,
                                     threads=threads,
                                     evalue=evalue,
                                     ).run(
                                           unpack.read_file,
                                           unpack.sequence_type(),
                                           daa_file_basename=output_search_file
                                           )
            search_result = [search_result]

        else:  # if the search_method isn't recognised
            raise Exception("Programming error: unexpected search_method %s" % search_method)

        orfm_regex = OrfM.regular_expression()
        
        if maximum_range:
            hits = self._get_read_names(
                                        search_result,  # define the span of hits
                                        maximum_range
                                        )
        else:   
            hits = self._get_sequence_directions(search_result)
    
        orf_hit_readnames = hits.keys() # Orf read hit names
        if unpack.sequence_type() == 'nucleotide':
            hits={(orfm_regex.match(key).groups(0)[0] if orfm_regex.match(key) else key): item for key, item in hits.iteritems()}
            hit_readnames = hits.keys() # Store raw read hit names 
        else:
            hit_readnames=orf_hit_readnames
        
        hit_reads_fasta, direction_information = self._extract_from_raw_reads(
                                                       hit_reads_fasta,
                                                       hit_readnames,
                                                       unpack.read_file,
                                                       unpack.format(),
                                                       hits
                                                       )
        

        if not hit_readnames:
            hit_read_counts = [0, len(hit_readnames)]
            result = DBSearchResult(None,
                                    search_result,
                                    hit_read_counts,
                                    None)

            return result, direction_information


        if unpack.sequence_type() == 'nucleotide':
            # Extract the orfs of these reads that hit the original search

            self._extract_orfs(
                               hit_reads_fasta,
                               extracting_orfm,
                               orf_hit_readnames,
                               hit_reads_orfs_fasta,
                               search_method,
                               list(search_result[0].each([SequenceSearchResult.QUERY_ID_FIELD,
                                                           SequenceSearchResult.ALIGNMENT_DIRECTION,
                                                           SequenceSearchResult.QUERY_FROM_FIELD,
                                                           SequenceSearchResult.QUERY_TO_FIELD])
                                    )
                               )

            hit_reads_fasta = hit_reads_orfs_fasta
        slash_endings=self._check_for_slash_endings(hit_readnames)
        result = DBSearchResult(hit_reads_fasta,
                                search_result,
                                [0, len([itertools.chain(*hits.values())])],  # array of hits [euk hits, true hits]. Euk hits alway 0 unless searching from 16S
                                slash_endings)  # Any reads that end in /1 or /2     

        if maximum_range:
            n_hits = sum([len(x["strand"]) for x in hits.values()])
        else:
            n_hits = len(hits.keys())
        logging.info("%s read(s) detected" % n_hits)
                
        return result, direction_information
    

    @T.timeit
    def nt_db_search(self, files, base, unpack, euk_check,
                     search_method, maximum_range, threads, evalue):
        '''
        Nucleotide database search pipeline - pipeline where reads are searched
        as nucleotides, and hits are identified using nhmmer searches

        Parameters
        ----------
        files : obj
            graftm_output_paths object.
        base : str
            The name of the input file, stripped of all suffixes, and paths.
            Used for creating file names with 'files' object.
        unpack : obj
            UnpackRawReads object, returns string command that will output
            sequences to stdout when called on command line
            (use: unpack.command_line())
        euk_check : bool
            True indicates the sample will be checked for eukaryotic reads,
            False indicates not.
        search_method : str
            The method for searching e.g. 'hmmsearch' or 'diamond'
        maximum_range : int
            Maximum range that a gene can extend within a contig. Any hits
            that extend beyond this length cannot be linked. max_range is defined
            as 1.5 X the average length of all full length genes used in the
            search database. This is defined in the CONTENTS.json file within a
            gpkg.
        threads : str
            Number of threads for hmmer to use
        evalue : str
            Evalue cutoff for hmmer to use
        Returns
        -------
        String path to amino acid fasta file of reads that hit
        '''

        # Define outputs
        hmmsearch_output_table = files.hmmsearch_output_path(base)
        hit_reads_fasta = files.fa_output_path(base)
        return \
            self.search_and_extract_nucleotides_matching_nucleotide_database(\
                                                      unpack,
                                                      euk_check,
                                                      search_method,
                                                      maximum_range,
                                                      threads,
                                                      evalue,
                                                      hmmsearch_output_table,
                                                      hit_reads_fasta)

    def search_and_extract_nucleotides_matching_nucleotide_database(self,
                                                      unpack,
                                                      euk_check,
                                                      search_method,
                                                      maximum_range,
                                                      threads,
                                                      evalue,
                                                      hmmsearch_output_table,
                                                      hit_reads_fasta):
        '''As per nt_db_search() except slightly lower level. Search an
        input read set (unpack) and then extract the sequences that hit.

        Parameters
        ----------
        hmmsearch_output_table: str
            path to hmmsearch output table
        hit_reads_fasta: str
            path to hit nucleotide sequences
        
        Returns
        -------
        direction_information: dict
                    
            {read_1: False
             ...
             read n: True}
            
            where True = Forward direction
            and False = Reverse direction
        
        result: DBSearchResult object containing file locations and hit 
        information
        '''

        if search_method == "hmmsearch":
            # First search the reads using the HMM
            search_result, table_list = self.nhmmer(
                                                    hmmsearch_output_table,
                                                    unpack,
                                                    threads,
                                                    evalue
                                                    )


        elif search_method == 'diamond':
            raise Exception("Diamond searches not supported for nucelotide databases yet")

        
        if maximum_range:  
            
            hits = self._get_read_names(
                                        search_result,  # define the span of hits
                                        maximum_range
                                        )
        else:   
            hits = self._get_sequence_directions(search_result)

        hit_readnames = hits.keys()
        
        if euk_check:
            euk_reads = self._check_euk_contamination(table_list)
            hit_readnames = set([read for read in hit_readnames if read not in euk_reads])
            hits = {key:item for key, item in  hits.iteritems() if key in hit_readnames}
            hit_read_count = [len(euk_reads), len(hit_readnames)]
        else:
            hit_read_count = [0, len(hit_readnames)]
        
        hit_reads_fasta, direction_information = self._extract_from_raw_reads(
                                                       hit_reads_fasta,
                                                       hit_readnames,
                                                       unpack.read_file,
                                                       unpack.format(),
                                                       hits
                                                       )
        
        if not hit_readnames:
            result = DBSearchResult(None,
                                  search_result,
                                  hit_read_count,
                                  None)
        else:
            slash_endings=self._check_for_slash_endings(hit_readnames)
            result = DBSearchResult(hit_reads_fasta,
                                    search_result,
                                    hit_read_count,
                                    slash_endings)

        if maximum_range:
            n_hits = sum([len(x["strand"]) for x in hits.values()])
        else:
            n_hits = len(hits)
        logging.info("%s read(s) detected" % n_hits)
        
        return result, direction_information
        
    @T.timeit
    def align(self, input_path, output_path, directions, pipeline, 
              filter_minimum):
        '''align - Takes input path to fasta of unaligned reads, aligns them to
        a HMM, and returns the aligned reads in the output path

        Parameters
        ----------
        input_path : str
        output_path : str
        reverse_direction : dict
            A dictionary of read names, with the entries being the complement
            strand of the read (True = forward, False = reverse)
        pipeline : str
            Either "P" or "D" corresponding to the protein and nucleotide (DNA)
            pipelines, respectively.


        Returns
        -------
        N/A - output alignment path known.
        '''

        # HMMalign the forward reads, and reverse complement reads.
        with tempfile.NamedTemporaryFile(prefix='for_conv_file', suffix='.fa') as fwd_fh:
            fwd_conv_file = fwd_fh.name
            with tempfile.NamedTemporaryFile(prefix='rev_conv_file', suffix='.fa') as rev_fh:
                rev_conv_file = rev_fh.name
                alignments = self._hmmalign(
                    input_path,
                    directions,
                    pipeline,
                    fwd_conv_file,
                    rev_conv_file)
                alignment_result = self.alignment_correcter(alignments,
                                                            output_path,
                                                            filter_minimum)
                return alignment_result
