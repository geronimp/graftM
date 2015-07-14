import subprocess
import os
import re
import itertools
import logging

import tempfile
import regex

from Bio import SeqIO
from collections import OrderedDict
from graftm.timeit import Timer
from graftm.housekeeping import HouseKeeping
from graftm.hmmsearcher import HmmSearcher, NhmmerSearcher
from graftm.orfm import OrfM, ZcatOrfM
from graftm.unpack_sequences import UnpackRawReads
from graftm.readHmmTable import HMMreader
from graftm.db_search_results import DBSearchResult

FORMAT_FASTA    = "FORMAT_FASTA"
FORMAT_FASTQ    = "FORMAT_FASTQ"
FORMAT_FASTQ_GZ = "FORMAT_FASTQ_GZ"
FORMAT_FASTA_GZ = "FORMAT_FASTA_GZ"
T=Timer()

class Hmmer:

    def __init__(self, search_hmm, aln_hmm=None):
        self.search_hmm = search_hmm
        self.aln_hmm = aln_hmm
        self.hk = HouseKeeping()
        

    def hmmalign(self, input_path, run_stats, for_file, rev_file, for_conv_file, rev_conv_file):
        '''
        hmmalign Align reads to the aln_hmm. Receives unaligned sequences and 
        aligns them.
        
        Parameters
        ----------
        input_path : str
            Filename of unaligned hits to be aligned
        run_stats : dict
            Dictionary containing run stats. Will have run time from alignment 
            added to it. Also tells us if there are reads in the reverse 
            compliment if in nucleotide pipeline
        for_file : str
            Output forward compliment path.
        rev_file : str
            Output reverse compliment path.
        for_conv_file : str
            Output alignment of forward reads with insertions removed
        rev_conv_file : str
            Output alignment of reverse reads with insertions removed
        
        Returns
        -------
        conv_files : array
            list containing the paths to files containing the forward and 
            reverse reads in aligned fasta format.
        
        '''

        
        for_file      = tempfile.NamedTemporaryFile(prefix='for_file', suffix='.fa').name
        rev_file      = tempfile.NamedTemporaryFile(prefix='rev_file', suffix='.fa').name
        for_conv_file = tempfile.NamedTemporaryFile(prefix='for_conv_file', suffix='.fa').name
        rev_conv_file = tempfile.NamedTemporaryFile(prefix='rev_conv_file', suffix='.fa').name
        
        # Align input reads to a specified hmm.
        if False in directions.values(): # Any that are in the reverse direction would be True
            reverse = []
            forward = []
            records = list(SeqIO.parse(open(input_path), 'fasta'))

            # Split the reads into reverse and forward lists
            for record in records:
                regex_match = orfm_regex.match(record.id)
                if regex_match:
                    read_id=regex_match.groups(0)[0]
                else:
                    read_id=record.id
                if directions[read_id] == True:
                    forward.append(record)
                elif directions[read_id] == False:
                    reverse.append(record)

                else:
                    raise Exception(logging.error('Programming error: hmmalign'))
                    exit(1)
            logging.debug("Found %i forward compliment reads" % len(forward))
            logging.debug("Found %i reverse compliment reads" % len(reverse))
            # Write reverse complement and forward reads to files
            
            
            with open(for_file, 'w') as for_aln:
                logging.debug("Writing forward compliment reads to %s" % for_file)
                for record in forward:
                    if record.id and record.seq: # Check that both the sequence and ID fields are there, HMMalign will segfault if not.
                        for_aln.write('>'+record.id+'\n')
                        for_aln.write(str(record.seq)+'\n')
                        # HMMalign and convert to fasta format
            if any(forward):
                cmd = 'hmmalign --trim %s %s | seqmagick convert --input-format stockholm - %s' % (self.aln_hmm,
                                                                                            for_file,
                                                                                            for_conv_file)
            else:
                cmd = 'touch %s' % (for_conv_file)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
            with open(rev_file, 'w') as rev_aln:
                logging.debug("Writing reverse compliment reads to %s" % rev_file)
                for record in reverse:
                    if record.id and record.seq:
                        rev_aln.write('>'+record.id+'\n')
                        rev_aln.write(str(record.seq.reverse_complement())+'\n')           
            cmd = 'hmmalign --trim %s %s | seqmagick convert --input-format stockholm - %s' % (self.aln_hmm,
                                                                                        rev_file,
                                                                                        rev_conv_file)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
            conv_files = [for_conv_file, rev_conv_file]
            return conv_files

        # If there are only forward reads, just hmmalign and be done with it.
        else:
            cmd = 'hmmalign --trim %s %s | seqmagick convert --input-format stockholm - %s' % (self.aln_hmm,
                                                                             input_path,
                                                                             for_conv_file)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
            conv_files = [for_conv_file]
            
            return conv_files
    
    def makeSequenceBinary(self, sequences, fm):
        cmd='makehmmerdb %s %s' % (sequences, fm)
        subprocess.check_call(cmd, shell=True)

    def hmmsearch(self, output_path, input_path, unpack, seq_type, threads, evalue, orfm):
        '''
        hmmsearch - Search raw reads for hits using search_hmm list
        
        Parameters
        ----------
        output_path : str
            path to output domtblout table
        input_path : str
            path to input sequences to search
        unpack : obj
            UnpackRawReads object, returns string command that will output
            sequences to stdout when called on command line 
            (use: unpack.command_line())
        seq_type : var
            variable containing a string, either 'nucleotide' or 'aminoacid'.
            Tells the pipeline whether or not to call ORFs on the sequence.
            If sequence is 'nucleotide', ORFs are called. If not, no ORFs.
        threads : str
            Number of threads to use. Passed to HMMsearch command. 
        evalue : str
            evalue cutoff for HMMsearch to use. Passed to HMMsearch command. 
        orfm : obj
            Object that builds the command chunch for calling ORFs on sequences
            coming through as stdin. Outputs to stdout. Calls command_line
            to construct final command line string.

        Returns
        -------
        output_table_list : array
            Includes the name of the output domtblout table given by hmmer          
        '''

        # Define the base hmmsearch command.
        logging.debug("Using %i HMMs to search" % (len(self.search_hmm)))
        output_table_list = []
        if len(self.search_hmm) > 1:
            for hmm in self.search_hmm:
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(output_path)[1])
                output_table_list.append(out)
        elif len(self.search_hmm) == 1:
            output_table_list.append(output_path)
        else:
            raise Exception("Programming error: expected 1 or more HMMs")
        
        # Choose an input to this base command based off the file format found.
        if seq_type == 'nucleotide': # If the input is nucleotide sequence
            input_cmd = orfm.command_line(input_path)
        elif seq_type == 'aminoacid': # If the input is amino acid sequence
            input_cmd=unpack.command_line()
        else:
            raise Exception('Programming Error: error guessing input sequence type')
        
        # Run the HMMsearches
        searcher = HmmSearcher(threads, '-E %s' % evalue)

        searcher.hmmsearch(input_cmd, self.search_hmm, output_table_list)
        
        return output_table_list

    def merge_forev_aln(self, aln_list, outputs): 
        '''
        merge_forev_aln - Merges forward and reverse alignments for a given run
        
        Parameters
        ----------
        aln_list : array
            List of the forward and reverse alignments for each of the runs 
            given to graftM. **MUST** be the following pattern: 
            [forward_run1, reverse_run1, forward_run2, reverse_run2 ...]
        outputs : array
            List of paths to output file to which the merged aligments from the
            aln_list will go into. Must be exactly half the size of the aln_list
            (i.e. one output file for every forward and reverse file provided)
        
        Returns
        -------
        Nothing - output files are known.
        '''
        
        
        orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')

        while len(aln_list)>0:
            forward_path=aln_list.pop(0)
            reverse_path=aln_list.pop(0)
            output_path=outputs.pop(0)
            logging.info('Merging pair %s, %s' % (os.path.basename(forward_path), os.path.basename(reverse_path)))
            forward_reads=SeqIO.parse(forward_path,'fasta')
            reverse_reads=SeqIO.to_dict(SeqIO.parse(reverse_path,'fasta'))
            reverse_reads={orfm_regex.match(key).groups(0)[0]: value for key, value in reverse_reads.iteritems() if orfm_regex.match(key)}
            with open(output_path, 'w') as out:
                for forward_record in forward_reads:
                    regex_match = orfm_regex.match(forward_record.id)
                    if regex_match:
                        id=regex_match.groups(0)[0]
                    else:
                        id=forward_record.id
                    forward_sequence=str(forward_record.seq)
                    try:
                        reverse_sequence=str(reverse_reads[id].seq)
                        new_seq=''
                        if len(forward_sequence)==len(reverse_sequence):
                            for f,r in zip(forward_sequence, reverse_sequence):
                                if f==r:
                                    new_seq+=f
                                elif f=='-' and r!='-':
                                    new_seq+=r
                                elif r=='-' and f!='-':
                                    new_seq+=f                               
                                elif f!='-' and r!='-':
                                    if f != r:
                                        new_seq+=f
                                else:
                                    new_seq+='-'
                        else:
                            logging.error('Alignments do not match')
                            raise Exception('Merging alignments failed: Alignments do not match')
                        out.write('>%s\n' % forward_record.id)
                        out.write('%s\n' % (new_seq))
                    except:
                        continue
                    

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
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(output_path)[1])
                output_table_list.append(out)
        elif len(self.search_hmm) == 1:
            output_table_list.append(output_path)
        else:
            raise Exception("Programming error: Expected 1 or more HMMs")
        input_pipe=unpack.command_line()
        searcher = NhmmerSearcher(threads, extra_args='-E %s' % evalue)
        searcher.hmmsearch(input_pipe, self.search_hmm, output_table_list)
        
        hmmtables=[HMMSearchResult.import_from_nhmmer_table(x) for x in output_table_list]
        
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
            list of all read names deemed to be eukaryotic

        '''
        euk_hit_table=HMMreader(hmm_hit_tables.pop(-1))
        other_hit_tables           = [HMMreader(x) for x in hmm_hit_tables]
        crossover                  = []
        reads_unique_to_eukaryotes = []   
        reads_with_better_euk_hit  = []     
        
        for hit in euk_hit_table.names():
            bits=[]
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
        
        run_stats['euk_uniq'] = len(reads_unique_to_eukaryotes)
        run_stats['euk_contamination'] = len(reads_with_better_euk_hit)
        euk_reads=reads_with_better_euk_hit + reads_unique_to_eukaryotes
        
        return euk_reads, run_stats

    def _extract_from_raw_reads(self, output_path, input_path, raw_sequences_path, input_file_format, hits):
        '''
        _extract_from_raw_reads - call fxtract to extract the hit sequences 
        of the hmm/diamond search from the raw sequences file. Output into 
        specified file
        
        Parameters
        ----------
        output_path : str
            Path of the desired output file
        input_path : str
            Path to a file containing read IDs, one per line.
        raw_sequences_path : str
            Path to the raw sequences
        input_file_format : var
            Variable, either FORMAT_FASTA_GZ, FORMAT_FASTQ_GZ, FORMAT_FASTQ or
            FORMAT_FASTA, denoting the format of the input sequence
        hits : dict
            A hash with the readnames as the keys and the spans as the values
        '''  

        def extractMultipleHits(stats, reads_path): 
            reads=SeqIO.to_dict(SeqIO.parse(reads_path, "fasta")) # open up reads as dictionary
            new_output_path = reads_path.split('.')[0] + '_split.fa' # define the output path
            with open(new_output_path, 'w') as out:
                for read_name, ranges in stats.iteritems(): # For each contig
                    index=1 
                    if len(ranges)>1: # if there are multiple hits in that contig
                        for r in ranges: # for each of those hits
                            new_record=reads[read_name][r[0]-1:r[1]] # subset the record by the span of that hit
                            new_record.id=new_record.id+'_split_%i' % index # give that subset record a new header
                            SeqIO.write(new_record, out, "fasta") # and write it to output
                            index+=1 # increment the split counter
                    else: # Otherwise, just write the read back to the file
                        SeqIO.write(reads[read_name], out, "fasta")
            return new_output_path # return the split reads file and celebrate.
                
        # Run fxtract to obtain reads form original sequence file
        fxtract_cmd = "fxtract -H -X -f %s " % input_path
        if input_file_format == FORMAT_FASTA:
            cmd = "%s %s > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "%s -z %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
        elif input_file_format == FORMAT_FASTA_GZ:
            cmd = "%s -z %s > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
        elif input_file_format == FORMAT_FASTQ:
            cmd = "%s %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
        else:
            raise Exception("Programming error")

        logging.debug("Running command: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
        
        # If any of the hits have more than one alignment span (i.e. >1 hit)
        if any([x for x in hits.values() if len(x)>1]):
            output_path = extractMultipleHits(hits, output_path) # split them into multiple reads
            return output_path
        else:
            return output_path
        
    def check_read_length(self, reads, pipe):
        lengths = []
        record_list = []
        # First check if the reverse pipe is happening, because the read names
        # are different.
        record_list += list(SeqIO.parse(open(reads, 'r'), 'fasta'))
        for record in record_list:
            lengths.append(len(record.seq))
        if pipe == "P":
            return (sum(lengths) / float(len(lengths)))/3
        elif pipe =="D":
            return sum(lengths) / float(len(lengths))
        
    def alignment_correcter(self, alignment_file_list, output_file_name):
        corrected_sequences = {}
        for alignment_file in alignment_file_list:
            insert_list = [] # Define list containing inserted positions to be removed (lower case characters)
            sequence_list = list(SeqIO.parse(open(alignment_file, 'r'), 'fasta'))
            for sequence in sequence_list: # For each sequence in the alignment
                for idx, nt in enumerate(list(sequence.seq)): # For each nucleotide in the sequence
                    if nt.islower(): # Check for lower case character
                        insert_list.append(idx) # Add to the insert list if it is
            insert_list = list(OrderedDict.fromkeys(sorted(insert_list, reverse = True))) # Reverse the list and remove duplicate positions
            for sequence in sequence_list: # For each sequence in the alignment
                new_seq = list(sequence.seq) # Define a list of sequences to be iterable list for writing
                for position in insert_list: # For each position in the removal list
                    del new_seq[position] # Delete that inserted position in every sequence
                corrected_sequences['>'+sequence.id+'\n'] = (''.join(new_seq)+'\n').replace('~','-')
        with open(output_file_name, 'w') as output_file: # Create an open file to write the new sequences to
                for fasta_id, fasta_seq in corrected_sequences.iteritems():
                    if any(c.isalpha() for c in fasta_seq):
                        output_file.write(fasta_id)
                        output_file.write(fasta_seq)

    def extract_orfs(self, input_path, raw_orf_path, hmmsearch_out_path, orf_titles_path, orfm, orf_out_path):
        '''
        extract_orfs - Return a path to a FASTA file containing ORFs that hit \
                       the HMM (self.search_hmm)
        
        Parameters
        ----------
        input_path : str
            Path to input sequences in fasta nucleotide format
        orfm: obj
            graftm.OrfM object with parameters already set
        orf_out_path
            Path to output fasta file, containing amino acid ORFs
        '''
        # Build the command
        output_table_list = []
        if len(self.search_hmm) > 1:
            for hmm in self.search_hmm:
                out = os.path.join(os.path.split(hmmsearch_out_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(hmmsearch_out_path)[1])
                output_table_list.append(out)
        elif len(self.search_hmm) == 1:
            output_table_list.append(hmmsearch_out_path)
        else:
            raise Exception("Programming error: expected 1 or more HMMs")
        
        # Call orfs on the sequences
        orfm_cmd = orfm.command_line()
        cmd = '%s %s > %s' % (orfm_cmd, input_path, raw_orf_path)
        
        logging.debug("Running command: %s" % cmd)
        subprocess.check_call(cmd, shell=True)

        cmd = 'cat %s' % raw_orf_path
        searcher = HmmSearcher(1)
        searcher.hmmsearch(cmd, self.search_hmm, output_table_list)
        
        with open(orf_titles_path, 'w') as output:

            reads= set(
                       list(
                            itertools.chain(
                                            *[HMMreader(x).names() for x \
                                              in output_table_list]
                                            )
                            )
                       )
            [output.write(x + '\n') for x in reads]
            
        # Extract the reads using the titles.
        cmd = 'fxtract -H -X -f %s %s > %s' % (orf_titles_path, raw_orf_path, orf_out_path)
        
        logging.debug("Running command: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
    
    def _get_read_names(self, search_result, r):
        '''
        _get_read_names - loops through hits and alingment spans of hmm hits
        and determines potentially linked hits (for example, if a gene in a 
        contig hits a hmm more than once, in two different conserved regions
        within a gene) and combines them into one 'hit'. The total span of the
        hits deemed to be linked is returned.
        
        Parameters
        ----------
        search_result : obj
            SequenceSearchResult object with all paramaters defined. Used here
            to create rows containing information on alignment direction and 
            alignment span.
        r : int
            Maximum range that a gene can extend within a contig. Any hits 
            that extend beyond this length cannot be linked. r is defined as
            1.5 X the average length of all full length genes used in the 
            database. This is defined in the CONTENTS.json file within a gpkg.
        Returns
        -------
            Dictionary where keys are the contig/read name. The value for each
            entry is an array lists, one per hit in each contig, each with the
            span (min and max) of the alignment. 
        '''
        
        splits = {} # Define an output dictionary to be filled
        for result in search_result: # Create a table (list of rows contain span, and compliment information
            spans = list(search_result[0].each([SequenceSearchResult.QUERY_ID_FIELD, 
                                                SequenceSearchResult.ALIGNMENT_DIRECTION,
                                                SequenceSearchResult.HIT_FROM_FIELD, 
                                                SequenceSearchResult.HIT_TO_FIELD]))
            for hit in spans: # For each of these rows (i.e. hits)
                i = hit[0] # set id to i
                c = hit[1] # set compliment to c
                ft = [min(hit[2:]), max(hit[2:])] # set span as ft (i.e. from - to)
                if ft[0] == ft[1]: continue # if the span covers none of the query, skip that entry (seen this before)
                if i not in splits: # If the hit hasnt been seen yet
                    splits[i]={'span':[ft],
                               'strand':[c]} # add the span and compliment as new entry
                else: # otherwise (if it has been seen)                   
                    for idx, entry in enumerate(splits[i]['span']): # for each previously added entry     
                        if splits[i]['strand'][idx] == c: # If the hit is on the same compliment strand
                            if min(entry) < min(ft): # if/else to determine which entry comes first (e.g. 1-5, 6-10 not 6-10, 1-5)
                                if max(ft)-min(entry) < r: # Check if they lie within range of eachother
                                    entry[1]=max(ft) # ammend the entry if they are
                                    break # And break the loop
                            else:
                                if max(entry)-min(ft) < r:  # Check if they lie within range of eachother
                                    entry[0] = min(ft) # ammend the entry if they are
                                    break # And break the loop
                    else: # if no break occured (no overlap)
                        splits[i]['span'].append(ft) # Add the new range to be split out in the future
                        splits[i]['strand'].append(c) # Add the complement strand as well
        
        return {key: entry['span'] for key, entry in splits.iteritems()} # return the dict, without strand information which isn't required.
    
    @T.timeit
    def aa_db_search(self, files, base, unpack, raw_reads, search_method, gpkg, 
                     threads, evalue, min_orf_length, restrict_read_length,
                     srch_aln_only):
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
        raw_reads : str
            The reads to be searched.
        search_method : str
            The method for searching, either 'hmmsearch' or 'diamond'
        gpkg : obj
            GraftM package object, created by graftm_package.py
        threads : str
            Number of threads for hmmer to use
        evalue : str
            evalue cutoff for hmmer to use
        min_orf_length : str
            minimum orf length for orfm to use
        restrict_read_length : str
            orf length to retrict orfm to.
        Returns
        -------
        hit_orfs : str
            The output fasta file of reads that hit 
        '''
        # Define outputs
        hmmsearch_output_table = files.hmmsearch_output_path(base)
        hit_reads_fasta        = files.fa_output_path(base)
        hit_reads_orfs_fasta   = files.orf_fasta_output_path(base)
        
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
                                           hmmsearch_output_table,
                                           raw_reads,
                                           unpack,
                                           unpack.sequence_type(),
                                           threads,
                                           evalue,
                                           orfm
                                           )

        elif search_method == 'diamond':
            #run diamond
            search_result =  Diamond(
                                     database=os.path.join(
                                                           gpkg.base_directory, \
                                                           gpkg.contents_hash[gpkg.DIAMOND_DATABASE_KEY]
                                                           ),
                                     threads=threads,
                                     evalue=evalue,
                                     ).run(
                                           raw_reads,
                                           unpack.sequence_type()
                                           )
            search_result=[search_result]

        else: # if the search_method isn't recognised
            raise Exception("Programming error: unexpected search_method %s" % search_method)
        
        with tempfile.NamedTemporaryFile(prefix='graftm_readnames') as readnames:
            # Write the names of hits to a tmpfile
            
            orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)') # to remove OrfM suffix from read names
            if srch_aln_only:
                hits = {key:[] for key in list(
                                               itertools.chain(
                                                               *itertools.chain(
                                                                                *[list(result.each([SequenceSearchResult.QUERY_ID_FIELD])) \
                                                                                  for result in search_result]
                                                                                )))}
            else:
                hits = self._get_read_names(search_result,  # define the span of hits
                                            gpkg.contents_hash[gpkg.RANGE_KEY])
            hit_readnames = [orfm_regex.match(name).groups(0)[0] if orfm_regex.match(name) else name for name in hits.keys()] # Define readnames
            for read in hit_readnames:
                readnames.write(read+'\n')

            readnames.flush()
            hit_reads_fasta = self._extract_from_raw_reads(
                                                           hit_reads_fasta,
                                                           readnames.name,
                                                           raw_reads,
                                                           unpack.format(),
                                                           hits
                                                           )

        if not hit_readnames:
            hit_read_counts=[0,len(hit_readnames)]
            result=DBSearchResult(None, 
                                  search_result, 
                                  hit_read_counts,
                                  None)
            return result
    
        if unpack.sequence_type() == 'nucleotide':
            # Extract the orfs of these reads that hit the original search
            self.extract_orfs(
                              hit_reads_fasta,
                              extracting_orfm,
                              hit_reads_orfs_fasta
                              )
            
            hit_reads_fasta=hit_reads_orfs_fasta
        
        result=DBSearchResult(hit_reads_fasta, 
                              search_result, 
                              [0,len(list(itertools.chain(*hits.values())))], # array of hits [euk hits, true hits]. Euk hits alway 0 unless searching from 16S
                              any([x for x in hit_readnames if x.endswith('\\1') or x.endswith('\\2')])) # Any reads that end in /1 or /2     
        return result
    
    @T.timeit
    def nt_db_search(self, files, base, unpack, raw_reads, euk_check, 
                     search_method, gpkg, threads, evalue, srch_aln_only):
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
        input_file_format : var
            The input format of the file, either FORMAT_FASTA or 
            FORMAT_FASTQ_GZ.
        raw_reads : str
            The reads to be searched.
        euk_check : bool
            True indicates the sample will be checked for eukaryotic reads, 
            False indicates not.
        search_method : str
            The method for searching e.g. 'hmmsearch' or 'diamond'
        gpkg : obj
            GraftM package object, created by graftm_package.py
        threads : str
            Number of threads for hmmer to use
        evalue : str
            Evalue cutoff for hmmer to use
        Returns
        -------
        hit_reads : str
            The output fasta file of reads that hit 
        '''
        
        # Define outputs
        hmmsearch_output_table = files.hmmsearch_output_path(base)
        hit_reads_fasta        = files.fa_output_path(base)
        
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

        with tempfile.NamedTemporaryFile(prefix='graftm_readnames') as readnames:
            if srch_aln_only:
                hits = {key:[[]] for key in list(
                                               itertools.chain(
                                                               *itertools.chain(
                                                                                *[list(result.each([SequenceSearchResult.QUERY_ID_FIELD])) \
                                                                                  for result in search_result]
                                                                                )))}
            else:
                hits = self._get_read_names(search_result,  # define the span of hits
                                            gpkg.contents_hash[gpkg.RANGE_KEY])
            hit_readnames = hits.keys()
            if euk_check:
                euk_reads  = self._check_euk_contamination(table_list)
                euk_reads  = set(euk_reads)
                prok_reads = set([hit_readnames for read in hit_readnames\
                                  if read not in euk_reads])
                for read in prok_reads: readnames.write(read+'\n')
                hit_read_count=[len(euk_reads), len(prok_reads)]
            else:
                [readnames.write(read+'\n') for read in hit_readnames]
                hit_read_count=[0,len(hit_readnames)]
            readnames.flush()
            
            hit_reads_fasta = self._extract_from_raw_reads(
                                                           hit_reads_fasta,
                                                           readnames.name,
                                                           raw_reads,
                                                           unpack.format(),
                                                           hits
                                                           )
        if not hit_readnames:
            result=DBSearchResult(None, 
                                  search_result, 
                                  hit_read_count,
                                  None)
        else:
            result=DBSearchResult(hit_reads_fasta, 
                                  search_result, 
                                  hit_read_count,
                                  any([x for x in hit_readnames if x.endswith('\1') or x.endswith('\2')]))
        
        return result
        
    @T.timeit
    def align(self, input_path, output_path, directions):
        '''align - Takes input path to fasta of unlaigned reads, aligns them to
        a HMM, and returns the aligned reads in the output path
        
        Parameters
        ----------
        input_path : str
        output_path : str
        reverse_direction : dict
            A dictionary of read names, with the entries being the complement 
            strand of the read (True = forward, False = reverse)
            
        Returns
        -------
        N/A - output alignment path known.

        '''
        start = timeit.default_timer() # Start search timer
        
        # First search the reads using the HMM
        hit_table = self.nhmmer(files.hmmsearch_output_path(base),
                                unpack,
                                args.threads,
                                args.eval)
        
        # Next, get a list of readnames
        run_stats, hit_readnames = self.csv_to_titles(files.readnames_output_path(base),
                                                      hit_table,
                                                      run_stats,
                                                      euk_check)
        if not hit_readnames:
            return False, run_stats

        # And extract them from the original sequence file
        run_stats['reads'], hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                                    hit_readnames,
                                                                    raw_reads,
                                                                    unpack.format(),
                                                                    run_stats['reads'])
        # Define the read length
        run_stats['read_length'] = self.check_read_length(hit_reads, "D")

        # Stop timing search and start timing euk check step.
        stop = timeit.default_timer()
        run_stats['search_t'] = str(int(round((stop - start), 0)) )
        start = timeit.default_timer()



        # Stop timing eukaryote check
        stop = timeit.default_timer()
        run_stats['euk_check_t'] = str(int(round((stop - start), 0)) )

        # Finally, return the hits
        return hit_reads, run_stats

    def align(self, files, args, run_stats, base, reads):

        # This pipeline takes unaligned reads, and aligns them agains a hmm,
        # regardless of their direction. Aligned reads with base insertions
        # removed are returned in the end. Times and commands are logged.

        start = timeit.default_timer()

        # HMMalign the forward reads, and reverse complement reads.
        self.hmmalign(reads,
                      run_stats,
                      files.output_for_path(base),
                      files.output_rev_path(base),
                      files.conv_output_for_path(base),
                      files.conv_output_rev_path(base))

        # Correct the alignment for base insertions.
        if run_stats['rev_true']:
            self.alignment_correcter([files.conv_output_for_path(base), files.conv_output_rev_path(base)],
                                     files.aligned_fasta_output_path(base))
        else:
            self.alignment_correcter([files.conv_output_for_path(base)],
                                     files.aligned_fasta_output_path(base))
        stop = timeit.default_timer()
        run_stats['aln_t'] = str(int(round((stop - start), 0)) )

        # Return
        return files.aligned_fasta_output_path(base), run_stats

