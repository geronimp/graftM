import subprocess
import os
import re
import timeit
import itertools
import logging

from Bio import SeqIO
from collections import OrderedDict

from graftm.housekeeping import HouseKeeping
from graftm.hmmsearcher import HmmSearcher, NhmmerSearcher
from graftm.orfm import OrfM, ZcatOrfM
from graftm.unpack_sequences import UnpackRawReads
from graftm.readHmmTable import HMMreader
from _struct import unpack
FORMAT_FASTA    = "FORMAT_FASTA"
FORMAT_FASTQ    = "FORMAT_FASTQ"
FORMAT_FASTQ_GZ = "FORMAT_FASTQ_GZ"
FORMAT_FASTA_GZ = "FORMAT_FASTA_GZ"

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
        Nothing - output files are known.
        
        '''
        # Align input reads to a specified hmm.
        if run_stats['rev_true']:
            read_info = run_stats['reads']
            reverse = []
            forward = []
            records = list(SeqIO.parse(open(input_path), 'fasta'))

            # Split the reads into reverse and forward lists
            for record in records:

                if read_info[record.id]['direction'] == '+':
                    forward.append(record)

                elif read_info[record.id]['direction'] == '-':
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

        # If there are only forward reads, just hmmalign and be done with it.
        else:
            cmd = 'hmmalign --trim %s %s | seqmagick convert --input-format stockholm - %s' % (self.aln_hmm,
                                                                             input_path,
                                                                             for_conv_file)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(cmd, shell=True)

    def makeSequenceBinary(self, sequences, fm):
        cmd='makehmmerdb %s %s' % (sequences, fm)
        subprocess.check_call(cmd, shell=True)

    def hmmsearch(self, output_path, input_path, unpack, seq_type, threads, eval, orfm):
        '''
        hmmsearch - Search raw reads for hits using search_hmm list
        
        Parameters
        ----------
        output_path : str
            path to output domtblout table
        input_path : str
            path to input sequences to search
        unpack : obj
            Object that builds the commad chunk for unpacking the raw sequences.
            First guesses file format, and unpacks appropriately. Calls 
            command_line to construct final command line string.
        seq_type : var
            variable containing a string, either 'nucleotide' or 'protein'.
            Tells the pipeline whether or not to call ORFs on the sequence.
            If sequence is 'nucleotide', ORFs are called. If not, no ORFs.
        threads : str
            Number of threads to use. Passed to HMMsearch command. 
        eval : str
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
        elif seq_type == 'protein': # If the input is amino acid sequence
            input_cmd=unpack.command_line()
        else:
            raise Exception('Programming Error: error guessing input sequence type')
        
        # Run the HMMsearches
        searcher = HmmSearcher(threads, eval)
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
        while len(aln_list)>0:
            forward_path=aln_list.pop(0)
            reverse_path=aln_list.pop(0)
            output_path=outputs.pop(0)
            logging.info('Merging pair %s, %s' % (os.path.basename(forward_path), os.path.basename(reverse_path)))
            forward_reads=SeqIO.parse(forward_path,'fasta')
            reverse_reads=SeqIO.to_dict(SeqIO.parse(reverse_path,'fasta'))
            
            with open(output_path, 'w') as out:
                for forward_record in forward_reads:
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
                    
    def nhmmer(self, output_path, unpack, threads, eval):
        ## Run a nhmmer search on input_path file and return the name of
        ## resultant output table. Keep log of command.
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
        
        searcher = NhmmerSearcher(threads, extra_args=eval)
        searcher.hmmsearch(input_pipe, self.search_hmm, output_table_list)
        
        return output_table_list

    def hmmtable_reader(self, hmmtable):
        hash = {}
        seen = {}
        
        def buildHash(hit, program):
            if program == 'hmmsearch':
                if float(hit[17]) - float(hit[18]) > 0:
                    len = float(hit[17]) - float(hit[18])
                elif float(hit[17]) - float(hit[18]) < 0:
                    len = float(hit[18]) - float(hit[17])
                read_hash= {'len': len,
                            'bit': float(hit[7]),
                            'hmmfrom':float(hit[16]),
                            'hmmto':float(hit[15]),
                            'alifrom':hit[17],
                            'alito':hit[18]}
            elif program == 'nhmmer':
                if float(hit[6]) - float(hit[7]) > 0:
                    len = float(hit[6]) - float(hit[7])
                elif float(hit[6]) - float(hit[7]) < 0:
                    len = float(hit[7]) - float(hit[6])
                read_hash = {'len':len,
                             'bit':float(hit[13]),
                             'direction':hit[11],
                             'hmmfrom':hit[4],
                             'hmmto':hit[5],
                             'alifrom':hit[6],
                             'alito':hit[7]}
            return read_hash
        for idx, table in enumerate(hmmtable):
            program = [line.rstrip().split()[2] for line in open(table).readlines() if line.startswith('# Program:')][0]
            for hit in [line.rstrip().split() for line in open(table).readlines() if not line.startswith('#')]:                
                read_name = hit[0]
                
                if read_name in seen: # If the read name has been seen before.. 
                    if seen[read_name]==idx:
                        hash[read_name].append(buildHash(hit, program))
                else:
                    hash[read_name]=[buildHash(hit, program)]
                seen[read_name]=idx

        logging.debug("Sequences found: %i" % (len(hash)))
        logging.debug("Sequenced found with >1 hit: %i" % (len([x for x in hash.values() if len(x)>1])))
        if any(hash.values()):
            logging.debug("Sequence length Average: %i Max: %i Min: %i" % (sum([x[0]['len'] for x in  hash.values()])/len(hash.values()),
                                                                           max([x[0]['len'] for x in  hash.values()]),
                                                                           min([x[0]['len'] for x in  hash.values()])))
            logging.debug("Bit score Average: %i Max: %i Min: %i" % (sum([x[0]['bit'] for x in  hash.values()])/len(hash.values()),
                                                                     max([x[0]['bit'] for x in  hash.values()]),
                                                                     min([x[0]['bit'] for x in  hash.values()])))            
        return hash
        
    def __check_euk_contamination(self, hmm_hit_tables, run_stats):
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
        run_stats : dict
            Updated form of the above.
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

    def csv_to_titles(self, output_path, input_path, run_stats, euk_check=False, euk_reads=[]):
        if euk_check:
            logging.info("Checking for Eukaryotic contamination")
            euk_reads, run_stats = self.__check_euk_contamination(
                                                    input_path,
                                                    run_stats
                                                    )
        
        # process hmmsearch/nhmmer results into a list of titles to 
        # <base_filename>_readnames.txt
        run_stats['reads'] = {key:item for key, item in \
                              self.hmmtable_reader(input_path).iteritems() \
                              if key not in euk_reads}
        count=sum([len(x) for x in run_stats['reads'].values()])
        
        # See if there are any reads in there reverse direction. 
        # Store True if so for later reference
        try:
            if any([x for x in sum(run_stats['reads'].values(), []) \
                    if x['direction'] =='-']):
                run_stats['rev_true'] = True
            else:
                run_stats['rev_true'] = False
        except KeyError:
            run_stats['rev_true'] = False
        if count > 0: # Return if there weren't any reads found
            logging.info('%s read(s) found' % (count))
        else: # Otherwise, report the number of reads
            logging.info('%s reads found, cannot continue with no information' % (len(run_stats['reads'].keys())))
            return run_stats, False
        
        # And write the read names to output
        orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
        with open(output_path, 'w') as output_file:
            for record in run_stats['reads'].keys():
                regex_match = orfm_regex.match(record)
                if regex_match is not None:
                    output_file.write(regex_match.groups(0)[0]+'\n')
                if regex_match is None:
                    output_file.write(record+'\n')
        return run_stats, output_path

    def extract_from_raw_reads(self, output_path, input_path, raw_sequences_path, input_file_format, read_stats):
        # Use the readnames specified to extract from the original sequence
        # file to a fasta formatted file.        
        def removeOverlaps(item):
            for a, b in itertools.combinations(item, 2):
                fromto_a=[int(a['alifrom']),int(a['alito'])]
                fromto_b=[int(b['alifrom']),int(b['alito'])]
                range_a=range(min(fromto_a), max(fromto_a))
                range_b=range(min(fromto_b), max(fromto_b))
                intersect_length=len(set(range_a).intersection(set(range_b)))
                if intersect_length > 0:
                    if range_a > range_b:
                        item.remove(b)
                    elif a in item:
                        item.remove(a)
                else:
                    continue
            return item
                    
        def extractMultipleHits(reads_path, stats):
            # Extra function that reads in hits and splits out the regions 
            # (usually in a contig) that hit the HMM as a distinct match.
            reads=SeqIO.to_dict(SeqIO.parse(reads_path, "fasta"))
            new_stats={}
            out_reads={}
            for key,item in stats.iteritems():
                item=removeOverlaps(item)

                if len(item)>1:
                    counter=0
                    for entry in item:
                        f=int(entry['alifrom'])-1
                        t=int(entry['alito'])-1
                        read_rename=key + '_%s' % str(counter)
                        out_reads[read_rename]=str(reads[key].seq)[f:t]
                        new_stats[read_rename]=entry
                        counter+=1
                else:
                    out_reads[key]=str(reads[key].seq)
                    new_stats[key]=item[0]
            out_path = reads_path[:-3]+'_split.fa'
            with open(out_path, 'w') as out:
                for key,item in out_reads.iteritems():
                    out.write(">%s\n" % (str(key)))
                    out.write("%s\n" % (str(item)))
            return new_stats, out_path
        
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
        
        # Check if there are reads that need splitting
        if any([x for x in read_stats if len(read_stats[x])>1]):
            read_stats, output_path=extractMultipleHits(output_path, read_stats) 
        else:
            new_stats={}
            for key, item in read_stats.iteritems():
                new_stats[key]=item[0]
            read_stats=new_stats
        return read_stats, output_path

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
        orfm: graftm.OrfM object with parameters already set
        '''
        'Extract only the orfs that hit the hmm, return sequence file with within.'        
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
            seen = []
            for table in output_table_list:
                for title in [x.split(' ')[0] for x in open(table).readlines() if not x.startswith('#')]:
                    if title not in seen:
                        output.write(str(title) + '\n')
                        seen.append(title)
                    else:
                        pass       
        
        # Extract the reads using the titles.
        cmd = 'fxtract -H -X -f %s %s > %s' % (orf_titles_path, raw_orf_path, orf_out_path)
        
        logging.debug("Running command: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
        
        # Return name of output file
        return orf_out_path

    def p_search(self, files, args, run_stats, base, unpack, raw_reads):
        '''Protein search pipeline - The searching step for the protein 
        pipeline, where hits are identified in the input reads through 
        HMMsearches
        
        Parameters
        ----------
        files : obj
            graftm_output_paths object.
        args : dict
            input arguments, including threads, evalue cutoffs, 
            min_orf_length, input_sequence_type, and restrict_read_length.
        base : str
            The name of the input file, stripped of all suffixes, and paths. 
            Used for creating file names with 'files' object.
        raw_reads : str
            The reads to be searched.

        Returns
        -------
        hit_orfs : str
            The output fasta file of reads that hit 
        run_stats : dict
            Updated run stats for the file being worked on (inluding run times 
            for each search and numbers of hits etc)
            
        Raises
        ------
        N/A
        
        Examples
        --------
        N/A
        '''
        start  = timeit.default_timer() # Start search timer
        
        unpack = UnpackRawReads(raw_reads)
        if unpack.is_zcattable():
            clazz = ZcatOrfM
        else:
            clazz = OrfM
        orfm = clazz(min_orf_length=args.min_orf_length,
                      restrict_read_length=args.restrict_read_length)
        extracting_orfm = OrfM(min_orf_length=args.min_orf_length,
                      restrict_read_length=args.restrict_read_length)

        hit_table = self.hmmsearch(files.hmmsearch_output_path(base),
                                   raw_reads,
                                   unpack,
                                   args.input_sequence_type,
                                   args.threads,
                                   args.eval,
                                   orfm)
        # Processing the output table to give you the readnames of the hits
        run_stats, hit_readnames = self.csv_to_titles(files.readnames_output_path(base),
                                                      hit_table,
                                                      run_stats)

        if not hit_readnames:
            return False, run_stats
        
        if args.input_sequence_type == 'nucleotide':
            # Only accept 1 HSP per protein sequence
            old_stats = run_stats['reads']
            read_stats = {}
            for orf_name, hsps in old_stats.iteritems():
                read_stats[orf_name] = [hsps[0]]
            run_stats['reads'] = read_stats
        
        # Extract the hits form the original raw read file
        run_stats['reads'], hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                                    hit_readnames,
                                                                    raw_reads,
                                                                    unpack.format(),
                                                                    run_stats['reads'])
        
        if args.input_sequence_type == 'nucleotide':
            # Extract the orfs of these reads that hit the original search
            hit_orfs = self.extract_orfs(hit_reads,
                                         files.orf_output_path(base),
                                         files.orf_hmmsearch_output_path(base),
                                         files.orf_titles_output_path(base),
                                         extracting_orfm,
                                         files.orf_fasta_output_path(base))
        elif args.input_sequence_type == 'protein':
            hit_orfs = hit_reads
        else:
            raise Exception('Programming Error')
        # Define the average read length of the hits
        run_stats['read_length'] = self.check_read_length(hit_orfs, "P")
        # Stop and log search timer
        stop = timeit.default_timer()
        run_stats['search_t'] = str(int(round((stop - start), 0)) )
        # Falsify some summary entries
        run_stats['euk_contamination'] = 'N/A'
        run_stats['euk_uniq'] = 'N/A'
        run_stats['euk_check_t'] = 'N/A'
        # Return hit reads, and summary hash
        return hit_orfs, run_stats

    def d_search(self, files, args, run_stats, base, unpack, raw_reads, euk_check):
        '''Nucleotide search pipeline - The searching step for the nucleotide
        pipeline, where hits are identified in the input reads through nhmmer 
        searches
        
        Parameters
        ----------
        files : obj
            graftm_output_paths object.
        args : dict
            input arguments, including threads, evalue cutoffs, 
            min_orf_length, input_sequence_type, and restrict_read_length.
        base : str
            The name of the input file, stripped of all suffixes, and paths. 
            Used for creating file names with 'files' object.
        input_file_format : var
            The input format of the file, either FORMAT_FASTA or 
            FORMAT_FASTQ_GZ.
        raw_reads : str
            The reads to be searched.
        euk_check : bool
            True False, whether to check the entire sample for euk reads.
        
        Returns
        -------
        hit_reads : str
            The output fasta file of reads that hit 
        run_stats : dict
            Updated run stats for the file being worked on (inluding run times 
            for each search and numbers of hits etc)
            
        Raises
        ------
        N/A
        
        Examples
        --------
        N/A
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

