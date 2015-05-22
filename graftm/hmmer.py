import subprocess
import os
import re
import timeit
import itertools


from Bio import SeqIO
from collections import OrderedDict

from graftm.messenger import Messenger
from graftm.housekeeping import HouseKeeping

FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class Hmmer:

    def __init__(self, search_hmm, aln_hmm=None):
        self.search_hmm = search_hmm
        self.aln_hmm = aln_hmm
        self.hk = HouseKeeping()

    def hmmalign(self, input_path, run_stats, cmd_log, for_file, rev_file, for_sto_file, rev_sto_file, for_conv_file, rev_conv_file):
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
                    raise Exception(Messenger().error_message('Programming error: hmmalign'))
                    exit(1)

            # Write reverse complement and forward reads to files
            with open(for_file, 'w') as for_aln:
                for record in forward:
                    if record.id and record.seq: # Check that both the sequence and ID fields are there, HMMalign will segfault if not.
                        for_aln.write('>'+record.id+'\n')
                        for_aln.write(str(record.seq)+'\n')

            with open(rev_file, 'w') as rev_aln:
                for record in reverse:
                    if record.id and record.seq:
                        rev_aln.write('>'+record.id+'\n')
                        rev_aln.write(str(record.seq.reverse_complement())+'\n')


            # HMMalign and convert to fasta format
            cmd = 'hmmalign --trim -o %s %s %s 2>/dev/null; seqmagick convert %s %s' % (for_sto_file,
                                                                                        self.aln_hmm,
                                                                                        for_file,
                                                                                        for_sto_file,
                                                                                        for_conv_file)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)
            cmd = 'hmmalign --trim -o %s %s %s 2>/dev/null; seqmagick convert %s %s' % (rev_sto_file,
                                                                                        self.aln_hmm,
                                                                                        rev_file,
                                                                                        rev_sto_file,
                                                                                        rev_conv_file)

            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)

        # If there are only forward reads, just hmmalign and be done with it.
        else:
            cmd = 'hmmalign --trim -o %s %s %s ; seqmagick convert %s %s' % (for_sto_file,
                                                                             self.aln_hmm,
                                                                             input_path,
                                                                             for_sto_file,
                                                                             for_conv_file)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)


    def hmmsearch(self, output_path, input_path, input_file_format, seq_type, threads, eval, min_orf_length, restrict_read_length, cmd_log):
        '''Run a hmmsearch on the input_path raw reads, and return the name
        of the output table. Keep a log of the commands.'''
        # Define the base hmmsearch command.
        output_table_list = []
        tee = ' | tee'
        hmm_number = len(self.search_hmm)
        if hmm_number > 1:
            for idx, hmm in enumerate(self.search_hmm):
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(output_path)[1])
                output_table_list.append(out)
                if idx + 1 == hmm_number:
                    tee += " | hmmsearch %s --cpu %s --domtblout %s %s - >/dev/null " % (eval, threads, out, hmm)
                elif idx + 1 < hmm_number:
                    tee += " >(hmmsearch %s --cpu %s --domtblout %s %s - >/dev/null) " % (eval, threads, out, hmm)
                else:
                    raise Exception("Programming Error.") 

        elif hmm_number == 1:
            tee = ' | hmmsearch %s --cpu %s --domtblout %s %s - >/dev/null' % (eval, threads, output_path, self.search_hmm[0])
            output_table_list.append(output_path)
        
        # Choose an input to this base command based off the file format found.
        if seq_type == 'nucleotide': # If the input is nucleotide sequence
            orfm_cmdline = self.orfm_command_line(min_orf_length, restrict_read_length)
            cmd = '%s %s %s ' % (orfm_cmdline, input_path, tee)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])
        
        elif seq_type == 'protein': # If the input is amino acid sequence
            if input_file_format == FORMAT_FASTQ_GZ: # If its gzipped
                cmd = "awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s) %s" % (input_path, tee) # Unzip it and feed it into the base command
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])
            elif input_file_format == FORMAT_FASTA: # If it is in fasta format
                cmd = "cat %s %s" % (input_path, tee) # It can be searched directly, no manpulation required
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])
            else:
                raise Exception('Programming Error: error guessing input file format')
        else:
            raise Exception('Programming Error: error guessing input sequence type')
        return output_table_list

    def nhmmer(self, output_path, input_path, input_file_format, threads, eval, cmd_log):
        ## Run a nhmmer search on input_path file and return the name of
        ## resultant output table. Keep log of command.
        output_table_list = []
        tee = ''
        hmm_number = len(self.search_hmm)
        if hmm_number > 1:
            for idx, hmm in enumerate(self.search_hmm):
                out = os.path.join(os.path.split(output_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(output_path)[1])
                output_table_list.append(out)
                if idx + 1 == hmm_number:
                    tee += " | nhmmer %s --cpu %s --tblout %s %s - >/dev/null " % (eval, threads, out, hmm)
                elif idx + 1 < hmm_number:
                    tee += " >(nhmmer %s --cpu %s --tblout %s %s - >/dev/null) " % (eval, threads, out, hmm)
                else:
                    raise Exception("Programming Error.")    
        elif hmm_number == 1:
            tee = ' | nhmmer %s --cpu %s --tblout %s %s - >/dev/null' % (eval, threads, output_path, self.search_hmm[0])
            output_table_list.append(output_path) 
        if input_file_format == FORMAT_FASTA:
            cmd = "cat %s | tee %s" % (input_path, tee)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])

        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s) | tee %s " % (input_path, tee)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])
        else:
            raise Exception(Messenger().message('ERROR: Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (input_path)))
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
                
                if read_name in seen.keys(): # If the read name has been seen before.. 
                    if seen[read_name]==idx:
                        hash[read_name].append(buildHash(hit, program))
                else:
                    hash[read_name]=[buildHash(hit, program)]
                seen[read_name]=idx
        return hash
        
    def check_euk_contamination(self, output_path, euk_free_output_path, input_path, run_stats, input_file_format, threads, evalue, raw_reads, base, cmd_log, euk_hmm):
        reads_with_better_euk_hit = []
        reads_unique_to_eukaryotes = []
        cutoff = float(0.9*run_stats['read_length'])
        # do a nhmmer using a Euk specific hmm
        nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s %s" % (threads, evalue, output_path, euk_hmm)

        if input_file_format == FORMAT_FASTA:
            cmd = "%s %s 2>&1 > /dev/null" % (nhmmer_cmd, raw_reads)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell = True)

        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s )) 2>&1 > /dev/null" %  (nhmmer_cmd, raw_reads)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])

        else:
            raise Exception(Messenger().error_message('Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (raw_reads)))


        # check for evalues that are lower, after eliminating hits with an
        # alignment length of < 90% the length of the whole read.
        euk_reads = self.hmmtable_reader([output_path])
        euk_crossover = [x for x in euk_reads.keys() if x in run_stats['reads'].keys()]
        reads_unique_to_eukaryotes = [x for x in euk_reads.keys() if x not in run_stats['reads'].keys()]
        
        for entry in euk_crossover: # for every cross match
            if euk_reads[entry][0]['bit'] >= float(run_stats['reads'][entry]['bit']):
                if euk_reads[entry][0]['len'] > cutoff:
                    reads_with_better_euk_hit.append(entry)
                elif euk_reads[entry][0]['len'] < cutoff:
                    continue
            else:
                continue

        # Return Euk contamination
        if len(reads_with_better_euk_hit) == 0:
            Messenger().message("No contaminating eukaryotic reads detected in %s" % (os.path.basename(raw_reads)))

        else:
            Messenger().message("Found %s read(s) that may be eukaryotic" % len(reads_with_better_euk_hit + reads_unique_to_eukaryotes))
        # Write a file with the Euk free reads.
        with open(euk_free_output_path, 'w') as euk_free_output:
            for record in list(SeqIO.parse(open(input_path, 'r'), 'fasta')):
                if record.id not in reads_with_better_euk_hit:
                    SeqIO.write(record, euk_free_output, "fasta")
        run_stats['euk_uniq'] = len(reads_unique_to_eukaryotes)
        run_stats['euk_contamination'] = len(reads_with_better_euk_hit)
        return run_stats, euk_free_output_path

    def filter_hmmsearch(self, output_hash, contents, args, input_file_format, cmd_log):
        for seq_file in sequence_file_list:
            hmmout_table_title = suffix[0]
            table_title_list.append(hmmout_table_title)
            hmmsearch_cmd = " hmmsearch --cpu %s %s -o /dev/null --domtblout %s %s " % (threads, eval, hmmout_table_title, self.hmm)
            # TODO: capture stderr and report if the check_call fails
            if input_file_format == FORMAT_FASTA or input_file_format == FORMAT_FASTQ_GZ:
                if contents.pipe == 'P':
                    cmd = 'orfm %s | %s /dev/stdin' % (seq_file, hmmsearch_cmd)
                    self.hk.add_cmd(cmd_log, cmd)
                    subprocess.check_call(["/bin/bash", "-c", cmd])
            else:
                Messenger().message('ERROR: Suffix on %s not recegnised\n' % (seq_file))
                exit(1)
            del suffix[0]
        return table_title_list

    def csv_to_titles(self, output_path, input_path, run_stats):
        ## process hmmsearch/nhmmer results into a list of titles to <base_filename>_readnames.txt
        run_stats['reads'] = self.hmmtable_reader(input_path)
        count=sum([len(x) for x in run_stats['reads'].values()])

        # See if there are any reads in there reverse direction. Store True if so for later reference

        try:
            if any([x for x in sum(run_stats['reads'].values(), []) if x['direction'] =='-']):
                run_stats['rev_true'] = True
            else:
                run_stats['rev_true'] = False
        except KeyError:
            run_stats['rev_true'] = False

        if count > 0: # Return if there weren't any reads found
            Messenger().message('%s read(s) found' % (count))
        else: # Otherwise, report the number of reads
            Messenger().message('%s reads found, cannot continue with no information' % (len(run_stats['reads'].keys())))
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

    def extract_from_raw_reads(self, output_path, input_path, raw_sequences_path, input_file_format, cmd_log, read_stats):
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
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)
        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "%s -z %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            self.hk.add_cmd(cmd_log, cmd)
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
                corrected_sequences['>'+sequence.id+'\n'] = ''.join(new_seq)+'\n'
        with open(output_file_name, 'w') as output_file: # Create an open file to write the new sequences to
                for fasta_id, fasta_seq in corrected_sequences.iteritems():
                    if any(c.isalpha() for c in fasta_seq):
                        output_file.write(fasta_id)
                        output_file.write(fasta_seq)
      
    def orfm_command_line(self, min_orf_length, restrict_read_length):
        '''Return a string to run OrfM with, assuming sequences are incoming on
        stdin'''
        if restrict_read_length:
            orfm_arg_l = " -l %d" % restrict_read_length
        else:
            orfm_arg_l = ''
        
        return 'orfm -m %d %s ' % (min_orf_length, orfm_arg_l)

    def extract_orfs(self, input_path, raw_orf_path, hmmsearch_out_path, orf_titles_path, min_orf_length, restrict_read_length, orf_out_path, cmd_log):
        'Extract only the orfs that hit the hmm, return sequence file with within.'        
        # Build the command
        output_table_list = []
        tee = ' | tee'
        hmm_number = len(self.search_hmm)
        if hmm_number > 1:
            for idx, hmm in enumerate(self.search_hmm):
                out = os.path.join(os.path.split(hmmsearch_out_path)[0], os.path.basename(hmm).split('.')[0] +'_'+ os.path.split(hmmsearch_out_path)[1])
                output_table_list.append(out)
                if idx + 1 == hmm_number:
                    tee += " | hmmsearch --domtblout %s %s - >/dev/null " % (out, hmm)
                elif idx + 1 < hmm_number:
                    tee += " >(hmmsearch --domtblout %s %s - >/dev/null) " % (out, hmm)
                else:
                    raise Exception("Programming Error.") 
        elif hmm_number == 1:
            tee = ' | hmmsearch --domtblout %s %s - >/dev/null' % (hmmsearch_out_path, self.search_hmm[0])
            output_table_list.append(hmmsearch_out_path)
        # Call orfs on the sequences

        orfm_cmd = self.orfm_command_line(min_orf_length, restrict_read_length)
        cmd = '%s %s > %s' % (orfm_cmd, input_path, raw_orf_path)
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)

        cmd = 'cat %s %s' % (raw_orf_path, tee)
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(['bash','-c',cmd])
        
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
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)
        
        # Return name of output file
        return orf_out_path

    def p_search(self, files, args, run_stats, base, input_file_format, raw_reads):
        # Main pipe of search step in protein pipeline:
        # recieves reads, and returns hits
        start = timeit.default_timer() # Start search timer
        # Searching raw reads with HMM
        hit_table = self.hmmsearch(files.hmmsearch_output_path(base),
                                   raw_reads,
                                   input_file_format,
                                   args.input_sequence_type,
                                   args.threads,
                                   args.eval,
                                   args.min_orf_length,
                                   args.restrict_read_length,
                                   files.command_log_path())
        # Processing the output table to give you the readnames of the hits
        run_stats, hit_readnames = self.csv_to_titles(files.readnames_output_path(base),
                                                      hit_table,
                                                      run_stats)

        if not hit_readnames:
            return False, run_stats
        # Extract the hits form the original raw read file
        run_stats['reads'], hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                                    hit_readnames,
                                                                    raw_reads,
                                                                    input_file_format,
                                                                    files.command_log_path(),
                                                                    run_stats['reads'])
        
        if args.input_sequence_type == 'nucleotide':
            # Extract the orfs of these reads that hit the original search
            hit_orfs = self.extract_orfs(hit_reads,
                                         files.orf_output_path(base),
                                         files.orf_hmmsearch_output_path(base),
                                         files.orf_titles_output_path(base),
                                         args.min_orf_length,
                                         args.restrict_read_length,
                                         files.orf_fasta_output_path(base),
                                         files.command_log_path())
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

    def d_search(self, files, args, run_stats, base, input_file_format, raw_reads, euk_check):
        # Main pipe of search step in nucleotide pipeline:
        # recieves reads, and returns hits
        start = timeit.default_timer() # Start search timer

        # First search the reads using the HMM
        hit_table = self.nhmmer(files.hmmsearch_output_path(base),
                                raw_reads,
                                input_file_format,
                                args.threads,
                                args.eval,
                                files.command_log_path())

        # Next, get a list of readnames
        run_stats, hit_readnames = self.csv_to_titles(files.readnames_output_path(base),
                                                      hit_table,
                                                      run_stats)
        if not hit_readnames:
            return False, run_stats

        # And extract them from the original sequence file
        run_stats['reads'], hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                                    hit_readnames,
                                                                    raw_reads,
                                                                    input_file_format,
                                                                    files.command_log_path(),
                                                                    run_stats['reads'])
        # Define the read length
        run_stats['read_length'] = self.check_read_length(hit_reads, "D")

        # Stop timing search and start timing euk check step.
        stop = timeit.default_timer()
        run_stats['search_t'] = str(int(round((stop - start), 0)) )
        start = timeit.default_timer()

        # Check for Eukarytoic contamination
        if euk_check:
            Messenger().message("Checking for Eukaryotic contamination")
            run_stats, hit_reads = self.check_euk_contamination(files.euk_contam_path(base),
                                                                files.euk_free_path(base),
                                                                hit_reads,
                                                                run_stats,
                                                                input_file_format,
                                                                args.threads,
                                                                args.eval,
                                                                raw_reads,
                                                                base,
                                                                files.command_log_path(),
                                                                args.euk_hmm_file)

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
                      files.command_log_path(),
                      files.output_for_path(base),
                      files.output_rev_path(base),
                      files.sto_for_output_path(base),
                      files.sto_rev_output_path(base),
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

