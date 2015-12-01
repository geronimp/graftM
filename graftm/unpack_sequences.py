import logging
import subprocess
from signal import signal, SIGPIPE, SIG_DFL
import os
import itertools
from string import lower

class UnpackRawReads:
    class UnexpectedFileFormatException(Exception): pass
    
    FORMAT_FASTA    = "FORMAT_FASTA"
    FORMAT_FASTQ    = "FORMAT_FASTQ"
    FORMAT_FASTQ_GZ = "FORMAT_FASTQ_GZ"
    FORMAT_FASTA_GZ = "FORMAT_FASTA_GZ"
    
    PROTEIN_SEQUENCE_TYPE = 'aminoacid'
    NUCLEOTIDE_SEQUENCE_TYPE = 'nucleotide'
    
    _EXTENSION_TO_FILE_TYPE = {'.fa': FORMAT_FASTA,
                               '.faa': FORMAT_FASTA,
                               '.fna': FORMAT_FASTA,
                               '.fasta': FORMAT_FASTA,
                               
                               '.fq': FORMAT_FASTQ,
                               '.fastq': FORMAT_FASTQ,
                               
                               '.fq.gz': FORMAT_FASTQ_GZ,
                               '.fastq.gz': FORMAT_FASTQ_GZ,
                               
                               '.fa.gz': FORMAT_FASTA_GZ,
                               '.faa.gz': FORMAT_FASTA_GZ,
                               '.fna.gz': FORMAT_FASTA_GZ,
                               '.fasta.gz': FORMAT_FASTA_GZ,
                               }
    
    def __init__(self, read_file):
        self.read_file   = read_file
    
    def _guess_sequence_type_from_string(self, seq):
        '''Return 'protein' if there is >10% amino acid residues in the 
        (string) seq parameter, else 'nucleotide'. Raise Exception if a
        non-standard character is encountered'''
        # Define expected residues for each sequence type
        aa_chars = ['P','V','L','I','M','F','Y','W','H','K','R','Q','N','E','D','S','X']
        aas = set(itertools.chain(aa_chars, [lower(a) for a in aa_chars]))
        na_chars = ['A', 'T', 'G', 'C', 'N', 'U']
        nas = set(itertools.chain(na_chars, [lower(a) for a in na_chars]))
        
        num_nucleotide = 0
        num_protein = 0
        count = 0
        for residue in seq:
            if residue in nas:
                num_nucleotide += 1
            elif residue in aas:
                num_protein += 1
            else:
                raise Exception('Encountered unexpected character when attempting to guess sequence type: %s' % (residue))
            count += 1
            if count >300: break
        if float(num_protein) / (num_protein+num_nucleotide) > 0.1:
            return self.PROTEIN_SEQUENCE_TYPE
        else:
            return self.NUCLEOTIDE_SEQUENCE_TYPE
    
    def sequence_type(self):
        '''Guess the type of input sequence provided to graftM (i.e. nucleotide
        or amino acid) and return'''
        if hasattr(self, "type"):
            return self.type
        else:
            # If its Gzipped and fastq make a small sample of the sequence to be 
            # read
            cmd='%s | head -n 2' % (self.command_line())
            first_seq = subprocess.check_output(
                                                cmd, 
                                                shell=True,
                                                preexec_fn=lambda:signal(SIGPIPE, SIG_DFL)
                                                )
            _, seq = tuple(first_seq.strip().split('\n'))
            self.type = self._guess_sequence_type_from_string(seq)
            logging.debug("Detected sequence type as %s" % self.type)
            return self.type

    def guess_sequence_input_file_format(self, sequence_file_path):
        '''Given a sequence file, guess the format and return. Raise an 
        exception if it cannot be guessed'''
        return self._EXTENSION_TO_FILE_TYPE[self._get_extension(sequence_file_path)]
    
    def format(self):
        return self.guess_sequence_input_file_format(self.read_file)
    
    def is_zcattable(self):
        return self.guess_sequence_input_file_format(self.read_file) in \
            (self.FORMAT_FASTA_GZ, self.FORMAT_FASTQ_GZ)
            
    def basename(self):
        '''Return the name of the file with the '.fasta' or 'fq.gz' etc 
        removed'''
        return os.path.basename(self.read_file)[:-len(self._get_extension(self.read_file))]
        
    def _get_extension(self, sequence_file_path):
        for ext in self._EXTENSION_TO_FILE_TYPE.keys():
            if sequence_file_path.endswith(ext): return ext
        raise self.UnexpectedFileFormatException("Unable to guess file format of sequence file: %s" % sequence_file_path)
    
    def command_line(self):
        '''Return a string to open read files with'''
        file_format=self.guess_sequence_input_file_format(self.read_file)
        logging.debug("Detected file format %s" % file_format)
        if file_format == self.FORMAT_FASTA:
            cmd="""cat %s""" % (self.read_file)
        elif file_format == self.FORMAT_FASTQ_GZ:
            cmd="""zcat %s | awk '{print ">" substr($0,2);getline;print;getline;getline}' -""" % (self.read_file)
        elif file_format == self.FORMAT_FASTA_GZ:
            cmd="""zcat %s""" % (self.read_file)
        elif file_format == self.FORMAT_FASTQ:
            cmd="""awk '{print ">" substr($0,2);getline;print;getline;getline}' %s""" % (self.read_file)
        logging.debug("raw read unpacking command chunk: %s" % cmd)
        return cmd 
