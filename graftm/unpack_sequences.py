import logging
import subprocess
import os
import itertools
import tempfile
import extern

from Bio import SeqIO
from signal import signal, SIGPIPE, SIG_DFL
from string import lower

from graftm.sequence_io import SequenceIO


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
        self.gapped_alignment = False
        self.file_type = self._guess_sequence_input_file_format(self.read_file)
        self.type = self._guess_sequence_type(self.read_file)
        
    def _ungap_sequences(self):
        logging.debug("Removing gaps from the input sequences")
        if self.format()==self.FORMAT_FASTA:
            reads=SeqIO.parse(open(self.read_file), "fasta")
        elif self.format()==self.FORMAT_FASTQ:
            reads=SeqIO.parse(open(self.read_file), "fastq")
        else:
            logging.warning("If this is a large sequence file this is going to take some time - Did I mention inputting alignments to GraftM is probably a bad idea?")
            with tempfile.NamedTemporaryFile(prefix='gapped_', suffix='.fa') as f_out:
                print self.format()
                cmd = '%s > %s ' % (self.command_line(), f_out.name)
                extern.run(cmd)
            reads=SeqIO.parse(open(input_fasta), "fasta")
        t_file = tempfile.NamedTemporaryFile(prefix='ungapped_', 
                                             suffix='.fa',
                                             delete=False)
        with open(t_file.name, 'w') as out:
            for record in reads:
                out.write(">%s\n%s\n" % (record.description,
                                         str(record.seq).translate(None, "*-")))
        self.read_file = t_file.name
                
    
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
            elif residue == '-':
                if self.gapped_alignment:
                    continue
                else:
                    logging.warning("Input sequences are gapped ('-' detected in sequences). These gaps will be removed, but the resultant reads may not be identical to the originals (if insertions and deletions have been removed, etc..)")
                    self.gapped_alignment = True
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
        return self.type

    def _guess_sequence_type(self, sequence_file_path):
        # If its Gzipped and fastq make a small sample of the sequence to be 
        # read
        cmd='%s | head -n 100' % (self.command_line())
        first_seq = subprocess.check_output(
                                            cmd, 
                                            shell=True,
                                            preexec_fn=lambda:signal(SIGPIPE, SIG_DFL)
                                            )
        
        seq = ''
        for line in first_seq.strip().split('\n'):
            if not line.startswith('>'):
                seq+=line
                
        type = self._guess_sequence_type_from_string(seq)
        logging.debug("Detected sequence type as %s" % type)    
        if self.gapped_alignment:
            self._ungap_sequences()
        
        return type
    
    def _guess_sequence_input_file_format(self, sequence_file_path):
        '''Given a sequence file, guess the format and return. Raise an 
        exception if it cannot be guessed'''
        return self._EXTENSION_TO_FILE_TYPE[self._get_extension(sequence_file_path)]
    
    def format(self):
        return self.file_type
    
    def is_zcattable(self):
        if self.file_type == self.FORMAT_FASTA_GZ or self.file_type == self.FORMAT_FASTQ_GZ:
            return True
        else:
            return False
            
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
        logging.debug("Detected file format %s" % self.file_type)
        if self.file_type == self.FORMAT_FASTA:
            cmd="""cat %s""" % (self.read_file)
        elif self.file_type == self.FORMAT_FASTQ_GZ:
            cmd="""zcat %s | awk '{print ">" substr($0,2);getline;print;getline;getline}' -""" % (self.read_file)
        elif self.file_type == self.FORMAT_FASTA_GZ:
            cmd="""zcat %s""" % (self.read_file)
        elif self.file_type == self.FORMAT_FASTQ:
            cmd="""awk '{print ">" substr($0,2);getline;print;getline;getline}' %s""" % (self.read_file)
        logging.debug("raw read unpacking command chunk: %s" % cmd)
        return cmd 
