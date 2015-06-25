import logging

FORMAT_FASTA    = "FORMAT_FASTA"
FORMAT_FASTQ    = "FORMAT_FASTQ"
FORMAT_FASTQ_GZ = "FORMAT_FASTQ_GZ"
FORMAT_FASTA_GZ = "FORMAT_FASTA_GZ"

class UnpackRawReads:
    def __init__(self, read_file):
        self.read_file   = read_file
                 
    def guess_sequence_input_file_format(self, sequence_file_path):
        '''Given a sequence file, guess the format and return. Raise an 
        exception if it cannot be guessed'''
        if sequence_file_path.endswith(('.fa', '.faa', '.fna', '.fasta')): 
            return FORMAT_FASTA
        if sequence_file_path.endswith(('.fq', '.fastq')):  
            return FORMAT_FASTQ
        elif sequence_file_path.endswith(('.fq.gz', '.fastq.gz')):
            return FORMAT_FASTQ_GZ
        elif sequence_file_path.endswith(('.fa.gz', '.faa.gz', '.fna.gz', '.fasta.gz')):
            return FORMAT_FASTA_GZ
        else:
            raise Exception("Unable to guess file format of sequence file: %s" % sequence_file_path)
    
    def format(self):
        return self.guess_sequence_input_file_format(self.read_file)
    
    def is_zcattable(self):
        return self.guess_sequence_input_file_format(self.read_file) in \
            (FORMAT_FASTA_GZ, FORMAT_FASTQ_GZ)
    
    def command_line(self):
        '''Return a string to open read files with'''
        file_format=self.guess_sequence_input_file_format(self.read_file)
        logging.debug("Detected file format %s" % file_format)
        if file_format == FORMAT_FASTA:
            cmd="""cat %s""" % (self.read_file)
        elif file_format == FORMAT_FASTQ_GZ:
            cmd="""zcat %s | awk '{print ">" substr($0,2);getline;print;getline;getline}' -""" % (self.read_file)
        elif file_format == FORMAT_FASTA_GZ:
            cmd="""zcat %s""" % (self.read_file)
        elif file_format == FORMAT_FASTQ:
            cmd="""awk '{print ">" substr($0,2);getline;print;getline;getline}' %s""" % (self.read_file)
        logging.debug("raw read unpacking command chunk: %s" % cmd)
        return cmd 
