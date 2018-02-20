import extern
from Bio import SeqIO
from StringIO import StringIO

class SequenceExtractor:
    def extract(self, reads_to_extract, database_fasta_file, output_file):
        '''Extract the reads_to_extract from the database_fasta_file and put them in
        output_file.

        Parameters
        ----------
        reads_to_extract: Iterable of str
            IDs of reads to be extracted
        database_fasta_file: str
            path the fasta file that containing the reads
        output_file: str
            path to the file where they are put

        Returns
        -------
        Nothing'''
        cmd = "fxtract -XH -f /dev/stdin '%s' > %s" % (
            database_fasta_file, output_file)

        extern.run(cmd, stdin='\n'.join(reads_to_extract))

    def extract_forward_and_reverse_complement(
            self, forward_reads_to_extract, reverse_reads_to_extract, database_fasta_file,
            output_file):
        '''As per extract except also reverse complement the sequences.'''
        self.extract(forward_reads_to_extract, database_fasta_file, output_file)
        cmd_rev = "fxtract -XH -f /dev/stdin '%s'" % database_fasta_file

        output = extern.run(cmd_rev, stdin='\n'.join(reverse_reads_to_extract))

        with open(output_file, 'a') as f:
            for record in SeqIO.parse(StringIO(output), 'fasta'):
                record.seq = record.reverse_complement().seq
                SeqIO.write(record, f, 'fasta')
