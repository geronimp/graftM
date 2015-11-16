import subprocess
import logging
import csv
import os

class SequenceSearchResult:
    QUERY_FROM_FIELD = 'query_from'
    QUERY_TO_FIELD = 'query_to'
    QUERY_LENGTH_FIELD = 'query_length'
    HIT_FROM_FIELD = 'hit_from'
    HIT_TO_FIELD = 'hit_to'
    ALIGNMENT_LENGTH_FIELD = 'alignment_length'
    ALIGNMENT_BIT_SCORE = 'alignment_bit_score'
    ALIGNMENT_DIRECTION = 'alignment_direction'
    HIT_ID_FIELD = 'hit_id'
    QUERY_ID_FIELD = 'query_id'
    HMM_NAME_FIELD = 'hmm_name'
    ACCESSION_ID_FIELD = 'accession_id'
    PERCENT_ID_FIELD = 'percent_id'
    MISMATCH_FIELD = "mismatch"
    EVALUE_FIELD = "evalue"
    
    
    
    def __init__(self):
        self.fields = []
        self.results = []
        
    def each(self, field_names):
        """Iterate over the results, yielding a list for each result, where
        each element corresponds to the field given in the field_name parameters
        
        Parameters
        ----------
        field_names: list of str
            The names of the fields to be returned during iteration
            
        Returns
        -------
        None
        
        Exceptions
        ----------
        raises something when a field name is not in self.fields
        """
        field_ids = []
        for f in field_names:
            # below raises error if the field name is not found, so
            # don't need to account for that.
            field_ids.append(self.fields.index(f))
        
        for r in self.results:
            yield([r[i] for i in field_ids])
        
class DiamondSearchResult(SequenceSearchResult):
    @staticmethod
    def import_from_daa_file(daa_filename):
        '''Generate new results object from the output of diamond blastx/p'''
        
        # blast m8 format is
        # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        res = DiamondSearchResult()
        res.fields = [
                       SequenceSearchResult.QUERY_ID_FIELD,
                       SequenceSearchResult.HIT_ID_FIELD,
                       SequenceSearchResult.PERCENT_ID_FIELD,
                       SequenceSearchResult.ALIGNMENT_LENGTH_FIELD,
                       SequenceSearchResult.MISMATCH_FIELD,
                       #skip
                       SequenceSearchResult.QUERY_FROM_FIELD,
                       SequenceSearchResult.QUERY_TO_FIELD,
                       SequenceSearchResult.HIT_FROM_FIELD,
                       SequenceSearchResult.HIT_TO_FIELD,
                       SequenceSearchResult.EVALUE_FIELD,
                       SequenceSearchResult.ALIGNMENT_BIT_SCORE,
                       # extras
                       SequenceSearchResult.ALIGNMENT_DIRECTION,
                       SequenceSearchResult.HMM_NAME_FIELD
                       ]
        
        cmd = "diamond view -a '%s'" % daa_filename
        logging.debug("Running cmd: %s" % cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        
        reader = csv.reader(stdout.decode('ascii').splitlines(),
                                delimiter='\t')
        if process.returncode != 0:
            raise Exception("Problem running diamond view with cmd: '%s'," 
                            "stderr was %s" % (cmd, stderr))

        for row in reader:
            # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            #    0       1     2      3        4        5      6     7    8      9    10     11
            query_start = int(row[6])
            query_end = int(row[7])
            res.results.append([row[0],
                                 row[1],
                                 row[2],
                                 row[3],
                                 row[4],
                                 query_start,
                                 query_end,
                                 int(row[8]),
                                 int(row[9]),
                                 row[10],
                                 row[11],
                                 query_start < query_end,
                                 os.path.basename(daa_filename)
                                 ])
        return res

class HMMSearchResult(SequenceSearchResult):
    @staticmethod
    def import_from_nhmmer_table(hmmout_path):
        '''Generate new results object from the output of nhmmer search'''
        # nhmmer format is
        # qseqid queryname hmmfrom hmmto alifrom alito envfrom envto sqlen strand evalue bitscore bias description
        #   0        2        4      5      6      7      8      9    10     11     12      13     14     15     
        res=HMMSearchResult()
        res.fields = [
                       SequenceSearchResult.QUERY_ID_FIELD,
                       SequenceSearchResult.HMM_NAME_FIELD,
                       SequenceSearchResult.ALIGNMENT_LENGTH_FIELD,
                       SequenceSearchResult.QUERY_FROM_FIELD,
                       SequenceSearchResult.QUERY_TO_FIELD,
                       SequenceSearchResult.HIT_FROM_FIELD,
                       SequenceSearchResult.HIT_TO_FIELD,
                       SequenceSearchResult.ALIGNMENT_BIT_SCORE,
                       SequenceSearchResult.ALIGNMENT_DIRECTION,
                       ]
        
        for row in [x.rstrip().split() for x in open(hmmout_path) if not x.startswith('#')]:
            alifrom    = int(row[6])
            alito      = int(row[7])
            aln_length = (alito-alifrom if alito-alifrom>0 else alifrom-alito)
            res.results.append([row[0],
                                row[2],
                                aln_length,
                                int(row[4]),
                                int(row[5]),
                                alifrom,
                                alito,
                                row[13],
                                alito > alifrom
                                ])
        return res
    @staticmethod
    def import_from_hmmsearch_table(hmmout_path):
        '''Generate new results object from the output of hmmsearch search'''
        # hmmsearch format is
        # qseqid tlen queryname qlen evalue bitscore bias hmmfrom hmmto alifrom alito envfrom envto acc
        #    0    2       3      5     6        7     8      15    16     17     18      19    20   21
        res=HMMSearchResult()
        res.fields = [
                       SequenceSearchResult.QUERY_ID_FIELD,
                       SequenceSearchResult.HMM_NAME_FIELD,
                       SequenceSearchResult.ACCESSION_ID_FIELD,
                       SequenceSearchResult.QUERY_LENGTH_FIELD,
                       SequenceSearchResult.ALIGNMENT_LENGTH_FIELD,
                       SequenceSearchResult.QUERY_FROM_FIELD,
                       SequenceSearchResult.QUERY_TO_FIELD,
                       SequenceSearchResult.HIT_FROM_FIELD,
                       SequenceSearchResult.HIT_TO_FIELD,
                       SequenceSearchResult.ALIGNMENT_BIT_SCORE,
                       SequenceSearchResult.ALIGNMENT_DIRECTION,
                       ]
        
        for row in [x.rstrip().split() for x in open(hmmout_path) if not x.startswith('#')]:
            alifrom    = int(row[17])
            alito      = int(row[18])
            aln_length = (alito-alifrom if alito-alifrom>0 else alifrom-alito)
            if alito != alifrom: #this actually happens..
                res.results.append([row[0],
                                    row[3],
                                    row[4],
                                    row[5],
                                    aln_length,
                                    int(row[15]),
                                    int(row[16]),
                                    alifrom,
                                    alito,
                                    row[7],
                                    True
                                    ])
        return res
        
        
        
        
        
        
        
        
         