import logging
from graftm.sequence_search_results import DiamondSearchResult
import tempfile
import subprocess
import os

class Diamond:
    def __init__(self, database, threads=None, evalue=None):
        self._database = database
        self._threads = threads
        self._evalue = evalue
        
    def run(self, input_sequence_file, input_sequence_type):
        '''Run input sequences in either blastp or blastx mode against the
        database specified in __init__.
            
        Parameters
        ----------
        input_sequence_file: str
            path to query sequences
        input_sequence_type: either 'nucleotide' or 'protein'
            the input_sequences are this kind of sequence
            
        Returns
        -------
        DiamondSearchResult
        '''
        
        cmd_list = ["diamond"]
        if input_sequence_type == 'protein':
            cmd_list.append('blastp')
        elif input_sequence_type == 'nucleotide':
            cmd_list.append('blastx')
        else:
            raise Exception("Programming error")
        
        with tempfile.NamedTemporaryFile(prefix='graftm_diamond') as t:
            for c in ['-k 1',
                      "-d",
                        self._database,
                        "-q",
                        "'%s'" % input_sequence_file,
                        "-a",
                        t.name]:
                cmd_list.append(c)
            if self._threads:
                cmd_list.append("--threads")
                cmd_list.append(str(self._threads))
            if self._evalue:
                cmd_list.append("--evalue")
                cmd_list.append(str(self._evalue))
                
            cmd = ' '.join(cmd_list)
            logging.debug("Running cmd: %s" % cmd)
            subprocess.check_call(cmd, shell=True)
            
            daa_name = "%s.daa" % t.name
            res = DiamondSearchResult.import_from_daa_file(daa_name)
            
            # Diamond makes an extra file, need to remove this specifically,
            # not just the tempfile
            os.remove(daa_name)
        return res
            
            
            
                