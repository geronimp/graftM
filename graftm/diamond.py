from graftm.sequence_search_results import DiamondSearchResult
import tempfile
import extern
import os
from graftm.unpack_sequences import UnpackRawReads

class Diamond:
    def __init__(self, database, threads=None, evalue=None):
        self._database = database
        self._threads = threads
        self._evalue = evalue
        
    def run(self, input_sequence_file, input_sequence_type, daa_file_basename=None):
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
        if input_sequence_type == UnpackRawReads.PROTEIN_SEQUENCE_TYPE:
            cmd_list.append('blastp')
        elif input_sequence_type == UnpackRawReads.NUCLEOTIDE_SEQUENCE_TYPE:
            cmd_list.append('blastx')
        else:
            raise Exception("Programming error")
        
        basename = daa_file_basename
        if basename is None:
            with tempfile.NamedTemporaryFile(prefix='graftm_diamond') as t:
                # we are just stealing the name, don't need the file itself
                basename = t.name
            
        for c in ['-k 1',
                  "-d",
                    self._database,
                    "-q",
                    "%s" % input_sequence_file,
                    "-a",
                    basename]:
            cmd_list.append(c)
        if self._threads:
            cmd_list.append("--threads")
            cmd_list.append(str(self._threads))
        if self._evalue:
            cmd_list.append("--evalue")
            cmd_list.append(str(self._evalue))

        cmd = ' '.join(cmd_list)
        extern.run(cmd)
        
        daa_name = "%s.daa" % basename
        res = DiamondSearchResult.import_from_daa_file(daa_name)
        
        if daa_file_basename is None:
            # Diamond makes an extra file, need to remove this
            os.remove(daa_name)
            
        return res
            
            
            
                