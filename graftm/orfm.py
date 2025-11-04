import logging
import re


class OrfM:
    def __init__(self, **kwargs):
        self.min_orf_length = kwargs.pop('min_orf_length',96)
        self.restrict_read_length = kwargs.pop('restrict_read_length',None)
        self.translation_table = kwargs.pop('translation_table')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
    def command_line(self, input_path=None):
        '''Return a string to run OrfM with, assuming sequences are incoming on
        stdin and printed to stdout
        
        Parameters
        ----------
        input_path: str
            path to the input path, or None for STDIN being the input
        '''
        
        if self.min_orf_length:
            orfm_arg_l = " -m %d" % self.min_orf_length
        else:
            orfm_arg_l = ''
        
        if self.restrict_read_length:
            orfm_arg_l += " -l %d" % self.restrict_read_length

        if self.translation_table:
            orfm_arg_l += " -c %d" % self.translation_table
            
        cmd = 'orfm %s ' % orfm_arg_l
        if input_path:
            cmd += input_path
        logging.debug("OrfM command chunk: %s" % cmd)
        return cmd
    
    @staticmethod
    def regular_expression():
        '''Return a compiled regular expression matching the OrfM output defline
        with the originalName, startPosition, frameNumber, orfNumber being
        the contents of the regular expression'''
        return re.compile(r'^(\S+)_(\d+)_(\d)_(\d+)')
