import logging
import re

class OrfM:
    def __init__(self, **kwargs):
        self.min_orf_length = kwargs.pop('min_orf_length',96)
        self.restrict_read_length = kwargs.pop('restrict_read_length',None)
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
            orfm_arg_l = " -l %d" % self.restrict_read_length
            
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
        return re.compile('^(\S+)_(\d+)_(\d)_(\d+)')

class ZcatOrfM(OrfM):
    '''Like OrfM, except separate the zcat and the orfm commands into two
    parts of a shell pipeline'''
    
    def command_line(self, input_path=None):
        if input_path is None:
            raise Exception("input_path required for ZcatOrfM")
        original = OrfM.command_line(self, input_path=None)
        return "zcat '%s' | %s" % (input_path, original)