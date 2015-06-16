class OrfM:
    def __init__(self, **kwargs):
        self.min_orf_length = kwargs.pop('min_orf_length',None)
        self.restrict_read_length = kwargs.pop('restrict_read_length',None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        self.min_orf_length = 96
        
    def command_line(self):
        '''Return a string to run OrfM with, assuming sequences are incoming on
        stdin and printed to stdout'''
        
        if self.min_orf_length:
            orfm_arg_l = " -m %d" % self.min_orf_length
        else:
            orfm_arg_l = ''
        
        if self.restrict_read_length:
            orfm_arg_l = " -l %d" % self.restrict_read_length
            
        return 'orfm %s ' % orfm_arg_l
