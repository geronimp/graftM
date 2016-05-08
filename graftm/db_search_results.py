
class DBSearchResult:
    '''
    DBSearchResult - Class for containing results from search pipeline in GraftM
    '''
    def __init__(self, output_reads, search_result, hit_read_count, has_slash_endings):
        self.output_reads   = output_reads
        self.search_result = search_result
        self.hit_count     = hit_read_count
        self.has_slash_endings = has_slash_endings
    
    def hit_fasta(self): # Return the path to the fasta file of hits
        return self.output_reads
    
    def search_objects(self): # Return SequenceSearchResult object with parameters defined
        return self.search_result
    
    def hit_count(self): # Return the number of hits found
        return self.hit_count  
    
