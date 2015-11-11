import os

class GraftMFiles:
    
    def __init__(self, old_title, outdir, direction):
        if direction == 'forward' or direction == 'reverse':
            self.basename = os.path.join(direction, old_title + '_' + direction)
        elif direction == False:
            self.basename = old_title
        else:
            raise Exception('Programming Error.')   
         
        self.outdir = outdir
    
    def search_otu_table(self):
        return os.path.join(self.outdir, "search_otu_table.txt")
    
    def hmmsearch_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s.hmmout.csv" % self.basename)
    
    def diamond_search_output_basename(self, out_path):
        '''Does not include the .daa part that diamond creates'''
        return os.path.join(self.outdir, out_path, "%s_diamond_search" % self.basename)
    
    def diamond_assignment_output_basename(self, out_path):
        '''Does not include the .daa part that diamond creates'''
        return os.path.join(self.outdir, out_path, "%s_diamond_assignment" % out_path)
    
    def read_tax_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_read_tax.tsv" % out_path)

    def jplace_output_path(self):
        return "placements.jplace"
    
    def euk_free_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_euk_free.fa" % self.basename)
    
    def euk_contam_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_euk_contam.txt" % self.basename)
    
    def summary_table_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_count_table.txt" % self.basename)
    
    def krona_output_path(self):
        return os.path.join(self.outdir, "krona.html")
    
    def aligned_fasta_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_hits.aln.fa" % self.basename)

    def orf_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_orf" % self.basename)

    def orf_titles_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_orf.titles" % self.basename)

    def orf_fasta_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_orf.fa" % self.basename)

    def conv_output_for_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_conv_for.faa" % self.basename)
    
    def output_for_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_for.faa" % self.basename)
    
    def output_rev_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_rev.faa" % self.basename)
    
    def conv_output_rev_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_conv_rev.faa" % self.basename)
    
    def comb_aln_fa(self):
        return os.path.join(self.outdir, "combined_alignment.aln.fa")

    def fa_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_hits.fa" % self.basename)     
        
    def readnames_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_readnames.txt" % self.basename)
    
    def sto_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s.sto" % self.basename)

    def basic_stats_path(self):
        return os.path.join(self.outdir, "basic_stats.txt")

    def for_aln_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_for_aln.fa" % self.basename)
        
    def rev_aln_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_rev_aln.fa" % self.basename)
    
    def combined_biom_output_path(self):
        return os.path.join(self.outdir, "graftm.biom")
    
    def combined_summary_table_output_path(self):
        return os.path.join(self.outdir, "combined_count_table.txt")
    
    def bootstrap_hmm_path(self):
        return os.path.join(self.outdir, "bootstrap.hmm")
      
    def base(self, out_path):
        return os.path.join(self.outdir, out_path, "%s" % self.basename)