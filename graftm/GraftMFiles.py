#!/usr/bin/env python
import os
import IPython

class GraftMFiles:

    def __init__(self, old_title, outdir):
        self.basename = os.path.basename(old_title.split(',')[0]).split('.')[0]
        self.outdir = outdir
    
    def forward_read_hmmsearch_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_for.hmmout.csv" % self.basename)

    def jplace_output_path(self):
        return "placements.jplace"
    
    def euk_free_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_euk_free.fa" % self.basename)
    
    def euk_contam_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_euk_contam.txt" % self.basename)

    def guppy_file_output_path(self):
        return "_placements.guppy"

    def summary_table_output_path(self):
        return "_count_table.txt"
    
    def krona_output_path(self):
        return "_krona.html"
    
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
    
    def conv_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_conv.faa" % self.basename)
    
    def comb_aln_fa(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_comb.aln.fa" % self.basename)
    
    def conv_output_rev_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_conv_rev.faa" % self.basename)

    def reverse_read_hmmsearch_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_rev.hmmout.csv" % self.basename)

    def fa_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_hits.fa" % self.basename)     
        
    def readnames_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_readnames.txt" % self.basename)

    def sto_for_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s.for.sto" % self.basename)
    
    def sto_rev_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s.rev.sto" % self.basename)
    
    def sto_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s.sto" % self.basename)

    def orf_hmmsearch_output_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_orf.hmmout.csv" % self.basename)
    
    def basic_stats_path(self):
        return os.path.join(self.outdir, "basic_stats.txt")

    def command_log_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_commands_log.txt" % self.basename)

    def for_aln_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_for_aln.fa" % self.basename)
        
    def rev_aln_path(self, out_path):
        return os.path.join(self.outdir, out_path, "%s_rev_aln.fa" % self.basename)
    
    def coverage_table_path(self):
        return '_coverage_table.txt'
      
    def base(self, out_path):
        return os.path.join(self.outdir, out_path, "%s" % self.basename)