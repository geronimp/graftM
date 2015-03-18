#!/usr/bin/env python
__author__ = "Joel Boyd, Ben Woodcroft"
__copyright__ = "Copyright 2014"
__credits__ = ["Joel Boyd", "Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd, Ben Woodcroft"
__email__ = "joel.boyd near uq.net.au, b.woodcroft near uq.edu.au"
__status__ = "Development"
__version__ = "0.4.0"

from graftm.Messenger import Messenger
from graftm.assembler import TaxoGroup
from Bio import SeqIO
import IPython
class Extract:
    
    def __init__(self):
        self.TG = TaxoGroup()
        
    def extract(self, args):
            if args.file:
                lineages = [x.rstrip() for x in open(args.seq, 'r').readlines()]
                Messenger().header("Acquiring reads for %s lineages" % len(lineages))
            else:
                lineages = [args.seq]
                Messenger().header("Acquiring reads for %s" % lineages[0])
        
            
            guppy = self.TG.guppy_splitter([args.profile + '/' + args.profile + '_placements.guppy'])
            Messenger().message("Parsing sequences")
            record_dict = SeqIO.to_dict(SeqIO.parse(args.profile + '/' + args.profile + '_hits.fa', "fasta"))            
            
            if not args.file:
            
                if args.non_cumil:
                    with open(args.profile+'_'+lineages[0]+'_Seqs.fa', 'w') as output:
                        for entry in guppy:
                            if guppy[entry][-1] == lineages[0]:
                                output.write('>' + entry[:-2] + '\n')
                                output.write(str(record_dict[entry[:-2]].seq) + '\n')
    
                
                elif not args.non_cumil:
                    with open(args.profile+'_'+lineages[0]+'_Seqs.fa', 'w') as output:
                        for entry in guppy:

                            if lineages[0] in guppy[entry]:
                                output.write('>' + entry[:-2] + '\n')
                                output.write(str(record_dict[entry[:-2]].seq) + '\n')
            elif args.file:
                for lin in lineages:
                    if args.non_cumil:
                        with open(args.profile+'_'+lin+'_Seqs.fa', 'w') as output:
                            for entry in guppy:
                                if guppy[entry][-1] == lin:
                                    output.write('>' + entry[:-2] + '\n')
                                    output.write(str(record_dict[entry[:-2]].seq) + '\n')
        
                    
                    elif not args.non_cumil:
                        with open(args.profile+'_'+lin+'_Seqs.fa', 'w') as output:
                            for entry in guppy:
        
                                if lin in guppy[entry]:
                                    output.write('>' + entry[:-2] + '\n')
                                    output.write(str(record_dict[entry[:-2]].seq) + '\n')