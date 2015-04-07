import os

from Bio import SeqIO

from graftm.messenger import Messenger
from graftm.assembler import TaxoGroup

class Extract:

    def __init__(self):
        self.tg = TaxoGroup()

    def extract(self, args):
            has_rev = False
            if args.file:
                lineages = [x.rstrip() for x in open(args.seq, 'r').readlines()]
                Messenger().header("Acquiring reads for %s lineages" % len(lineages))
            else:
                lineages = [args.seq]
                Messenger().header("Acquiring reads for %s" % lineages[0])
            if os.path.isdir( os.path.join(args.profile, 'forward') ) and os.path.isdir( os.path.join(args.profile, 'reverse') ):
                has_rev = True
                # Forward read guppy and sequences
                f_guppy = self.tg.guppy_splitter(os.path.join(args.profile, 'forward', 'placements.guppy'), '0.75')
                f_record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(args.profile, 'forward' , os.path.basename(args.profile)+'_forward_hits.fa'), "fasta"))
                # Reverse read guppy and sequences
                r_guppy = self.tg.guppy_splitter(os.path.join(args.profile, 'reverse', 'placements.guppy'), '0.75')
                r_record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(args.profile, 'reverse' , os.path.basename(args.profile)+'_reverse_hits.fa'), "fasta"))
            elif not os.path.isdir( os.path.join(args.profile, 'forward') ) and not os.path.isdir( os.path.join(args.profile, 'reverse') ):  
                guppy = self.tg.guppy_splitter(os.path.join(args.profile, 'placements.guppy'), '0.75')
                record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(args.profile, os.path.basename(args.profile)+'_hits.fa'), "fasta"))
            else:
                raise Exception("Programming Error: Incomplete or misformatted GraftM profile.")
            import IPython
            IPython.embed()

            
            Messenger().message("Parsing sequences")
            record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(args.profile, '*_hits.fa'), "fasta"))

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
