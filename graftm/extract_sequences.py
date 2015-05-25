import os
import json
import logging

from Bio import SeqIO

from graftm.assembler import TaxoGroup
from graftm.classify import Classify

class Extract:

    def __init__(self):
        self.tg = TaxoGroup()

    def extract(self, args):
            if not os.path.isdir(args.profile):
                logging.error("Profile doesn't exist")
                exit(0)
            runs = [os.path.join(args.profile, f) for f in os.listdir(args.profile) if os.path.isdir(os.path.join(args.profile, f))]
            d = {}
            for r in runs:
                if len([f for f in os.listdir(r) if f == 'reverse' or f == 'forward']) > 0:
                    files = [os.path.join(f, 'forward') for f in runs] + [os.path.join(f, 'reverse') for f in runs]
                    prefix = ['forward', 'reverse']
                    for f in files:
                        p = prefix.pop(0)
                        d[os.path.basename(r) + '_' + p] = {'jplace': os.path.join(f, 'placements.jplace'),
                                                  'reads': os.path.join(f, os.path.basename(r) + '_%s_hits.fa' % p)}
                elif len([f for f in os.listdir(r) if f == 'reverse' or f == 'forward']) == 0:
                    d[os.path.basename(r)] = {'jplace': os.path.join(r, 'placements.jplace'),
                                              'reads': os.path.join(r, os.path.basename(r) + '_hits.fa')}           
            for run in d.keys():
                try:
                    reads = SeqIO.to_dict(SeqIO.parse(d[run]["reads"], "fasta"))
                    guppy=Classify('/srv/db/graftm/0/SILVA_p94_e91.gpkg/SILVA_p94_e91.refpkg/p94_e91_taxinfo.txt').assignPlacement(d[run]["jplace"], args.cutoff, 'reads')

                except:
                    logging.info("Cannot use unfinished run: %s" % run)
                    continue
                out_path = os.path.join(os.path.split(d[run]["reads"])[0], 'reads')
                if not os.path.isdir(out_path):
                    os.mkdir(out_path)
                
                if not args.file:                    
                    for lin in args.seq.split(','):
                        logging.info("Finding sequences for %s in %s" % (lin, run))
                        if args.non_cumil:
                            with open(os.path.join(out_path, run+'_'+lin+'_non_cumil_seqs.fa'), 'w') as output:
                                for entry in guppy:

                                    if guppy[entry]['placement'][-1] == lin:
                                        output.write('>' + entry + '\n')
                                        output.write(str(reads[entry].seq) + '\n')
                        elif not args.non_cumil:
                            with open(os.path.join(out_path, run+'_'+lin+'_seqs.fa'), 'w') as output:
                                for entry in guppy:
                                    import IPython
                                    IPython.embed()
                                    if lin in guppy[entry]['placement']:
                                        output.write('>' + entry + '\n')
                                        output.write(str(reads[entry].seq) + '\n')
                elif args.file:
                    for lin in [x.rstrip() for x in open(args.seq, 'r').readlines()]:
                        logging.info("Finding sequences for %s in %s" % (lin, run))
                        if args.non_cumil:
                            with open(os.path.join(out_path, run+'_'+lin+'_seqs.fa'), 'w') as output:
                                for entry in guppy:
                                    if guppy[entry][-1] == lin:
                                        output.write('>' + entry[:-2] + '\n')
                                        output.write(str(record_dict[entry[:-2]].seq) + '\n')
                        elif not args.non_cumil:
                            with open(os.path.join(out_path, run+'_'+lin+'_seqs.fa'), 'w') as output:
                                for entry in guppy:
    
                                    if lin in guppy[entry]:
                                        output.write('>' + entry[:-2] + '\n')
                                        output.write(str(record_dict[entry[:-2]].seq) + '\n')
