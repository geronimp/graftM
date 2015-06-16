import random
import os
import shutil
import subprocess
import json
import tempfile
import logging
from timeit import itertools

# Constants - don't change them evar.
FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class UninstalledProgramError(Exception):
    pass

class HouseKeeping:
    ### Functions for setting up the graftM pipeline to run correctly, and 
    ### general housekeeping things like adding and removing directories and 
    ### files.
    
    def __init__(self): pass

    def contents(self, path):
        contents = os.path.join(path, 'CONTENTS.json')
        if os.path.isfile(contents):
            return json.load(open(contents, 'r'))
        else:
            return None
        
    def guess_sequence_input_file_format(self, sequence_file_path):
        '''Given a sequence file, guess the format and return. Raise an 
        exception if it cannot be guessed'''
        if sequence_file_path.endswith(('.fa', '.faa', '.fna', '.fasta')):  # Check the file type
            return FORMAT_FASTA
        elif sequence_file_path.endswith(('.fq.gz', '.fastq.gz')):
            return FORMAT_FASTQ_GZ
        else:
            raise Exception("Unable to guess file format of sequence file: %s" % sequence_file_path)
        
    def guess_sequence_type(self, input_file, file_format):
        '''Guess the type of input sequence provided to graftM (i.e. nucleotide
        or amino acid) and return'''
        
        # Define expected residues for each sequence type
        aas = set(['P','V','L','I','M','F','Y','W','H','K','R','Q','N','E','D','S'])
        nas = set(['A', 'T', 'G', 'C', 'N', 'U'])
        
        # If its Gzipped and fastq make a small sample of the sequence to be 
        # read
        if file_format == FORMAT_FASTQ_GZ:
            return 'nucleotide'
        else:
            with open(input_file) as in_file:
                head = [next(in_file).rstrip() for x in xrange(2)]
                for nucl in set(head[1]):
                    if nucl not in nas and nucl in aas:
                        return 'protein'
                    elif nucl not in nas and nucl not in aas:
                        raise Exception(logging.error('Encountered unexpected character when attempting to guess sequence type: %s' % (nucl)))
                    else:
                        continue
                return 'nucleotide'

    def set_euk_hmm(self, args):
        'Set the hmm used by graftM to cross check for euks.'
        if hasattr(args, 'euk_hmm_file'):
            pass
        elif not hasattr(args, 'euk_hmm_file'):
            setattr(args, 'euk_hmm_file', os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'constants', '18S.hmm'))
        else:
            raise Exception('Programming Error: setting the euk HMM')    

    
    def setpipe(self, hmm):
        ## Read in the hmm file and return the evalue trusted cutoff and the
        ## type of HMM
        hmm_tc = None
        hmm_read = [x.rstrip() for x in open(hmm, 'r').readlines() if x.startswith('ALPH') or x.startswith('TC')]
        if len(hmm_read) > 2:
            raise Exception('Possibly Programming Error: Misread HMM library.')
        for item in hmm_read:
            if item.startswith('ALPH'):
                if item.split()[1] == 'DNA' or item.split()[1] == 'RNA':
                    hmm_type = 'D'
                elif item.split()[1] == 'amino':
                    hmm_type = 'P'
            elif item.startswith('TC'):
                hmm_tc = item.split()[1]
            else:
                raise Exception('Possibly Programming Error: Misread HMM library.')
        return hmm_type, hmm_tc
        
    def delete(self, delete_list):
        for item in delete_list:
            try:
                os.remove(item)
            except:
                pass

    def make_working_directory(self, directory_path, force):
        if force:
            shutil.rmtree(directory_path, ignore_errors=True)
            os.mkdir(directory_path)
        else:
            try:
                os.mkdir(directory_path)
            except:
                logging.error('Directory %s already exists. Exiting to prevent over-writing' % directory_path)
                raise Exception('Directory %s already exists. Exiting to prevent over-writing'% directory_path)
    
    def parameter_checks(self, args):
        ## Check that the necessary files are in place
        if args.subparser_name == 'graft':
            # Check that the placement cutoff is between 0.5 and 1
            if float(args.placements_cutoff) < float(0.5) or float(args.placements_cutoff) > float(1.0):
                logging.info('Please specify a confidence level (-d) between 0.5 and 1.0! Found: %s' % args.placements_cutoff)
                exit(1)

            # Set string for hmmsearch evalue
            args.eval = '-E %s' % args.eval
            
            self._check_file_existence(args.forward)
            if hasattr(args, 'reverse'):
                self._check_file_existence(args.reverse)
                    
            # Determine the File format based on the suffix
            input_file_format = self.guess_sequence_input_file_format(args.forward[0])
            logging.debug("Detected file format %s" % input_file_format)
            sequence_file_list = []
            if hasattr(args, 'reverse'):
                if len(args.forward) != len(args.reverse):
                    logging.error('Confusing input. There appears to be different numbers of forward and reverse files specified')
                for i, forward_file in enumerate(args.forward):
                    sequence_file_list.append([forward_file, args.reverse[i]])
            else:
                sequence_file_list = [[f] for f in args.forward]
                
            return sequence_file_list, input_file_format
        else:
            return
        
    def _check_file_existence(self, files):
        '''Iterate through files and exit(1) if any do not pass
        os.path.isfile'''
        for f in files:
            if not os.path.isfile(f):
                logging.error("The file '%s' does not appear to exist, stopping" % f)
                exit(1)
            
    
    def checkCreatePrerequisites(self):
        uninstalled_programs = []
        prerequisites = {'taxit': 'https://github.com/fhcrc/taxtastic',
                         'FastTree': 'http://www.microbesonline.org/fasttree/',
                         'seqmagick': 'https://github.com/fhcrc/seqmagick',
                         'hmmalign': 'http://hmmer.janelia.org/'}
        for program in prerequisites.keys():
            if self.which(program):
                pass
            else:
                uninstalled_programs.append(program)
        if uninstalled_programs:
            msg = "The following programs must be installed to run GraftM create\n"
            logging.info(msg)
            for program in uninstalled_programs:
                l = '\t%s\t%s' % (program, prerequisites[program])
                print l
                msg += l+"\n"
            exit(0)
            
    def set_attributes(self, args):
        # Check the presence of all prerequisite programs needed for GraftM
        uninstalled_programs = []
        prerequisites = {'orfm': 'https://github.com/wwood/OrfM',
                        'nhmmer': 'http://hmmer.janelia.org/',
                        'hmmsearch': 'http://hmmer.janelia.org/',
                        'fxtract': 'https://github.com/ctSkennerton/fxtract',
                        'pplacer': 'http://matsen.fhcrc.org/pplacer/',
                        'seqmagick': 'https://github.com/fhcrc/seqmagick',
                        'ktImportText': 'http://sourceforge.net/p/krona/home/krona/'}
        for program in prerequisites.keys():
            if self.which(program):
                pass
            else:
                uninstalled_programs.append(program)
        if uninstalled_programs:
            msg = "The following programs must be installed to run GraftM\n"
            logging.info(msg)
            for program in uninstalled_programs:
                l = '\t%s\t%s' % (program, prerequisites[program])
                print l
                msg += l+"\n"
            raise UninstalledProgramError(msg)
        else:
            pass
        
        # Read graftM package and assign HMM and refpkg file
        
        if hasattr(args, 'graftm_package'):
            if not os.path.isdir(args.graftm_package):
                raise Exception("%s does not exist. Are you sure you provided the correct path?" % args.graftm_package)
            if self.contents(args.graftm_package) is None:
                raise Exception('Misformatted GraftM package')
            elif self.contents(args.graftm_package) is not None:
                c = self.contents(args.graftm_package)
                if hasattr(args, 'search_hmm_files'): # If a hmm is specified, overwrite the one graftM package
                    setattr(args, 'aln_hmm_file', os.path.join(args.graftm_package, c['aln_hmm']))
                    setattr(args, 'reference_package', os.path.join(args.graftm_package, c['rfpkg']))
                else:
                    setattr(args, 'search_hmm_files', [])
                    for hmm in c['search_hmm']:
                        args.search_hmm_files.append(os.path.join(args.graftm_package, hmm))
                    setattr(args, 'aln_hmm_file', os.path.join(args.graftm_package, c['aln_hmm']))
                    setattr(args, 'reference_package', os.path.join(args.graftm_package, c['rfpkg']))
        elif hasattr(args, 'search_hmm_files'):
            if not hasattr(args, 'aln_hmm_file'):
                if len(args.search_hmm_files) == 1:
                    setattr(args, 'aln_hmm_file', args.search_hmm_files[0])
                else:
                    raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")
        elif hasattr(args, 'search_hmm_list_file'):
            setattr(args, 'search_hmm_files', [x.rstrip() for x in open(args.search_hmm_list_file).readlines()])

            if not hasattr(args, 'aln_hmm_file'):
                raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")
        else:
            raise Exception('No refpkg or HMM specified: Do not know what to search with.')
        return
    
    def which(self, program):
        '''Credits to BamM and http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python'''
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
        return None

