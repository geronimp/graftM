import random
import os
import shutil
import subprocess

from graftm.messenger import Messenger
from test.test_MimeWriter import OUTPUT

# Constants - don't change them evar.
FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class HouseKeeping:

    def __init__(self): pass

    def make_filter_directories(self, filter_list, input_list, force):

        output_hash = {}

        for sequence_file in input_list:

            base_sequence = os.path.basename(sequence_file).split('.')[0]

            if force:
                shutil.rmtree(base_sequence, ignore_errors=True)

            try:
                os.mkdir(base_sequence)
            except:
                Messenger().header('Directory %s already exists. Exiting to prevent over-writing\n' % base_sequence)
                exit(1)

            output_hash[sequence_file] = {'dir': base_sequence}

            for f in filter_list:
                base_filter = os.path.basename(f.split('.')[0])
                filter_output = os.path.join(base_sequence, base_filter)
                os.mkdir(filter_output)
                output_hash[sequence_file][base_filter+ '.hmm'] = filter_output

        return output_hash


    def read_contents(self, gpkg_path):

        c = os.path.join(gpkg_path, "contents.txt")

        contents_hash = {}

        f =  [x.rstrip() for x in open(c, 'r').readlines()]

        for entry in f:
            s = entry.split('=')
            contents_hash[s[0]] = s[1].split(',')

        contents_hash['HMM'] = [os.path.join(gpkg_path, x) for x in contents_hash['HMM']]
        return contents_hash


    # Given a Return the guessed file format, or raise an Exception if
    def guess_sequence_input_file_format(self, sequence_file_path):
        if sequence_file_path.endswith(('.fa', '.faa', '.fna')):  # Check the file type
            return FORMAT_FASTA
        elif sequence_file_path.endswith(('.fq.gz', '.fastq.gz')):
            return FORMAT_FASTQ_GZ
        else:
            raise Exception("Unable to guess file format of sequence file: %s" % sequence_file_path)

    def guess_sequence_type(self, input_file, file_format):

        aas = set(['P','V','L','I','M','F','Y','W','H','K','R','Q','N','E','D','S'])
        nas = set(['A', 'T', 'G', 'C', 'N', 'U'])

        if file_format == 'FORMAT_FASTQ_GZ':
            filename = "/tmp/graftM_sample_"+str(random.randint(1, 100))+".fa"
            #unzip, and head into tmp_file, read tmp file, and return sequence type
            cmd = "head -n2 <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s 2> /dev/null) 2> /dev/null ) > %s " % (input_file, filename)
            subprocess.check_call(["/bin/bash", "-c", cmd])
            input_file = filename

        with open(input_file) as in_file:
            head = [next(in_file).rstrip() for x in xrange(2)]
            for nucl in set(head[1]):
                if nucl not in nas and nucl in aas:
                    return 'protein'
                elif nucl not in nas and nucl not in aas:
                    Messenger().error_message('Encountered unexpected character when attempting to guess sequence type: %s' % (nucl))
                    exit(1)
                else:
                    continue

            return 'nucleotide'

    def set_euk_hmm(self, args):
        if hasattr(args, 'euk_hmm_file'):
            pass
        elif not hasattr(args, 'euk_hmm_file'):
            setattr(args, 'euk_hmm_file', os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'constants', '18S.hmm'))
        else:
            raise Exception('Programming Error: setting the euk HMM')    
        return
    
    def setpipe(self, hmm):
        t = [x.rstrip() for x in open(hmm, 'r').readlines() if x.startswith('ALPH')][0].split()[1]

        if t == 'DNA':
            return 'D'
        elif t == 'amino':
            return 'P'
        else:
            Messenger().error_message('Unrecognised HMM library. Not sure which pipeline to use.')
            exit(1)

    def delete(self, delete_list):
        for item in delete_list:
            try:
                os.remove(item)

            except:
                pass

    def add_cmd(self, cmd_log, command):

        if not os.path.isfile(cmd_log):
            open(cmd_log, 'a').close()

        with open(cmd_log, 'a') as log:
            log.write(command + '\n')


    def make_working_directory(self, directory_path, force):

        if force:
            shutil.rmtree(directory_path, ignore_errors=True)
            os.mkdir(directory_path)

        else:
            try:
                os.mkdir(directory_path)
            except:
                Messenger().header('Directory %s already exists. Exiting to prevent over-writing\n' % directory_path)
                exit(1)

    def parameter_checks(self, args):
        ## Check that the necessary files are in place
        if args.subparser_name == 'graft':
            # Check that the placement cutoff is between 0.5 and 1
            if float(args.placements_cutoff) < float(0.5) or float(args.placements_cutoff) > float(1.0):
                Messenger().message('Please specify a confidence level (-d) between 0.5 and 1.0! Found: %s' % args.placements_cutoff)
                exit(1)

            # Set string for hmmsearch evalue
            args.eval = '-E %s' % args.eval

            # Determine the File format based on the suffix
            input_file_format = self.guess_sequence_input_file_format(args.forward)
            sequence_file_list = []
            if hasattr(args, 'reverse'):
                fors = args.forward.split(',')
                revs = args.reverse.split(',')

                for i in range(0, len(fors)):
                    sequence_file_list.append([fors[i], revs[i]])

            elif not hasattr(args, 'reverse'):
                for f in args.forward.split(','):
                    f = [f]
                    sequence_file_list.append(f)
            else:
                Messenger().error_message('Confusing input. Did you specify the same amount of reverse and forward read files?')
            return sequence_file_list, input_file_format
        else:
            return
        
    def set_attributes(self, args):
        # Read graftM package and assign HMM and refpkg file
        if hasattr(args, 'graftm_package'):
            if hasattr(args, 'hmm_file'): # If a hmm is specified, overwrite the one graftM package
                setattr(args, 'reference_package', os.path.join(args.graftm_package, [x for x in os.listdir(args.graftm_package) if x.endswith('.refpkg')][0]))
            elif not hasattr(args, 'hmm_file'): 
                setattr(args, 'hmm_file', os.path.join(args.graftm_package, [x for x in os.listdir(args.graftm_package) if x.endswith('.hmm')][0]))
                setattr(args, 'reference_package', os.path.join(args.graftm_package, [x for x in os.listdir(args.graftm_package) if x.endswith('.refpkg')][0]))
        elif not hasattr(args, 'graftm_package'): # Or if no graftM package and hmm is specified
            if hasattr(args, 'hmm_file') and args.search_only or args.search_and_align_only: # and the placement step is skipped
                pass # That's okay
            elif hasattr(args, 'hmm_file') and not args.search_only or not args.search_and_align_only: # But if a hmm is specified
                raise Exception('No refpkg specified.')
        else:
            raise Exception("Programming Error: Assigning hmm and refpkg attributes to args")       
        return

