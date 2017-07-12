import os
import shutil
import logging
import inspect
from graftm.graftm_package import GraftMPackage

class InvalidFileExtensionError(Exception):
    pass

class UninstalledProgramError(Exception):
    pass

class HouseKeeping:
    ### Functions for setting up the graftM pipeline to run correctly, and
    ### general housekeeping things like adding and removing directories and
    ### files.

    HMMSEARCH_AND_DIAMOND_SEARCH_METHOD = 'hmmsearch+diamond'
    
    DIAMOND_SEARCH_METHOD   = 'diamond'
    HMMSEARCH_SEARCH_METHOD = 'hmmsearch'


    def file_basename(self, file):
        '''
        Strips the path and last extension from the file variable. If the
        extension is found to be valid, the basename will be returned. Otherwise
        an error will be raise and graftM will exit
        '''
        valid_extensions = set(".tree",
                               ".tre")
        split_file = os.path.basename(file).split('.')
        base, suffix = '.'.join(split_file[:-1]), split_file[-1]

        if suffix in valid_extensions:
            return base
        else:
            logging.error("Invalid file extension found on file: %s" % file)
            logging.error("For trees, please provide a file with one of the \
following extensions: %s" % ' '.join(valid_extensions.keys()))
            raise InvalidFileExtensionError


    def set_euk_hmm(self, args):
        'Set the hmm used by graftM to cross check for euks.'
        if hasattr(args, 'euk_hmm_file'):
            pass
        elif not hasattr(args, 'euk_hmm_file'):
            # set to path based on the location of bin/graftM, which has
            # a more stable relative path to the HMM when installed through
            # pip.
            setattr(args, 'euk_hmm_file', os.path.join(os.path.dirname(inspect.stack()[-1][1]),'..','share', '18S.hmm'))
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

            self._check_file_existence(args.forward)
            if args.reverse:
                self._check_file_existence(args.reverse)

            # Determine the File format based on the suffix

            sequence_file_list = []
            if args.reverse:
                if len(args.forward) != len(args.reverse):
                    logging.error('Confusing input. There appears to be different numbers of forward and reverse files specified')
                for i, forward_file in enumerate(args.forward):
                    sequence_file_list.append([forward_file, args.reverse[i]])
            else:
                sequence_file_list = [[f] for f in args.forward]

            return sequence_file_list

    def _check_file_existence(self, files):
        '''Iterate through files and exit(1) if any do not pass
        os.path.isfile'''
        for f in files:
            if not os.path.isfile(f):
                logging.error("The file '%s' does not appear to exist, stopping" % f)
                exit(1)



    def get_maximum_range(self, hmm):
        '''
        If no maximum range has been specified, and if using a hmm search, a
        maximum range can be determined by using the length of the HMM

        Parameters
        ----------
        hmm : str
            path to hmm profile

        Returns
        -------
        Length to search to when linking hits on a single contig
        '''
        length=int([x for x in open(hmm) if x.startswith("LENG")][0].split()[1])
        max_length=round(length*1.5, 0)

        return max_length

    def set_attributes(self, args):


        # Read graftM package and assign HMM and refpkg file
        if args.no_merge_reads:
            setattr(args, 'merge_reads', False)
        else:
            if args.reverse:
                setattr(args, 'merge_reads', True)
            else:
                setattr(args, 'merge_reads', False)
                
        if args.graftm_package:
            if not os.path.isdir(args.graftm_package):
                raise Exception("%s does not exist. Are you sure you provided the correct path?" % args.graftm_package)
            else:
                gpkg = GraftMPackage.acquire(args.graftm_package)
                if hasattr(args, 'search_hmm_files'): # If a hmm is specified, overwrite the one graftM package
                    setattr(args, 'aln_hmm_file', gpkg.alignment_hmm_path())
                    setattr(args, 'reference_package', gpkg.reference_package_path())
                else:
                    setattr(args, 'search_hmm_files', [])
                    for hmm in gpkg.search_hmm_paths():
                        args.search_hmm_files.append(hmm)
                    setattr(args, 'aln_hmm_file', gpkg.alignment_hmm_path())
                    setattr(args, 'reference_package', gpkg.reference_package_path())

        elif hasattr(args, 'search_diamond_files'):
            if args.search_method == self.DIAMOND_SEARCH_METHOD:
                if hasattr(args, 'aln_hmm_file'):
                    pass
                else:
                    raise Exception("aln_hmm_file not specified")
            else:
                raise Exception("Specified DIAMOND databases when not using the diamond search pipeline. Using: %s" % (args.search_method))

        elif hasattr(args, 'search_hmm_files'):
            if args.search_method == self.HMMSEARCH_SEARCH_METHOD:
                if not hasattr(args, 'aln_hmm_file'):
                    if len(args.search_hmm_files) == 1:
                        if not args.search_only:
                            setattr(args, 'aln_hmm_file', args.search_hmm_files[0])
                    else:
                        raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")

            else:
                raise Exception("Specified HMM search_hmm_files when not using the hmmsearch pipeline. Using: %s" % (args.search_method))

        elif hasattr(args, 'search_hmm_list_file'):
            if args.search_method == self.HMMSEARCH_SEARCH_METHOD:
                setattr(args, 'search_hmm_files', [x.rstrip() for x in open(args.search_hmm_list_file).readlines()])
                if not hasattr(args, 'aln_hmm_file'):
                    if not args.search_only:
                        raise Exception("Multiple search HMMs specified, but aln_hmm_file not specified")
            else:
                raise Exception("Specified HMM search_hmm_files when not using the hmmsearch pipeline. Using: %s" % (args.search_method))

        else:
            if args.search_only:
                if args.search_diamond_file:
                    args.search_method = self.DIAMOND_SEARCH_METHOD
                    args.search_hmm_files = None
            else:
                raise Exception('No gpkg, HMM, or DIAMOND database was specified, so there is no reference database to search with.')

