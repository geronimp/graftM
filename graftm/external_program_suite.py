import extern
import logging

class ExternalProgramSuite:

    _FASTTREE="FastTreeMP"

    package_urls = {'taxit': 'https://github.com/fhcrc/taxtastic',
             'FastTree': 'http://www.microbesonline.org/fasttree/',
             'orfm': 'https://github.com/wwood/OrfM',
            'fxtract': 'https://github.com/ctSkennerton/fxtract',
            'pplacer': 'http://matsen.fhcrc.org/pplacer/',
            'krona': 'http://sourceforge.net/p/krona/home/krona/',
            'mafft': 'http://mafft.cbrc.jp/alignment/software/',
            'diamond': 'https://github.com/bbuchfink/diamond',
            'hmmer': 'http://hmmer.janelia.org/'}

    programs_to_packages = {'hmmalign': 'hmmer',
                            'hmmsearch': 'hmmer',
                            'nhmmer': 'hmmer',
                            'ktImportText': 'krona',
                            'FastTreeMP': 'FastTree'}

    programs_to_possibilities = {
                                 _FASTTREE: ["FastTreeMP",
                                              "fasttree",
                                              "fasttreeMP",
                                              "FastTree"]
                                 }

    def __init__(self, program_list):
        '''Given a list of executable names, check that they are available
        on the PATH, raising an exception otherwise

        Parameters
        ----------
        program_list: iterable of str
            list of program names to check
        '''

        uninstalled_programs = []
        for program in program_list:
            if program in ExternalProgramSuite.programs_to_possibilities:
                check = [p for p in ExternalProgramSuite\
                                        .programs_to_possibilities[program]
                         if extern.which(p)]
                if len(check)>1:
                    logging.warning("Program found with multiple commands: %s. \
Arbitrarily selecting %s" % (' '.join(program_list), check[0]))

                for key, item in ExternalProgramSuite\
                                        .programs_to_possibilities\
                                        .iteritems():
                    if program in item:
                        program = key

                if len(check) > 0:
                    if key == self._FASTTREE:
                        self.fasttree = check[0]
                else:
                    uninstalled_programs.append(program)
            else:
                if not extern.which(program):
                    uninstalled_programs.append(program)

        if any(uninstalled_programs):
            msg = "The following programs appear to be missing, and need to be updated or installed before GraftM can continue:"
            logging.error(msg)
            for program in sorted(uninstalled_programs):
                try:
                    package_name = ExternalProgramSuite.programs_to_packages[program]
                    print_name = '/'.join((program, package_name))
                except KeyError:
                    package_name = program
                    print_name = program

                logging.error('\t%25s\t%s' % (print_name, ExternalProgramSuite.package_urls[package_name]))
            exit(1)
