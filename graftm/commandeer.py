import subprocess
import logging
import shlex


class CommandeerCalledProcessError(subprocess.CalledProcessError):
    def __init__(self, command, returncode, stderr, stdout):
        self.command = command
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = stdout
        
    def __str__(self):
        return "Command %s returned non-zero exit status %i.\n"\
            "STDERR was: %sSTDOUT was: %s" % (self.command,
                                                self.returncode,
                                                self.stderr,
                                                self.stdout)


class Commandeer:
    @staticmethod
    def run(command):
        '''
        Run a subprocess.check_output() with the given command with shell=True,
        returning the stdout. If the command fails, raise a 
        SubprocessFailedException that includes the $stderr as part of the
        error message
        
        Parameters
        ----------
        command: str
            command to run
        '''
        logging.debug("Running commandeer cmd: %s" % command)

        process = subprocess.Popen(["bash","-c", command],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
            
        if process.returncode != 0:
            raise CommandeerCalledProcessError(command,
                                               process.returncode,
                                               stderr,
                                               stdout)
        return stdout
