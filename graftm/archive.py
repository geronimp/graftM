import tarfile
import logging
import itertools
import os
import tempdir
import shutil

from graftm_package import GraftMPackage

class ArchiveDefaultOptions:
    force=False

class Archive:
    def create(self, input_package_path, output_package_path, **kwargs):
        """Create an archived GraftM package

        Parameters
        ----------
        input_package_path: str
            path to gpkg to be archived
        output_pacakge_path: str
            output package path
        kwargs:
            force: bool
                overwrite an existing directory
        """
        force = kwargs.pop('force',ArchiveDefaultOptions.force)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Archiving GraftM package '%s' as '%s'" % (input_package_path, output_package_path))
        gpkg = GraftMPackage.acquire(input_package_path)
        if gpkg.version != 3:
            raise Exception("Archiving GraftM packages only works with format 3 packages")

        self._setup_output(output_package_path, force)

        logging.debug("Compressing contents for archive")
        with tarfile.open(output_package_path, 'w:gz') as tar:
            for path in itertools.chain([gpkg.contents_file_path(),
                                         gpkg.alignment_hmm_path(),
                                         gpkg.reference_package_path(),
                                         gpkg.unaligned_sequence_database_path()],
                                        [hmm for hmm in gpkg.search_hmm_paths() if hmm != gpkg.alignment_hmm_path()]):
                logging.debug("Compressing '%s'" % path)
                # Put the gpkg folder itself in the archive so as not to tar bomb.
                tar.add(path,
                        os.path.join(os.path.basename(os.path.abspath(gpkg._base_directory)),
                                     os.path.basename(path)))
            
        logging.info("Archive successfully created")

    def extract(self, archive_path, output_package_path, **kwargs):
        '''Extract an archived GraftM package.

        Parameters
        ----------
        archive_path: str
            path to archive
        output_package_path: str
            path to where to put the extracted file
        kwargs:
            force: bool
                overwrite an existing directory
        '''
        force = kwargs.pop('force', ArchiveDefaultOptions.force)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        logging.info("Un-archiving GraftM package '%s' from '%s'" % (output_package_path, archive_path))
        archive = os.path.abspath(archive_path)
        output = os.path.abspath(output_package_path)
        
        self._setup_output(output, force)
        
        with tempdir.in_tempdir():
            with tarfile.open(archive) as tar:
                tar.extractall()
                for tarinfo in tar:
                    folder = os.path.dirname(tarinfo.name)
                    break

            # recreate diamond file
            gpkg = GraftMPackage.acquire(folder)
            if gpkg.version != 3:
                raise Exception("Encountered an archived gpkg of unexpected version: %d" % gpkg.version)
            if gpkg.diamond_database_path():
                logging.debug("Creating diamond DB")
                gpkg.create_diamond_db()

            shutil.move(folder, output)

        logging.info("Archive successfully extracted")

    def _setup_output(self, path, force):
        '''Clear the way for an output to be placed at path'''
        # Allow for special case of output being a pipe
        if os.path.isdir(path) or os.path.isfile(path):
            if force:
                logging.warn("Deleting previous file/directory '%s'" % path)
                if os.path.isfile(path):
                    os.remove(path)
                else:
                    shutil.rmtree(path)
            else:
                raise Exception("Cowardly refusing to overwrite already existing path at %s" % path)
        

