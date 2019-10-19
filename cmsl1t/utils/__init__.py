import os
import logging
from rootpy import ROOT, ROOTError

logger = logging.getLogger(__name__)


def load_ROOT_library(library, lib_path='build'):
    PROJECT_ROOT = os.environ.get('PROJECT_ROOT', os.getcwd())
    if library in ROOT.gSystem.GetLibraries():
        # already loaded
        return
    # try adding from lib_path
    path_to_lib = os.path.join(PROJECT_ROOT, lib_path, library)
    if os.path.exists(path_to_lib):
        try:
            ROOT.gSystem.Load(path_to_lib)
        except ROOTError as e:
            msg = 'Could not load {0}: {1}'.format(path_to_lib, e)
            logger.error(msg)
    else:
        logger.error('Library {0} does not exist'.format(library))
