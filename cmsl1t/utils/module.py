from importlib import import_module
import os
import logging
import sys
from .. import PROJECT_ROOT

logger = logging.getLogger(__name__)


def exists(module_name):
    try:
        import_module(module_name)
    except ImportError as e1:
        logger.debug(
            'Error importing a module "{}": {}'.format(module_name, e1))
        if '.' not in module_name:
            return False
        # check if it is a member of a module instead
        tokens = module_name.split('.')
        module_name = '.'.join(tokens[:-1])
        try:
            m = import_module(module_name)
        except ImportError as e2:
            logger.debug(
                'Error importing a module "{}": {}'.format(module_name, e2))
            return False
        return hasattr(m, tokens[-1])
    else:
        return True


def load_L1TNTupleLibrary(lib_name='L1TAnalysisDataformats.so'):
    import ROOT
    external_includes = os.path.join(PROJECT_ROOT, 'external')
    if external_includes not in ROOT.gSystem.GetIncludePath():
        ROOT.gSystem.AddIncludePath('-I"{}"'.format(external_includes))
    if lib_name not in ROOT.gSystem.GetLibraries():
        path_to_lib = os.path.join(PROJECT_ROOT, 'build', lib_name)
        ROOT.gSystem.Load(path_to_lib)
        if lib_name not in ROOT.gSystem.GetLibraries():
            logger.error('Could not load ROOT library {0}'.format(lib_name))
            sys.exit(-1)
