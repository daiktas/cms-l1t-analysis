import logging
import six

import ROOT
from rootpy.tree import TreeChain

from cmsl1t.utils.root_glob import glob
from cmsl1t.utils.module import load_L1TNTupleLibrary

load_L1TNTupleLibrary()
logger = logging.getLogger(__name__)
sumTypes = ROOT.l1t.EtSum


def _get_input_files(paths):
    # TODO: move this replacement into cmsl1t.config
    input_files = []
    for p in paths:
        if '*' in p:
            input_files.extend(glob(p))
        else:
            input_files.append(p)
    return input_files


def _create_alias_map(ntuple_map):
    '''
        Translates the ntuple_map into a format that can more easily processed
        for creating event attributes, e.g.
        content:
            <tree path>:
                branches:
                  <branchname>.<leafname>:
                    aliases:
                      - event.myAlias
        to
        {'myAlias': (<tree path>, <branchname>.<leafname>)}

    '''
    # TODO: move this replacement into cmsl1t.config
    aliasMap = {}
    trees = ntuple_map['content']
    for treeName, content in trees.items():
        for branchName, branchContent in content['branches'].items():
            for alias in branchContent['aliases']:
                tidyAlias = alias.replace('event.', '')
                aliasMap[tidyAlias] = (treeName, branchName)
    return aliasMap


class EventReader(object):

    def __init__(self, input_files, ntuple_map, nevents=-1):
        '''
            Reads ntuple_info as defined by bin/create-map-file
        '''
        self._treeNames = ntuple_map['content'].keys()
        self._aliasMap = _create_alias_map(ntuple_map)
        self.input_files = _get_input_files(input_files)
        self.nevents = nevents
        self._trees = {}

        self._load_trees()

    def _load_trees(self):
        for treeName in self._treeNames:
            try:
                self._trees[treeName] = TreeChain(
                    treeName,
                    self.input_files,
                    cache=True,
                    events=self.nevents,
                )
            except RuntimeError:
                logger.warn(
                    "Cannot find tree: {0} in input file".format(treeName))
                continue

    def __contains__(self, name):
        return name in self._aliasMap.keys()

    def __iter__(self):
        # event loop
        for trees in six.moves.zip(*self._trees.itervalues()):
            yield Event(self._trees, self._aliasMap)


class Event(object):

    def __init__(self, trees, mapping):
        self._map = mapping
        self._trees = trees
        self._cache = {}

    def __getattr__(self, name):
        if name in object.__getattribute__(self, '_cache'):
            return object.__getattribute__(self, '_cache')[name]

        if name not in object.__getattribute__(self, '_map'):
            return object.__getattribute__(self, name)
        treeName, treeAttr = object.__getattribute__(self, '_map')[name]
        tree = object.__getattribute__(self, '_trees')[treeName]
        obj = tree
        for attr in treeAttr.split('.'):
            obj = getattr(obj, attr)
        object.__getattribute__(self, '_cache')[name] = obj
        return obj

    def __getitem__(self, name):
        return object.__getattribute__(self, '__getattr__')(name)
