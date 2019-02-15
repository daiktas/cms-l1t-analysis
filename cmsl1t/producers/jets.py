from importlib import import_module
from ..jet import Jet, L1Jet, CaloJet, PFJet
from .base import BaseProducer

# TODO:
# 1. good reco jets (use jet filters)
# 2. matched L1 jets


def _load_filter_module(filterPath):
    print(filterPath)
    if filterPath is None:
        return None
    tokens = filterPath.split('.')
    module_path = '.'.join(tokens[:-1])
    function = tokens[-1]
    mod = import_module(module_path)
    return getattr(mod, function)


class Producer(BaseProducer):

    def __init__(self, inputs, outputs, **kwargs):
        self.__set_jetType(kwargs)
        super(Producer, self).__init__(inputs, outputs, **kwargs)

    def __set_jetType(self, params):
        if params and 'jetType' in params:
            self._jetType = params['jetType']
        else:
            self._jetType = 'Calo'
        self._init_jet_class()

        if params and 'filter' in params:
            self._jetFilter = _load_filter_module(params['filter'])
        else:
            self._jetFilter = None

    def _init_jet_class(self):
        if self._jetType == 'PF':
            self._expected_input_order = [
                'et', 'eta', 'phi', 'etCorr', 'cemef', 'chef',
                'cMult', 'mef', 'muMult', 'nemef', 'nhef', 'nMult',
            ]
            self._jetClass = PFJet
        elif self._jetType == 'Calo':
            self._expected_input_order = ['et', 'eta', 'phi', 'etCorr', ]
            self._jetClass = CaloJet
        elif self._jetType == 'Gen':
            self._expected_input_order = ['pt', 'eta', 'phi']
            self._jetClass = Jet
        elif self._jetType == 'L1':
            self._expected_input_order = ['et', 'eta', 'phi', 'bx']
            self._jetClass = L1Jet
        else:
            self._expected_input_order = []
            self._jetClass = None

    def produce(self, event):
        variables = [event[i] for i in self._inputs]
        jets = [self._jetClass(*args) for args in zip(*variables)]
        if 'L1' in self._jetType:
            jets = [jet for jet in jets if jet.bx == 0]
        if self._jetFilter:
            jets = self._jetFilter(jets)

        # sort by ET, largest first
        sorted_jets = sorted(
            jets,
            key=lambda jet: jet.etCorr if jet.etCorr else jet.et,
            reverse=True,
        )

        setattr(event, self._outputs[0], sorted_jets)
        return True
