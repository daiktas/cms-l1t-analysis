from __future__ import absolute_import
import six
import logging
import numpy as np

import cmsl1t.geometry as geo
from cmsl1t.utils.iterators import pairwise
from . import HistogramsByPileUpCollection

logger = logging.getLogger(__name__)


class ResolutionCollection(HistogramsByPileUpCollection):
    '''
        Specialisation of BaseHistCollection for making resolution plots.

        :Example:
            >>> res = ResolutionCollection()
            >>> res.add_variable('jetEt')
            >>> ...
            >>> res.set_pileup(pu)
            >>> res.set_region_by_eta(jet.eta())
            >>> res.fill('jetEt', jet.et())
            or for direct access
            >>> res[pu]['jetEt'][region].fill(jet.et())
    '''
    BINS = {
        'position': np.arange(-0.3, 0.3, 0.005),
        'energy': np.arange(-1, 1.5, 0.05),
        'default': np.arange(-1, 1, 0.05)
    }

    def __init__(self, pileupBins=[0, 13, 20, 999],
                 regions=geo.eta_regions):
        r'''

            :param dict regions: A dictionary of region_name: <function>.
                The function should take 1 value, \eta, and return a boolean.
        '''
        self._dimensions = 3
        HistogramsByPileUpCollection.__init__(
            self, pileupBins=pileupBins, dimensions=self._dimensions)
        self._pileUp = 0
        self._currentRegions = []
        self._pileUpBins = pileupBins
        self._regions = regions

    def set_region_by_eta(self, eta):
        self._currentRegions = [name for name in six.iterkeys(
            self._regions) if geo.is_in_region(name, eta, self._regions)]

    def fill(self, hist_name, x, w=1.0):
        h = self[self._pileUp][hist_name]
        if not h:
            logger.error('Histogram {0} does not exist'.format(hist_name))
            return
        if not self._currentRegions:
            logger.warn(
                'No valid current regions. Did you set_region_by_eta()?')
        for region in self._currentRegions:
            h[region].fill(x, w)

    def add_variable(self, variable, bins=[]):
        from rootpy.plotting import Hist
        if variable in self.keys():
            logger.warn('Variable {0} already exists!')
            return
        hist_names = []
        add_name = hist_names.append

        for puBinLower, puBinUpper in pairwise(self._pileUpBins):
            for region in six.iterkeys(self._regions):
                name = '{0}_{1}_pu{2}To{3}'.format(
                    variable, region, puBinLower, puBinUpper)
                if not self[puBinLower][variable][region]:
                    add_name(name)
                    self[puBinLower][variable][region] = Hist(bins, name=name)
        logger.debug('Created {0} histograms: {1}'.format(
            len(hist_names), ', '.join(hist_names)))
