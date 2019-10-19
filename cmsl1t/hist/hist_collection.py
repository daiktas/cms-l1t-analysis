import collections
import bisect
from exceptions import RuntimeError, KeyError, NotImplementedError
from copy import deepcopy
import logging
from cmsl1t.hist.factory import HistFactory
from cmsl1t.hist.binning import Base as BinningBase


__all__ = ["HistCollectionView", "HistogramCollection"]


logger = logging.getLogger(__name__)


class HistCollectionView(object):

    def __init__(self, bin_indices, hist_list):
        self.histograms = hist_list
        self.bin_indices = bin_indices

    def __getattr__(self, attr):
        return [getattr(hist, attr) for hist in self.histograms]

    def __method(self, method_name, *vargs, **kwargs):
        for hist in self.histograms:
            getattr(hist, method_name)(*vargs, **kwargs)

    def fill(self, *vargs, **kwargs):
        self.__method("fill", *vargs, **kwargs)

    def __iter__(self):
        for hist in self.histograms:
            yield hist

    def __len__(self):
        return len(self.histograms)

    def items(self):
        for bin_hist_pair in zip(self.bin_indices, self.histograms):
            yield bin_hist_pair


class HistogramCollection(object):
    '''
    The histogram collection needs a few things:
     - it needs to be able to essentially have binned maps of histograms
     - needs to know how to create new histograms
    '''

    def __init__(self, dimensions, histogram_factory, *vargs, **kwargs):
        '''
            Should dimensions include or exclude histogram names?
        '''
        if not isinstance(dimensions, list):
            dimensions = [dimensions]
        for dim in dimensions:
            if not isinstance(dim, BinningBase):
                raise RuntimeError("non-Dimension object given to histogram")
        self.__dimensions = dimensions
        self.shape = tuple([len(dim) for dim in dimensions])

        if isinstance(histogram_factory, str):
            histogram_factory = HistFactory(histogram_factory,
                                            *vargs,
                                            **kwargs)
        self.values = self._prepare_collection(dimensions, histogram_factory)

    def _prepare_collection(self, dimensions, histogram_factory,
                            bin_indices=[], depth=0):
        this_dim = deepcopy(dimensions[0])
        remaining_dims = dimensions[1:]
        for bin in this_dim.iter_all():
            if remaining_dims:
                value = self._prepare_collection(remaining_dims,
                                                 histogram_factory,
                                                 bin_indices + [bin],
                                                 depth + 1)
                this_dim.set_value(bin, value)
            else:
                indices = bin_indices + [bin]
                labels = {d.label: index for d, index in
                          zip(self.__dimensions, indices)}
                # TODO: Should fill proper bin labels here and pass through
                hist = histogram_factory(labels=labels)
                this_dim.set_value(bin, hist)
        return this_dim

    @classmethod
    def _flatten_bins(self, bins):
        flattened_bins = []
        for dimension in bins:
            if len(flattened_bins) == 0:
                for index in dimension:
                    flattened_bins.append([index])
            else:
                new_bins = []
                for previous in flattened_bins:
                    new_bins += [previous + [index] for index in dimension]
                flattened_bins = new_bins
        output_bin_list = []
        for bin in flattened_bins:
            output_bin_list.append(tuple(bin))
        return output_bin_list

    def _find_bins(self, keys):
        # In python 3.3, this becomes collections.abc.Sequence
        if not isinstance(keys, collections.Sequence):
            keys = [keys]

        n_keys = len(keys)

        # Check every dimension if it contains these values
        bins = []
        for key, dimension in zip(keys, self.__dimensions[:n_keys]):
            bins.append(dimension.find_all_bins(key))

        # Some dimensions might return multiple values, flatten returned arrays
        bins = self._flatten_bins(bins)

        return bins

    def get_bin_contents(self, bin_list):
        if isinstance(bin_list, (str, int)):
            bin_list = [bin_list]
        value = self.values
        for index in bin_list:
            value = value.get_bin_contents(index)
        return value

    def __getitem__(self, keys):
        '''
            Supposed to handle
                coll[x]
            and
                coll[x, y, z]
        '''
        bin_indices = self._find_bins(keys)
        objects = [self.get_bin_contents(bins) for bins in bin_indices]
        return HistCollectionView(bin_indices, objects)

    def shape(self):
        return self.shape

    def __len__(self):
        return len(self.__dimensions[0])

    def __iter__(self):
        # # In python >3.3 we should do
        # yield from self.__dimensions[0]
        for bin in self.__dimensions[0]:
            yield bin

    def iter_all(self):
        # # In python >3.3 we should do
        # yield from self.__dimensions[0]
        for bin in self.__dimensions[0].iter_all():
            yield bin

    def flat_items(self):
        all_bin_indices = [list(dim) for dim in self.__dimensions]
        flat_bins = self._flatten_bins(all_bin_indices)
        for bin in flat_bins:
            yield bin, self.get_bin_contents(bin)

    def flat_items_all(self):
        all_bin_indices = [list(dim.iter_all()) for dim in self.__dimensions]
        flat_bins = self._flatten_bins(all_bin_indices)
        for bin in flat_bins:
            yield bin, self.get_bin_contents(bin)

    def __iadd__(self, other):
        for bin, hist in self.flat_items_all():
            hist += other.get_bin_contents(bin)
        return self
