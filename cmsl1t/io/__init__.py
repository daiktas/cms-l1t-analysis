from rootpy.io.pickler import dump, load
from rootpy.plotting.hist import Hist, _HistBase
from cmsl1t.hist.hist_collection import HistogramCollection
# no pickles without dill
import dill


def to_root(obj, output_file):
    '''
        Saves the obj into a ROOT file
    '''
    if isinstance(output_file, str) and not output_file.endswith('.root'):
        output_file += '.root'
    dump(obj, output_file)


def from_root(input_file):
    '''
        Loads the obj from a ROOT file
    '''
    reloaded = load(input_file, use_proxy=False)

    # Histogram objects are "owned" by the current TDirectory.  Need to unhook
    # them, so that they're not free-ed when we close that the directory
    for obj in reloaded.__dict__.values():
        if isinstance(obj, (Hist, _HistBase)):
            obj.SetDirectory(None)
        elif isinstance(obj, HistogramCollection):
            for bins, hist in obj.flat_items_all():
                if isinstance(hist, (Hist, _HistBase)):
                    hist.SetDirectory(None)

    return reloaded
