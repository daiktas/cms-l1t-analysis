from rootpy.ROOT import gPad
from cmsl1t.io import to_root, from_root
from copy import deepcopy
import logging
logger = logging.getLogger(__name__)


class BasePlotter(object):
    """
    A Base class to be used by the various plotters
    """

    def __init__(self, directory_name):
        self.directory_name = directory_name
        self._is_built = False

    def set_plot_output_cfg(self, outdir, fmt):
        """
        Store the output directory and image formats
        """
        self.output_dir = outdir
        self.output_format = fmt

    def get_output_image(self, name):
        return "{}/{}.{}".format(self.output_dir, name, self.output_format)

    def save_canvas(self, canvas=None, name=None):
        if not canvas:
            canvas = gPad.func()
        if not name:
            name = canvas.GetName()
        out_path = "{}/{}.{}".format(self.output_dir, name, self.output_format)
        canvas.SaveAs(out_path)

    def build(self, *vargs, **kwargs):
        """
        We're reading in tree data so build the histograms etc from scratch.
        Not used when we're restarting from previously generated histograms.
        Calls the create_histograms() method which the derived class should implement.
        """
        self.create_histograms(*vargs, **kwargs)
        self._is_built = True

    def create_histograms(self):
        """
        Create the histograms that this plotter needs.

        This is not done in an init function so that we can by-pass this in the
        case where we reload things from disk
        """
        raise NotImplementedError("create_histograms() needs to be implemented")

    def to_root(self, outfile):
        """
        Write this plotter to a root file.
        ROOT (and rootpy) take care of persisting the ROOT objects held by this class
        """
        to_root(self, outfile)
        return True

    def from_root(self, filename):
        """
        Reload histograms from existing files on disk.

        Might be called multiple times, in which case histograms should be
        merged together with existing ones
        """
        reloaded = from_root(filename)
        return self.merge_in(reloaded)

    def merge_in(self, reloaded):
        # Have already been initialised, so merge this in
        if self._is_built:
            if self._is_consistent(reloaded):
                return self._merge(reloaded)
            else:
                logger.error("Trying to reload an inconsistent histogram")
                return False

        # Have not yet been initialised so need copy over values
        final = reloaded.__dict__
        final.update(self.__dict__)
        self.__dict__ = final
        self._is_built = True
        return True

    def fill(self):
        """
        Has to be overloaded by users code.

        Fill the histograms with event data. The actual arguments passed to
        this function is up to the derived code.
        """
        raise NotImplementedError("fill() needs to be implemented")

    def draw(self):
        """
        Has to be overloaded by users code.

        Draw the histograms on a canvas.
        Can produce multiple output plots.
        Use self.save_canvas() to actually save the plot
        """
        raise NotImplementedError("draw() needs to be implemented")

    def get_stats(self, summary_bins=[], summary_label=''):
        """
        Has to be overloaded by users code.

        Return the statistics in each histogram on a canvas.
        Should return a Dataframe

        :param summary_bins bins to summarise the stats into
        :param summary_label label to assign to summary columns
        """
        raise NotImplementedError("get_stats() needs to be implemented")
