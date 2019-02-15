from __future__ import print_function
import pandas as pd

from cmsl1t.plotting.base import BasePlotter
from cmsl1t.hist.hist_collection import HistogramCollection
import cmsl1t.hist.binning as bn
from cmsl1t.utils.draw import draw, label_canvas
from cmsl1t.recalc.resolution import get_resolution_function

from rootpy.context import preserve_current_style
from rootpy.plotting import Legend


class RateVsPileupPlot(BasePlotter):

    def __init__(self, online_name):
        name = ["rate_vs_pileup", online_name]
        super(RateVsPileupPlot, self).__init__("__".join(name))
        self.online_name = online_name

    def create_histograms(self,
                          online_title,
                          thresholds, n_bins, low, high, legend_title=""):
        """ This is not in an init function so that we can by-pass this in the
        case where we reload things from disk """
        self.online_title = online_title
        self.thresholds = thresholds
        self.thresholds = bn.GreaterThan(thresholds, "threshold", True)
        self.legend_title = legend_title
        name = ["rate_vs_pileup", self.online_name]
        name += ["thresh_{threshold}"]
        name = "__".join(name)
        title = " ".join([self.online_name, " rate vs pileup",
                          "passing threshold: {threshold}"])
        self.plots = HistogramCollection([self.thresholds],
                                         "Hist1D", n_bins, low, high,
                                         name=name, title=title)

        self.filename_format = name

    def fill(self, pileup, online):
        self.plots[online].fill(pileup)

    def draw(self, with_fits=False):

        for (threshold, ), hist in self.plots.flat_items_all():
            hists = []
            labels = []
            fits = []
            thresholds = []
            if not isinstance(threshold, int):
                continue
            label_template = '{online_title} > {threshold} GeV'
            label = label_template.format(
                online_title=self.online_title,
                threshold=self.thresholds.bins[threshold],
            )
            hist.Divide(self.plots.get_bin_contents([bn.Base.everything]))
            hist.drawstyle = "EP"
            hist.SetMarkerSize(0.5)
            hist.SetMarkerColor(2)
            # if with_fits:
            #    fit = self.fits.get_bin_contents([threshold])
            #    fits.append(fit)
            hists.append(hist)
            labels.append(label)
            thresholds.append(threshold)

            self.__make_overlay(hists, fits, labels, thresholds)

    def overlay_with_emu(self, emu_plotter, with_fits=False):

        for (threshold, ), hist in self.plots.flat_items_all():
            hists = []
            labels = []
            fits = []
            thresholds = []
            if not isinstance(threshold, int):
                continue
            label_template = '{online_title} > {threshold} GeV'
            label = label_template.format(
                online_title=self.online_title,
                threshold=self.thresholds.bins[threshold],
            )
            hist.Divide(self.plots.get_bin_contents([bn.Base.everything]))
            hist.drawstyle = "EP"
            hist.SetMarkerSize(0.5)
            hist.SetMarkerColor(1)
            # if with_fits:
            #    fit = self.fits.get_bin_contents([threshold])
            #    fits.append(fit)
            hists.append(hist)
            labels.append(label)
            thresholds.append(threshold)
            for (emu_threshold, ), emu_hist in emu_plotter.plots.flat_items_all():
                if emu_threshold == threshold:
                    label_template = '{online_title} > {threshold} GeV'
                    emu_label = label_template.format(
                        online_title=emu_plotter.online_title,
                        threshold=emu_plotter.thresholds.bins[threshold],
                    )
                    emu_hist.Divide(emu_plotter.plots.get_bin_contents([bn.Base.everything]))
                    emu_hist.drawstyle = "EP"
                    emu_hist.SetMarkerSize(0.5)
                    emu_hist.SetMarkerColor(2)
                    # if with_fits:
                    #    emu_fit = self.fits.get_bin_contents([threshold])
                    #    fits.append(emu_fit)
                    hists.append(emu_hist)
                    labels.append(emu_label)
                    thresholds.append(emu_threshold)

            self.__make_overlay(hists, fits, labels, thresholds)

    def __make_overlay(self, hists, fits, labels, thresholds, suffix=""):
        with preserve_current_style():
            # Draw each rate vs pileup (with fit)
            xtitle = "# reco vertices"
            ytitle = "a.u."
            canvas = draw(hists, draw_args={"xtitle": xtitle, "ytitle": ytitle})
            if fits:
                for fit, hist in zip(fits, hists):
                    fit["asymmetric"].linecolor = hist.GetLineColor()
                    fit["asymmetric"].Draw("same")

            # Add labels
            label_canvas()

            # Add a legend
            legend = Legend(
                len(hists),
                header=self.legend_title,
                topmargin=0.35,
                rightmargin=0.7,
                leftmargin=0.3,
                textsize=0.025,
                entryheight=0.028,
            )
            for hist, label in zip(hists, labels):
                legend.AddEntry(hist, label)
                legend.SetBorderSize(0)
                legend.Draw()

            # Save canvas to file
            name = self.filename_format.format(threshold=thresholds[0])
            self.save_canvas(canvas, name + suffix)

    def _is_consistent(self, new):
        """
        Check the two plotters are the consistent, so same binning and same axis names
        """
        return all([self.thresholds.bins == new.thresholds.bins,
                    self.online_name == new.online_name,
                    ])

    def _merge(self, other):
        """
        Merge another plotter into this one
        """
        self.plots += other.plots
        return self.plots

    def get_stats(self, summary_bins=[], summary_label=''):
        summary_columns = list(self._summary_columns(summary_bins, summary_label))
        stats = list(self._collect_stats(summary_bins, summary_label))
        df = pd.DataFrame(stats)
        return df[['identifier', 'total', 'overflow'] + summary_columns]

    def _summary_columns(self, summary_bins, summary_label):
        for lower, upper in zip(summary_bins[:-1], summary_bins[1:]):
            yield '{} {}-{}'.format(summary_label, lower, upper)

    def _collect_stats(self, summary_bins, summary_label):
        normalisation = self.plots.get_bin_contents([bn.Base.everything]).integral(overflow=True)
        for (threshold, ), hist in self.plots.flat_items():
            human_readable_threshold = '{0} > {1} GeV'.format(self.online_title, self.thresholds.bins[threshold])
            rhist = hist.rebinned(summary_bins)
            stats = {}
            summary_columns = self._summary_columns(summary_bins, summary_label)
            for summary_column, y in zip(summary_columns, rhist.y()):
                stats[summary_column] = y * normalisation
            total = sum(stats.values())
            overflow = rhist.integral(overflow=True) * normalisation - total
            header = dict(identifier=human_readable_threshold, total=total, overflow=overflow)
            header.update(stats)
            yield header
