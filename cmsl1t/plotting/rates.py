from __future__ import division
import numpy as np
import pandas as pd
from cmsl1t.plotting.base import BasePlotter
from cmsl1t.hist.hist_collection import HistogramCollection
from cmsl1t.hist.factory import HistFactory
import cmsl1t.hist.binning as bn
from cmsl1t.utils.draw import draw, label_canvas
from cmsl1t.utils.hist import cumulative_hist, normalise_to_collision_rate
from cmsl1t.utils.hist import normalise_to_unit_area


from rootpy.context import preserve_current_style
from rootpy.plotting import Legend


class RatesPlot(BasePlotter):

    def __init__(self, online_name):
        name = ["rates", online_name]
        super(RatesPlot, self).__init__("__".join(name))
        self.online_name = online_name

    def create_histograms(self,
                          online_title,
                          pileup_bins, n_bins, low, high, legend_title=""):
        """ This is not in an init function so that we can by-pass this in the
        case where we reload things from disk """
        self.online_title = online_title
        self.pileup_bins = bn.Sorted(pileup_bins, "pileup",
                                     use_everything_bin=True)
        self.legend_title = legend_title

        name = ["rate_vs_threshold", self.online_name, "pu_{pileup}"]
        name = "__".join(name)
        title = " ".join([self.online_name, "vs.", "in PU bin: {pileup}"])
        title = ";".join([title, self.online_title])
        self.plots = HistogramCollection([self.pileup_bins],
                                         "Hist1D", n_bins, low, high,
                                         name=name, title=title)
        self.filename_format = name

    def fill(self, pileup, online):
        self.plots[pileup].fill(online)

    def draw(self, with_fits=False):
        hists = []
        labels = []
        fits = []
        for (pile_up, ), hist in self.plots.flat_items_all():
            h = cumulative_hist(hist)
            h = normalise_to_collision_rate(h)
            if pile_up == bn.Base.everything:
                h.linestyle = "dashed"
                label = "Everything"
            elif isinstance(pile_up, int):
                h.drawstyle = "EP"
                label = "~ {:.0f}".format(
                    self.pileup_bins.get_bin_center(pile_up))
            else:
                continue
            h.SetMarkerSize(0.5)
            hists.append(h)
            labels.append(label)
            # if with_fits:
            #     fits.append(self.fits.get_bin_contents([pile_up]))
        self.__make_overlay(hists, fits, labels, "Number of events")

        normed_hists = list(normalise_to_unit_area(hists))
        self.__make_overlay(normed_hists, fits, labels,
                            "Fraction of events", "__shapes")

    def overlay_with_emu(self, emu_plotter, with_fits=False):

        hists = []
        labels = []
        fits = []

        hist = self.plots.get_bin_contents([bn.Base.everything])
        hist = cumulative_hist(hist)

        hist.drawstyle = "EP"
        hist.SetMarkerSize(0.5)
        hist.SetMarkerColor(1)
        # if with_fits:
        #    fit = self.fits.get_bin_contents([threshold])
        #    fits.append(fit)
        hists.append(hist)
        labels.append("HW")

        emu_hist = emu_plotter.plots.get_bin_contents([bn.Base.everything])
        emu_hist = cumulative_hist(emu_hist)

        emu_hist.drawstyle = "EP"
        emu_hist.SetMarkerSize(0.5)
        emu_hist.SetMarkerColor(2)
        # if with_fits:
        #    emu_fit = self.fits.get_bin_contents([threshold])
        #    fits.append(emu_fit)
        hists.append(emu_hist)
        labels.append("Emu")

        self.__make_overlay(hists, fits, labels, "# Events", setlogy=True)

    def __make_overlay(self, hists, fits, labels, ytitle, suffix="", setlogy=False):
        with preserve_current_style():
            # Draw each resolution (with fit)
            xtitle = self.online_title
            canvas = draw(hists, draw_args={
                          "xtitle": xtitle, "ytitle": ytitle, "logy": setlogy})
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
                rightmargin=0.2,
                leftmargin=0.8,
                entryheight=0.028
            )
            for hist, label in zip(hists, labels):
                legend.AddEntry(hist, label)
            legend.SetBorderSize(0)
            legend.Draw()

            # Save canvas to file
            name = self.filename_format.format(pileup="all")
            self.save_canvas(canvas, name + suffix)

    def _is_consistent(self, new):
        """
        Check the two plotters are the consistent, so same binning and same axis names
        """
        return all([self.pileup_bins.bins == new.pileup_bins.bins,
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
        for (pileup, ), hist in self.plots.flat_items():
            human_readable_threshold = '{0}; PU > {1}'.format(self.online_title, self.pileup_bins.bins[pileup])
            rhist = hist.rebinned(summary_bins)
            stats = {}
            summary_columns = self._summary_columns(summary_bins, summary_label)
            for summary_column, y in zip(summary_columns, rhist.y()):
                stats[summary_column] = y
            total = sum(stats.values())
            overflow = rhist.integral(overflow=True) - total
            header = dict(identifier=human_readable_threshold, total=total, overflow=overflow)
            header.update(stats)
            yield header
