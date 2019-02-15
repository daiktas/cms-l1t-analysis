from BaseAnalyzer import BaseAnalyzer
from cmsl1t.producers import jets
from cmsl1t.producers.match import get_matched_l1_jet
import numpy as np
from cmsl1t.energySums import EnergySum, Met

import cmsl1t.hist.binning as bn
from cmsl1t.plotting.resolution import ResolutionPlot
from cmsl1t.plotting.efficiency import EfficiencyPlot
from cmsl1t.hist.hist_collection import HistogramCollection
from cmsl1t.utils.hist import normalise_to_unit_area
from cmsl1t.utils.draw import draw, label_canvas

from rootpy.context import preserve_current_style
from rootpy.plotting import Legend
from rootpy import asrootpy, ROOT

import random

#inherit from some cmsl1t class ideally

class Resolution(ResolutionPlot):

    def create_histograms(self,
                          online_title, offline_title,
                          pileup_bins, n_bins, low, high, legend_title=''):
        """ This is not in an init function so that we can by-pass this in the
        case where we reload things from disk """
        self.online_title = online_title
        self.offline_title = offline_title
        self.legend_title = legend_title
        self.pileup_bins = bn.Sorted(pileup_bins, "pileup",
                                     use_everything_bin=True)

        name = ["resolution", self.online_name,
                self.offline_name, "pu_{pileup}"]

        name = "__".join(name)
        title = " ".join(
            [self.online_name, "vs.", self.offline_name, "in PU bin: {pileup}"])
        title = ";".join([title, self.offline_title, self.online_title])
        print self.pileup_bins

        self.plots = HistogramCollection([self.pileup_bins],
                                         "Hist1D", n_bins, low, high,
                                         name=name, title=title)

        name = name.join("__MC")
        title = title.join("__MC")

        self.plots_mc = HistogramCollection([self.pileup_bins],
                                         "Hist1D", n_bins, low, high,
                                         name=name, title=title)

        #self.filename_format = name.replace('/','_').replace('#','').replace('pT', 'p_{T}').replace('HT', 'H_{T}')
        self.filename_format = name.replace('/', '_')

    def __init__(self, offline_name, online_name, low=-1, high=2, n_bins=100, resolution_type="energy", pu_bins=[0]):
        super(Resolution, self).__init__(resolution_type, online_name, offline_name)
        print dir(__builtins__)
        self.low = low
        self.high = high
        self.create_histograms(online_name, offline_name, pu_bins, n_bins, low, high, legend_title='')

    def __make_overlay(self, hists, fits, labels, ytitle, suffix=""):
        with preserve_current_style():
            # Draw each resolution (with fit)
            # TODO: this feels like it does not belong here
            for hist in hists:
                hist.GetYaxis().SetRangeUser(0, 0.1)
                hist.GetYaxis().SetTitleOffset(1.4)

            xtitle = self.resolution_method.label.format(
                on=self.online_title, off=self.offline_title)
            canvas = draw(hists, draw_args={
                          "xtitle": xtitle, "ytitle": ytitle})
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
                rightmargin=0.3,
                leftmargin=0.7,
                textsize=0.03,
                entryheight=0.03,
            )

            for hist, label in zip(hists, labels):
                legend.AddEntry(hist, label)
            legend.SetBorderSize(0)
            legend.Draw()

            ymax = 1.2 * max([hist.GetMaximum() for hist in hists])
            line = ROOT.TLine(0., 0., 0., ymax)
            line.SetLineColor(15)

            # Save canvas to file
            #name = name.format(pileup="all").replace('pT', 'p_{T}').replace('HT', 'H_{T}')

        #self.filename_format = name.replace('/','_').replace('#','').replace('pT', 'p_{T}').replace('HT', 'H_{T}')
            name = self.filename_format.format(pileup="all")
            self.save_canvas(canvas, name.replace('pT', 'p_{T}').replace('HT', 'H_{T}') + suffix)


    def output(self,output_folder):
        self.set_plot_output_cfg(output_folder, "png")

    def fill(self, label, pileup, offline, online):
        difference = self.resolution_method(online, offline)
        if label == 1:
            self.plots_mc[pileup].fill(difference)
        else:
            self.plots[pileup].fill(difference)

    def draw(self, with_fits=False):
        hists = []
        labels = []
        fits = []
        for (pile_up, ), hist in self.plots.flat_items_all():
            if pile_up == bn.Base.everything:
                # hist.linestyle = "dashed"
                hist.drawstyle = ResolutionPlot.drawstyle
                label = "data"
            elif isinstance(pile_up, int):
                hist.drawstyle = ResolutionPlot.drawstyle
                if self.pileup_bins.get_bin_upper(pile_up) < 500:
                    label = "{:.0f} \\leq PU < {:.0f}".format(
                        self.pileup_bins.get_bin_lower(pile_up),
                        self.pileup_bins.get_bin_upper(pile_up),
                    )
                else:
                    label = "{:.0f} < PU".format(
                        self.pileup_bins.get_bin_lower(pile_up))
            else:
                continue
            hist.SetMarkerSize(0.5)
            hist.SetLineWidth(1)
            hists.append(hist)
            labels.append(label)

        normed_hists = list(normalise_to_unit_area(hists))
        for hist in normed_hists:
            hist.GetYaxis().SetRangeUser(-0.1, 1.1)

        for (pile_up, ), hist in self.plots_mc.flat_items_all():
            if pile_up == bn.Base.everything:
                # hist.linestyle = "dashed"
                hist.drawstyle = ResolutionPlot.drawstyle
                label = "MC"
            elif isinstance(pile_up, int):
                hist.drawstyle = ResolutionPlot.drawstyle
                if self.pileup_bins.get_bin_upper(pile_up) < 500:
                    label = "{:.0f} \\leq PU < {:.0f}".format(
                        self.pileup_bins.get_bin_lower(pile_up),
                        self.pileup_bins.get_bin_upper(pile_up),
                    )
                else:
                    label = "{:.0f} < PU".format(
                        self.pileup_bins.get_bin_lower(pile_up))
            else:
                continue
            hist.SetMarkerSize(0.5)
            hist.SetLineWidth(1)
            hists.append(hist)
            labels.append(label)

        normed_hists = list(normalise_to_unit_area(hists))
        for hist in normed_hists:
            hist.GetYaxis().SetRangeUser(-0.1, 1.1)

        self.__make_overlay(normed_hists, fits, labels, "a.u.")

class Efficiency(EfficiencyPlot):

    def create_histograms(
            self, online_title, offline_title, pileup_bins, thresholds,
            n_bins, low, high=400, legend_title=""):
        """ This is not in an init function so that we can by-pass this in the
        case where we reload things from disk """
        self.online_title = online_title
        self.offline_title = offline_title
        self.pileup_bins = bn.Sorted(pileup_bins, "pileup",
                                     use_everything_bin=True)
        self.thresholds = bn.GreaterThan(thresholds, "threshold")
        self.legend_title = legend_title

        name = ["efficiency", self.online_name, self.offline_name]
        name += ["thresh_{threshold}", "pu_{pileup}"]
        name = "__".join(name)
        title = " ".join([self.online_name, " in PU bin: {pileup}",
                          "and passing threshold: {threshold}"])

        self.filename_format = name.replace('/', '_')

        def make_efficiency(labels):
            this_name = "efficiency" + name.format(**labels)
            this_title = title.format(**labels)
            '''Checking type of 'low' to see whether it's int (x-range minimum)
                    or array (bin edges) for constructing TEfficiency'''
            if isinstance(low, np.ndarray):
                eff = asrootpy(
                    ROOT.TEfficiency(this_name, this_title, n_bins, low)
                )
            else:
                eff = asrootpy(
                    ROOT.TEfficiency(this_name, this_title, n_bins, low, high)
                )
            eff.drawstyle = EfficiencyPlot.drawstyle
            return eff
        self.efficiencies = HistogramCollection(
            [self.pileup_bins, self.thresholds],
            make_efficiency
        )

        name = name.join("__MC")
        title = title.join("__MC")

        self.efficiencies_mc = HistogramCollection(
            [self.pileup_bins, self.thresholds],
            make_efficiency
        )


    def __init__(self, offline_name, online_name, thresholds, low=0, high=400, n_bins=100, pu_bins=[0]):
        super(Efficiency, self).__init__(online_name, offline_name)
        global output_folder
        self.low = low
        self.high = high
        self.create_histograms(online_name, offline_name, pu_bins, thresholds, n_bins, low, high)

    def fill(self, label, pileup, offline, online):

        if label == 1:
            efficiencies = {(pu, thresh): eff
                            for (pu,), thresholds in self.efficiencies_mc[pileup].items()
                            for thresh, eff in thresholds.items()}

        else:
            efficiencies = {(pu, thresh): eff
                            for (pu,), thresholds in self.efficiencies[pileup].items()
                            for thresh, eff in thresholds.items()}

       
        for (pu_bin, threshold_bin), efficiency in efficiencies.items():
            threshold = self.thresholds.get_bin_center(threshold_bin)
            passed = False
            if isinstance(threshold, str):
                continue
            elif online > threshold:
                passed = True
            efficiency.fill(passed, offline)

    def __make_overlay(self, pileup, threshold, hists, fits, labels, header):
        with preserve_current_style():
            name = self.filename_format.format(pileup=pileup,
                                               threshold=threshold)

            xmin = hists[0].GetTotalHistogram().GetBinLowEdge(1)
            xmax = hists[0].GetTotalHistogram().GetBinLowEdge(hists[0].GetTotalHistogram().GetNbinsX() + 1)

            # Draw each efficiency (with fit)
            draw_args = {"xtitle": self.offline_title, "ytitle": "Efficiency", "xlimits": [xmin, xmax]}

            canvas = draw(hists, draw_args=draw_args)
            if len(fits) > 0:
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
                rightmargin=0.3,
                leftmargin=0.7,
                textsize=0.025,
                entryheight=0.028,
            )

            for hist, label in zip(hists, labels):
                legend.AddEntry(hist, label)

            legend.SetBorderSize(0)
            legend.Draw()

            for val in [0.25, 0.5, 0.75, 0.95, 1.]:
                line = ROOT.TLine(xmin, val, xmax, val)
                line.SetLineStyle("dashed")
                line.SetLineColor(15)
                line.Draw()

            for thresh in self.thresholds.bins:
                line = ROOT.TLine(thresh, 0., thresh, 1.)
                line.SetLineStyle("dashed")
                line.SetLineColor(15)
                line.Draw()

            # Save canvas to file
            self.save_canvas(canvas, name.replace('pT', 'p_{T}').replace('HT', 'H_{T}'))


    def draw(self, with_fits=False):
        # Fit the efficiencies if requested
        if with_fits:
            self.__fit_efficiencies()

        # Overlay the "all" piloe-up bin for each threshold

        all_pileup_effs = self.efficiencies.get_bin_contents(
            [bn.Base.everything])

        all_pileup_mc_effs = self.efficiencies_mc.get_bin_contents(
            [bn.Base.everything])

        hists = []
        labels = []
        fits = []
        label_template = '{label}: {online_title} > {threshold} GeV'

        for threshold in all_pileup_effs.iter_all():

            if not isinstance(threshold, int):
                continue
            hist = all_pileup_effs.get_bin_contents(threshold)
            hist.drawstyle = EfficiencyPlot.drawstyle_data
            self._dynamic_bin(hist)
            hists.append(hist)

            label = label_template.format(
                label="data",
                online_title=self.online_title,
                threshold=self.thresholds.bins[threshold],
            )
            labels.append(label)
            if with_fits:
                fits.append(self.fits.get_bin_contents(
                    [bn.Base.everything, threshold]))

        for threshold in all_pileup_mc_effs.iter_all():
            if not isinstance(threshold, int):
                continue
            hist = all_pileup_mc_effs.get_bin_contents(threshold)
            hist.drawstyle = EfficiencyPlot.drawstyle_data
            self._dynamic_bin(hist)
            hists.append(hist)

            label = label_template.format(
                label="mc",
                online_title=self.online_title,
                threshold=self.thresholds.bins[threshold],
            )
            labels.append(label)
            if with_fits:
                fits.append(self.fits.get_bin_contents(
                    [bn.Base.everything, threshold]))


        self.__make_overlay(
            "all", "all", hists, fits, labels, self.online_title,
        )

        if 'HiRange' not in self.filename_format:
            # Overlay individual pile-up bins for each threshold
            for threshold in self.thresholds:
                hists = []
                labels = []
                fits = []
                for pileup in self.pileup_bins:
                    if not isinstance(pileup, int):
                        continue
                    hist = self.efficiencies.get_bin_contents([pileup, threshold])
                    hist.drawstyle = EfficiencyPlot.drawstyle_data
                    self._dynamic_bin(hist)
                    hists.append(hist)
                    if with_fits:
                        fits.append(self.fits.get_bin_contents(
                            [pileup, threshold]))
                    labels.append(str(self.pileup_bins.bins[pileup]))
                    self.__make_overlay(pileup, threshold, hists,
                                    fits, labels, "PU bin")

        # Produce the fit summary plot
        if with_fits:
            self.__summarize_fits()

    def output(self, output_folder):
        self.set_plot_output_cfg(output_folder, "png")

ALL_THRESHOLDS = dict(
    HT=[120, 200],
    METBE=[80, 120],
    METHF=[80, 120],
    JetET=[30, 60]
)

class Analyzer(BaseAnalyzer):

    def __init__(self, **kwargs):
        super(Analyzer, self).__init__(**kwargs)
        self.output_folder = kwargs["output_folder"]

        self.object_dict = {}

        self.object_dict["METBE_res"] = Resolution("calo MET BE", "L1 MET")
        self.object_dict["METBE_eff"] = Efficiency("calo MET BE / GeV", "L1 MET", ALL_THRESHOLDS["METBE"])

        self.object_dict["METHF_res"] = Resolution("calo MET HF", "L1 MET")
        self.object_dict["METHF_eff"] = Efficiency("calo MET HF / GeV", "L1 MET", ALL_THRESHOLDS["METHF"])

        self.object_dict["caloHT_res"] = Resolution("calo HT", "L1 HT", low=-1, high=3)
        self.object_dict["caloHT_eff"] = Efficiency("calo HT / GeV", "L1 HT", ALL_THRESHOLDS["HT"])

        self.object_dict["pfHT_res"] = Resolution("PF HT", "L1 HT", low=-1, high=3)
        self.object_dict["pfHT_eff"] = Efficiency("PF HT / GeV", "L1 HT", ALL_THRESHOLDS["HT"])

        self.object_dict["caloJetEtBE_res"] = Resolution("calo jet pT BE", "L1 jet pT")
        self.object_dict["caloJetEtBE_eff"] = Efficiency("calo jet pT BE / GeV", "L1 jet pT", ALL_THRESHOLDS["JetET"])

        self.object_dict["caloJetEtHF_res"] = Resolution("calo jet pT HF", "L1 jet pT")
        self.object_dict["caloJetEtHF_eff"] = Efficiency("calo jet pT HF / GeV", "L1 jet pT", ALL_THRESHOLDS["JetET"])

        self.object_dict["pfJetEtBE_res"] = Resolution("PF jet pT BE", "L1 jet pT")
        self.object_dict["pfJetEtBE_eff"] = Efficiency("PF jet pT BE / GeV", "L1 jet pT", ALL_THRESHOLDS["JetET"])

        self.object_dict["pfJetEtHF_res"] = Resolution("PF jet pT HF", "L1 jet pT")
        self.object_dict["pfJetEtHF_eff"] = Efficiency("PF jet pT HF / GeV", "L1 jet pT", ALL_THRESHOLDS["JetET"])

        self.object_dict["METBEPHI_res"] = Resolution("calo MET BE #phi", "L1 MET #phi", resolution_type="phi")
        self.object_dict["METHFPHI_res"] = Resolution("calo MET HF #phi", "L1 MET #phi", resolution_type="phi")

        self.object_dict["caloJetPhiBE_res"] = Resolution("calo jet BE #phi", "L1 jet #phi", resolution_type="phi", low=-0.5, high=0.5)
        self.object_dict["pfJetPhiBE_res"] = Resolution("PF jet BE #phi", "L1 jet #phi", resolution_type="phi", low=-0.5, high=0.5)

        self.object_dict["caloJetPhiHF_res"] = Resolution("calo jet HF #phi", "L1 jet #phi", resolution_type="phi", low=-0.5, high=0.5)
        self.object_dict["pfJetPhiHF_res"] = Resolution("PF jet HF #phi", "L1 jet #phi", resolution_type="phi", low=-0.5, high=0.5)

        self.object_dict["caloJetEtaBE_res"] = Resolution("calo jet BE #eta", "L1 jet #eta", resolution_type="eta", low=-0.5, high=0.5)
        self.object_dict["pfJetEtaBE_res"] = Resolution("PF jet BE #eta", "L1 jet #eta", resolution_type="eta", low=-0.5, high=0.5)

        self.object_dict["caloJetEtaHF_res"] = Resolution("calo jet HF #eta", "L1 jet #eta", resolution_type="eta", low=-0.5, high=0.5)
        self.object_dict["pfJetEtaHF_res"] = Resolution("PF jet HF #eta", "L1 jet #eta", resolution_type="eta", low=-0.5, high=0.5)

    def prepare_for_events(self, reader):
        return True

    def reload_histograms(self, input_file):
        return True

    def fill_histograms(self, entry, event):

        # cut on trigger
        # log scale
        # cut on muon in MC

        pileup = event['Vertex_nVtx']

        run = event['run']

        hlts = event['hlt']

        hlt_passed = False

        for hlt in hlts:
            if "HLT_IsoMu24" in hlt.Data():
                hlt_passed = True
                break

        if not hlt_passed:
            return True
            
        #if run == 1 and event["L1Upgrade_nMuons"] == 0:
            #return True

        if run != 1 and not event.MetFilters_hbheNoiseFilter:
            print "filters"
            return True

        caloHT = EnergySum(event['Sums_caloHt'])
        pfHT = EnergySum(event['Sums_Ht'])

        caloMETBE = Met(event['Sums_caloMetBE'], event['Sums_caloMetPhiBE'])
        caloMETHF = Met(event['Sums_caloMet'], event['Sums_caloMetPhi'])
        pfMET_NoMu = Met(event['Sums_pfMetNoMu'], event['Sums_pfMetNoMuPhi'])

        l1HT = event['l1Sums_Htt']
        l1MET = event['l1Sums_Met']
        l1METHF = event['l1Sums_MetHF']
        l1MET_NoMu = event['l1Sums_MetHF']

        goodPFJets = event['goodPFJets']
        caloJets = event['caloJets']
        l1Jets = event['l1Jets']

        if pfHT < 120. or pfMET_NoMu < 40.:
            return True

        label = run

        self.object_dict["METBE_res"].fill(label, pileup, caloMETBE.et, l1MET.et)
        self.object_dict["METBE_eff"].fill(label, pileup, caloMETBE.et, l1MET.et)
        self.object_dict["METHF_res"].fill(label, pileup, caloMETHF.et, l1METHF.et)
        self.object_dict["METHF_eff"].fill(label, pileup, caloMETHF.et, l1METHF.et)
        self.object_dict["METBEPHI_res"].fill(label, pileup, caloMETBE.phi, l1MET.phi)
        self.object_dict["METHFPHI_res"].fill(label, pileup, caloMETHF.phi, l1METHF.phi)

        if caloHT.et != 0 and pfHT.et != 0 and l1HT.et != 0:
            self.object_dict["caloHT_res"].fill(label, pileup, caloHT.et, l1HT.et)
            self.object_dict["caloHT_eff"].fill(label, pileup, caloHT.et, l1HT.et)
            self.object_dict["pfHT_res"].fill(label, pileup, caloHT.et, l1HT.et)
            self.object_dict["pfHT_eff"].fill(label, pileup, caloHT.et, l1HT.et)

        if caloJets:
            jet = caloJets[0]
            matched_l1_jet = get_matched_l1_jet(jet, l1Jets, deltaR=0.3)
            if matched_l1_jet and jet.etCorr > 30.:
                if abs(jet.eta) < 3.:
                    self.object_dict["caloJetEtBE_eff"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["caloJetEtBE_res"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["caloJetPhiBE_res"].fill(label, pileup, jet['phi'], matched_l1_jet['phi'])
                    self.object_dict["caloJetEtaBE_res"].fill(label, pileup, jet['eta'], matched_l1_jet['eta'])
                else:
                    self.object_dict["caloJetEtHF_eff"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["caloJetEtHF_res"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["caloJetPhiHF_res"].fill(label, pileup, jet['phi'], matched_l1_jet['phi'])
                    self.object_dict["caloJetEtaHF_res"].fill(label, pileup, jet['eta'], matched_l1_jet['eta'])


        if goodPFJets:
            jet = goodPFJets[0]
            matched_l1_jet = get_matched_l1_jet(jet, l1Jets, deltaR=0.3)
            if matched_l1_jet and jet.etCorr > 30.:
                if abs(jet.eta) < 3.:
                    self.object_dict["pfJetEtBE_eff"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["pfJetEtBE_res"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["pfJetPhiBE_res"].fill(label, pileup, jet['phi'], matched_l1_jet['phi'])
                    self.object_dict["pfJetEtaBE_res"].fill(label, pileup, jet['eta'], matched_l1_jet['eta'])
                else: 
                    self.object_dict["pfJetEtHF_eff"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["pfJetEtHF_res"].fill(label, pileup, jet.etCorr, matched_l1_jet.etCorr)
                    self.object_dict["pfJetPhiHF_res"].fill(label, pileup, jet['phi'], matched_l1_jet['phi'])
                    self.object_dict["pfJetEtaHF_res"].fill(label, pileup, jet['eta'], matched_l1_jet['eta'])

        return True

    def write_histograms(self):
        #self.efficiencies.to_root(self.get_histogram_filename())
        return True

    def make_plots(self):

        for key, obj in self.object_dict.items():
            obj.output(self.output_folder)
            obj.draw()

        # Something like this needs to be implemented still
        # self.efficiencies.draw_plots(self.output_folder, "png")
        return True
