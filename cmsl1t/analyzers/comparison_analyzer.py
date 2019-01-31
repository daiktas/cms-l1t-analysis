from BaseAnalyzer import BaseAnalyzer
from cmsl1t.producers import jets
from cmsl1t.producers.match import get_matched_l1_jet
import numpy as np
from cmsl1t.energySums import EnergySum, Met

from cmsl1t.plotting.resolution import ResolutionPlot

#inherit from some cmsl1t class ideally
class Resolution():
    def __init__(self, online, offline, low=-1, high=2, n_bins=100, resolution_type="energy", pu_bins=[0]):
        self.low = low
        self.high = high
        self.res = ResolutionPlot(resolution_type, online, offline)
        self.res.create_histograms(online, offline, pu_bins, n_bins, low, high, legend_title='legend_title')
        self.res.set_plot_output_cfg(".", "pdf")
    def draw(self):
        self.res.draw()

    def fill(self, pileup, online, offline):
        self.res.fill(pileup, offline, online) 

resolutions = [ 
                Resolution("calo MET", "L1 MET"),
                Resolution("calo MET BE", "L1 MET"),
                Resolution("calo HT", "L1 HT"),
                Resolution("PF HT", "L1 HT"),
                Resolution("calo jet pT", "L1 jet pT"),
                Resolution("PF jet pT", "L1 jet pT"),
                ]


class Analyzer(BaseAnalyzer):

    def __init__(self, **kwargs):
        super(Analyzer, self).__init__(**kwargs)

    def prepare_for_events(self, reader):
        return True

    def reload_histograms(self, input_file):
        return True

    def fill_histograms(self, entry, event):

        if not event.MetFilters_hbheNoiseFilter:
            return True

        pileup = event['Vertex_nVtx']

        caloHT = EnergySum(event['Sums_caloHt'])
        pfHT = EnergySum(event['Sums_Ht'])
        caloMETBE = EnergySum(event['Sums_caloMetBE'], event['Sums_caloMetPhiBE'])
        # check this
        caloMETHF = Met(event['Sums_caloMet'], event['Sums_caloMetPhi'])
        pfMET_NoMu = Met(event['Sums_pfMetNoMu'], event['Sums_pfMetNoMuPhi'])

        l1HT = event['l1Sums_Htt']
        l1MET = event['l1Sums_Met']
        l1METHF = event['l1Sums_MetHF']
        l1MET_NoMu = event['l1Sums_MetHF']

        goodPFJets = event['goodPFJets']
        caloJets = event['caloJets']
        l1Jets = event['l1Jets']

        print caloMet.et

        resolutions[0].fill(pileup, caloMet, l1Met)
        resolutions[1].fill(pileup, caloMetBE, l1Met)
        resolutions[2].fill(pileup, caloHt, l1Ht)
        resolutions[3].fill(pileup, Ht, l1Ht)

        for jet in caloJets:
            matched_l1_jet = get_matched_l1_jet(jet, l1Jets, deltaR=0.4)
            if matched_l1_jet:
                resolutions[4].fill(pileup, jet['et'], matched_l1_jet['et'])

        for jet in goodPFJets:
            matched_l1_jet = get_matched_l1_jet(jet, l1Jets, deltaR=0.4)
            if matched_l1_jet:
                resolutions[5].fill(pileup, jet['et'], matched_l1_jet['et'])

        return True

    def write_histograms(self):
        #self.efficiencies.to_root(self.get_histogram_filename())
        return True

    def make_plots(self):

        [res.draw() for res in resolutions]
        # Something like this needs to be implemented still
        # self.efficiencies.draw_plots(self.output_folder, "png")
        return True
