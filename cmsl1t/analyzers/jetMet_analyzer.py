import os
from cmsl1t.analyzers.BaseAnalyzer import BaseAnalyzer
from cmsl1t.plotting.efficiency import EfficiencyPlot
from cmsl1t.collections import EfficiencyCollection
from cmsl1t.plotting.onlineVsOffline import OnlineVsOffline
from cmsl1t.plotting.resolution import ResolutionPlot
from cmsl1t.plotting.resolution_vs_X import ResolutionVsXPlot
from cmsl1t.playground.jetfilters import pfJetFilter
# from cmsl1t.playground.metfilters import pfMetFilter
from cmsl1t.filters import pfMetFilter
from cmsl1t.filters import LuminosityFilter
import cmsl1t.recalc.met as recalc
from cmsl1t.energySums import EnergySum, Met
from math import pi
import pprint
from collections import namedtuple
import numpy as np
import csv
import cmsl1t
from cmsl1t.jet import match


def types(doEmu, doReco, doGen):

    sum_types = []
    jet_types = []
    if doReco:
        sum_types.extend(["caloHT", "pfHT", "caloMETHF",
                          "caloMETBE", "pfMET_NoMu"])
        jet_types.extend(["caloJetET_BE", "caloJetET_HF",
                          "pfJetET_BE", "pfJetET_HF"])
    if doGen:
        sum_types.extend(["genHT", "genMETHF", "genMETBE"])
        jet_types.extend(["genJetET_BE", "genJetET_HF"])
    if doEmu:
        sum_types += [t + '_Emu' for t in sum_types]
        jet_types += [t + '_Emu' for t in jet_types]

    return sum_types, jet_types


# Eta ranges so we can put |\eta| < val as the legend header on the
# efficiency plots.
ETA_RANGES = dict(
    HT="|\\eta| < 2.4",
    METBE="|\\eta| < 3.0",
    METHF="|\\eta| < 5.0, 50 \\leq PU < 60",
    JetET_BE="|\\eta| < 3.0",
    JetET_HF="3.0 < |\\eta| < 5.0",
)

ALL_THRESHOLDS = dict(
    HT=[120, 200, 320],
    METBE=[80, 100, 120],
    METHF=[80, 100, 120],
    JetET=[35, 90, 120]
)


def extractSums(event, doEmu, doReco, doGen):
    offline = dict()
    online = dict()
    if doReco:
        offline.update(dict(
            caloHT=EnergySum(event.Sums_caloHt),
            pfHT=EnergySum(event.Sums_Ht),
            caloMETBE=Met(event.Sums_caloMetBE, event.Sums_caloMetPhiBE),
            caloMETHF=Met(event.Sums_caloMet, event.Sums_caloMetPhi),
            pfMET_NoMu=Met(event.Sums_pfMetNoMu, event.Sums_pfMetNoMuPhi),
        ))
        online.update(dict(
            caloHT=event.l1Sums_Htt,
            pfHT=event.l1Sums_Htt,
            caloMETBE=event.l1Sums_Met,
            caloMETHF=event.l1Sums_MetHF,
            pfMET_NoMu=event.l1Sums_MetHF,
        ))
        if doEmu:
            offline.update(dict(
                caloHT_Emu=EnergySum(event.Sums_caloHt),
                pfHT_Emu=EnergySum(event.Sums_Ht),
                caloMETBE_Emu=Met(event.Sums_caloMetBE,
                                  event.Sums_caloMetPhiBE),
                caloMETHF_Emu=Met(event.Sums_caloMet, event.Sums_caloMetPhi),
                pfMET_NoMu_Emu=Met(event.Sums_pfMetNoMu,
                                   event.Sums_pfMetNoMuPhi),
            ))
            online.update(dict(
                caloHT_Emu=event.l1EmuSums_Htt,
                pfHT_Emu=event.l1EmuSums_Htt,
                caloMETBE_Emu=event.l1EmuSums_Met,
                caloMETHF_Emu=event.l1EmuSums_MetHF,
                pfMET_NoMu_Emu=event.l1EmuSums_MetHF,
            ))

    if doGen:
        offline.update(dict(
            genHT=event.genSums_HT,
            genMETHF=event.genSums_MetHF,
            genMETBE=event.genSums_MetBE,
        ))
        online.update(dict(
            genHT=event.l1Sums_Htt,
            genMETHF=event.l1Sums_MetHF,
            genMETBE=event.l1Sums_Met,
        ))
        if doEmu:
            offline.update(dict(
                genHT_Emu=event.genSums_HT,
                genMETHF_Emu=event.genSums_MetHF,
                genMETBE_Emu=event.genSums_MetBE,
            ))
            online.update(dict(
                genHT_Emu=event.l1EmuSums_Htt,
                genMETHF_Emu=event.l1EmuSums_MetHF,
                genMETBE_Emu=event.l1EmuSums_Met,
            ))
    return offline, online


class Analyzer(BaseAnalyzer):

    def __init__(self, **kwargs):
        super(Analyzer, self).__init__(**kwargs)

        lumiMuDict = dict()
	run_lumi_csv = os.path.join(cmsl1t.PROJECT_ROOT, 'run_lumi.csv')
        with open(run_lumi_csv, 'rb') as runLumiFile:
            reader = csv.reader(runLumiFile, delimiter=',')
            for line in reader:
                lumiMuDict[(int(line[1]),int(line[2]))] = float(line[3])
        self._lumiMu = lumiMuDict

        self._lumiFilter = None
        self._lumiJson = self.params['lumiJson']
        if self._lumiJson:
            self._lumiFilter = LuminosityFilter(self._lumiJson)

        self._lastRunAndLumi = (-1, -1)
        self._processLumi = True

        # TODO: this needs changing, these should be analyser parameters
        # or even move out into separate calls of the same analyzer
        loaded_trees = self.params['load_trees']
        self._doVertex = 'recoTree' in loaded_trees
        self._doEmu = 'emuUpgrade' in loaded_trees
        self._doReco = 'recoTree' in loaded_trees
        self._doGen = 'genTree' in loaded_trees

        self._sumTypes, self._jetTypes = types(
            self._doEmu, self._doReco, self._doGen)

        for name in self._sumTypes:
            eff_plot = EfficiencyPlot("L1", "offline_" + name)
            res_plot = ResolutionPlot("energy", "L1", "offline_" + name)
            twoD_plot = OnlineVsOffline("L1", "offline_" + name)
            self.register_plotter(eff_plot)
            self.register_plotter(res_plot)
            self.register_plotter(twoD_plot)
            setattr(self, name + "_eff", eff_plot)
            setattr(self, name + "_res", res_plot)
            setattr(self, name + "_2D", twoD_plot)

            eff_plot_HR = EfficiencyPlot("L1", "offline_" + name + "_HiRange")
            twoD_plot_HR = OnlineVsOffline(
                "L1", "offline_" + name + "_HiRange")
            self.register_plotter(eff_plot_HR)
            self.register_plotter(twoD_plot_HR)
            setattr(self, name + "_eff_HR", eff_plot_HR)
            setattr(self, name + "_2D_HR", twoD_plot_HR)

        for angle in self._sumTypes:
            name = angle + "_phi"
            if 'HT' in angle:
                continue
            res_plot = ResolutionPlot("phi", "L1", "offline_" + name)
            twoD_plot = OnlineVsOffline("L1", "offline_" + name)
            self.register_plotter(res_plot)
            self.register_plotter(twoD_plot)
            setattr(self, name + "_res", res_plot)
            setattr(self, name + "_2D", twoD_plot)

        for name in self._jetTypes:
            eff_plot = EfficiencyPlot("L1", "offline_" + name)
            res_plot = ResolutionPlot("energy", "L1", "offline_" + name)
            twoD_plot = OnlineVsOffline("L1", "offline_" + name)
            self.register_plotter(eff_plot)
            self.register_plotter(res_plot)
            self.register_plotter(twoD_plot)
            setattr(self, name + "_eff", eff_plot)
            setattr(self, name + "_res", res_plot)
            setattr(self, name + "_2D", twoD_plot)

            eff_plot_HR = EfficiencyPlot("L1", "offline_" + name + "_HiRange")
            twoD_plot_HR = OnlineVsOffline(
                "L1", "offline_" + name + "_HiRange")
            self.register_plotter(eff_plot_HR)
            self.register_plotter(twoD_plot_HR)
            setattr(self, name + "_eff_HR", eff_plot_HR)
            setattr(self, name + "_2D_HR", twoD_plot_HR)

        self.res_vs_eta_CentralJets = ResolutionVsXPlot(
            "energy", "", "", "pfJet_eta")
        self.register_plotter(self.res_vs_eta_CentralJets)
        self.res_vs_eta_CentralGenJets = ResolutionVsXPlot(
            "energy", "", "", "genJet_eta")
        self.register_plotter(self.res_vs_eta_CentralGenJets)

    def prepare_for_events(self, reader):
        puBins = self.params['pu_bins']
        puBins_HR = [0, 999]

        Config = namedtuple(
            "Config",
            "name off_title on_title min max",
        )
        cfgs = []
        if self._doReco:
            cfgs.extend([
                Config("caloHT", "Offline Calo HT", "L1 HT", 30, 830),
                Config("pfHT", "Offline PF HT", "L1 HT", 30, 830),
                Config("caloMETHF", "Offline Calo MET HF", "L1 MET HF", 0, 400),
                Config("caloMETBE", "Offline Calo MET BE", "L1 MET BE", 0, 400),
                Config("pfMET_NoMu", "Offline PF MET NoMu",
                       "L1 MET HF", 0, 400),
                Config("caloJetET_BE", "Offline Central Calo Jet ET",
                       "L1 Jet ET", 20, 420),
                Config("caloJetET_HF", "Offline Forward Calo Jet ET",
                       "L1 Jet ET", 20, 420),
                Config("pfJetET_BE", "Offline Central PF Jet ET",
                       "L1 Jet ET", 20, 420),
                Config("pfJetET_HF", "Offline Forward PF Jet ET",
                       "L1 Jet ET", 20, 420),
            ])
        if self._doGen:
            cfgs.extend([
                Config("genHT", "Gen HT", "L1 HT", 30, 830),
                Config("genMETHF", "Gen MET HF", "L1 MET HF", 0, 400),
                Config("genMETBE", "Gen MET BE", "L1 MET BE", 0, 400),
                Config("genJetET_BE", "Central Gen Jet ET",
                       "L1 Jet ET", 20, 420),
                Config("genJetET_HF", "Forward Gen Jet ET",
                       "L1 Jet ET", 20, 420),
            ])

        self._plots_from_cfgs(cfgs, puBins)
        if self._doEmu:
            self._plots_from_cfgs(cfgs, puBins, emulator=True)
        self._plots_from_cfgs(cfgs, puBins_HR, high_range=True)
        if self._doEmu:
            self._plots_from_cfgs(
                cfgs, puBins_HR, emulator=True, high_range=True)

        self.res_vs_eta_CentralJets.build(
            "Online Jet energy (GeV)",
            "Offline Jet energy (GeV)",
            "Offline Jet Eta (rad)",
            puBins,
            50, -0.5, 3.5, 50, -5.0, 5.0,
        )
        self.res_vs_eta_CentralGenJets.build(
            "Online Jet energy (GeV)",
            "Offline Jet energy (GeV)",
            "Offline Jet Eta (rad)",
            puBins,
            50, -0.5, 3.5, 50, -5.0, 5.0,
        )
        return True

    def _plots_from_cfgs(self, cfgs, puBins, emulator=False, high_range=False):
        suffix = ""
        prefix = ""
        if high_range:
            suffix = '_HR'
        if emulator:
            prefix = '_Emu'

        if 'thresholds' in self.params:
            allThresholds = self.params['thresholds']
        else:
            allThresholds = ALL_THRESHOLDS

        for cfg in cfgs:

            eff_plot = getattr(self, cfg.name + prefix + "_eff" + suffix)
            twoD_plot = getattr(self, cfg.name + prefix + "_2D" + suffix)

            thresholds = []
            for l1trig, thresh in allThresholds.items():
                if emulator and "Emu" not in l1trig or not emulator and "Emu" in l1trig:
                    continue
                if l1trig.replace("_Emu", "") in cfg.name:
                    thresholds = thresh
                    break
            if 'pfMET' in cfg.name:
                if emulator:
                    thresholds = allThresholds.get("METHF_Emu")
                else:
                    thresholds = allThresholds.get("METHF")

            etaRange = ""
            for l1trig, etaRange in ETA_RANGES.items():
                if l1trig in cfg.name:
                    etaRange = etaRange
                    break
            if 'pfMET' in cfg.name:
                etaRange = ETA_RANGES.get("METHF")

            params = [
                cfg.on_title, cfg.off_title + " (GeV)", puBins, thresholds,
                80, cfg.min, cfg.max,
            ]
            if high_range:
                params = [
                    cfg.on_title, cfg.off_title + " (GeV)", puBins, thresholds,
                    105, 0, 2100,
                ]
                if "HT" in cfg.name:
                    params = [
                        cfg.on_title, cfg.off_title +
                        " (GeV)", puBins, thresholds,
                        105, 30, 2130,
                    ]
                if "JetET" in cfg.name:
                    params = [
                        cfg.on_title, cfg.off_title +
                        " (GeV)", puBins, thresholds,
                        105, 20, 2120,
                    ]
                if "MET" in cfg.name:
                    params = [
                        cfg.on_title, cfg.off_title +
                        " (GeV)", puBins, thresholds,
                        100, 0, 2000,
                    ]

            eff_plot.build(*params, legend_title=etaRange)
            params.remove(thresholds)
            twoD_plot.build(*params)

            if high_range:
                continue
            res_plot = getattr(self, cfg.name + prefix + "_res" + suffix)
            res_plot.build(cfg.on_title, cfg.off_title,
                           puBins, 150, -3, 3, legend_title=etaRange)

            if not hasattr(self, cfg.name + prefix + "_phi_res"):
                continue
            res_plot = getattr(self, cfg.name + prefix + "_phi_res")
            twoD_plot = getattr(self, cfg.name + prefix + "_phi_2D")
            twoD_plot.build(
                cfg.on_title + " Phi (rad)",
                cfg.off_title + " Phi (rad)",
                puBins, 100,
                -pi,
                2 * pi,
            )
            res_plot.build(
                cfg.on_title + " Phi",
                cfg.off_title + " Phi",
                puBins,
                100,
                -2 * pi,
                2 * pi,
                legend_title=etaRange,
            )

    def fill_histograms(self, entry, event):

        if not self._passesLumiFilter(event.run, event.lumi):
            return True

        offline, online = extractSums(
            event, self._doEmu, self._doReco, self._doGen)

        recoNVtx = 1
        genNVtx = 1

        if self._doReco or self._doVertex:
            recoNVtx = event.Vertex_nVtx
        if self._doGen:
            genNVtx = event.Generator_nVtx

        pileup = self._lumiMu[(event['run'],event['lumi'])]
        #print pileup
        #if pileup >= 60 or pileup < 50:
        #    return True


        for name in self._sumTypes:
            if 'pfMET' in name and not pfMetFilter(event):
                continue
            on = online[name]
            off = offline[name]
            for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                if '_res' in suffix and (on.et < 0.01 or off.et < 30):
                    continue
                pileup = recoNVtx
                if 'gen' in name:
                    pileup = genNVtx
                getattr(self, name + suffix).fill(pileup, off.et, on.et)
            if hasattr(self, name + "_phi_res"):
                getattr(self, name + "_phi_res").fill(pileup, off.phi, on.phi)
                getattr(self, name + "_phi_2D").fill(pileup, off.phi, on.phi)

        if self._doReco and self._doEmu:
            goodRefJets = event.goodPFJets

            for refJet in goodRefJets:
                l1Jet = match(refJet, event.l1EmuJets)
                if not l1Jet:
                    continue
                if refJet.etCorr > 30.:
                    self.res_vs_eta_CentralJets.fill(
                        recoNVtx, refJet.eta, refJet.etCorr, l1Jet.et)

        if self._doGen and self._doEmu:
            goodRefJets = event.goodGenJets

            for refJet in goodRefJets:
                l1Jet = match(refJet, event.l1EmuJets)
                if not l1Jet:
                    continue
                if refJet.etCorr > 30.:
                    self.res_vs_eta_CentralGenJets.fill(
                        genNVtx, refJet.eta, refJet.etCorr, l1Jet.et)

        if self._doGen:
            leadingGenJet = None
            if event.goodGenJets:
                leadingGenJet = event.goodGenJets[0]

            if leadingGenJet:

                genFillRegions = []
                if abs(leadingGenJet.eta) < 3.0:
                    genFillRegions = ['BE']
                else:
                    genFillRegions = ['HF']

                if self._doEmu:
                    genL1EmuJet = match(leadingGenJet, event.l1EmuJets)
                    if genL1EmuJet:
                        genL1EmuJetEt = genL1EmuJet.et
                    else:
                        genL1EmuJetEt = 0.

                    for region in genFillRegions:
                        for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                            if '_res' in suffix and (genL1EmuJetEt == 0 or leadingGenJet.etCorr < 30):
                                continue
                            name = 'genJetET_{0}_Emu{1}'.format(region, suffix)
                            getattr(self, name).fill(
                                genNVtx, leadingGenJet.etCorr, genL1EmuJetEt,
                            )

                genL1Jet = match(leadingGenJet, event.l1Jets)
                if genL1Jet:
                    genL1JetEt = genL1Jet.et
                else:
                    genL1JetEt = 0.

                for region in genFillRegions:
                    for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                        if '_res' in suffix and (genL1JetEt == 0 or leadingGenJet.etCorr < 30):
                            continue
                        name = 'genJetET_{0}{1}'.format(region, suffix)
                        getattr(self, name).fill(
                            genNVtx, leadingGenJet.etCorr, genL1JetEt,
                        )

        if self._doReco:
            leadingPFJet, leadingCaloJet = None, None
            if event.goodPFJets:
                leadingPFJet = event.goodPFJets[0]
            if event.caloJets:
                leadingCaloJet = event.caloJets[0]

            if leadingPFJet and leadingPFJet.etCorr > 20:
                pfFillRegions = []
                if abs(leadingPFJet.eta) < 3.0:
                    pfFillRegions = ['BE']
                else:
                    pfFillRegions = ['HF']

                if self._doEmu:
                    pfL1EmuJet = match(leadingPFJet, event.l1EmuJets)
                    if pfL1EmuJet:
                        pfL1EmuJetEt = pfL1EmuJet.et
                    else:
                        pfL1EmuJetEt = 0.

                    for region in pfFillRegions:
                        for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                            if '_res' in suffix and (pfL1EmuJetEt == 0 or leadingPFJet.etCorr < 30):
                                continue
                            name = 'pfJetET_{0}_Emu{1}'.format(region, suffix)
                            getattr(self, name).fill(
                                recoNVtx, leadingPFJet.etCorr, pfL1EmuJetEt,
                            )

                pfL1Jet = match(leadingPFJet, event.l1Jets)
                if pfL1Jet:
                    pfL1JetEt = pfL1Jet.et
                else:
                    pfL1JetEt = 0.

                for region in pfFillRegions:
                    for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                        if '_res' in suffix and (pfL1JetEt == 0 or leadingPFJet.etCorr < 30):
                            continue
                        name = 'pfJetET_{0}{1}'.format(region, suffix)
                        getattr(self, name).fill(
                            recoNVtx, leadingPFJet.etCorr, pfL1JetEt,
                        )

            if leadingCaloJet and leadingCaloJet.etCorr > 20:

                caloFillRegions = []
                if abs(leadingCaloJet.eta) < 3.0:
                    caloFillRegions = ['BE']
                else:
                    caloFillRegions = ['HF']

                if self._doEmu:
                    caloL1EmuJet = match(leadingCaloJet, event.l1EmuJets)
                    if caloL1EmuJet:
                        caloL1EmuJetEt = caloL1EmuJet.et
                    else:
                        caloL1EmuJetEt = 0.

                    for region in caloFillRegions:
                        for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                            if '_res' in suffix and (caloL1EmuJetEt == 0 or leadingCaloJet.etCorr < 30):
                                continue
                            name = 'caloJetET_{0}_Emu{1}'.format(
                                region, suffix)
                            getattr(self, name).fill(
                                recoNVtx, leadingCaloJet.etCorr, caloL1EmuJetEt,
                            )

                caloL1Jet = match(leadingCaloJet, event.l1Jets)
                if caloL1Jet:
                    caloL1JetEt = caloL1Jet.et
                else:
                    caloL1JetEt = 0.

                for region in caloFillRegions:
                    for suffix in ['_eff', '_res', '_2D', '_eff_HR', '_2D_HR']:
                        if '_res' in suffix and (caloL1JetEt == 0 or leadingCaloJet.etCorr < 30):
                            continue
                        name = 'caloJetET_{0}{1}'.format(region, suffix)
                        getattr(self, name).fill(
                            recoNVtx, leadingCaloJet.etCorr, caloL1JetEt,
                        )

        return True

    def _passesLumiFilter(self, run, lumi):
        if self._lumiFilter is None:
            return True
        if (run, lumi) == self._lastRunAndLumi:
            return self._processLumi

        self._lastRunAndLumi = (run, lumi)
        self._processLumi = self._lumiFilter(run, lumi)

        return self._processLumi

    def make_plots(self):
        """
        Custom version, does what the normal one does but also overlays whatever you like.
        """
        # for plot in self.all_plots:
        #    plot.draw()

        if self._doReco:
            getattr(self, 'caloHT_eff').draw()
            getattr(self, 'pfHT_eff').draw()
            getattr(self, 'caloMETBE_eff').draw()
            getattr(self, 'caloMETHF_eff').draw()
            getattr(self, 'pfMET_NoMu_eff').draw()
            getattr(self, 'caloJetET_BE_eff').draw()
            getattr(self, 'caloJetET_HF_eff').draw()
            getattr(self, 'pfJetET_BE_eff').draw()
            getattr(self, 'pfJetET_HF_eff').draw()

            getattr(self, 'caloHT_eff_HR').draw()
            getattr(self, 'pfHT_eff_HR').draw()
            getattr(self, 'caloMETBE_eff_HR').draw()
            getattr(self, 'caloMETHF_eff_HR').draw()
            getattr(self, 'pfMET_NoMu_eff_HR').draw()
            getattr(self, 'caloJetET_BE_eff_HR').draw()
            getattr(self, 'caloJetET_HF_eff_HR').draw()
            getattr(self, 'pfJetET_BE_eff_HR').draw()
            getattr(self, 'pfJetET_HF_eff_HR').draw()

            getattr(self, 'caloHT_res').draw()
            getattr(self, 'pfHT_res').draw()
            getattr(self, 'caloMETBE_res').draw()
            getattr(self, 'caloMETHF_res').draw()
            getattr(self, 'pfMET_NoMu_res').draw()
            getattr(self, 'caloJetET_BE_res').draw()
            getattr(self, 'caloJetET_HF_res').draw()
            getattr(self, 'pfJetET_BE_res').draw()
            getattr(self, 'pfJetET_HF_res').draw()

            if self._doEmu:
                getattr(self, 'caloHT_eff').overlay_with_emu(
                    getattr(self, 'caloHT_Emu_eff'))
                getattr(self, 'pfHT_eff').overlay_with_emu(
                    getattr(self, 'pfHT_Emu_eff'))
                getattr(self, 'caloMETBE_eff').overlay_with_emu(
                    getattr(self, 'caloMETBE_Emu_eff'))
                getattr(self, 'caloMETHF_eff').overlay_with_emu(
                    getattr(self, 'caloMETHF_Emu_eff'))
                getattr(self, 'pfMET_NoMu_eff').overlay_with_emu(
                    getattr(self, 'pfMET_NoMu_Emu_eff'))
                getattr(self, 'caloJetET_BE_eff').overlay_with_emu(
                    getattr(self, 'caloJetET_BE_Emu_eff'))
                getattr(self, 'caloJetET_HF_eff').overlay_with_emu(
                    getattr(self, 'caloJetET_HF_Emu_eff'))
                getattr(self, 'pfJetET_BE_eff').overlay_with_emu(
                    getattr(self, 'pfJetET_BE_Emu_eff'))
                getattr(self, 'pfJetET_HF_eff').overlay_with_emu(
                    getattr(self, 'pfJetET_HF_Emu_eff'))

                getattr(self, 'caloHT_res').overlay_with_emu(
                    getattr(self, 'caloHT_Emu_res'))
                getattr(self, 'pfHT_res').overlay_with_emu(
                    getattr(self, 'pfHT_Emu_res'))
                getattr(self, 'caloMETBE_res').overlay_with_emu(
                    getattr(self, 'caloMETBE_Emu_res'))
                getattr(self, 'caloMETHF_res').overlay_with_emu(
                    getattr(self, 'caloMETHF_Emu_res'))
                getattr(self, 'pfMET_NoMu_res').overlay_with_emu(
                    getattr(self, 'pfMET_NoMu_Emu_res'))
                getattr(self, 'caloJetET_BE_res').overlay_with_emu(
                    getattr(self, 'caloJetET_BE_Emu_res'))
                getattr(self, 'caloJetET_HF_res').overlay_with_emu(
                    getattr(self, 'caloJetET_HF_Emu_res'))
                getattr(self, 'pfJetET_BE_res').overlay_with_emu(
                    getattr(self, 'pfJetET_BE_Emu_res'))
                getattr(self, 'pfJetET_HF_res').overlay_with_emu(
                    getattr(self, 'pfJetET_HF_Emu_res'))

        if self._doGen:
            getattr(self, 'genHT_eff').draw()
            getattr(self, 'genMETBE_eff').draw()
            getattr(self, 'genMETHF_eff').draw()
            getattr(self, 'genJetET_BE_eff').draw()
            getattr(self, 'genJetET_HF_eff').draw()

            getattr(self, 'genHT_eff_HR').draw()
            getattr(self, 'genMETBE_eff_HR').draw()
            getattr(self, 'genMETHF_eff_HR').draw()
            getattr(self, 'genJetET_BE_eff_HR').draw()
            getattr(self, 'genJetET_HF_eff_HR').draw()

            getattr(self, 'genHT_res').draw()
            getattr(self, 'genMETBE_res').draw()
            getattr(self, 'genMETHF_res').draw()
            getattr(self, 'genJetET_BE_res').draw()
            getattr(self, 'genJetET_HF_res').draw()

            if self._doEmu:
                getattr(self, 'genHT_eff').overlay_with_emu(
                    getattr(self, 'genHT_Emu_eff'))
                getattr(self, 'genMETBE_eff').overlay_with_emu(
                    getattr(self, 'genMETBE_Emu_eff'))
                getattr(self, 'genMETHF_eff').overlay_with_emu(
                    getattr(self, 'genMETHF_Emu_eff'))
                getattr(self, 'genJetET_BE_eff').overlay_with_emu(
                    getattr(self, 'genJetET_BE_Emu_eff'))
                getattr(self, 'genJetET_HF_eff').overlay_with_emu(
                    getattr(self, 'genJetET_HF_Emu_eff'))

                getattr(self, 'genHT_res').overlay_with_emu(
                    getattr(self, 'genHT_Emu_res'))
                getattr(self, 'genMETBE_res').overlay_with_emu(
                    getattr(self, 'genMETBE_Emu_res'))
                getattr(self, 'genMETHF_res').overlay_with_emu(
                    getattr(self, 'genMETHF_Emu_res'))
                getattr(self, 'genJetET_BE_Emu_res').overlay_with_emu(
                    getattr(self, 'genJetET_BE_Emu_res'))
                getattr(self, 'genJetET_HF_Emu_res').overlay_with_emu(
                    getattr(self, 'genJetET_HF_Emu_res'))

        return True
