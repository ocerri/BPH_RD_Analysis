import numpy as np
import ROOT as rt

class pileupReweighter(object):
    def __init__(self, mcSkimFile, cat, dataDate='200515'):
        loc = '../data/cmsRD/ParkingBPH{}/'+'Run2018D-05May2019promptD-v1_RDntuplizer_PrescaleVertices_{}_CAND.root'.format(dataDate)
        fAuxPileupRD = []

        hPileupTarget = None

        for i in range(1, 6):
            fAuxPileupRD.append(rt.TFile.Open(loc.format(i), 'READ'))
            if hPileupTarget is None:
                hPileupTarget = fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg).Clone()
            else:
                hPileupTarget.Add(fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg))

        hPileupTarget.Scale(1./hPileupTarget.Integral())

        fAuxPileupMC = rt.TFile.Open(mcSkimFile, 'READ')
        hPileupGen = fAuxPileupMC.Get('hAllNvtx')
        hPileupGen.Scale(1./hPileupGen.Integral())

        weights = np.ones(hPileupGen.GetNbinsX())
        s = 0
        for i in range(weights.shape[0]):
            if hPileupGen.GetBinContent(i+1) == 0:
                continue
            weights[i] = hPileupTarget.GetBinContent(i+1)/hPileupGen.GetBinContent(i+1)
            s += hPileupGen.GetBinContent(i+1) * weights[i]

        self.weightsPileupMC = weights/s

        for f in fAuxPileupRD + [fAuxPileupMC]:
            f.Close()

    def getPileupWeights(arrNvtx, selection=None):
        x = arrNvtx
        if not selection is None:
            x = x[selection]
        return self.weightsPileupMC[x.astype(np.int)]
