import numpy as np
import ROOT as rt

class pileupReweighter(object):
    def __init__(self, mcFile, mcHistoName, trg=None, dataPileupFile=None, dataHistoName='pileup'):
        if dataPileupFile is None:
            if trg is None:
                print 'Must specify a trigger to load default data pileup'
                raise
            dataPileupFile = '/storage/af/group/rdst_analysis/BPhysics/data/PileupHistograms/' + trg + '_allParts.root'

        hPileupTarget = None

        fData = rt.TFile.Open(dataPileupFile, 'READ')
        hPileupTarget = fData.Get(dataHistoName).Clone('hPileupData')
        hPileupTarget_norm = float(hPileupTarget.Integral())

        fMC = rt.TFile.Open(mcFile, 'READ')
        hPileupGen = fMC.Get(mcHistoName).Clone('hPileupMC')
        hPileupGen_norm = float(hPileupGen.Integral())

        weights = np.ones(hPileupGen.GetNbinsX())
        s = 0
        for i in range(weights.shape[0]):
            ih = i + 1
            if hPileupGen.GetBinContent(ih) == 0:
                continue
            x = hPileupGen.GetBinCenter(ih)
            yT = hPileupTarget.GetBinContent(hPileupTarget.FindBin(x)) / hPileupTarget_norm
            weights[i] = yT / (hPileupGen.GetBinContent(ih)/hPileupGen_norm)
            s += hPileupGen.GetBinContent(ih) * weights[i]

        self.weightsPileupMC = hPileupGen_norm*weights/s
        # print ['{:.3f}'.format(x) for x in self.weightsPileupMC]
        print '[pileup reweighting] wMax={:.1f}, wMin={:.1f}'.format(np.max(self.weightsPileupMC), np.min(self.weightsPileupMC))

        fData.Close()
        fMC.Close()

    def getPileupWeights(self, arrNvtx, selection=None, clip_range=(0,10)):
        x = np.rint(arrNvtx).astype(np.int)
        if not selection is None:
            x = x[selection]
        return np.clip(self.weightsPileupMC[x.astype(np.int)], clip_range[0], clip_range[1])

# class pileupReweighter(object):
#     def __init__(self, mcSkimFile, cat, histoName='hAllNvtx', dataDate='200515'):
#         loc = '../data/cmsRD/ParkingBPH{}/'+'Run2018D-05May2019promptD-v1_RDntuplizer_PrescaleVertices_{}_CAND.root'.format(dataDate)
#         fAuxPileupRD = []
#
#         hPileupTarget = None
#
#         for i in range(1, 6):
#             fAuxPileupRD.append(rt.TFile.Open(loc.format(i), 'READ'))
#             if hPileupTarget is None:
#                 hPileupTarget = fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg).Clone()
#             else:
#                 hPileupTarget.Add(fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg))
#
#         hPileupTarget.Scale(1./hPileupTarget.Integral())
#
#         fAuxPileupMC = rt.TFile.Open(mcSkimFile, 'READ')
#         hPileupGen = fAuxPileupMC.Get(histoName)
#         hPileupGen.Scale(1./hPileupGen.Integral())
#
#         weights = np.ones(hPileupGen.GetNbinsX())
#         s = 0
#         for i in range(weights.shape[0]):
#             if hPileupGen.GetBinContent(i+1) == 0:
#                 continue
#             weights[i] = hPileupTarget.GetBinContent(i+1)/hPileupGen.GetBinContent(i+1)
#             s += hPileupGen.GetBinContent(i+1) * weights[i]
#
#         self.weightsPileupMC = weights/s
#
#         for f in fAuxPileupRD + [fAuxPileupMC]:
#             f.Close()
#
#     def getPileupWeights(self, arrNvtx, selection=None):
#         x = arrNvtx
#         if not selection is None:
#             x = x[selection]
#         return self.weightsPileupMC[x.astype(np.int)]
