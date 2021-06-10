import numpy as np
import re

def exclusiveTrigger(j, ev, trgAcc, trgNegate = []):
    if not hasattr(ev, 'trgMu_'+trgAcc):
        return False
    if getattr(ev, 'trgMu_'+trgAcc)[j] == 0:
        return False
    for t in trgNegate:
        if hasattr(ev, t):
            if getattr(ev, 'trgMu_'+t)[j] == 1:
                return False
    return True

def trigger_selection(iMu, ev, cat, muPt, muEta):
    if ev.trgMu_L1_dR[iMu] > 0.5:
        return False
    ptThr = float(re.search('Mu[0-9]+_', cat.trg).group(0)[2:-1])
    if np.abs(ev.trgMu_L1_pt[iMu]) < ptThr:
        return False
    if np.abs(ev.trgMu_L1_eta[iMu]) > 1.5:
        return False

    if not exclusiveTrigger(iMu, ev, 'HLT_' + cat.trg):
        return False

    if muPt < cat.min_pt or muPt > cat.max_pt:
        return False
    if not ev.trgMu_sigdxy_BS[iMu] > cat.minIP:
        return False
    if not abs(muEta) < 1.5:
        return False
    return True

def category_selection(j, ev, evEx, cat, saveTrgMu=False):
    idxTrigger = []
    if cat == 'probe':
        raise

    passed = [False, False]
    # print '-'*30
    if ev.mup_isTrg[j] >= 0:
        # print 'mup'
        # print ev.mup_isTrg[j]
        passed[0] = trigger_selection(int(ev.mup_isTrg[j]), ev, cat, evEx.mup_pt, evEx.mup_eta)
    if ev.mum_isTrg[j] >= 0:
        # print 'mum'
        # print ev.mum_isTrg[j]
        passed[1] = trigger_selection(int(ev.mum_isTrg[j]), ev, cat, evEx.mum_pt, evEx.mum_eta)
    if passed[0] and passed[1]: passed[np.random.randint(2)] = False
    if saveTrgMu:
        if passed[0]:
            evEx.trgMu_pt = evEx.mup_pt
            evEx.trgMu_eta = evEx.mup_eta
            evEx.trgMu_sigdxy = ev.trgMu_sigdxy_BS[int(ev.mup_isTrg[j])]
            evEx.otherMu_pt = evEx.mum_pt
        elif passed[1]:
            evEx.trgMu_pt = evEx.mum_pt
            evEx.trgMu_eta = evEx.mum_eta
            evEx.trgMu_sigdxy = ev.trgMu_sigdxy_BS[int(ev.mum_isTrg[j])]
            evEx.otherMu_pt = evEx.mup_pt
        else:
            evEx.trgMu_pt = -1
            evEx.trgMu_eta = -999999999
            evEx.trgMu_sigdxy = -1
            evEx.otherMu_pt = -1

    return np.sum(passed) > 0

def candidate_selection(j, ev, evEx, skipCut=None):
    if not abs(evEx.mup_eta) < 2.2:
        return False
    if not evEx.mup_pt > 3.5:
        return False
    if not ev.mup_dxy_PV[j] < 3:
        return False

    if not abs(evEx.mum_eta) < 2.2:
        return False
    if not evEx.mum_pt > 3.5:
        return False
    if not ev.mum_dxy_PV[j] < 3:
        return False

    if not ev.pval_mumu[j] > 0.1:
        return False
    if not abs(evEx.mass_mumu - 3.096916) < 0.08:
        return False
    if not evEx.Jpsi_pt > 4.5:
        return False
    if not ev.cosT_Jpsi_PV[j] > 0.95:
        return False

    if not evEx.K_pt > 0.8:
        return False
    if not abs(evEx.K_eta) < 2.4:
        return False
    if not abs(evEx.pi_eta) < 2.4:
        return False
    if not evEx.pi_pt > 0.8:
        return False

    if not ev.pval_piK[j] > 0.1:
        return False
    if not (7 in skipCut):
        if not abs(evEx.mass_candKst - 0.8955) <  0.07:
            return False
    if not evEx.mass_KK > 1.035:
        return False
    if not ev.sigdxy_vtxKst_PV[j] > 5:
        return False

    if not ev.pval_mumupiK[j] > 0.1:
        return False
    if not (7 in skipCut):
        if not abs(evEx.mass_mumupiK - 5.27963) < 0.275:
            return False

    return True

    # aux = abs(ev.mass_piK[j] - 0.8955) <  0.07
    # aux &= abs(ev.mass_piK[j] - 0.8955) < abs(ev.mass_piK_CPconj[j] - 0.8955)
    # aux &= ev.mass_KK[j] > 1.035
    # aux &= abs(ev.mass_mumupiK[j] - 5.27963) < 0.275
    # return aux

candidateSelection_stringList = [
    'abs(mass_mumu - 3.0969) < 0.08',
    'abs(mass_piK - 0.8955) <  0.07',
    'mum_pt > 3.5',
    'mup_pt > 3.5',
    'Jpsi_pt > 4.5',
    'pval_mumu > 0.1',
    'abs(mum_eta) < 2.2',
    'abs(mup_eta) < 2.2',
    'cosT_Jpsi_PV > 0.95',
    'mum_dxy < 3',
    'mup_dxy < 3',
    'pval_piK > 0.1',
    'fabs(mass_piK - 0.895) < fabs(mass_piK_CPconj - 0.895)',
    'mass_KK > 1.035',
    'K_sigdxy_PV > 2',
    'pi_sigdxy_PV > 2',
    'sigdxy_vtxKst_PV > 5',
    'K_pt > 0.8',
    'pval_mumupiK > 0.1',
    'pi_pt > 0.8',
    'abs(mass_mumupiK - 5.27963) < 0.275',
]

candidateSelection_nameList = [
    'mass_mumu',
    'mass_piK',
    'mum_pt',
    'mup_pt',
    'Jpsi_pt',
    'pval_mumu',
    '|mum_eta|',
    '|mup_eta|',
    'cosT_Jpsi_PV',
    'mum_dxy',
    'mup_dxy',
    'pval_piK',
    'piK VS CPconj',
    'mass_KK',
    'K_IP',
    'pi_IP',
    'sigdxy_vtxKst',
    'K_pt',
    'pval_mumupiK',
    'pi_pt',
    'mass_mumupiK',
]
