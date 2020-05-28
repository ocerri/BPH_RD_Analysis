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

def trigger_selection(j, ev, cat):
    aux = exclusiveTrigger(j, ev, 'HLT_' + cat.trg)
    aux &= ev.trgMu_pt[j] > cat.min_pt
    aux &= ev.trgMu_pt[j] < cat.max_pt
    aux &= ev.trgMu_sigdxy[j] > cat.minIP
    aux &= abs(ev.trgMu_eta[j]) < 1.5
    return aux

def candidate_selection(j, ev, skipCut=None):
    aux = ev.mum_pt[j] > 3.5
    aux &= abs(ev.mum_eta[j]) < 2.2
    aux &= ev.mum_dxy[j] < 3
    aux &= ev.mup_pt[j] > 3.5
    aux &= abs(ev.mup_eta[j]) < 2.2
    aux &= ev.mup_dxy[j] < 3
    aux &= ev.pval_mumu[j] > 0.1
    aux &= abs(ev.mass_mumu[j] - 3.0969) < 0.08
    aux &= ev.Jpsi_pt[j] > 4.5
    aux &= ev.cosT_Jpsi_PV[j] > 0.95
    aux &= ev.K_pt[j] > 0.8
    aux &= ev.K_sigdxy_PV[j] > 2
    aux &= ev.pi_pt[j] > 0.8
    aux &= ev.pi_sigdxy_PV[j] > 2
    aux &= ev.pval_piK[j] > 0.1
    aux &= abs(ev.mass_piK[j] - 0.8955) <  0.07
    aux &= abs(ev.mass_piK[j] - 0.8955) < abs(ev.mass_piK_CPconj[j] - 0.8955)
    aux &= ev.mass_KK[j] > 1.035
    aux &= ev.sigdxy_vtxKst_PV[j] > 5
    aux &= ev.pval_mumupiK[j] > 0.1
    aux &= abs(ev.mass_mumupiK[j] - 5.27963) < 0.275
    return aux

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
