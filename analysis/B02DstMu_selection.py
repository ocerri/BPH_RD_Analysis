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
    aux = True
    if not skipCut == 1:
        aux &= ev.pval_piK[j] > 0.1
    if not skipCut == 2:
        aux &= ev.K_pt[j] > 0.8
    if not skipCut == 3:
        aux &= abs(ev.K_eta[j]) < 2.4
    if not skipCut == 4:
        aux &= ev.sigdxy_K_PV[j] > 2
    if not skipCut == 5:
        aux &= ev.pi_pt[j] > 0.8
    if not skipCut == 6:
        aux &= abs(ev.pi_eta[j]) < 2.4
    if not skipCut == 7:
        aux &= ev.sigdxy_pi_PV[j] > 2
    if not skipCut == 8:
        aux &= abs(ev.mass_piK[j] - 1.864) < 0.05
    if not skipCut == 9:
        aux &= ev.sigdxy_vtxD0_PV[j] > 2
    if not aux:
        return False

    if not skipCut == 10:
        aux = ev.pis_pt[j] > 0.4
    if not skipCut == 11:
        aux &= abs(ev.pis_eta[j]) < 2.4
    if not skipCut == 12:
        aux &= ev.sigdxy_pis_PV[j] > 2
    if not skipCut == 13:
        aux &= ev.pval_D0pis[j] > 0.1

    if not skipCut == 14:
        aux &= abs(ev.mass_D0pis[j] - 2.01) < 0.03
    if not skipCut == 15:
        aux &= 1e3*abs(ev.mass_D0pis[j] - ev.mass_piK[j] - 0.14543) < 2.5

    if not skipCut == 16:
        aux &= ev.pval_D0pismu[j] > 0.1
    if not skipCut == 17:
        aux &= ev.cos_D0pismu_PV[j] > 0.99
    if not skipCut == 18:
        aux &= ev.q2_D0pismu[j] > -2.
    if not skipCut == 19:
        aux &= ev.q2_D0pismu[j] < 12
    if not skipCut == 20:
        aux &= ev.mass_D0pismu[j] < 7.
    if not aux:
        return False

    if not skipCut == 21 and not ev.nTksAdd[j] == 0:
        idx_st = 0
        for jjj in range(j):
            idx_st += int(ev.nTksAdd[jjj])

        idx_stop = int(idx_st + ev.nTksAdd[j])
        for jj in range(idx_st, idx_stop):
            if ev.tksAdd_massVis[jj] < 5.28:
                return False

    return True

candidateSelection_stringList = [
    'pval_piK > 0.1',
    'K_pt > 0.8',
    'abs(K_eta) < 2.4',
    'K_IP > 2',
    'pi_pt > 0.8',
    'abs(pi_eta) < 2.4',
    'pi_IP > 2',
    'abs(mass_piK - 1.864) < 0.05',
    'sigdxy_vtxD0_PV > 2',
    'pis_pt > 0.4',
    'abs(pis_eta) < 2.4',
    'pis_IP > 2',
    'pval_D0pis > 0.1',
    'abs(mass_D0pis - 2.01) < 0.03',
    '1e3*abs(mass_D0pis - mass_piK - 0.14543) < 2.5',
    'pval_D0pismu > 0.1',
    'cos_D0pismu_PV > 0.99',
    'q2 > -2.',
    'q2 < 12',
    'mass_D0pismu < 7.',
    'N_lowMassAddTks == 0',
]

candidateSelection_nameList = [
    'pval_piK',
    'K_pt',
    '|K_eta|',
    'K_IP',
    'pi_pt',
    '|pi_eta|',
    'pi_IP',
    'mass_piK',
    'sigdxy_vtxD0_PV',
    'pis_pt',
    '|pis_eta|',
    'pis_IP',
    'pval_D0pis',
    'mass_D0pis',
    '|m_D0pis - m_piK|',
    'pval_D0pismu',
    'cos_D0pismu_PV',
    'min q2',
    'max q2',
    'm D0pismu',
    'N tks',
]
