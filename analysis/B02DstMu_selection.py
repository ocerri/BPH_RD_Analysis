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

def candidate_selection(j, ev, skipCut=[]):
    aux = True
    if not (1 in skipCut):
        aux &= ev.pval_piK[j] > 0.1
    if not (2 in skipCut):
        aux &= ev.K_pt[j] > 0.8
    if not (3 in skipCut):
        aux &= abs(ev.K_eta[j]) < 2.4
    if not (4 in skipCut):
        aux &= ev.sigdxy_K_PV[j] > 2
    if not (5 in skipCut):
        aux &= ev.pi_pt[j] > 0.8
    if not (6 in skipCut):
        aux &= abs(ev.pi_eta[j]) < 2.4
    if not (7 in skipCut):
        aux &= ev.sigdxy_pi_PV[j] > 2
    if not (8 in skipCut):
        aux &= abs(ev.mass_piK[j] - 1.864) < 0.05
        aux &= ev.mass_piK_hKK[j] > 1.91 and ev.mass_piK_hpipi[j] < 1.83
    if not (9 in skipCut):
        aux &= ev.sigdxy_vtxD0_PV[j] > 2
    if not aux:
        return False

    if not (10 in skipCut):
        aux = ev.pis_pt[j] > 0.4
    if not (11 in skipCut):
        aux &= abs(ev.pis_eta[j]) < 2.4
    if not (12 in skipCut):
        aux &= ev.sigdxy_pis_PV[j] > 2
    if not (13 in skipCut):
        aux &= ev.pval_D0pis[j] > 0.1

    if not (14 in skipCut):
        aux &= abs(ev.mass_D0pis[j] - 2.01026) < 0.03
    if not (15 in skipCut):
        aux &= 1e3*abs(ev.mass_D0pis[j] - ev.mass_piK[j] - 0.14543) < 2.

    if not (16 in skipCut):
        aux &= ev.pval_D0pismu[j] > 0.1
    if not (17 in skipCut):
        aux &= ev.cos_D0pismu_PV[j] > 0.99
    if not (18 in skipCut):
        aux &= ev.q2_D0pismu[j] > -2.
    if not (19 in skipCut):
        aux &= ev.q2_D0pismu[j] < 12
    if not (20 in skipCut):
        aux &= ev.mass_D0pismu[j] < 7.
    if not aux:
        return False

    if not (21 in skipCut) and not ev.nTksAdd[j] == 0:
        idx_st = 0
        for jjj in range(j):
            idx_st += int(ev.nTksAdd[jjj])

        idx_stop = int(idx_st + ev.nTksAdd[j])
        for jj in range(idx_st, idx_stop):
            if ev.tksAdd_massVis[jj] < 5.28 and ev.tksAdd_cos_PV[jj]>0.95:
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
