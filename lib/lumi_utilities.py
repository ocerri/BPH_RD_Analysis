import os

def getLumiReport(datasets_loc):
    prescale = {'A':1./6., 'B':1./6., 'C':1./5., 'D':1./5.}
    lumi_max = 40.7

    lumi_dic = {'A':{}, 'B':{}, 'C':{}, 'D':{}}
    lumi_tot = 0
    for n in datasets_loc:
        aux = n.split('/')[-2:]
        lumi_rep = '../data/cmsRD/lumiReport/'
        lumi_rep += aux[0] + '_' + aux[1].replace('CAND.root', 'bricalc.csv')
        if os.path.isfile(lumi_rep):
            lines = open(lumi_rep, 'r').readlines()
            lumi_raw = float(lines[-3][:-1].split(',')[-1])

            aux = lumi_rep.split('/')[-1]
            ds = aux.split('_')[0]
            era = aux.split('_')[1][:8]
            lumi = lumi_raw * prescale[era[-1]]
            lumi_dic[era[-1]][ds] = lumi
            lumi_tot += lumi
    #         print ds, era, '{:.1f} fb^-1'.format(lumi)
        else:
            print 'Lumi brilcalc not found for', n

    print 'Lumi tot: {:.1f} fb^-1'.format(lumi_tot)
    lumi_analyzed = 100*lumi_tot/lumi_max
    print 'Lumi analyzed: {:.1f}%'.format(lumi_analyzed)
    return lumi_tot, lumi_dic
