import sys
sys.path.append('../lib')
import numpy as np
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

import tdrstyle
tdrstyle.setTDRStyle()

def plot_gridVarQ2(CMS_lumi, binning, histo, scale_dic):
    canvas = rt.TCanvas('c_out', 'c_out', 50, 50, 2*600, 400*binning['q2'][0])
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.Divide(2, binning['q2'][0], 0.005, 0.005)
    
    canvas.dnd = []

    vars_to_plot = ['M2_miss', 'Est_mu']

    xAx_title = {'M2_miss':'m^{2}_{miss} [GeV^{2}]', 'Est_mu':'E_{#mu}* [GeV]'}
    label_dic = {'data' : 'Data',
                 'mu' : 'B#rightarrow D*#mu#nu',
                 'tau' : 'B#rightarrow D*#tau#nu'
                }
    
    rt.TGaxis.SetMaxDigits(3)

    max_entries = dict(zip(vars_to_plot, [0]*len(vars_to_plot)))
    for k, h_dic in histo.iteritems():
        if 'M2' in k:
            max_entries['M2_miss'] = max(h_dic['data'].GetMaximum(), max_entries['M2_miss'])
        elif 'Est' in k:
            max_entries['Est_mu'] = max(h_dic['data'].GetMaximum(), max_entries['Est_mu'])

    for i_q2 in range(binning['q2'][0]):
        w_q2 = (binning['q2'][2] - binning['q2'][1])/binning['q2'][0]
        q2_l = binning['q2'][1] + w_q2 * i_q2
        q2_h = binning['q2'][1] + w_q2 * (i_q2+1)
        q2_txt = '{:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h)

        for i_v, vark in enumerate(vars_to_plot):
            cat_name = vark+'_q2-'+str(i_q2)
            h_dic = histo[cat_name]

            h = h_dic['data'].Clone('h_aux_data_'+cat_name)
            h.SetLineColor(1)
            h.SetLineWidth(1)
            h.SetMarkerColor(1)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xAx_title[vark])
            h.GetXaxis().SetTitleSize(0.07)
            h.GetXaxis().SetLabelSize(0.07)
            h.GetYaxis().SetTitleOffset(1.13)
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleSize(0.06)
            h.GetYaxis().SetLabelSize(0.07)
            iunits = xAx_title[vark].find('[') + 1
            h.GetYaxis().SetTitle('Candidates / {:.2f} '.format(h.GetBinWidth(1)) + xAx_title[vark][iunits:-1])
            h.GetYaxis().SetRangeUser(0, max_entries[vark]*1.2)

            h_tau = h_dic['tau'].Clone('h_aux_tau_'+cat_name)
            h_tau.Scale(scale_dic['tau'])
            h_tau.SetLineWidth(0)
            h_tau.SetFillColor(rt.kRed-4)
            h_tau.SetFillStyle(1)
            h_tau.Sumw2(0)

            h_mu = h_dic['mu'].Clone('h_aux_mu_'+cat_name)
            h_mu.Scale(scale_dic['mu'])
            h_mu.Add(h_tau)
            h_mu.SetLineWidth(0)
            h_mu.SetFillColor(rt.kAzure+1)
            h_mu.SetFillStyle(1) 
            h_mu.Sumw2(0)

            i_pad = i_q2*2 + i_v + 1
            pad = canvas.cd(i_pad)
            pad.SetBottomMargin(0.2)
            pad.SetTopMargin(0.07)
            pad.SetRightMargin(0.05)
            pad.SetLeftMargin(0.14)

            h.Draw('E1')
            h_mu.Draw('SAME')
            h_tau.Draw('SAME')
            h.Draw('SAMEE1')

            l = rt.TLatex()
            l.SetTextAlign(11)
            l.SetTextSize(0.06)
            l.SetTextFont(42)
            l.DrawLatexNDC(0.18, 0.85, q2_txt)

            CMS_lumi.CMS_lumi(pad, -1, 33, 0.75*1.5, 0.6*1.5)

            if i_pad == 1:
                leg = rt.TLegend(0.65, 0.4, 0.9, 0.7)
                leg.SetTextFont(42)
                leg.SetTextAlign(12)
                leg.SetLineWidth(0)
                leg.SetBorderSize(0)
                leg.AddEntry(h, label_dic['data'], 'lep')
                leg.AddEntry(h_mu, label_dic['mu'], 'f')
                leg.AddEntry(h_tau, label_dic['tau'], 'f')
                leg.Draw()
                canvas.dnd.append(leg)

            canvas.dnd.append([h, h_tau, h_mu])

    canvas.Draw()
    return canvas