import sys
sys.path.append('../lib')
import numpy as np
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

import tdrstyle
tdrstyle.setTDRStyle()

col_dic = {'mu': rt.kAzure+1, 'tau': rt.kRed-4, 'Hc':rt.kGreen+1, 'Dstst': rt.kViolet-7}

def plot_gridVarQ2(CMS_lumi, binning, histo, scale_dic, min_y=1e-4, draw_pulls=False, pulls_ylim=[0.8, 1.2], logy=False, iPad_legend=1):
    col_dic = {'mu': rt.kAzure+1, 'tau': rt.kRed-4, 'Hc':rt.kGreen+1, 'Dstst': rt.kViolet-7}

    canvas = rt.TCanvas('c_out', 'c_out', 50, 50, 2*600, 400*len(binning['q2'])-1)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.Divide(2, len(binning['q2'])-1, 0.001, 0.001)

    canvas.dnd = []

    vars_to_plot = ['M2_miss', 'Est_mu']

    xAx_title = {'M2_miss':'m^{2}_{miss} [GeV^{2}]', 'Est_mu':'E_{#mu}* [GeV]'}
    label_dic = {'data' : 'Data',
                 'mu'   : 'B#rightarrow D*#mu#nu',
                 'tau'  : 'B#rightarrow D*#tau#nu',
                 'Hc'   : 'B#rightarrow D*H_{c}',
                 'Dstst': 'B#rightarrow D**#mu#nu'
                }

    rt.TGaxis.SetMaxDigits(3)

    max_entries = dict(zip(vars_to_plot, [0]*len(vars_to_plot)))
    for k, h_dic in histo.iteritems():
        if 'M2' in k:
            max_entries['M2_miss'] = max(h_dic['data'].GetMaximum(), max_entries['M2_miss'])
        elif 'Est' in k:
            max_entries['Est_mu'] = max(h_dic['data'].GetMaximum(), max_entries['Est_mu'])

    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        q2_txt = '{:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h)

        for i_v, vark in enumerate(vars_to_plot):
            cat_name = vark+'_q2bin'+str(i_q2)
            h_dic = histo[cat_name]

            i_pad = i_q2*2 + i_v + 1
            pad_master = canvas.cd(i_pad)

            if not draw_pulls:
                pad = rt.TPad('pmain_'+cat_name, 'pmain_'+cat_name, 0, 0, 1, 1)
                pad.SetBottomMargin(0.2)
            else:
                pad = rt.TPad('pmain_'+cat_name, 'pmain_'+cat_name, 0, 0.25, 1, 1)
                pad.SetBottomMargin(0.)
            pad.SetTopMargin(0.07)
            pad.SetRightMargin(0.05)
            pad.SetLeftMargin(0.13)
            pad.Draw()
            pad.cd()

            h = h_dic['data'].Clone('h_aux_data_'+cat_name)
            if 'data' in scale_dic.keys(): h.Scale(scale_dic['data'])
            h.SetLineColor(1)
            h.SetLineWidth(1)
            h.SetMarkerColor(1)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xAx_title[vark])
            h.GetXaxis().SetTitleSize(0.07)
            h.GetXaxis().SetLabelSize(0.07)
            h.GetYaxis().SetTitleOffset(1.16)
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleSize(0.06)
            h.GetYaxis().SetLabelSize(0.07)
            iunits = xAx_title[vark].find('[') + 1
            h.GetYaxis().SetTitle('Candidates / {:.2f} '.format(h.GetBinWidth(3)) + xAx_title[vark][iunits:-1])
            max_y = max_entries[vark]*1.3
            if 'data' in scale_dic.keys():
                max_y *= scale_dic['data']
            h.GetYaxis().SetRangeUser(min_y, max_y)
            h_list = [h]

            if 'Dstst' in h_dic.keys():
                h_Dstst = h_dic['Dstst'].Clone('h_aux_Dstst_'+cat_name)
                if 'Dstst' in scale_dic.keys(): h_Dstst.Scale(scale_dic['Dstst'])
                h_Dstst.SetLineWidth(0)
                h_Dstst.SetFillColor(col_dic['Dstst'])
                h_Dstst.SetFillStyle(1)
                h_Dstst.Sumw2(0)
                h_list.append(h_Dstst)

            if 'Hc' in h_dic.keys():
                h_Hc = h_dic['Hc'].Clone('h_aux_Hc_'+cat_name)
                if 'Hc' in scale_dic.keys(): h_Hc.Scale(scale_dic['Hc'])
                if 'Dstst' in h_dic.keys():
                    h_Hc.Add(h_Dstst)
                h_Hc.SetLineWidth(0)
                h_Hc.SetFillColor(col_dic['Hc'])
                h_Hc.SetFillStyle(1)
                h_Hc.Sumw2(0)
                h_list.append(h_Hc)

            h_tau = h_dic['tau'].Clone('h_aux_tau_'+cat_name)
            if 'tau' in scale_dic.keys(): h_tau.Scale(scale_dic['tau'])
            if 'Hc' in h_dic.keys():
                h_tau.Add(h_Hc)
            elif 'Dstst' in h_dic.keys():
                h_tau.Add(h_Dstst)
            h_tau.SetLineWidth(0)
            h_tau.SetFillColor(col_dic['tau'])
            h_tau.SetFillStyle(1)
            h_tau.Sumw2(0)
            h_list.append(h_tau)

            h_mu = h_dic['mu'].Clone('h_aux_mu_'+cat_name)
            if 'mu' in scale_dic.keys(): h_mu.Scale(scale_dic['mu'])
            h_mu.Add(h_tau)
            h_mu.SetLineWidth(0)
            h_mu.SetFillColor(col_dic['mu'])
            h_mu.SetFillStyle(1)
            h_mu.Sumw2(0)
            h_list.append(h_mu)


            h.Draw('E1')
            h_mu.DrawCopy('SAME')
            h_tau.Draw('SAME')
            if 'Hc' in h_dic.keys():
                h_Hc.Draw('SAME')
            if 'Dstst' in h_dic.keys():
                h_Dstst.Draw('SAME')
            h.Draw('SAMEE1') #Draw it a second time to bring it in foreground

            l = rt.TLatex()
            l.SetTextAlign(11)
            l.SetTextSize(0.06)
            l.SetTextFont(42)
            l.DrawLatexNDC(0.18, 0.85, q2_txt)

            CMS_lumi.CMS_lumi(pad, -1, 33, cmsTextSize=0.75*1.2, lumiTextSize=0.6*1.2)
            if logy:
                pad.SetLogy()

            if i_pad == iPad_legend:
                leg = rt.TLegend(0.65, 0.4, 0.9, 0.7)
                leg.SetTextFont(42)
                leg.SetTextAlign(12)
                leg.SetLineWidth(0)
                leg.SetBorderSize(0)
                leg.AddEntry(h, label_dic['data'], 'lep')
                leg.AddEntry(h_mu, label_dic['mu'], 'f')
                leg.AddEntry(h_tau, label_dic['tau'], 'f')
                if 'Hc' in h_dic.keys():
                    leg.AddEntry(h_Hc, label_dic['Hc'], 'f')
                if 'Dstst' in h_dic.keys():
                    leg.AddEntry(h_Dstst, label_dic['Dstst'], 'f')
                leg.Draw()
                canvas.dnd.append(leg)

            canvas.dnd.append([pad, h_list])

            if draw_pulls:
                pad_master.cd()

                pad = rt.TPad('ppull_'+cat_name, 'ppull_'+cat_name, 0, 0, 1, 0.25)
                pad.SetBottomMargin(0.5)
                pad.SetTopMargin(0)
                pad.SetRightMargin(0.05)
                pad.SetLeftMargin(0.13)
                pad.Draw()
                pad.cd()

                h_dr = h.Clone('h_aux_dataratio_'+cat_name)
                h_dr.GetYaxis().SetTitle('RD/MC')
                h_tot = h_dic['total'].Clone('h_aux_total_'+cat_name)
                g_up = rt.TGraph()
                g_up.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), 1)
                g_down = rt.TGraph()
                g_down.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), 1)
                for i in range(1, h_dr.GetNbinsX()+1):
                    c = h_dr.GetBinContent(i)
                    e = h_dr.GetBinError(i)
                    c_MC = h_tot.GetBinContent(i)
                    e_MC = h_tot.GetBinError(i)
                    h_dr.SetBinContent(i, c/c_MC)
                    h_dr.SetBinError(i, e/c_MC)
                    x_low = h_tot.GetBinCenter(i) - 0.5*h.GetBinWidth(i)
                    x_up = h_tot.GetBinCenter(i) + 0.5*h.GetBinWidth(i)
                    g_up.SetPoint(2*i-1, x_low, (c_MC+e_MC)/c_MC)
                    g_up.SetPoint(2*i, x_up, (c_MC+e_MC)/c_MC)
                    g_down.SetPoint(2*i-1, x_low, (c_MC-e_MC)/c_MC)
                    g_down.SetPoint(2*i, x_up, (c_MC-e_MC)/c_MC)
                g_up.SetPoint(2*i+1, x_up, 1)
                g_down.SetPoint(2*i+1, x_up, 1)

                h_dr.GetYaxis().SetRangeUser(pulls_ylim[0], pulls_ylim[1])
                h_dr.GetYaxis().SetTitleOffset(0.35)
                h_dr.GetYaxis().SetTitleSize(0.2)
                h_dr.GetYaxis().SetLabelSize(0.2)
                h_dr.GetYaxis().SetNdivisions(402)
                h_dr.GetXaxis().SetTitleOffset(0.95)
                h_dr.GetXaxis().SetTitleSize(0.22)
                h_dr.GetXaxis().SetLabelSize(0.22)
                h_dr.GetXaxis().SetTickSize(0.07)

                h_dr.Draw('E1')

                g_up.SetFillColor(rt.kGray)
                g_up.SetFillStyle(1)
                g_up.Draw('F')
                g_down.SetFillColor(rt.kGray)
                g_down.SetFillStyle(1)
                g_down.Draw('F')
                l = rt.TLine()
                l.SetLineColor(rt.kGray+1)
                l.SetLineWidth(1)
                l.SetLineStyle(9)
                x_low = h_tot.GetBinCenter(1)-0.5*h.GetBinWidth(1)
                x_high = h_tot.GetBinCenter(i)+0.5*h.GetBinWidth(i)
                l.DrawLine(x_low, 1, x_high, 1)
                h_dr.Draw('sameE') #redraw it to bring it to front

                canvas.dnd.append([pad, h_dr, h_tot, g_up, g_down])

    canvas.Draw()
    return canvas
