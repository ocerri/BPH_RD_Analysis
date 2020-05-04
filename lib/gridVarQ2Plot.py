import sys
sys.path.append('../lib')
import numpy as np
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

import tdrstyle
tdrstyle.setTDRStyle()

col_dic = {'mu': rt.kAzure+1, 'tau': rt.kRed-4, 'Hc':rt.kGreen+1, 'BpDstst': rt.kOrange-3, 'B0Dstst': rt.kViolet-7}

label_dic = {'data' : 'Data',
             'mu'   : 'B#rightarrow D*#mu#nu',
             'tau'  : 'B#rightarrow D*#tau#nu',
             'Hc'   : 'B#rightarrow D*H_{c}',
             'BpDstst': 'B^{+}#rightarrow D**#mu#nu',
             'B0Dstst': 'B_{0}#rightarrow D**#mu#nu'
            }

fillStyleVar = [1, 3345, 3354]
sampleDstst = {
'BpDstst': ['DstPip', 'DstPipPi0'],
'B0Dstst': ['DstPi0', 'DstPipPim', 'DstPi0Pi0']
}


def createLegend(h_list, h_dic, canvas, loc=[0.65, 0.4, 0.9, 0.7], cat_name=''):
    leg = rt.TLegend(loc[0], loc[1], loc[2], loc[3])
    leg.SetTextFont(42)
    leg.SetTextAlign(12)
    leg.SetLineWidth(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_list[0], label_dic['data'], 'lep')
    for n in ['mu', 'tau', 'Hc']:
        for h in h_list:
            if 'h_aux_'+n+cat_name == h.GetName():
                leg.AddEntry(h, label_dic[n], 'f')
                break
    for k, sampleList in sampleDstst.iteritems():
        present = np.sum([n in h_dic.keys() for n in sampleList])
        if not present: continue
        h = rt.TH1D('hAuxLeg_'+k, 'hAuxLeg_'+k, 1, 0, 1)
        h.SetLineWidth(0)
        h.SetFillColor(col_dic[k])
        h.SetFillStyle(1)
        canvas.dnd.append(h)
        leg.AddEntry(h, label_dic[k], 'f')
    return leg

def plot_gridVarQ2(CMS_lumi, binning, histo, scale_dic={}, min_y=1e-4, draw_pulls=False, pulls_ylim=[0.8, 1.2], logy=False, iPad_legend=1, max_y_shared=False, mergeDstst=True, pullsRatio=False):
    canvas = rt.TCanvas('c_out', 'c_out', 50, 50, 2*600, 400*len(binning['q2'])-1)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.Divide(2, len(binning['q2'])-1, 0.001, 0.001)

    canvas.dnd = []

    vars_to_plot = ['M2_miss', 'Est_mu']

    xAx_title = {'M2_miss':'m^{2}_{miss} [GeV^{2}]', 'Est_mu':'E_{#mu}* [GeV]'}

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
            max_y = 1.3*h.GetMaximum()
            if max_y_shared:
                max_y = max_entries[vark]*1.3
                if 'data' in scale_dic.keys():
                    max_y *= scale_dic['data']
            h.GetYaxis().SetRangeUser(min_y, max_y)
            h_list = [h]

            for k, sampleList in sampleDstst.iteritems():
                for iproc, procName in enumerate(sampleList):
                    if not procName in h_dic.keys(): continue
                    h = h_dic[procName].Clone('h_aux_'+procName+'_'+cat_name)
                    if procName in scale_dic.keys(): h.Scale(scale_dic[procName])
                    h.SetLineWidth(0)
                    h.SetFillColor(col_dic[k])
                    h.SetFillStyle(1)
                    if not mergeDstst: h.SetFillStyle(fillStyleVar[iproc])
                    h.Sumw2(0)
                    hh = h_list[-1]
                    if not hh.GetName() == 'h_aux_data_'+cat_name:
                        h.Add(hh)
                    h_list.append(h)


            procOrder = ['Hc', 'tau', 'mu']
            for iProc, procName in enumerate(procOrder):
                if not procName in h_dic.keys(): continue
                h = h_dic[procName].Clone('h_aux_'+procName+'_'+cat_name)
                if procName in scale_dic.keys(): h.Scale(scale_dic[procName])
                h.SetLineWidth(0)
                h.SetFillColor(col_dic[procName])
                h.SetFillStyle(1)
                h.Sumw2(0)
                hh = h_list[-1]
                if not hh.GetName() == 'h_aux_data_'+cat_name:
                    h.Add(hh)
                h_list.append(h)

            h_list[0].Draw('E1')
            for h in reversed(h_list[1:]):
                h.Draw('SAME')
            h_list[0].Draw('SAMEE1') #Draw it a second time to bring it in foreground

            l = rt.TLatex()
            l.SetTextAlign(11)
            l.SetTextSize(0.06)
            l.SetTextFont(42)
            l.DrawLatexNDC(0.18, 0.85, q2_txt)

            CMS_lumi.CMS_lumi(pad, -1, 33, cmsTextSize=0.75*1.2, lumiTextSize=0.6*1.2)
            if logy:
                pad.SetLogy()

            if i_pad == iPad_legend:
                loc=[0.6, 0.35, 0.92, 0.7]
                if draw_pulls:
                    loc=[0.6, 0.25, 0.92, 0.7]
                leg = createLegend(h_list, h_dic, canvas, loc=loc, cat_name='_'+cat_name)
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

                h_dr = h_list[0].Clone('h_aux_dataratio_'+cat_name)
                h_dr.GetYaxis().SetTitle('Pull')
                if pullsRatio:
                    h_dr.GetYaxis().SetTitle('RD/MC')
                h_tot = h_dic['total'].Clone('h_aux_total')
                g_up = rt.TGraph()
                yExp = 1 if pullsRatio else 0
                g_up.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), yExp)
                g_down = rt.TGraph()
                g_down.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), yExp)
                for i in range(1, h_dr.GetNbinsX()+1):
                    c = h_dr.GetBinContent(i)
                    e = h_dr.GetBinError(i)
                    c_MC = h_tot.GetBinContent(i)
                    e_MC = h_tot.GetBinError(i)
                    if pullsRatio:
                        h_dr.SetBinContent(i, c/c_MC)
                        h_dr.SetBinError(i, e/c_MC)
                    else:
                        if e > 0:
                            h_dr.SetBinContent(i, (c-c_MC)/e)
                        else:
                            h_dr.SetBinContent(i, (c-c_MC)/e_MC)
                        h_dr.SetBinError(i, 0)

                    x_low = h_tot.GetBinCenter(i) - 0.5*h_tot.GetBinWidth(i)
                    x_up = h_tot.GetBinCenter(i) + 0.5*h_tot.GetBinWidth(i)
                    if pullsRatio:
                        y_up = (c_MC+e_MC)/c_MC
                        y_down = (c_MC-e_MC)/c_MC
                    else:
                        y_up = e_MC/e if e > 0 else 0
                        y_down = -y_up
                    g_up.SetPoint(2*i-1, x_low, y_up)
                    g_up.SetPoint(2*i, x_up, y_up)
                    g_down.SetPoint(2*i-1, x_low, y_down)
                    g_down.SetPoint(2*i, x_up, y_down)
                g_up.SetPoint(2*i+1, x_up, yExp)
                g_down.SetPoint(2*i+1, x_up, yExp)

                h_dr.GetYaxis().SetTitleOffset(0.35)
                h_dr.GetYaxis().SetTitleSize(0.2)
                h_dr.GetYaxis().SetLabelSize(0.2)
                h_dr.GetXaxis().SetTitleOffset(0.95)
                h_dr.GetXaxis().SetTitleSize(0.22)
                h_dr.GetXaxis().SetLabelSize(0.22)
                h_dr.GetXaxis().SetTickSize(0.07)
                if pullsRatio:
                    h_dr.GetYaxis().SetRangeUser(pulls_ylim[0], pulls_ylim[1])
                    h_dr.GetYaxis().SetNdivisions(402)
                else:
                    h_dr.SetLineWidth(2)
                    ax = h_dr.GetYaxis()
                    ax.SetRangeUser(-4, 4)
                    ax.SetNdivisions(-8)
                    for i in range(1,10):
                        if i == 2: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'-3#sigma')
                        elif i == 8: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'3#sigma')
                        elif i == 5: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'0')
                        else: ax.ChangeLabel(i,-1,0,-1,-1,-1,'')

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
                x_high = h_tot.GetBinCenter(h_dr.GetNbinsX())+0.5*h.GetBinWidth(h_dr.GetNbinsX())
                l.DrawLine(x_low, yExp, x_high, yExp)
                h_dr.Draw('sameE1') #redraw it to bring it to front
                h_dr.Draw('sameP')

                canvas.dnd.append([pad, h_dr, h_tot, g_up, g_down])

    canvas.Draw()
    return canvas

def plot_SingleCategory(CMS_lumi,
                        h_dic,
                        scale_dic={},
                        min_y=1e-4,
                        draw_pulls=False,
                        pulls_ylim=[0.8, 1.2],
                        pullsRatio=False,
                        logy=False,
                        mergeDstst=True,
                        tag='',
                        xtitle='',
                        addText='',
                        addTextPos=[0.18, 0.83],
                        legLoc=[0.65, 0.4, 0.9, 0.7]
                        ):
    if len(tag) > 0 and not tag[0] == '_':
        tag = '_' + tag
    canvas = rt.TCanvas('c_'+tag, 'c_'+tag, 50, 50, 600, 450)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.dnd = []
    pad_master = canvas.cd()

    if not draw_pulls:
        pad = rt.TPad('pmain'+tag, 'pmain'+tag, 0, 0, 1, 1)
        pad.SetBottomMargin(0.2)
    else:
        pad = rt.TPad('pmain'+tag, 'pmain'+tag, 0, 0.25, 1, 1)
        pad.SetBottomMargin(0.)
    pad.SetTopMargin(0.07)
    pad.SetRightMargin(0.05)
    pad.SetLeftMargin(0.13)
    pad.Draw()
    pad.cd()

    h = h_dic['data'].Clone('h_aux_data'+tag)
    if 'data' in scale_dic.keys(): h.Scale(scale_dic['data'])
    h.SetLineColor(1)
    h.SetLineWidth(1)
    h.SetMarkerColor(1)
    h.SetMarkerStyle(20)
    h.GetXaxis().SetTitle(xtitle+' [GeV]')
    h.GetXaxis().SetTitleSize(0.07)
    h.GetXaxis().SetLabelSize(0.07)
    h.GetYaxis().SetTitleOffset(1.16)
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetLabelSize(0.07)
    h.GetYaxis().SetTitle('Candidates / {:.2f} '.format(h.GetBinWidth(3)) + 'GeV')
    max_y = np.max([hhh.GetMaximum() for hhh in h_dic.values()])*1.3
    h.GetYaxis().SetRangeUser(min_y, max_y)
    h_list = [h]

    procOrder = ['tau', 'Hc', 'mu']
    for iProc, procName in enumerate(procOrder):
        if not procName in h_dic.keys(): continue
        h = h_dic[procName].Clone('h_aux'+'_'+procName+tag)
        if procName in scale_dic.keys(): h.Scale(scale_dic[procName])
        h.SetLineWidth(0)
        h.SetFillColor(col_dic[procName])
        h.SetFillStyle(1)
        h.Sumw2(0)
        hh = h_list[-1]
        if not hh.GetName() == 'h_aux_data'+tag:
            h.Add(hh)
        h_list.append(h)

    for k, sampleList in sampleDstst.iteritems():
        for iproc, procName in enumerate(sampleList):
            if not procName in h_dic.keys(): continue
            h = h_dic[procName].Clone('h_aux'+'_'+procName+tag)
            if procName in scale_dic.keys(): h.Scale(scale_dic[procName])
            h.SetLineWidth(0)
            h.SetFillColor(col_dic[k])
            h.SetFillStyle(1)
            if not mergeDstst: h.SetFillStyle(fillStyleVar[iproc])
            h.Sumw2(0)
            hh = h_list[-1]
            if not hh.GetName() == 'h_aux_data'+tag:
                h.Add(hh)
            h_list.append(h)


    h_list[0].Draw('E1')
    for h in reversed(h_list[1:]):
        h.Draw('SAME')
    h_list[0].Draw('SAMEE1') #Draw it a second time to bring it in foreground

    if addText:
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.DrawLatexNDC(addTextPos[0], addTextPos[1], addText)

    CMS_lumi.CMS_lumi(pad, -1, 33, cmsTextSize=0.75*1.2, lumiTextSize=0.6*1.2)
    if logy:
        pad.SetLogy()

    leg = createLegend(h_list, h_dic, canvas, loc=legLoc, cat_name=tag)
    leg.Draw()
    canvas.dnd.append([pad, h_list, leg])

    if draw_pulls:
        pad_master.cd()

        pad = rt.TPad('ppull'+tag, 'ppull'+tag, 0, 0, 1, 0.25)
        pad.SetBottomMargin(0.5)
        pad.SetTopMargin(0)
        pad.SetRightMargin(0.05)
        pad.SetLeftMargin(0.13)
        pad.Draw()
        pad.cd()

        h_dr = h_list[0].Clone('h_aux_dataratio'+tag)
        h_dr.GetYaxis().SetTitle('Pull')
        if pullsRatio:
            h_dr.GetYaxis().SetTitle('RD/MC')
        h_tot = h_dic['total'].Clone('h_aux_total'+tag)
        g_up = rt.TGraph()
        yExp = 1 if pullsRatio else 0
        g_up.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), yExp)
        g_down = rt.TGraph()
        g_down.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), yExp)
        for i in range(1, h_dr.GetNbinsX()+1):
            c = h_dr.GetBinContent(i)
            e = h_dr.GetBinError(i)
            c_MC = h_tot.GetBinContent(i)
            e_MC = h_tot.GetBinError(i)
            if pullsRatio:
                h_dr.SetBinContent(i, c/c_MC)
                h_dr.SetBinError(i, e/c_MC)
            else:
                if e > 0: h_dr.SetBinContent(i, (c-c_MC)/e)
                elif e_MC > 0: h_dr.SetBinContent(i, (c-c_MC)/e_MC)
                else: h_dr.SetBinContent(i, 0)
                h_dr.SetBinError(i, 0)

            x_low = h_tot.GetBinCenter(i) - 0.5*h_tot.GetBinWidth(i)
            x_up = h_tot.GetBinCenter(i) + 0.5*h_tot.GetBinWidth(i)
            if pullsRatio:
                y_up = (c_MC+e_MC)/c_MC
                y_down = (c_MC-e_MC)/c_MC
            else:
                y_up = e_MC/e if e > 0 else 0
                y_down = -y_up
            g_up.SetPoint(2*i-1, x_low, y_up)
            g_up.SetPoint(2*i, x_up, y_up)
            g_down.SetPoint(2*i-1, x_low, y_down)
            g_down.SetPoint(2*i, x_up, y_down)
        g_up.SetPoint(2*i+1, x_up, yExp)
        g_down.SetPoint(2*i+1, x_up, yExp)

        h_dr.GetYaxis().SetTitleOffset(0.35)
        h_dr.GetYaxis().SetTitleSize(0.2)
        h_dr.GetYaxis().SetLabelSize(0.2)
        h_dr.GetXaxis().SetTitleOffset(0.95)
        h_dr.GetXaxis().SetTitleSize(0.22)
        h_dr.GetXaxis().SetLabelSize(0.22)
        h_dr.GetXaxis().SetTickSize(0.07)
        if pullsRatio:
            h_dr.GetYaxis().SetRangeUser(pulls_ylim[0], pulls_ylim[1])
            h_dr.GetYaxis().SetNdivisions(402)
        else:
            h_dr.SetLineWidth(2)
            ax = h_dr.GetYaxis()
            ax.SetRangeUser(-4, 4)
            ax.SetNdivisions(-8)
            for i in range(1,10):
                if i == 2: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'-3#sigma')
                elif i == 8: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'3#sigma')
                elif i == 5: ax.ChangeLabel(i,-1,-1,-1,-1,-1,'0')
                else: ax.ChangeLabel(i,-1,0,-1,-1,-1,'')

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
        x_low = h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1)
        x_high = h_tot.GetBinCenter(h_dr.GetNbinsX())+0.5*h_tot.GetBinWidth(h_dr.GetNbinsX())
        l.DrawLine(x_low, yExp, x_high, yExp)
        h_dr.Draw('sameE1') #redraw it to bring it to front
        h_dr.Draw('sameP')

        canvas.dnd.append([pad, h_dr, h_tot, g_up, g_down])

    canvas.Draw()
    return canvas
