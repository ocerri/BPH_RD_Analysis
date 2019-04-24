import numpy as np
import ROOT as rt
import root_numpy as rtnp

std_color_list = [1, 2, 4, 8, 6, 28, 43, 7, 25]

def create_TH1D(x, name='h', title=None, binning=[None, None, None], weights=None, h2clone=None):
    if title is None:
        title = name
    if h2clone == None:
        if binning[1] is None:
            binning[1] = min(x)
        if binning[2] is None:
            if ((np.percentile(x, 95)-np.percentile(x, 50))<0.2*(max(x)-np.percentile(x, 95))):
                binning[2] = np.percentile(x, 90)
            else:
                binning[2] = max(x)
        if binning[0] is None:
            bin_w = 4*(np.percentile(x,75) - np.percentile(x,25))/(len(x))**(1./3.)
            binning[0] = int((binning[2] - binning[1])/bin_w)

        h = rt.TH1D(name, title, binning[0], binning[1], binning[2])
    else:
        h = h2clone.Clone(name)
        h.SetTitle(title)
        h.Reset()

    rtnp.fill_hist(h, x, weights=weights)
    h.binning = binning
    return h

def create_TH2D(sample, name='h', title=None, binning=[None, None, None, None, None, None], weights=None, axis_title = ['','']):
    if title is None:
        title = name
    if binning[1] is None:
        binning[1] = min(sample[:,0])
    if binning[2] is None:
        binning[2] = max(sample[:,0])
    if binning[0] is None:
        bin_w = 4*(np.percentile(sample[:,0],75) - np.percentile(sample[:,0],25))/(len(sample[:,0]))**(1./3.)
        binning[0] = int((binning[2] - binning[1])/bin_w)

    if binning[4] is None:
        binning[4] = min(sample[:,1])
    if binning[5] == None:
        binning[5] = max(sample[:,1])
    if binning[3] == None:
        bin_w = 4*(np.percentile(sample[:,1],75) - np.percentile(sample[:,1],25))/(len(sample[:,1]))**(1./3.)
        binning[3] = int((binning[5] - binning[4])/bin_w)

    h = rt.TH2D(name, title, binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
    rtnp.fill_hist(h, sample, weights=weights)
    h.SetXTitle(axis_title[0])
    h.SetYTitle(axis_title[1])
    h.binning = binning
    return h


def make_ratio_plot(h_list_in, title = "", label = "", in_tags = None, ratio_bounds = [0.1, 4], draw_opt = 'E1'):
    h_list = []
    if in_tags == None:
        tag = []
    else:
        tag = in_tags
    for i, h in enumerate(h_list_in):
        h_list.append(h.Clone('h{}aux{}'.format(i, label)))
        if in_tags == None:
            tag.append(h.GetTitle())

    c_out = rt.TCanvas("c_out_ratio"+label, "c_out_ratio"+label, 600, 800)
    pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.SetLeftMargin(0.15)
    # pad1.SetGridx()
    pad1.Draw()
    pad1.cd()

    leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    c_out.cd(1)

    for i, h in enumerate(h_list):
        if i == 0:
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitle("")
            h.GetYaxis().SetRangeUser(0, 1.05*max(map(lambda x: x.GetMaximum(), h_list)))
            h.GetYaxis().SetTitleOffset(1.5)
            h.GetYaxis().SetTitleSize(0.05)
            h.GetYaxis().SetLabelSize(0.05)
            h.SetTitle(title)
            h.DrawCopy(draw_opt)
        else:
            h.DrawCopy(draw_opt+"same")

        leg.AddEntry(h, tag[i], "lep")

    leg.Draw("same")

    c_out.cd()
    pad2 = rt.TPad("pad2", "pad2", 0, 0, 1, 0.3)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.15)
    # pad2.SetGrid()
    pad2.Draw()
    pad2.cd()

    for i, h in enumerate(h_list):
        if i == 0:
            continue
        elif i == 1:
            h.Divide(h_list[0])
            h.GetYaxis().SetTitleOffset(0.6)
            h.GetYaxis().SetRangeUser(ratio_bounds[0], ratio_bounds[1])
            h.GetYaxis().SetTitleSize(0.12)
            h.GetYaxis().SetLabelSize(0.12)
            h.GetYaxis().SetNdivisions(506)
            h.GetXaxis().SetTitleOffset(0.95)
            h.GetXaxis().SetTitleSize(0.12)
            h.GetXaxis().SetLabelSize(0.12)
            h.GetXaxis().SetTickSize(0.07)
            h.SetYTitle('Ratio with {}'.format(tag[0]))
            h.SetTitle("")
            h.DrawCopy(draw_opt)

        else:
            h.Divide(h_list[0])
            h.DrawCopy("same"+draw_opt)

        ln = rt.TLine(h.GetXaxis().GetXmin(), 1, h.GetXaxis().GetXmax(), 1)
        ln.SetLineWidth(3)
        ln.SetLineColor(h_list_in[0].GetLineColor())
        ln.DrawLine(h.GetXaxis().GetXmin(), 1, h.GetXaxis().GetXmax(), 1)


    pad2.Update()

    c_out.pad1 = pad1
    c_out.pad2 = pad2
    c_out.h_list = h_list
    c_out.leg = leg

    return c_out
