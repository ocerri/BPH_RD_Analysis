import sys,os
from optparse import OptionParser
from subprocess import call,check_output

usage=''' Compute and compare ratio with combine and C-P intervals.
        * python compare.py
        * python compare.py -n 3000 -d 5000 --muD='4500,5500' --rmin=0.55 --rmax=0.75
        '''
parser=OptionParser(usage=usage)
parser.add_option("-n","--numerator",help="Numerator [%default]",type='int',default=3)
parser.add_option("-d","--denominator",help="Denominator [%default]",type='int',default=5)

parser.add_option("","--format",help="Format string for numbers [%default]",type='string',default="%.2f")

parser.add_option("","--muD",help="Range for muD: [-xx,+yy] [%default]",type='string',default="")
parser.add_option("","--rmin",help="min boundary for rMin [%default]",type='string',default="0")
parser.add_option("","--rmax",help="min boundary for rMax [%default]",type='string',default="4")
opts,args=parser.parse_args()

opts.poi="r"

def mycall(cmd):
    status = call(cmd,shell=True)
    if status != 0:
        print "<*> ERROR: unable to execute '"+cmd+"'"
        raise IOError
    
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

## UTILITIES FUNCTIONS
def findQuantile(pts,cl):
	#gr is a list of r,nll
	# start by walking along the variable and check if crosses a CL point
	if cl<=0:  
	 min=pts[0][0]
	 mincl=pts[0][1]
	 for pt in pts: 
		if pt[1]<mincl: 
			mincl=pt[1]
			min = pt[0]
     
	 return min,min

	crossbound = [ pt[1]<=cl for pt in pts ]
	rcrossbound = crossbound[:]
	rcrossbound.reverse()

	minci = 0
	maxci = len(crossbound)-1
	min = pts[0][0]
	max = pts[maxci][0]

	for c_i,c in enumerate(crossbound): 
		if c : 
			minci=c_i
			break
	
	for c_i,c in enumerate(rcrossbound): 
		if c : 
			maxci=len(rcrossbound)-c_i-1
			break

	if minci>0: 
		y0,x0 = pts[minci-1][0],pts[minci-1][1]
		y1,x1 = pts[minci][0],pts[minci][1]
		min = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)
		
	if maxci<len(crossbound)-1: 
		y0,x0 = pts[maxci][0],pts[maxci][1]
		y1,x1 = pts[maxci+1][0],pts[maxci+1][1]
		max = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)

	return min,max

def smoothNLL_v2(gr,res,w=0.3):

  print "-> Smooth v2"
  minVal = min([re[0] for re in res])
  maxVal = max([re[0] for re in res])

  gr2=ROOT.TGraph()

  myfunc = ROOT.TF1("myfunc","pol2",-100,100)
  for p in range(100):
    x = minVal+p*((maxVal-minVal)/100.)
    fitmin = x-w
    fitmax = x+w
    myfunc.SetRange(fitmin,fitmax)
    gr.Fit("myfunc","RQN")
    gr.Fit("myfunc","RQNM")
    y = myfunc.Eval(x)
    if y<5:
        gr2.SetPoint(gr2.GetN(),x,y)
  return gr2

def smoothNLL_v3(gr,res,w=0.3,delta=0.05):
    print "-> Smooth v3 <-"
    gr2 = smoothNLL_v2(gr,res,w)
    gr3 = ROOT.TGraph()
    for x,y in res:
        y2 = gr2.Eval(x)
        if abs(y-y2) < delta:
            gr3.SetPoint(gr3.GetN(), x,y ) 
    gr4 = smoothNLL_v2(gr3,res,w)
    return gr4

def cleanSpikes1D(rfix):

 # cindex is where deltaNLL = 0 (pre anything)
 MAXDER = 1.0
 for i,r in enumerate(rfix):
   if abs(r[1]) <0.001: cindex = i

 lhs = rfix[0:cindex]; lhs.reverse()
 rhs= rfix[cindex:-1]
 keeplhs = []
 keeprhs = []

 for i,lr in enumerate(lhs): 
   if i==0: 
        prev = lr[1]
        idiff = 1
   if abs(lr[1]-prev) > MAXDER :
        idiff+=1
        continue 
   keeplhs.append(lr)
   prev = lr[1]
   idiff=1
 keeplhs.reverse()

 for i,rr in enumerate(rhs):
   if i==0: 
    prev = rr[1]
    idiff = 1
   if abs(rr[1]-prev) > MAXDER : 
   	idiff+=1
   	continue 
   keeprhs.append(rr)
   prev = rr[1]
   idiff=1
 
 rfix = keeplhs+keeprhs
 
 rkeep = []
 #now try to remove small jagged spikes
 for i,r in enumerate(rfix):
   if i==0 or i==len(rfix)-1: 
   	rkeep.append(r)
   	continue
   tres = [rfix[i-1][1],r[1],rfix[i+1][1]]
   mean = float(sum(tres))/3.
   mdiff = abs(max(tres)-min(tres))
   if abs(tres[1] - mean) > 0.6*mdiff :continue
   rkeep.append(r)
 return rkeep

### PARSE DATACARD AND RUN COMBINE
mycall("cp datacard2.txt datacard_tmp.txt")
mycall("sed -i'' 's/$NUM/%d/g' datacard_tmp.txt"%(opts.numerator))
mycall("sed -i'' 's/$DEN/%d/g' datacard_tmp.txt"%(opts.denominator))
mycall("text2workspace.py --X-allow-no-background -o tmp.root datacard_tmp.txt")
extra=""
if opts.muD !="":
    #extra += "--setPhysicsModelParameterRanges muD="+opts.muD
    extra += "--setParameterRanges muD="+opts.muD
mycall("combine -M MultiDimFit --algo=grid --points=1000 --rMin=%s --rMax=%s %s tmp.root"%(opts.rmin,opts.rmax,extra))

fIn=ROOT.TFile.Open("higgsCombineTest.MultiDimFit.mH120.root")
tree = fIn.Get('limit')

res=[]
for i in range(tree.GetEntries()):
  tree.GetEntry(i)
  xv = getattr(tree,opts.poi)
  if tree.deltaNLL<0 : print "Warning, found -ve deltaNLL = ",  tree.deltaNLL, " at ", xv 
  if 2*tree.deltaNLL < 100:
    res.append([xv,2*tree.deltaNLL])
res.sort()

obs=ROOT.TGraph()
for re, nll in res: 
    if nll>=0. and nll<5:
        obs.SetPoint(obs.GetN(),re,nll)

m,m1 = findQuantile(res,0);
l,h  = findQuantile(res,1); ## 1sigma
l2,h2  = findQuantile(res,4); ## 2sigma

xmin = m
eplus = h-m
eminus = m-l

c=ROOT.TCanvas('c','c',800,800)
c.SetBottomMargin(0.15)
c.SetLeftMargin(0.15)
c.SetTopMargin(0.05)
c.SetRightMargin(0.05)

dummy= ROOT.TH1D("dummy","dummy",100,0,10)
dummy.GetXaxis().SetTitle("Ratio")
dummy.GetYaxis().SetTitle("-2 * #Delta logL")
dummy.GetXaxis().SetTitleOffset(1.3)
dummy.GetYaxis().SetTitleOffset(1.5)
dummy.Draw("AXIS")
dummy.Draw("AXIS X+ Y+ SAME")
dummy.GetYaxis().SetRangeUser(0,5)
dummy.GetXaxis().SetRangeUser(float(opts.rmin),float(opts.rmax))


w=0.2
y=1 ## all centered at 1
oneSigma_LR = ROOT . TPave (l,y-w,h,y+w)
twoSigma_LR = ROOT . TPave (l2,y-w,h2,y+w)

oneSigma_LR.SetFillColor(ROOT.kYellow)
oneSigma_LR.SetLineColor(ROOT.kYellow)
twoSigma_LR.SetFillColor(ROOT.kGreen)
twoSigma_LR.SetLineColor(ROOT.kGreen)

oneSigma_LR.SetBorderSize(0)
twoSigma_LR.SetBorderSize(0)

twoSigma_LR.Draw()
oneSigma_LR.Draw()

bf = ROOT.TGraph()
bf.SetPoint(0,m,y)
bf.SetMarkerStyle(28)
bf.SetMarkerSize(1.5)
bf.SetMarkerColor(ROOT.kBlack)
bf.Draw("P SAME")

obs.Draw("L SAME")

## DERIVE C-P intervals

t = opts.numerator
s = opts.denominator + opts.numerator

CL1=1.-2.*ROOT.RooStats.SignificanceToPValue(1)
CL2=1.-2.*ROOT.RooStats.SignificanceToPValue(2)

cp1 = ROOT.TEfficiency.ClopperPearson(s,t, CL1,True)
h_CP = cp1/(1.-cp1)
cp1 = ROOT.TEfficiency.ClopperPearson(s,t, CL1,False)
l_CP = cp1/(1.-cp1)

cp2 = ROOT.TEfficiency.ClopperPearson(s,t, CL2,True)
h2_CP = cp2/(1.-cp2)
cp2 = ROOT.TEfficiency.ClopperPearson(s,t, CL2,False)
l2_CP = cp2/(1.-cp2)

##
midp1 = ROOT.TEfficiency.MidPInterval(s,t, CL1,True)
h_MIDP = midp1/(1.-midp1)
midp1 = ROOT.TEfficiency.MidPInterval(s,t, CL1,False)
l_MIDP = midp1/(1.-midp1)

midp2 = ROOT.TEfficiency.MidPInterval(s,t, CL2,True)
h2_MIDP = midp2/(1.-midp2)
midp2 = ROOT.TEfficiency.MidPInterval(s,t, CL2,False)
l2_MIDP = midp2/(1.-midp2)

nominal = float(opts.numerator)/float(opts.denominator)
print "BestFit (Likelihood): %4.4f [%4.4g, %4.4g]" % ( xmin, l,h )
print "BestFit (CP): %4.4f [%4.4g, %4.4g]" % ( nominal, l_CP,h_CP )
print "BestFit (MIDP): %4.4f [%4.4g, %4.4g]" % ( nominal, l_MIDP, h_MIDP )

ROOT.gStyle.SetEndErrorSize(5)
##
y=1
cp = ROOT.TGraphAsymmErrors()
cp.SetPoint(0,nominal,y)
cp.SetPointError(0,nominal-l_CP,h_CP-nominal, 0,0)
cp.SetPoint(1,nominal,y)
cp.SetPointError(1,nominal-l2_CP,h2_CP-nominal, 0,0)
cp.SetMarkerStyle(20)
cp.SetMarkerColor(ROOT.kBlack)
cp.SetLineWidth(4)
cp.Draw("PE SAME")

y=1.07
midp = ROOT.TGraphAsymmErrors()
midp.SetPoint(0,nominal,y)
midp.SetPointError(0,nominal-l_MIDP,h_MIDP-nominal, 0,0)
midp.SetPoint(1,nominal,y)
midp.SetPointError(1,nominal-l2_MIDP,h2_MIDP-nominal, 0,0)
midp.SetMarkerStyle(4)
midp.SetMarkerColor(ROOT.kRed)
midp.SetLineColor(ROOT.kRed)
midp.SetLineWidth(2)
midp.Draw("PE SAME")

y=0.93
err = nominal * ROOT.TMath.Sqrt( 1./opts.numerator + 1./opts.denominator)
std = ROOT.TGraphAsymmErrors()
std.SetPoint(0,nominal,y)
std.SetPointError(0,err,err, 0,0)
std.SetPoint(1,nominal,y)
std.SetPointError(1,2*err,2*err, 0,0)
std.SetMarkerStyle(0)
std.SetMarkerColor(ROOT.kBlue+2)
std.SetLineColor(ROOT.kBlue+2)
std.SetLineWidth(2)
std.Draw("PE SAME")

leg = ROOT.TLegend(.23,.66,.54,.91)
leg.AddEntry(obs,"LR","L")
leg.AddEntry(oneSigma_LR,"1#sigma LR","F")
leg.AddEntry(twoSigma_LR,"1#sigma LR","F")
leg.AddEntry(cp,"C-P","PL")
leg.AddEntry(midp,"MidP","PL")
leg.AddEntry(std,"std error","PL")
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw()

txt= ROOT.TLatex()
txt.SetNDC()
txt.SetTextFont(42)
txt.SetTextSize(0.035)
txt.SetTextAlign(11)

xtxt = 0.55
ytxt = 0.5
wtxt = 0.04
if opts.format != "":
    txt.DrawLatex(xtxt,ytxt+0*wtxt,"#bf{Ratio of %d/%d:}"%(opts.numerator,opts.denominator))
    txt.DrawLatex(xtxt,ytxt-1*wtxt,("LR: ["+opts.format+","+opts.format+"] ["+opts.format+","+opts.format+"]")%(l,h,l2,h2))
    txt.DrawLatex(xtxt,ytxt-2*wtxt,("CP: ["+opts.format+","+opts.format+"] ["+opts.format+","+opts.format+"]")%(l_CP,h_CP,l2_CP,h2_CP))
    txt.DrawLatex(xtxt,ytxt-3*wtxt,("MP: ["+opts.format+","+opts.format+"] ["+opts.format+","+opts.format+"]")%(l_MIDP,h_MIDP,l2_MIDP,h2_MIDP))
    txt.DrawLatex(xtxt,ytxt-4*wtxt,("STD: ["+opts.format+","+opts.format+"] ["+opts.format+","+opts.format+"]")%(nominal-err,nominal+err,nominal-2*err,nominal+2*err))

c.Update()
raw_input("ok?")

