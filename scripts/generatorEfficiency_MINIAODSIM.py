import sys

import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Lumis
from DataFormats.FWLite import Handle

handle = {}
handle['genFilter'] = [Handle('GenFilterInfo'), ('genFilterEfficiencyProducer', '', 'SIM')]
handle['genProduct'] = [Handle('GenLumiInfoProduct'), ('generator', '', 'SIM')]

for lumi in Lumis(sys.argv[1]):
    prods = {}
    for k,v in handle.iteritems():
        lumi.getByLabel(v[1], v[0])
        prods[k] = v[0].product()
    N_cuts = prods['genFilter'].numEventsPassed()
    N_gen = prods['genFilter'].numEventsTotal()
    xs = prods['genProduct'].getProcessInfos()[0].lheXSec()

    print N_gen, N_cuts, xs.value(), xs.error()
