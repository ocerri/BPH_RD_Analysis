import time, sys

class ProgressBar():
    def __init__(self, maxEntry, percentPrecision=5, showEvery=2, headLabel=''):
        self.maxEntry = maxEntry
        self.percentPrecision = percentPrecision

        nStep = int(100/percentPrecision)
        self.nStep = nStep if nStep <= maxEntry else maxEntry
        self.setpSize = int(maxEntry/self.nStep)

        self.headLabel = headLabel

        self.maxPrintoutLen = 0
        self.showEvery = showEvery
        self.lastShown = showEvery+1

    def show(self, entry, tail_label = ''):
        if entry%self.setpSize==0 or time.time() - self.lastShown > self.showEvery:
            self.lastShown = time.time()
            if entry>0:
                sys.stdout.write('\r')
            else:
                self.startTime = time.time()

            Progress = float(entry)/self.maxEntry
            nStepDone = int(Progress*self.nStep)

            outLine = self.headLabel
            outLine += '['+'#'*nStepDone + '-'*(self.nStep-nStepDone) +']'+'  {}%'.format(int(100*Progress))

            if entry>0:
                timeleft = (self.maxEntry - float(entry))*(time.time() - self.startTime)/float(entry)
                if timeleft<181:
                    outLine += " - ETA:{:5.0f} s   ".format(timeleft)
                elif timeleft<10801:
                    timeleft/=60
                    outLine += " - ETA:{:5.1f} min ".format(timeleft)
                else:
                    timeleft/=3600
                    outLine += " - ETA:{:5.1f} h   ".format(timeleft)

            outLine += tail_label
            if len(outLine) > self.maxPrintoutLen:
                self.maxPrintoutLen = len(outLine)
            sys.stdout.write(outLine)
            sys.stdout.flush()

        if entry==self.maxEntry-1:
            outLine = '\r' + self.headLabel
            outLine +='['+ '#'*self.nStep +']  100% - Tot. time: {:.1f} s'.format(time.time() - self.startTime)
            if self.maxPrintoutLen > len(outLine):
                outLine += ' '*(self.maxPrintoutLen - len(outLine) + 2)
            outLine += '\n'
            sys.stdout.write(outLine)
            sys.stdout.flush()
