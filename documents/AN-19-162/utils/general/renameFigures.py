#!/usr/bin/env python

"""Script to replace original file names with numbered names, based on output from makeManifest pass on TeX output
    """

__version__ = "$Revision:  $"
#$HeadURL: $
#$Id: $
import re
import os
import string
import glob

def checkDuplicate(name, namesList):
    """Check to see if there is a match (up to, and possibly through, the extention) for the [oldName, newName] pairs in namesList.
        Returns the match (including) or empty and sets the match[0] to name.
        """

    returnVal = ''
    for testName in namesList:
        if (name == testName[0][0:len(name)]):
            testName[0] = name
            returnVal = name
            break
    return returnVal

def main(argv):
    import sys
    from optparse import OptionParser

    usage = "Usage: %prog [options]  [TeX filename] [file list file]"
    pat = re.compile("\$Revision:\s+(\d+)\s+\$")
    version = pat.search(__version__)
    global opts
    parser = OptionParser(usage=usage, version=version)
    (opts, args) = parser.parse_args()
    with open(args[1],"r",encoding="utf-8") as f:
        lines = f.readlines()
        if (lines[0] == ''):
            lines = lines[1:] # skip blank leading line?
        lines = [l.rstrip('\n') for l in lines]
        print(lines)
        names = [l.rsplit(' ') for l in lines] # no embedded blanks allowed!
    with open(args[0],"r+") as f:
        # look for filenames in the TeX that did not include an extension
        p = re.compile("""
            \\\\includegraphics(  # pull in all included graphics files: assumed file has been trimmed of comments
            \[.*\])?            # ignore any options
            {(.*?)              # pull in the filename up to any extension (non-greedy): second group
            ((\.)(.*))?}     # and look for any extension: fifth group; full filename is second+third; fouth is '.'
            """, re.VERBOSE)
        t = f.read()
        f.seek(0)
        m = p.findall(t)
        badnames = [i[1]+i[2] for i in m if not(i[2])]
        #badnames = []
        #for i in m:
        #    if (i[2]):
        #        badnames.append(i[1]+i[2])
        #    else:
        #        badnames.append(i[1]+'.pdf')
        # look for a match in names
        for i in badnames:
            if not(checkDuplicate(i,names)):
                print('>>> Failed to match {0}'.format(i))
        # Now we look for a name clash
        pName = re.compile('(.*/)?(.*)')
        inputNames = [pName.match(i[0]).group(2) for i in names]
        outputNames = [pName.match(i[1]).group(2) for i in names]
        ok = True
        for name in inputNames:
            if name in outputNames:
                print('>>> Name clash! Check/fix {0}'.format(name))
                ok = False
        if not(ok):
            return
        for name in names:
            if (name[0] != name[1]):
                t = t.replace('{'+name[0]+'}','{'+name[1]+'}',1) # avoid non-unique shorter versions by requiring bracketing {}
        f.write(t)
        f.truncate()

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])