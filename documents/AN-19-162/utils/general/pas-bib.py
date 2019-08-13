#!/usr/bin/env python

"""Script to generate a BibTeX file(s) of the CMS PASs/Papers.

    Author: G. Alverson, Northeastern University
    Creation Date: 22-SEP-2010
    """

__version__ = "$Revision: 44192 $"
#$HeadURL: svn+ssh://svn.cern.ch/reps/tdr2/utils/trunk/general/pas-bib.py $
#$Id: pas-bib.py 44192 2011-03-06 17:00:00Z alverson $

import re
import shutil
import ldap
import socket
import os
import shelve
from xml.dom import Node

def extractBalanced(text, delim):
    """ Extract a delimited section of text: available opening delimiters are '{', '"', and  '<' """
    delims = {"{":"}", '"':'"', "<":">"} # matching closing delims
    pin = text.find(delim) + 1
    nbraces = 1;
    pout = pin
    while nbraces > 0:
        if pout > len(text): 
            print "extractBalanced >>> Error parsing text"
            return [0, None] # probably unmatched } inside TeX comment string
        if text[pout:pout+2] == '\\'+delim: # look for escaped delim
            pout += 2
        else:
            if text[pout:pout+2] == '\\'+delims[delim]:
                pout += 2
            else:
                if text[pout:pout+1] == delims[delim]:
                    nbraces -= 1
                elif text[pout:pout+1] == delim:
                    nbraces += 1
                pout += 1
    return [pout, text[pin:pout-1]]

class DocListException(Exception):
    """Base class for exceptions in this module."""
    pass

class DocListBadTag(DocListException):
    """Doc name does not parse correctly."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
       
class DocListBadXML(DocListException):
    """Doc XML does not parse correctly."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
        

class DocList:
    """This class is for creating a BibTeX file from the CDS database. """
    
    def __init__(self, verbosity, remoteUser, noteType, overwrite):
        #self._bib = shelve.open('pas-bib-shelf',writeback=True) # basic information on documents
        self._bib = dict()
        self._dlist = dict() # list of all documents of current type in repository (svn)
        self._verbosity = verbosity
        self._noteRepo = "svn.cern.ch/reps/tdr2" # svn repository location
        self._overwrite = False # don't overwrite current bib info with that from CDS or svn
        if remoteUser:
            self._noteRepo = remoteUser+"@"+self._noteRepo # used to override username for remote svn access (if necessary). Best w/ key agents.
        self._noteType = noteType # notes (really pas); alternative is 
        self._overwrite = overwrite
        if self._noteType == "papers":
            self._bibFilename = "CMSPapersBib.bib"
        else:
            self._bibFilename ="pasBib.bib"

        
    def close(self):
        """ Cleans up and writes the persistent copy of the information """
        #self._bib.close()
        
    def getDocInfoFromBib(self):
        """ Open an existing bib file and load it into memory """
        if (self._verbosity > 3): print "+> getDocInfoFromBib"
        if not os.path.exists(self._bibFilename):
            if self._verbosity > 0:
                print "! No existing bib file\n"
            return
        f = open(self._bibFilename,"r")
        text = f.read()
        mtag = re.compile("\s*([\w|-]*)\s*,") # looking for match to document tag
        ftag = re.compile(",{0,1}\s*(\w*)\s*=\s*") # looking for field inside document listing; optional comma from end of previous entry
        urltag = re.compile("http://cdsweb.cern.ch/record/(\d*)")
        pin = 0

        while pin < len(text):
            # find the next article (always ARTICLE)
            nextArticle = text[pin:].find('@ARTICLE')
            article = dict()
            if nextArticle > -1:
                pin += nextArticle
                [pout, entry] = extractBalanced(text[pin:],"{")
                nextp = 0
                # process document tag
                taginfo = mtag.match(entry[nextp:])
                tag = taginfo.group(1)
                if (self._verbosity > 4): print "Tag: ",tag
                nextp += taginfo.span()[1]
                # process document fields
                while nextp < len(entry):
                    fieldinfo = ftag.match(entry[nextp:])
                    if fieldinfo:                
                        fieldname = fieldinfo.group(1).upper() # uppercase as we must pick convention for dict labels
                        nextp += fieldinfo.span()[1]
                        [n, field] = extractBalanced(entry[nextp:],entry[nextp:nextp+1])
                        if (self._verbosity > 4): print "    ",fieldname, ":", entry[nextp+1:nextp+n-1]
                        nextp += n
                        if fieldname == "URL":
                            article["cdsID"] =  urltag.match(field).group(1) 
                        elif fieldname == "TITLE":
                            article["svnTitle"] = field
                        elif fieldname == "YEAR":
                            article["cdsDate"] = field
                    else:
                        nextp = len(entry) # jump out on trailing junk at end of entries
                pin += nextp
                self._bib[tag[-10:]] = article # uses last 10 chars, example: CMS_PAS_BPH-10-002 
            else:
                pin = len(text)
                if (self._verbosity > 4): print "Done with importing bib entries"
        f.close()        
            
    def getDatafieldValue(self, record, xmltag, subfield):
        """ Parse a CDS XML record and extract subfield values."""

        if (self._verbosity > 3): print "+> getDatafieldValue"
        data = record.getElementsByTagName("datafield")
        for datum in data:
            if (datum.getAttribute("tag") == xmltag):
                subdata = datum.getElementsByTagName("subfield")
                for subdatum in subdata:
                    if (subdatum.getAttribute("code")==subfield) :
                        text = subdatum.firstChild
                        if (not text):
                            return ('')
                        else:
                            val = text.data
                        return (val)
        return ('')
    
    def getDatafieldValueList(self, record, xmltag, subfield):
        """ Parse a CDS XML record and extract a list of values for the subfield."""

        if (self._verbosity > 3): print "+> getDatafieldValueList"
        data = record.getElementsByTagName("datafield")
        val = []
        for datum in data:
            if (datum.getAttribute("tag") == xmltag):
                subdata = datum.getElementsByTagName("subfield")
                for subdatum in subdata:
                    if (subdatum.getAttribute("code")==subfield) :
                        text = subdatum.firstChild
                        if text:
                            val.append(text.data)
        if val:
            return (val)
        else:
            return ('')

    def getDocInfoFromCDS(self, tag):
        """ Get the list of PASs/Papers from CDS and parse. If tag is present, only do tag."""
        from urllib2 import urlopen
        from xml.dom.minidom import parse

        if (self._verbosity > 2): print "+> getDocInfoFromCDS"
        tagparse = re.compile('(?:CMS-){0,1}([A-Za-z]{3})-(\d{2})-(\d{3})') # parse XXX-YY-NNN
        dateparse = re.compile('(\d{4})|.*(\d{4})') # either YYYYMMDD or YYYY-MM-DD, but only need year. Also 'DD Mmm YYYY'

        url = "http://cdsweb.cern.ch/search?cc=CMS+Physics+Analysis+Summaries&of=xm"
        if self._noteType == "papers":
                    url = "http://cdsweb.cern.ch/search?cc=CMS%20Papers&of=xm"
        if (tag):
            m = tagparse.match(tag)
            if ( not m):
                raise DocListException ( "Note tag %s is not of the form XXX-YY-NNN" % tag )
            url = url+"&p=037__a:CMS-PAS-"+m.group(1)+"-"+m.group(2)+"-"+m.group(3)

        xobj = urlopen(url)
        dom = parse(xobj)
        commentnode = dom.firstChild
        if (commentnode.nodeType != Node.COMMENT_NODE):
            raise  DocListException ( "Did not get total record count as XML comment" )
            # expect: " Search-Engine-Total-Number-Of-Results: 149"
        mtotal = re.match('[^\d]+(\d+)\s*$',commentnode.data)
        if (not mtotal):
            raise  DocListException ( "Did not get total record count as XML comment: %s" % commentnode.data )
        totalrecs = int(mtotal.group(1))
        records = dom.getElementsByTagName("record")
        processed = 0
        while ( processed < totalrecs ):
            if (opts.verbose): print "Retrieved ", records.length, "records"
            for record in records:
                controls =  record.getElementsByTagName("controlfield")
                for control in controls:
                    if (control.getAttribute("tag") == "001"):
                        text = control.firstChild
                        if (not text): raise DocListBadXML
                        cdsID = text.data
                        break
                if (not cdsID): raise DocListBadXML
                cdsTag = None
                if self._noteType == "papers":
                    cdsTagList = self.getDatafieldValueList(record,"088","a") # example: CMS-PAS-EXO-10-005 or CERN-PH-EP-...
                    for t in cdsTagList:
                        m = tagparse.search(t)
                        if m:
                            cdsTag = m.group(1)+'-'+m.group(2)+'-'+m.group(3)
                            break
                else:
                    cdsTag = self.getDatafieldValue(record,"037","a") # example: CMS-PAS-EXO-10-005
                    m = tagparse.search(cdsTag)
                    if (not m):
                        raise DocListError ( "Note tag returned from CDS %s is not of the form XXX-YY-NNN" % tag )
                if cdsTag:
                    cdsTitle = self.getDatafieldValue(record,"245","a")
                    cdsDate = self.getDatafieldValue(record,"269","c")
                    if self._noteType == "papers": 
                        cdsJournal = self.getDatafieldValue(record,"773","p")
                        cdsDoi = self.getDatafieldValue(record,"773","a")
                        cdsVolume = self.getDatafieldValue(record,"773","v")
                        cdsYear = self.getDatafieldValue(record,"773","y")
                        cdsPages = self.getDatafieldValue(record,"773","c")
                        cdsArXiv = self.getDatafieldValue(record,"037","a")
                    if (not cdsTitle): 
                        cdsTitle="++> FIX ME: Title not found! <++"
                    if (not cdsDoi):
                        if (self._verbosity > 0) : print " >> CDS Paper missing doi entry (will skip as unpublished). Tag: ",cdsTag
                    else:                    
                        if (self._verbosity > 1) : 
                            try:
                                print " >> CDS Info-->Tag: ",cdsTag, "cdsId: ", cdsID, "cdsTitle: ",cdsTitle
                            except UnicodeEncodeError:
                                print " >> CDS Info-->Tag: ",cdsTag, "cdsId: ", cdsID, "cdsTitle: -- contains unprintable unicode chars -- "
                            
                        tag = m.group(1)+"-"+m.group(2)+"-"+m.group(3)
                        tag = tag.encode('latin-1') # necessary for pre python 3 shelve
                        n = dateparse.match(cdsDate)
                        if (not n):
                            cdsDate = "20"+m.group(2) # rough guess at year if not properly placed in CDS record 
                        else:
                            cdsDate = max(n.groups()) # can be in either 1st or 2nd group
                        bibentry = { "tag": cdsTag, "cdsTitle": cdsTitle, "cdsID": cdsID, "cdsDate": cdsDate}
                        if self._noteType == "papers":
                            bibentry.update({"cdsJournal": cdsJournal, "cdsDoi": cdsDoi, "cdsVolume": cdsVolume, "cdsYear": cdsYear, "cdsPages": cdsPages, "cdsArXiv": cdsArXiv})
                        if (not tag in self._bib) or self._overwrite:
                            if (self._verbosity > 1) : print " >> New/overwritten entry: ", tag
                            self._bib[tag] = bibentry
            processed += records.length
            if ( processed < totalrecs ):   # should think of some way to avoid checking twice
                xobj = urlopen(url+"&jrec={0}".format(processed))
                dom = parse(xobj)
                records = dom.getElementsByTagName("record")

    def titleFromSVNtoBib(self, tag):
        """ Extract the title directly from the svn repository"""
        import tempfile
        import subprocess
        
        if (self._verbosity > 2): print "+> titleFromSVNtoBib"
        if  tag and "svnTitle" in self._bib[tag] and not self._overwrite:
            return
        wd = os.getcwd()
        tmpd = tempfile.mkdtemp(prefix="qsvn-")
        if (self._verbosity > 1): print " >> Using temp directory "+tmpd
        os.chdir(tmpd)
        dump = open(os.devnull,"w")
        if (self._verbosity > 1): print " >> Fetching repository info from svn: this takes a while"
        subprocess.call("svn co -N svn+ssh://"+self._noteRepo, stdout=dump ) #lxplus
        subprocess.call("svn update -N tdr2/"+self._noteType, stdout=dump)
        if not type == "papers" : subprocess.call("svn update -N tdr2/"+self._noteType, stdout=dump)
        dump.close
        if (self._verbosity > 1): print " >> Getting list of all " + self._noteType + " from the repository"
        os.chdir("tdr2/"+self._noteType)
        proc = subprocess.Popen("svn list", stdout=subprocess.PIPE)
        svn_list = proc.communicate()[0] # get list of possible notes/papers
        self._dlist = re.split("/*\r*\n",svn_list) # remove eol indicator and possible directory indicator (/)
        doc_match = re.compile("[A-Z]{3}-\d{2}-\d{3}}") # only PAS/Papers
        self._dlist = filter( doc_match.search, self._dlist ) 
        if (self._verbosity > 1): print " >> Getting LaTeX titles from the repository"
        if tag:
            if True: # tag in self._bib:
                title = self.extractTitleFromSVN(tag)
                if (self._verbosity > 1): print ">>> Inserting ",tag
                self._bib[tag]["svnTitle"] = title
            else:
                if (self._verbosity > 1): print "+> Tag %s not found in CDS" % tag
        else:
            for doc, data in self._bib.iteritems():
                if (not "svnTitle" in self._bib[doc]) or self._overwrite:
                    title = self.extractTitleFromSVN(doc)
                    if (self._verbosity > 1): print ">>> Inserting ",tag
                    self._bib[doc]["svnTitle"] = title
        # delete tmp directory... to do
        os.chdir(wd)
                
        
    def extractTitleFromSVN(self, tag):
        """ export the basic doc directory and extract the title from the head tex file (must be "tag".tex) """
        
        import subprocess
        
        if (self._verbosity > 2): print "+> extractTitleFromSVN"
        fname = self._noteType+"/"+tag+"/trunk/"+tag+".tex"
        dump = open(os.devnull,"w")
        subprocess.call("svn export svn+ssh://"+self._noteRepo+"/"+fname, stdout=dump)
        dump.close
        f = open(tag+".tex","r")
        text = f.read()
        # extract the title looking for balanced braces: there is no extract_bracketed in standard python
        pin = text.find('\\title')
        if  pin > -1:
            [x, title] = extractBalanced(text[pin:],"{")
            if not title:
                print "extractTitleFromSVN >>> Empty title for "+tag
            return title
        else:
            return ""
        
            
     
            
    def generateBib(self):
    
        from datetime import datetime 
        
        f = open(self._bibFilename,"w") # wipes out old version
        f.write("BibFile generated by pas-bib version {0}, ".format(version)+datetime.utcnow().strftime("%Y-%m-%d %H:%m UTC")+"\n")
        if self._noteType == "papers":
            keybase = "CMS-PAPERS-"
        else:
            keybase = "CMS-PAS-"
        for doc in sorted(self._bib, key=lambda z: z[4:6]+z[0:2]+z[7:10], reverse=True ): # sort descending by (year, group, number)
            vals = self._bib[doc]
            if (self._noteType != "papers" or "cdsDoi" in vals):  # test missing cdsDoi to tag unpublished papers
                entry  = '@ARTICLE{'+keybase+'{0},\n'.format(doc)
                entry += '      AUTHOR      = "{CMS Collaboration}",\n'
                entry += '      COLLABORATION = {CMS},\n'
                if "svnTitle" in vals :
                    stripped = re.sub(r"\\\\"," ",vals["svnTitle"]) # remove TeX line breaks
                    entry += '      TITLE       = "{0}",\n'.format(stripped)
                else:
                    entry += '      TITLE       = "{0}",\n'.format(vals["cdsTitle"])
                if self._noteType == "papers":
                    entry += '      JOURNAL     = "{0}",\n'.format(vals["cdsJournal"])
                    entry += '      VOLUME      = "{0}",\n'.format(vals["cdsVolume"])
                    entry += '      DOI         = "{0}",\n'.format(vals["cdsDoi"])          
                    entry += '      YEAR        = "{0}"\n'.format(vals["cdsDate"][0:4]) # no comma for last entry  
                else:
                    entry += '      URL         = "http://cdsweb.cern.ch/record/{0}",\n'.format(vals["cdsID"])
                    entry += '      JOURNAL     = "CMS Physics Analysis Summary",\n'
                    #entry += '      VOLUME      = "\href{http://cdsweb.cern.ch/record/'+'{0}'.format(vals["cdsID"])+'}{'+'CMS-PAS-{0}'.format(doc)+'}",\n'
                    entry += '      VOLUME      =  "CMS-PAS-{0}",\n'.format(doc)
                    entry += '      YEAR        = "{0}"\n'.format(vals["cdsDate"][0:4]) # no comma for last entry
                entry += '}\n'
                f.write(entry)
        f.close()
            
                        
         
def main(argv):
    import sys
    from optparse import OptionParser

    usage = "Usage: %prog [options]  [tag]\n\t This program will update "
    pat = re.compile("\$Revision:\s+(\d+)\s+\$")
    global version
    version = pat.search(__version__).group(1)
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-v", "--verbosity", action="count", dest="verbose", default=False,
                        help="trace script execution")
    parser.add_option(  "--remote-user", action="store", dest="remoteUser",
                        help="remote user name for svn+ssh when Kerberos authentication is unavailable")
    parser.add_option("-s",  "--svn", action="store_true", dest="svn", help="get titles from svn repository in preference to CDS")
    parser.add_option("-t", "--type", action="store", dest="type", default="notes", choices=("papers","notes"),
                        help="note type: notes (PAS) [default], or papers")
    parser.add_option("-o","--overwrite", action="store_false", dest="overwrite",
                        help="normally existing bib entries are not overwritten with CDS or svn information")
    global opts
    (opts, args) = parser.parse_args()
    if opts.verbose:
        print "\tVerbosity = %s" % opts.verbose
    tag = ""
    if len(args) > 0:
        tag = args[len(args)-1]
        
    pas = DocList(opts.verbose, opts.remoteUser, opts.type, opts.overwrite)
    pas.getDocInfoFromBib()
    pas.getDocInfoFromCDS(tag)
    if opts.svn:
        pas.titleFromSVNtoBib(tag)
    pas.generateBib()
    
    pas.close()
        
if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
    print opts
