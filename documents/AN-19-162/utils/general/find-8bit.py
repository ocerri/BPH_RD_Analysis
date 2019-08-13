#!/usr/bin/env python
from __future__ import print_function

"""Script look for 8 bit characters
    """

__version__ = "$Revision: 5556 $"
#$HeadURL: svn+ssh://alverson@svn.cern.ch/reps/admin/tdr2/conf/new-note.py $
#$Id: new-note.py 5556 2011-01-14 13:30:27Z alverson $

import re
import os
import string
import glob

def main(argv):
	import sys
	from optparse import OptionParser

	usage = "Usage: %prog [options]  [filenames]"
	pat = re.compile("\$Revision:\s+(\d+)\s+\$")
	version = pat.search(__version__)
	global opts
	parser = OptionParser(usage=usage, version=version)
	(opts, args) = parser.parse_args()

	# check if console understands Unicode: unix or Windows; latter can be set with chcp 65001; chcp returns current code page
    # for PowerShell, $OutputEncoding = New-Object -typename System.Text.UTF8Encoding, but some commands don't recognize UTF 8
	# unicodeConsole = (sys.stdout.encoding == 'UTF-8' or sys.stdout.encoding == 'cp65001' or sys.stdout.encoding == 'utf8') 
	unicodeConsole = 'utf8'

	for file in glob.glob(args[0]):
		with open(file,"r",encoding="utf-8") as f:
			print("Scanning {} for Unicode".format(file))
			try:
				body = f.read()
			except UnicodeDecodeError as E:
				print(">Problem with Unicode while reading in the file")
				if (E.reason == 'invalid start byte' or E.reason == 'invalid continuation byte'):
					print('>> looks like {} at location {}'.format(E.reason,E.start))
					print('>> The hundred prior characters:')
					print(E.object[E.start-101:E.start-1])
					return
				else:
					raise
			p0 = re.compile(r"[^\x00-\x7F]")
			pm0 = p0.search(body)
			if pm0:
				p = re.compile("(.*?)([^\\x00-\\x7F]+)(.*?)")
				start = 0
				pm = p.search(body,start)
				while pm:
					x = pm.group(1)
					y = pm.group(2)
					z = pm.group(3)
					print(">non-ASCII characters at positions {}-{}".format(pm.start(2),pm.end(2)))
					if unicodeConsole:
						print(">{}>>...<<{}".format(pm.group(1),pm.group(3)))
						print(">{}>>{}<<{}".format(pm.group(1),pm.group(2),pm.group(3)))
					else:
						try:
							print(">{}^{}".format(pm.group(1),pm.group(3)))
						except UnicodeEncodeError:
							print(">{}... + additional Unicode in the same line".format(pm.group(1)))
					pm = p.search(body,pm.end(2))
			else:
				print("No Unicode characters (x80-xFF) found")


if __name__ == "__main__":
	import sys
	main(sys.argv[1:])
