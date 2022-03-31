#!/usr/bin/env python

import sys, os, pickle, time, json, yaml, re
from glob import glob
from prettytable import PrettyTable
import argparse


def make_results_table(version=None, cat=None, includeDirs=None, nPulls=4, skipDirs=None, allowedTags=None, verbose=False):
    dirs = []

    if version is not None and cat is not None:
        results_dir = os.environ['HOME'] + '/public_html/BPH_RDst/Combine/'
        matchingDirs = glob(results_dir + '*'+version+'*'+cat+'*')
        auxDirs = []
        for d in sorted(matchingDirs, key=lambda path: os.stat(path).st_ctime):
            if (not os.path.basename(d).startswith('v'+version+'_')) and (not os.path.basename(d).startswith(version+'_')):
                continue
            if skipDirs is not None:
                if d in skipDirs:
                    continue
            auxDirs.append(d)
        dirs += auxDirs

    if includeDirs is not None:
        auxDirs = []
        for d in includeDirs:
            if cat is not None:
                if cat in d:
                    pass
                elif '{}' in d:
                    d = d.format(cat)
                elif 'low' in d:
                    d = d.replace('low', cat)
                elif 'mid' in d:
                    d = d.replace('mid', cat)
                elif 'high' in d:
                    d = d.replace('high', cat)
                elif 'comb' in d:
                    d = d.replace('comb', cat)
                else:
                    if verbose:
                        print 'No match for category', cat, 'and input', d

            if os.path.isdir(d):
                if d[-1] == '/':
                    d = d[:-1]
                auxDirs.append(d)
        dirs += auxDirs

    if isinstance(dirs, list):
        if len(dirs) > 0:
            if verbose:
                print 'Fetching directories:'
                for d in dirs: print d
        else:
            print 'No existing directories'
            return
    else:
        print '[ERROR] dirs:', dirs
        return

    table = PrettyTable()
    table.field_names = ['Version', 'Sat. GoF', 'Scan [%]', 'Cat comp', 'Top pulls','[sig]']
    table.align['Top pulls'] = 'l'
    table.align['[sigma]'] = 'r'

    for dd in dirs:
        tag = os.path.basename(dd)
        if version is not None:
            tag = tag[tag.find(version+'_')+len(version+'_'):]
        if cat is not None:
            tag = tag[:tag.find(cat)]
        if tag.endswith('_'):
            tag = tag[:-1]
        if not tag:
            tag = 'Baseline'

        if allowedTags is not None:
            match = False
            for tagPattern in allowedTags:
                if re.match(tagPattern, tag):
                    match = True
            if not match:
                continue

        if len(tag) > 30:
            tag = tag[:30] + '...'
        # print 'Tag:', tag


        GoF_file = dd + '/GoF_results.txt'
        pval_sat = '-'
        if os.path.isfile(GoF_file):
            with open(GoF_file) as f:
                for line in reversed(f.readlines()):
                    if 'algoSat' in line:
                        break
                data = [x for x in line.split(' ') if x]
                data[2] = data[2][1:-2]
                pval_sat = data[1][:-4] + ' ('
                if float(data[2]) > 0:
                    pval_sat += data[2] + '%' + ')'
                else:
                    pval_sat += 'x{:.2f}'.format(float(data[1])/float(data[-2])) + ')'


        scan_file = dd + '/scan_results.txt'
        scan_result = '-'
        if os.path.isfile(scan_file):
            with open(scan_file) as f:
                line = f.readlines()[-1][:-1]
                data = [x for x in line.split(' ') if x]
                errUp = float(data[2][1:])
                errDw = float(data[4][1:])
                scan_result = ['{:.1f} + {:.1f}'.format(100*float(data[1]), 100*errUp)]
                scan_result.append(' '*scan_result[0].find('+') + '- {:.1f}'.format(100*errDw))
                if 'Upper lims' in line:
                    scan_result.append('< ' + data[8][:-1])

        catComp_file = dd + '/categoriesCompatibility.txt'
        catComp = '-'
        if os.path.isfile(catComp_file):
            with open(catComp_file) as f:
                content = f.readlines()
                if len(content):
                    line = content[-1]
                    catComp = line.split(' ')[-1][1:-5] + '%'
    #     print catComp

        pulls_file = dd + '/scanNuisanceOut_Base.txt'
        pulls = [['-', '-']]*nPulls
        if os.path.isfile(pulls_file):
            with open(pulls_file) as f:
                for i, line in enumerate(f.readlines()[3:nPulls+3]):
                    pulls[i] = [x for x in line[:-1].replace('|', '').split(' ') if x]
                    if pulls[i][0].startswith('prop_bin'):
                        pulls[i][0] = '~'+pulls[i][0][len('prop_bin'):]

        if dd != dirs[0]:
            table.add_row(len(table.field_names)*[''])
        table.add_row([tag, pval_sat, scan_result[0], catComp, pulls[0][0], pulls[0][1]])
        for i in range(1, max(len(pulls), len(scan_result))):
            auxScan = ''
            if i < len(scan_result):
                auxScan = scan_result[i]
            auxPull = ['', '']
            if i < len(pulls):
                auxPull = pulls[i]
            table.add_row(['', '', auxScan, '', auxPull[0], auxPull[1]])

    return table


def print_with_title(table, title='', show=True):
    if not title:
        print table
        return

    tw = len(table.get_string().split('\n')[0])

    if len(title) > tw-7:
        title = title[:tw-7] + '...'

    nb = tw - len(title) - 2
    nbb, nba = [nb/2, nb/2] if nb%2==0 else [(nb-1)/2, 1+(nb-1)/2]

    out = '+' + '-'*(tw-2) + '+' + '\n'
    out += '|' + nbb*' ' + title + nba*' ' + '|\n'
    out += table.get_string()
    out += '\n'
    if show:
        print out
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize fit results in a table',
                                     epilog='Test example: ./versionSummaryTable.py -v v15_0',
                                     add_help=True)
    parser.add_argument ('--version', '-v', default=None, help='Fit version, automatically fetch matching directories.')
    parser.add_argument ('--category', '-c', default=['low', 'mid', 'high', 'comb'], choices=['low', 'mid', 'high', 'comb'], nargs='*', help='Category.')
    parser.add_argument ('--directories', '-d', default=None, nargs='*', help='Directories to include in the table.')
    parser.add_argument ('--skipDir', default=None, nargs='*', help='Directories to skip in the table.')
    parser.add_argument ('--allowedTags', default=None, nargs='+', help='Regular expression for allowed tags.')
    parser.add_argument ('--nPulls', default=4, type=int, help='Number of parameters to display in pulls.')
    parser.add_argument ('--output', '-o', default=None, help='Output destination.')
    parser.add_argument ('--verbose', default=False, action='store_true', help='Verbose.')
    args = parser.parse_args()

    fullOutput = ''
    if len(args.category):
        for ccc in args.category:
            table = make_results_table(version=args.version, cat=ccc, includeDirs=args.directories,
                                       nPulls=args.nPulls, skipDirs=args.skipDir, allowedTags=args.allowedTags,
                                       verbose=args.verbose)
            fullOutput += print_with_title(table, ccc)
    else:
        table = make_results_table(version=args.version, includeDirs=args.directories,
                                   nPulls=args.nPulls, skipDirs=args.skipDir, allowedTags=args.allowedTags,
                                   verbose=args.verbose)
        print table
        fullOutput += table.get_string()

    if args.output is not None:
        # outfile = results_dir + 'summary_table_{}_{}.txt'.format(version, cat)
        outfile = args.output
        os.system('date > '+outfile)
        with open(outfile, 'a') as f:
            f.write(fullOutput)
