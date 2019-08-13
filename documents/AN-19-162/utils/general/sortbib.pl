#!/usr/bin/env perl

############
#
# perl script to sort bib file entries on the author key
#
# sortbib_all.pl file_list
#
# returns the original bib entries sorted on the author
#
# Original: 2005 by G. Alverson, Northeastern University HEP Group
# Last Modified: 2006/05 for full sorting (Thanks to L. Tuura)
#
############


use Text::Balanced qw/extract_bracketed
                      extract_delimited/; # for matching {...}, {{...}}, etc
use File::Glob 'bsd_glob'; # non-os dependent globbing                      
use strict;  # safety mode                    
use Getopt::Long qw(:config require_order); # the require order means that for
                                            #  "scriptName command -switch", -switch is passed in ARGV

# Basic options
my $bibkey = "author"; # can change this to sort on other keys
my $verbose = '';
my $help = '';

GetOptions ('verbose!'    => \$verbose,    # negatable: --noverbose; turn on verbose output
            'help|?'      => \$help,       # echo basic operations and options
            'bibkey=s'    => \$bibkey,     # key to sort on, eg, author, year
            );

if ($help || substr($ARGV[0],0,1) eq 'h')
{
    &print_usage();
    exit;
}

my @author;
my @entries;
if ($^O eq "MSWin32" && $ARGV[0] =~ /\*/ )
{
    @ARGV = bsd_glob($ARGV[0]);
}    
foreach my $file (@ARGV) # grab a list of files from the command line
{
    open (my $fh, $file) or die "Can't open file $file\n";
    $_ .= do { local( $/ ); <$fh> };
}


while (/^(@.*?){/gm) # must start with e.g., @Article{...}
{
    my $entryLabel = $1; # e.g., @Article
    pos = pos() - 1; # back up over the {
    my @result = extract_bracketed($_,"{}"); # look for the entry itself: {...}
    $_ = $result[0]; # entire entry
    my $entry = $_; # keep a copy 
    /\{\s*(.*?),/;     # get key: @Article{ keystring,
    my $tag = $1;
    my $author = ""; # default to none
    if ($bibkey =~ 'tag') 
    {
        $author = $tag;
    }
    else
    { 
        if (/$bibkey\s*=(.*)/sig)
        {
            $_ = $1;
            if ($author = extract_delimited) # pick up author = "xxx"
            {
                #print "Author = $author\n"; 
            }
            elsif ($author = extract_bracketed($_,"{}"))
            {
                #print "Author = $author\n"; # pick up author = {xxx}
            }
            else
            {
                warn "Mis-formatted entry?: $tag\n";
                # may be a macro which doesn't require quotes or braces in 
                # a journal field
                if ($bibkey =~ 'journal')
                { 
                    # pick up everything until the next comma
                    /(.*?),/;
                    $author = $1;
                }
            }
        }
    }
    # implement as array of hashes
    push(@entries,{'DATA'=>$entryLabel.$entry, 'AUTHORS'=>$author}); 
    $_ = $result[1]; #rest of file
}

my @sorted_entries = sort{ lc($$a{'AUTHORS'}) cmp lc($$b{'AUTHORS'}) } @entries;

foreach my $string (@sorted_entries)
{
    print $$string{'DATA'},"\n";
}

#############################################################################
# Subroutines:
#############################################################################


#############################################################################
sub print_usage {
    print "$0 Basic usage: $0 inputBibFile.bib > sortedBibFile.bib \n";
    print "$0 Options: --bibkey=key sets the sortkey to key. \n";
    print "            Use key=tag to sort on the citation label.\n";
    print "            Default = author. \n";
}
