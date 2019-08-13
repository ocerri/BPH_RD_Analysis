#!/usr/bin/env perl

use Getopt::Long;
use File::Copy;
use File::Spec::Functions; # catfile
use File::Basename;
use Cwd;              # like Unix pwd but more portable
use Cwd 'abs_path';
use File::Glob 'bsd_glob';
use HTML::Entities;

my $verbose = '';
my $html_directory = 'pages'; 
my $figure_directory = 'Figures';
my $execute = 1;
GetOptions ('verbose!' => \$verbose,
            'help|?' => \$help,
            'html=s' => \$html_directory,
            'figures=s' => \$figure_directory,
            'execute!' => \$execute);

if ($help)
{
    print
    "\n >>> patch_html2: Patches the html output of the Photoshop Web Gallery option to contain \n
    links to the original figures. Run this script in the directory with the ThumbnailFrame.htm file.
    
            Options:\n
              - html: directory containing the html file to patch. Defaults to $html_directory.
              - figures: directory containing the figure files. Defaults to $figure_directory.
              - verbose: produce diagnostic output. Defaults to $verbose.
              - execute: bool to turn off actual execution for debugging. Defaults to $execute.
              - help: produce this message.
              
            Example:\n
            patch_html.pl 
            \n\n";
    exit;
}
##

    my @HTML_Files = bsd_glob("$html_directory/Figure_*.htm");
    foreach $file (@HTML_Files)
    {
       my $stem = basename($file,".htm");
       my @fig_File = bsd_glob("$figure_directory/$stem.*");
       
       if ($verbose)
       {
         print "Processing $file and $fig_File[0]\n";
       }  
       # using the $^I option with <> allows for in-place editing (a la perl -p -i.bak xxxx)
       {local ($^I,@ARGV) = ('.bak',$file);
         if ($execute)
         {
           while (<>)
           {
             my $out = sprintf("../%s",encode_entities($fig_File[0],' '));
#             print $out;
             s|789MARKER|$out|;
             print;
           }
         }
       }  
     }  
         
       
