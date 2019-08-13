#!/usr/bin/env perl

use Getopt::Long;
use File::Copy;

my $verbose = '';
my $baseDir = 'C:/Documents and Settings/George Alverson/My Documents/CMS/TDRs/ptdr2/ptdr2/fig';
my $outDir = "Figures";
my $logFile = "ptdr2_temp.log";
my $execute = '1';
GetOptions ('verbose!' => \$verbose,
            'help|?' => \$help,
            'baseDir=s' => \$baseDir,
            'outDir=s' => \$outDir,
            'logFile=s' => \$logFile,
            'execute!' => \$execute);

if ($help)
{
    print
    "\n >>> figstrip: pulls references to figures from a TDR log file and copies the actual figures
            to a directory with names equal to the figure numbers.\n

            Note that the log file must be produced using a modified \\includegraphics command.
            
            Since the logo is in the \"general\" directory and has no path information in the log file,
            this script should be run there.

            Options:\n
              - baseDir: base of input figure directory. Defaults to $baseDir
              - outDir: output directory for files. Defaults to $outDir.
              - logFile: the TDR log file. Defaults to $logFile.
              - verbose: produce diagnostic output. Defaults to $verbose.
              - help: produce this message.\n\n
              
            Example:\n
              figstrip3.pl -logFile=\"..\\tmp\\ptdr2_temp.log\"
            \n\n";
              
            
    exit;
}

## see more documentation at the end



      my $figNum; # one of following categories: 0, 1.29, CP 2
      my $figNam; # filename less extension of included figure
      my $figPath; # filename with associated path

      my $lastFigPart = undef; #a, b, c, etc.
      my $lastFigNum = undef;
      my $lastFigPath = undef;
      my $lastNewFig = undef; # filename for new file
      my $lastExt = undef;  # .jpg or .pdf

#--------------------------------------------

      open (LOGFILE, $logFile) || die ("can't open the log file: $!");
      while (<LOGFILE> )
      {
        chomp;
# look for start of figure blocks
        if (/^<789FIG (.*)$/)
        {
#---------- pull out the name/number
            $figPath = undef;
            if (substr($_,-1,1) ne '>') # line wrapped
            {
              my $line = $_;
              $_ = <LOGFILE>;
              $_ = $line.$_;   #concatenate both lines
              # try again
            }
            if (!/^<789FIG (\S*)\s(.*)>$/) {die "Unsuccessfully tried to wrap 789 line:\n$_\n";}
            $figNam = $1;
            $figNum = $2;
#---------- now get the path to the file
            while (!$figPath && ($_ = <LOGFILE>))
            {
              if (/^<use (.*)/)
              {
                my $tmp = $1;
                if ($tmp =~ /(.*)>/)
                {
                   $figPath = $1;
                }
                else
                {
                   $_ = <LOGFILE>; # get next line
                   if (/(.*)>/)
                   {
                     $figPath = $tmp.$1;
                   }
                   else
                   {
                     die "Was looking for second line of path and got lost!\n";
                   }
                }
#                print "FIGPATH: ",$figPath,"\n";
              }
            }
            my $ext = substr($figPath,-4,4);
            if ($ext !~ '.jpg' && $ext !~ '.pdf')
            {
               die "Error extracting extension: ",$ext," found from ",$figPath,"\n";
            }
#---------  block processing complete. Have figure info. Now format it.
# Generate figure name
            my $newFig = "Figure_".&formatNumber($figNum);
            my $cmdString;
            my $outFile;
            if (!$lastFigPath) # first pass
            {
            #fall through
            }
            elsif ($figNum ne $lastFigNum)
            {
              if ($lastFigPart)
              {
                $outFile = "$outDir/$lastNewFig-$lastFigPart$lastExt";
                $lastFigPart = undef;
              }
              else
              {
                $outFile = "$outDir/$lastNewFig$lastExt";
              }
            }
            else
            {
              if (!$lastFigPart)
              {
                $lastFigPart = 'a';
              }
              $outFile = "$outDir/$lastNewFig-$lastFigPart$lastExt";
              $cmdString = "copy $lastFigPath $outFile";
              $lastFigPart = chr(ord($lastFigPart)+1);
            }
            if ($outFile)
            {
               if (!$execute)
               {
                 print "$baseDir/$lastFigPath --> $outFile\n";
               }
               else
               {
                 if ($verbose) {print "$baseDir$lastFigPath --> $outFile\n"};
                 copy ("$baseDir/$lastFigPath",$outFile) or die "Copy of $lastFigPath failed\n";
               }
            }
            $lastFigNum = $figNum;
            $lastNewFig = $newFig;
            $lastFigPath = $figPath;
            $lastExt = $ext;

         }
      }

sub formatNumber
{
#   format string  of form nn.mm as 0nn-0mm. Ignore strings with non-numeric data.
#
    my $arg = shift;
    $arg =~ s/\s/-/g; # remove any blanks
    if ($arg =~ /[^\d\.]/) 
    { 
      if ($arg =~ /^CP-(\d*)/)
      {
        my $out = sprintf(' CP-%03d',$1+1);
        return ($out);
      }
      elsif ($arg =~ /^([A-Z])\.(\d*)/)
      {
        my $out = sprintf('%s-%03d',$1,$2+1);
        return ($out);
      }
      else
      {    
        print "formatNumber: bailout: $arg\n";
        return ($arg);
      } 
    }#non-numeric, return as is
    $arg =~ /(\d*)\.*(\d*)/;
    my $num = $1;
    my $frac = $2;
    my $out;
    if ($frac eq '')
    {
      $out = sprintf('%03d',$num);
    }
    else
    {
      $out = sprintf('%03d',$num)."-".sprintf('%03d',$frac+1);
    }
    return ($out);
}

#---------------------------
#
#    figstrip.pl
#
#    Purpose: to create a directory of figures as used in a TDR, numbered as in the document
#
#    Requires: Modified cms-tdr.cls
#
# tracinggraphics must be turned on and the definition of \Ginclude@graphics must be modified with the following insertion:
#{0 0 .9}{#1}\typeout{<789FIG  #1 \thefigure>}\OldGinclude@graphics{#1}}
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
#
# Here's what we then look for in the logfile:
#----
#<789FIG introduction/BunchStructure-323_jpg 1.0>
#<introduction/BunchStructure-323_jpg.jpg, id=3354, 484.9317pt x 289.44135pt>
#File: introduction/BunchStructure-323_jpg.jpg Graphic file (type jpg)
#<use introduction/BunchStructure-323_jpg.jpg> [3] [4 <C:\Documents and Settings
#\George Alverson\My Documents\CMS\PTDR\TDR\ptdr1\tex\..\fig\introduction/BunchS
#tructure-323_jpg.jpg>]
#----
# whenever there is a warning, will see this form:
#----
#<789FIG ecal/Performance/Resolution_Comp_Hodo_3x3_704_tdr 1.6>
#
#Warning: pdflatex (file ecal/Performance/Resolution_Comp_Hodo_3x3_704_tdr.pdf):
# pdf inclusion: found pdf version <1.6>, but at most version <1.5> allowed
#
#<ecal/Performance/Resolution_Comp_Hodo_3x3_704_tdr.pdf, id=3716, 569.12625pt x
#532.99126pt>
#File: ecal/Performance/Resolution_Comp_Hodo_3x3_704_tdr.pdf Graphic file (type
#pdf)
#<use ecal/Performance/Resolution_Comp_Hodo_3x3_704_tdr.pdf> [16 <C:\Documents a
#nd Settings\George Alverson\My Documents\CMS\PTDR\TDR\ptdr1\tex\..\fig\ecal/Per
#formance/Resolution_Comp_Hodo_3x3_704_tdr.pdf>]
#----
# with multiple files in a figure (a),(b),(c), get:
#----
#<789FIG tracker/TrackerReco/muon_all_pt-col 1.10>
#<tracker/TrackerReco/muon_all_pt-col.pdf, id=3838, 569.12625pt x 546.04pt>
#File: tracker/TrackerReco/muon_all_pt-col.pdf Graphic file (type pdf)
#<use tracker/TrackerReco/muon_all_pt-col.pdf>
#<789FIG tracker/TrackerReco/muon_all_tip-col 1.10>
#<tracker/TrackerReco/muon_all_tip-col.pdf, id=3840, 569.12625pt x 546.04pt>
#File: tracker/TrackerReco/muon_all_tip-col.pdf Graphic file (type pdf)
#<use tracker/TrackerReco/muon_all_tip-col.pdf>
#<789FIG tracker/TrackerReco/muon_all_lip-col 1.10>
#<tracker/TrackerReco/muon_all_lip-col.pdf, id=3842, 569.12625pt x 546.04pt>
#File: tracker/TrackerReco/muon_all_lip-col.pdf Graphic file (type pdf)
#<use tracker/TrackerReco/muon_all_lip-col.pdf>
#---
# normal appendix
#---
#<789FIG appendix2/JetSystematicsFig B.1>
#---
# colour plates
#---
#<789FIG introduction/Hgg_m120.pdf CP 0>
#---
# Additional gotcha's:
# 1) first figure (logo on front) is figure 0, no .xx
# 2) figures for colour plates are "CP  0", etc, and are off by 1 (CP 0 is really CP 1).
# 3) a last single figure has to be dealt with by hand.
# 4) if there is are figures from appendices D, E, ... the Colour Plate figures will be out of order.
# 5) mixing cygwin and DOS shell may cause EOL confusion.
# 6) becasue the logo is in general and doesn't have full path information, this should be run from general
#-------------------------------------------------------------------------
