#!/usr/bin/env perl

use strict;
use warnings;


use Getopt::Long;
use File::Copy;
use File::Basename;
use File::Spec;
use Cwd;
use Cwd qw(abs_path);
use IO::File;
use XML::Twig;
use Text::Balanced qw (extract_bracketed);
#use PDF::API2; #not installed on lxplus


#my $VERSION = sprintf "%d.%03d", q$Revision: 436283 $ =~ /(\d+)/g; # not with git
my $verbose;
my $texFile = 'D:/tdr2/papers/XXX-08-000/XXX-08-000.tex';
my $doc = 'XXX-08-000_temp.pdf';
my $style = 'paper';
my $baseDir = 'D:/tdr2/papers/XXX-08-000';
my $outDir = "D:/tdr2/utils/tmp";
my $logFile = "D:/tdr2/utils/tmp/XXX-08-000_temp.log";
my $execute = '1';
my $reportTitle = '1';
my $contactAddress;
my $help;
my $updateRecord;
my $dataNotMC = 1;
my $thumbnails = 1; # generate thumbnails of the figures
my $checkRefs = 1; # currently needs Python >= 2.6.4
my $artType = 3;
my $gitDate;
#Note is a special article type: goes into a 980_a
my @artTypes = ("Particle Physics - Experiment", "Nuclear Physics - Experiment", "Detectors and Experimental Techniques", "Note");
GetOptions ('verbose!' => \$verbose,
            'help|?' => \$help,
            'texFile=s' => \$texFile,
            'doc=s' => \$doc,
            'style=s' => \$style,
            'baseDir=s' => \$baseDir,
            'outDir=s' => \$outDir,
            'logFile=s' => \$logFile,
            'contactAddress=s' => \$contactAddress,
            'execute!' => \$execute,
            'updateRecord=i' => \$updateRecord,
            'thumbs!' => \$thumbnails,
            'data!' => \$dataNotMC,
            'checkRefs!' => \$checkRefs,
            'artType=i' => \$artType,
            'date=s' => \$gitDate
             );
#$verbose=1;
if ($verbose)
{
    print "$0 called with texFile = $texFile, doc= $doc, style = $style, baseDir = $baseDir, outDir = $outDir, logFile = $logFile, contact = $contactAddress, article type = $artType\n";
}
if ($help)
{
    print
    "\n >>> makeManifest: pulls references to figures from a TDR log file and copies the actual figures
            to a directory with names equal to the figure numbers.\n

            Note that the log file must be produced using a modified \\includegraphics command.

            Since the logo is in the \"general\" directory and has no path information in the log file,
            this script should be run there.

            Options/Parameters:\n
              - texFile: tex file with pdftitle, pdfauthor, etc. Defaults to $texFile.\n
              - doc: base document (in pdf). Defaults to $doc.\n
              - style: note style. Defaults to $style.\n
              - baseDir: base of input directory. Defaults to $baseDir.\n
              - outDir: output directory for files. Defaults to $outDir.\n
              - logFile: the TDR log file. Defaults to $logFile.\n
              - updateRecord: update the given record number at CDS. Defaults to $updateRecord.\n
              - contactAddress: e-mail contact address. Defaults to $contactAddress.\n
              - verbose: produce diagnostic output. Defaults to $verbose.\n
              - thumbs: create thumbnails for the figures. Requires ImageMagick and Ghostscript (can be done with Photoshop + script, but not implemented). Defaults to $thumbnails.\n
              - checkRefs: check the references against the standard CMS recommendations
              - artType: set the article description to (0) 'Particle Physics - Experiment', (1) 'Nuclear Physics - Experiment', or (2) 'Detectors and Experimental Techniques'. Defaults to $artType, which is generic note.\n
              - date: the date the document was closed
              - help: produce this message.\n\n

            Example:\n
              $0 -logFile=\"..\\tmp\\ptdr2_temp.log\"
            \n\n";


    exit;
}

## see more documentation at the end


      my $baseURL = $outDir;
      my $figNum; # one of following categories: 0, 1.29, CP 2
      my $figNam; # filename less extension of included figure
      my $figPath; # filename with associated path
      my $figFileList = dirname($logFile)."/FigList.txt"; # file to store included figure names and ordered file names

      my $lastFigPart = undef; #a, b, c, etc.
      my $lastFigNum = undef;
      my $lastFigPath = undef;
      my $lastNewFig = undef; # filename for new file
      my $lastExt = undef;  # .jpg or .pdf
      my $skipfirst = 1; # first includegraphics is _almost_ always the experiment logo
      my $figOffset = 1;

      my $newFig;
      my $cmdString;
      my $outFile;
      my $ext;

      my $preface =<<EOD; # for checking abstract LaTeX
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<!DOCTYPE html>
<html>
<head>
<title>MathJax TeX Test Page</title>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({tex2jax: {inlineMath: [['\$','\$'], ['\\\\(','\\\\)']]}});
</script>
<script type="text/javascript"
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script> </head> <body> <p>Below are the title and abstract after
passing through the same formatter [MathJaX] as is used on CDS. You
should be able to spot those instances where comments (%) or macros
have sneaked into the abstract. Unrecognized macros in math mode
will be in red, but those in text mode will not stand out.</p> <hr>
EOD

#
#--------------------------------------------
#     extract Author/Title info from the TeX file
#
   my $title = '';
   my $author = '';

   open(FILE, $texFile) || die("can't open TeX file $texFile: $!");
   $_ .= do { local( $/ ); <FILE> }; #grab entire content!
   close(FILE);
   my $gitYear = substr($gitDate,0,4);
      # extract metadata
   m/\\hypersetup(.*)/s; # find hypersetup
   my $xtract = $1;
   my $substring = extract_bracketed($xtract,'{}');
   # pull out everything after the pdftitle mark
   if ($substring =~ m/pdftitle\s*=\s*(.*)/s)
   {
       $xtract = $1;
       my $subsubstring = extract_bracketed($xtract,'{}');
       $title = substr($subsubstring,1,-1); #strip off the bracketing character
   }
   else
   {
      print "Cannot find title in metadata\n";
   }
   if ($substring =~ m/pdfauthor\s*=\s*\{(.*?)\}/s)
   { $author = $1}
   else
   {
       print "Cannot find author in metadata\n";
   }
   if ($reportTitle)
   {
       print "\n", ">>> --- Metadata report on $doc --- <<<";
       print "\n", ">>>";
       print "\n", ">>> PDF TITLE: ",$title;
       print "\n", ">>> PDF AUTH : ",$author,"\n\n";
   }
   # extract the LaTeX title
   my $ltitle;
   my $pos = index($_,'\title');
   if ($pos)
   {
       my $stext = substr($_,$pos+6); #skip \title
       $ltitle = extract_bracketed($stext,'{}');
       $ltitle = substr($ltitle,1,-1); # remove braces
       # now clean up  whitespace
       $ltitle =~ s/^\s*//; #begin
       $ltitle =~ s/\s*$//; #end
       $ltitle =~ s/\s{2,}/ /mg; #multiple interior...
       if ($reportTitle) { print ">>> LaTeX TITLE: ",$ltitle,"\n";}
   }
   else
   {
       print "$0: no LaTeX title found\nExiting...\n";
       return;
   }
   # extract the LaTeX abstract
   my $abstract;
   $pos = index($_,'\abstract');
   if ($pos)
   {
       my $stext = substr($_,$pos+9); #skip \abstract
       $abstract = extract_bracketed($stext,'{}');
       $abstract = substr($abstract,1,-1); # remove braces
       # now clean up  whitespace
       $abstract =~ s/^\s*//; #begin
       $abstract =~ s/\s*$//; #end
       $abstract =~ s/\s{2,}/ /mg; #multiple interior...
       if ($reportTitle)
       {
           print ">>> LaTeX ABSTRACT: ",$abstract,"\n";
           print ">>> ---- <<<\n";
           # write out html for MathJax test

           my $htmlFile = File::Spec->catfile(dirname($logFile),$doc);
           $htmlFile =~ s/\.pdf//; #remove the .pdf extension
           $htmlFile = $htmlFile.".html";
           open(FILE, "> $htmlFile") || die("can't open file $htmlFile $! for HTML test");
           print FILE $preface, "\n", $ltitle,'<br/><br/>',$abstract, "\n</body></html>";
           close(FILE);
           print "\n\n>>> View the file below in a web browser to check the format of the abstract on CDS. <<<\n>>> HTML Check file: ",$htmlFile,"\n\n\n";



       }
       if ($abstract =~ m/(?<!\\)\%/m) # look for % without preceeding \
       {
           print ">>> -------------------------------------------------------------------------- <<<\n";
           print ">>>!!!    The abstract appears to have an unescaped TeX comment character (%). <<<\n";
           print ">>>       No comments are allowed within the abstract!                         <<<\n";
           print ">>> -------------------------------------------------------------------------- <<<\n";
       }
   }
   else
   {
       print "$0: no abstract found\nExiting...\n";
       return;
   }

    # extract from PDF. This requires the PDF::API2 module, which is not
    # available on lxplus
    #    my $pdf = PDF::API2->open($doc);
    #    my  %h = $pdf->info;
    ##    print %h;
    #   if ($verbose)
    #   {
    #       print "\n", ">>> TITLE: ",$title,"\n";
    #       print "\n", ">>> AUTH : ",$author,"\n";
    #   }
    #
    #--------------------------------------------

# Generate the XML
      my $noteCode = $baseDir; # get canonical name from directory name
      my @parts = split m+/+, $noteCode;
      $noteCode = $parts[$#parts];
#     don't want to include the collection parent: from J-Y LM, 2008/07/18
#      my $collection = XML::Twig::Elt->new(collection=>{xmlns=>"http://www.loc.gov/MARC21/slim"});
      my $record = XML::Twig::Elt->new('record');
      # see http://www.loc.gov/marc/bibliographic/ for tag mappings, but FFT is CERN specific.
      # currently documentation is available at http://invenio-demo.cern.ch/help/admin/bibupload-admin-guide
      # tag 100 is only for changes: from J-Y LM, 2008/07/18
      #my $control = XML::Twig::Elt->new(controlfield=>{ tag=>"100", ind1=>" ", ind2=>" "});
      #$control->paste('last_child',$record);
      if ($updateRecord)
      {
          my $control = XML::Twig::Elt->new(controlfield=>{ tag=>"001", ind1=>" ", ind2=>" "},$updateRecord);
          $control->paste('last_child',$record);
      }
      # CMS note number (MARC Source of Acquisition)
      my $xml_style = $style; # remove redundant CMS
      if (lc($style) eq 'cmspaper') {$xml_style = 'paper';}
      &insert_datafield("037",$record,"CMS-".uc($xml_style)."-".$noteCode);

      # Title (MARC Title Statement)
      &insert_datafield("245",$record,$title);

      # Other standard MARC data as will be implemented by CDS (as per message 2013-03-05)
      # [See also https://cds.cern.ch/help/admin/howto-marc]
      my $marcCollab = XML::Twig::Elt->new(datafield=>{ tag=>"710", ind1=>" ", ind2=>" "});
      &insert_subfield("g",$marcCollab,"CMS Collaboration");
      $marcCollab->paste('last_child',$record);
      my $marcAuth = XML::Twig::Elt->new(datafield=>{ tag=>"110", ind1=>" ", ind2=>" "});
      &insert_subfield("a",$marcAuth,"CMS Collaboration");
      $marcAuth->paste('last_child',$record);
      # Accelerator field (as per Annette Holtkamp request of 2013-03-28
      my $marcAccel = XML::Twig::Elt->new(datafield=>{ tag=>"693", ind1=>" ", ind2=>" "});
      &insert_subfield("a",$marcAccel,"CERN LHC");
      &insert_subfield("e",$marcAccel,"CMS");
      $marcAccel->paste('last_child',$record);

      # Year (MARC Publication)
      my $pubyear = XML::Twig::Elt->new(datafield=>{ tag=>"260", ind1=>" ", ind2=>" "});
      &insert_subfield("c",$pubyear,$gitYear);
      $pubyear->paste('last_child',$record);

      # MARC/CERN pubdata:
      my $pubdat = XML::Twig::Elt->new(datafield=>{ tag=>"269", ind1=>" ", ind2=>" "});
      &insert_subfield("a",$pubdat,"Geneva");
      &insert_subfield("b",$pubdat,"CERN");
      # &insert_subfield("c",$pubdat,$gitYear); needs YYYYMMDD
      $pubdat->paste('last_child',$record);
      my $cernPub = XML::Twig::Elt->new(datafield=>{ tag=>"690", ind1=>"C", ind2=>" "});
      &insert_subfield("a",$cernPub,"CERN");
      $cernPub->paste('last_child',$record);

      # Data/MC indicator: tag 653, ind1:1, subfield 9: CMS, subfield a: Data/Monte-Carlo: from J-Y LM, 2010/02/10
      my $tagDataMC = XML::Twig::Elt->new(datafield=>{tag=>"653", ind1=>"1", ind2=>" "});
      &insert_subfield("9",$tagDataMC,"CMS");
      if ($dataNotMC)
      {
          &insert_subfield("a",$tagDataMC,"Data");
      }
      else
      {
          &insert_subfield("a",$tagDataMC,"Monte-Carlo");
      }
      $tagDataMC->paste('last_child',$record);

#      # CDS
#      &insert_datafield("859",$record,'cds.support@cern.ch');

      # 520 tags: (MARC Electronic Location and Access: ind1 = 3 =>abstract, subfield a=> Summary)
      # Abstract
      my $tag520 = XML::Twig::Elt->new(datafield=>{tag=>"520", ind1=>" ", ind2=>" "});
      &insert_subfield("a",$tag520,$abstract);
      $tag520->paste('last_child',$record);

      # 856 tags: (MARC Electronic Location and Access: ind1 = 4 =>HTTP, subfield y=> Link text)
      # CMS analysis page url
      &insert_856datafield($record,"http://cms.cern.ch/iCMS/analysisadmin/cadi?ancode=".$noteCode,"Additional information for the analysis (restricted access)");


      # Authorlist: need url to supply authorlist to replace "The CMS Collaboration"
      #&insert_856datafield($record,"The CMS Collaboration","CERN");

      # 980 tag: (MARC Equivalence or Cross-Reference Series Personal Name/Title)"
      if ($style eq 'pas')
      {
          &insert_datafield("980",$record,"NOTE");
          &insert_datafield("980",$record,"CMS-PHYSICS-ANALYSIS-SUMMARIES");
      }
      elsif ($style eq 'cmspaper')
      {
          &insert_datafield("980",$record,"CMS_Papers"); # should be ?CMSPUBDRAFTFINAL? ?CMSPUBDRAFT?
      }

      # 859:8560 tags: (email:  subfield f)
      my $tag8560 = XML::Twig::Elt->new(datafield=>{tag=>"859", ind1=>" ", ind2=>" "});
      &insert_subfield("f",$tag8560,$contactAddress);
      $tag8560->paste('last_child',$record);
      # 65017 tag: CERN and document type. Options are "Particle Physics - Experiment", "Nuclear Physics - Experiment", and "Detectors and Experimental Techniques"
      if ($artType == 3)
      {
          my $tag980 = XML::Twig::Elt->new(datafield=>{tag=>"980", ind1=>" ", ind2=>" "});
          &insert_subfield("a",$tag980,$artTypes[$artType]);
          $tag980->paste('last_child',$record);
      }
      else
      {
          my $tag650 = XML::Twig::Elt->new(datafield=>{tag=>"650", ind1=>"1", ind2=>"7"});
          &insert_subfield("2",$tag650,"SzGeCERN");
          &insert_subfield("a",$tag650,$artTypes[$artType]);
          $tag650->paste('last_child',$record);
      }

      # FFT
      my $fft_tag = XML::Twig::Elt->new(datafield=>{tag=>"FFT", ind1=>" ", ind2=>" "});
      my $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"a"},$baseURL."/".$doc);
      $doc_tag->paste('last_child',$fft_tag);
      $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"t"},"Main");
      $doc_tag->paste('last_child',$fft_tag);
      $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"d"},"Fulltext in PDF");
      $doc_tag->paste('last_child',$fft_tag);
      $fft_tag->paste('last_child',$record);

      my $start = time;
      # add each figure
      my $convertCmd = "";
      if ($thumbnails)
      {
          if ($^O eq "MSWin32")
          {
              if (-d $ENV{'ProgramFiles'}."\\ImageMagick-7.0.7-Q16")
              {
                $ENV{'PATH'} = $ENV{'ProgramFiles'}."\\ImageMagick-7.0.7-Q16;".$ENV{'ProgramFiles'}."\\gs\\gs9.22\\bin;$ENV{'PATH'}";
                $convertCmd = "magick.exe";
              }
              elsif (-d $ENV{'ProgramFiles'}."\\ImageMagick-6.9.2-Q16")
              {
                $ENV{'PATH'} = $ENV{'ProgramFiles'}."\\ImageMagick-7.0.7-Q16;".$ENV{'ProgramFiles'}."\\gs\\gs9.18\\bin;$ENV{'PATH'}";
                $convertCmd = "convert.exe";
              }
              else {print('***No image conversion program found***');}
          }
          else
          {
              $convertCmd = "convert";
              if (-f '/afs/cern.ch/cms/external/gs') # temporary patch on lxplus for problems with central gs. See http://bugs.ghostscript.com/show_bug.cgi?id=690676.
              {
                  $ENV{'PATH'} = '/afs/cern.ch/cms/external'.':'.$ENV{'PATH'};
              }
          }
          if ($verbose) {print "> Using convert version...\n"; system($convertCmd,"-version"); print "\n";}
      }
      if ($verbose) {print "> Opening $logFile to look for figures\n";}
      open (LOGFILE, $logFile) || die ("can't open the log file $logFile: $!");
      if ($verbose) {print ("opening $figFileList for writing\n");}
      open (FIGLIST, ">", $figFileList) || die ("can't open file for figure list, $figFileList: $!");
local $/ = "\012"; # necessary when called, for some unknown reason
      while (<LOGFILE> )
      {
        chomp;
# look for start of figure blocks
        if (/^<789FIG (.*)$/)
        {
#---------- pull out the name/number
            if ($verbose) {print "> ", $_ ,"\n"};
            $figPath = undef;
            my $line = $_;
            $line =~ s/\n+$//; # now occasionally has trailing white space
            while (substr($line,-1,1) ne '>')
            {
              # line wrapped
              $_ = <LOGFILE>;
              $line = $line.$_;   #concatenate both lines
              $line =~ s/\n+$//;
            }
            if (not $line =~ /^<789FIG (\S*)\s(.*)>$/) {die "Unsuccessfully tried to wrap 789 line:\n$_\n";}
            $figNam = $1;
            $figNum = $2;
            if ($figNum =~ /(\d+)([[:alpha:]])$/) # subfloat pre-increments the figNum, sometimes!
            {
                if ($1 == 0) {$figNum = 0;}
                else {$figNum = $1-1;}
             }

#---------- now get the path to the file
            while (!$figPath && ($_ = <LOGFILE>))
            {
              if (/^\s*<use (.*)/)
              {
                $line = $1;
                if (not $line =~ /(.*)>/)
                {
                  while (not $line =~ /(.*)>/)
                  {
                      $_ = <LOGFILE>; # get next line
                      $line = $line.$_;
                      $line =~ s/\n+$//; # now occasionally has trailing white space
                  }
                  $line =~ /(.*)>/;
                  $figPath = $1;
                }
                else
                {
                  $figPath = $1;
                }
                if ($verbose) {print "> FIGPATH: ",$figPath,"\n";}
              }
            }
            $ext = substr($figPath,-4,4);
            if ($ext ne '.jpg' && $ext ne '.pdf' && $ext ne '.png' && $ext ne '.jpeg')
            {
               die "Error extracting extension: ",$ext," found from ",$figPath,"\n";
            }
#---------  block processing complete. Have figure info. Now format it.
            if ( $skipfirst or $figPath =~ /BigDraft/ )
            {
                $skipfirst = undef;
            }
            else
            {
                my ($outFile2, $startingFig) = &genFigName;
                if ($outFile2)
                {
                    print FIGLIST ("$startingFig $outFile2\n");
                    if ($verbose) {print ">...Processing $outFile2\n";}
                    $fft_tag = XML::Twig::Elt->new(datafield=>{tag=>"FFT", ind1=>" ", ind2=>" "});
                    $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"a"},$outDir."/".$outFile2); # need to put back URL if go to server from afs
                    $doc_tag->paste('last_child',$fft_tag);
                    $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"t"},"Additional");
                    $doc_tag->paste('last_child',$fft_tag);
                    $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"z"},"Figure");
                    $doc_tag->paste('last_child',$fft_tag);
                    $outFile2 =~ m/^(\S+)\.\S{3,4}$/s;
                    $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"d"},$1);
                    $doc_tag->paste('last_child',$fft_tag);
                    my @convertArgs = ("$outDir/$outFile2", "-trim", "-thumbnail", "240", "-define", "pdf:use-cropbox=true",  "$outDir/$1-thumb.png");
                    if ($thumbnails)
                    {
                      if ($verbose) {print("@convertArgs");}
                      system($convertCmd, @convertArgs)==0 or die "system call to create thumbnails with args @convertArgs failed: $?";
                      $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"x"},"$outDir/$1-thumb.png");
                      $doc_tag->paste('last_child',$fft_tag);
                    }
                    if (0)
                    {
                      @convertArgs = ("$outDir/$outFile2", "-trim", "-density", "600", "-quality", "100", "-thumbnail", "240", "-define", "pdf:use-cropbox=true",  "$outDir/$1-thumb.png");
                      system($convertCmd, @convertArgs)==0 or die "system call to create thumbnails with args @convertArgs failed: $?";
                      $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"x"},"$outDir/$1-thumb.png");
                      $doc_tag->paste('last_child',$fft_tag);
                    }
                    $fft_tag->paste('last_child',$record);
                }
            }
      }
      }
      close(LOGFILE);
      # kick out last figure
      $figNum = $figNum + 1; # otherwise sets figPart
      my ($outFile3, $startingFig) = &genFigName;
      print FIGLIST ("$startingFig $outFile3\n");
      if ($outFile3)
      {
          $fft_tag = XML::Twig::Elt->new(datafield=>{tag=>"FFT", ind1=>" ", ind2=>" "});
          $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"a"},$outDir."/".$outFile3); #ditto here for URL
          $doc_tag->paste('last_child',$fft_tag);
          $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"t"},"Additional");
          $doc_tag->paste('last_child',$fft_tag);
          $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"z"},"Figure");
          $doc_tag-> paste('last_child',$fft_tag);
          $outFile3 =~ m/^(\S+)\.\S{3,4}$/s;
          $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"d"},$1);
          $doc_tag->paste('last_child',$fft_tag);
          my @convertArgs = ("$outDir/$outFile3", "-trim", "-thumbnail", "240", "-define", "pdf:use-cropbox=true", "$outDir/$1-thumb.png");
          if ($thumbnails)
          {
              system($convertCmd, @convertArgs)==0 or die "system call to create thumbnails with args @convertArgs failed: $?";
              $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"x"},"$outDir/$1-thumb.png");
              $doc_tag->paste('last_child',$fft_tag);
          }
          if (0)
          {
            @convertArgs = ("$outDir/$outFile3", "-trim", "-density", "600", "-quality", "100", "-thumbnail", "240", "-define", "pdf:use-cropbox=true",  "$outDir/$1-thumb.png");
            system($convertCmd, @convertArgs)==0 or die "system call to create thumbnails with args @convertArgs failed: $?";
            $doc_tag = XML::Twig::Elt->new(subfield=>{code=>"x"},"$outDir/$1-thumb.png");
            $doc_tag->paste('last_child',$fft_tag);
          }
          $fft_tag->paste('last_child',$record);
      }
#      $record->paste('last_child',$collection);  #: from J-Y LM, 2008/07/18
#debug      print("Thumbnail generation time: ",time - $start, "\n");

      my $MANIFEST = IO::File->new("> ".$outDir."/manifest.xml") or die $!;
      if ($verbose) {print "> Going to print out XML Manifest to $outDir/manifest.xml\n";}
#      open ($MANIFEST, "> ".$outDir."/manifest.xml") or die "Can't open manifest file\n";
#     print the <?xml...> line to make good xml; for xml fragment, leave it out
#      print $MANIFEST qq|<?xml version="1.0" encoding="UTF-8"?>\n|;
      $record->print($MANIFEST);
      close $MANIFEST;
      if ($verbose) {print "> Closed manifest file: $outDir/manifest.xml\n";}
      close(FIGLIST);

sub genFigName
{
    # Generate figure name
    
    my $outFile = undef;
    $newFig = "Figure_".&formatNumber($figNum+$figOffset);
    if (!$lastFigPath) # first pass
    {
    #fall through
    }
    elsif ($figNum ne $lastFigNum)
    {
      if ($lastFigPart)
      {
        $outFile = "$lastNewFig-$lastFigPart$lastExt";
        $lastFigPart = undef;
      }
      else
      {
        $outFile = "$lastNewFig$lastExt";
      }
    }
    else
    {
      if (!$lastFigPart)
      {
        $lastFigPart = 'a';
      }
      $outFile = "$lastNewFig-$lastFigPart$lastExt";
      $cmdString = "copy $lastFigPath $outDir/$outFile";
      $lastFigPart = chr(ord($lastFigPart)+1);
    }
    if ($outFile)
    {
       if (!$execute)
       {
         print ">> $baseDir/$lastFigPath --> $outDir/$outFile\n";
       }
       else
       {
         if ($verbose) {print ">> $baseDir/$lastFigPath --> $outDir/$outFile\n"};
         if (!copy ("$baseDir/$lastFigPath",$outDir."/".$outFile))
         {
             print ">> Copy of >>>> $lastFigPath failed\n";
             $outFile = undef;
         }
       }
    }
    my $oldFigPath = $lastFigPath;
    $lastFigNum = $figNum;
    $lastNewFig = $newFig;
    $lastFigPath = $figPath;
    $lastExt = $ext;
    return ($outFile, $oldFigPath);

}
sub formatNumber
{
    #
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
#############################################################################
sub insert_subfield {
    #
    # arguments:
    # code (ex: "a" or )
    # field: xml twig to insert into (as last child)
    # value (ex: the title)
    #
    
    my $code = shift;
    my $field = shift;
    my $value = shift;
    my $subid = XML::Twig::Elt->new(subfield=>{code=>$code},$value);
    $subid->paste('last_child',$field);
}
#############################################################################
sub insert_856datafield {
    # 856 records are urls: http://cdsdev.cern.ch/help/admin/bibupload-admin-guide
    
    my $report = shift;
    my $url = shift;
    my $descrip = shift;
    my $id = XML::Twig::Elt->new(datafield=>{tag=>"856", ind1=>"4", ind2=>" "});
    my $subid = XML::Twig::Elt->new(subfield=>{code=>"u"},$url);
    $subid->paste('last_child',$id);
    $subid = XML::Twig::Elt->new(subfield=>{code=>"y"},$descrip);
    $subid->paste('last_child',$id);
    $id->paste('last_child',$record);
}
#############################################################################
sub insert_datafield {
    #
    # arguments:
    # tag (ex: "037")
    # record xml twig to insert into (as last child)
    # value (ex: the title)
    #
    my $tag = shift;
    my $record = shift;
    my $value = shift;
    my $id = XML::Twig::Elt->new(datafield=>{tag=>$tag, ind1=>" ", ind2=>" "});
    my $subid = XML::Twig::Elt->new(subfield=>{code=>"a"},$value);
    $subid->paste('last_child',$id);
    $id->paste('last_child',$record);
    }
    #---------------------------
    #
    #    (was) figstrip.pl
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
    # x 1) first figure (logo on front) is figure 0, no .xx
    # 2) figures for colour plates are "CP  0", etc, and are off by 1 (CP 0 is really CP 1).
    # x 3) a last single figure has to be dealt with by hand.
    # 4) if there is are figures from appendices D, E, ... the Colour Plate figures will be out of order.
    # 5) mixing cygwin and DOS shell may cause EOL confusion.
    # x 6) becasue the logo is in general and doesn't have full path information, this should be run from general
    #-------------------------------------------------------------------------


