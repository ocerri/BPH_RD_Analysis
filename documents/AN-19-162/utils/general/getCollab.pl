#!/usr/bin/env perl

#  getCollab
# shim to replace old myCollab with standalone routine to hide the SSL routines,
# which are not present on many machines. Now, if not present, just no authorlist.
sub myCollabF 
{

  # Adds the appendix with the collaboration list
  # default is append (see second arg in the open command below)
  # args are
  #  1. the tex filename
  #  2. the type of collaboration list desired: tex or xml
  #  3. the tex directory to copy from

    use strict;
    use warnings;
    use LWP::UserAgent qw( SSL_VERIFY_NONE SSL_CA_FILE ); #for URL handling
    use IO::Socket::SSL;
    use File::Copy;
    use File::Basename;
    use File::Spec qw( catfile );
    use Encode;


    my (@args) = @ARGV;
    $_ =$args[0];
    /^.*[\\\/](.*)_temp\.tex/; #extract the tag
    my $tag = $1;
    my $type = $args[1];
    my $tex_dir = $args[2];
    my $body = '';
    my $local_dir = dirname(__FILE__);
    my $certs = File::Spec->catfile($local_dir,'CERN-bundle.pem') ;

    if ($tag =~ /CFT-09-.*/) {$tag = 'CRAFT-V2-09'}; # all CFT-09's have a common authorlist??
    # name on server; xml file is now authorlistN

    my $tagfile;
    if ($type eq 'xml')
      {$tagfile = $tag.'-authorlistN.'.$type;}
    else
      {$tagfile = $tag.'-authorlist.'.$type;}


    my $tagfileout = $tag.'-authorlist.'.$type; # our name; could just use 'authorlist.'.$type

    if (-e $tex_dir.'/'.$tagfile) {
      copy($tex_dir.'/'.$tagfile,$tagfileout); 
      print ">>>  Used local copy of author list file ", $tagfileout,"\n";
      }
    else {
      #my $url = 'https://cms-secr.web.cern.ch/cms-secr/Documents/authorListCMS/'.$tagfile;
      my $url = 'https://cmsdoc.cern.ch/cms/cmssecr/www/Documents/authorListCMS/'.$tagfile;
      #  For new format AL generation
      my $Uri = "https://icms.cern.ch/tools/api/alFile";
      my $ua;
      if (-e $certs)
      {
        $ua =  LWP::UserAgent->new(ssl_opts=>{verify_hostname=>0, SSL_CA_FILE=>$certs});
      }
      else 
      {
        $ua =  LWP::UserAgent->new(ssl_opts=>{verify_hostname=>0, verify_hostname=>SSL_VERIFY_NONE});
      }
      my $form = qq({"fName":"$tagfile"}); # plain POST with form data doesn't work; have to hand-code it this way
      my $response = $ua->post($Uri, 'Content-type'=>'application/x-www-form-urlencoded', 'Content'=>$form);
      # Try new format first
      if (!$response->is_success) 
      {
        print "::: Problem getting author list: ",$tagfile, " :::\n";
        print("    Status: ",$response->status_line);
        print("    Trying old author list site...\n");
        # Try old format
        $response = $ua->get($url);
        if (!$response->is_success) 
        {
          print "::: $url not found :::\n";
          print("    Status: ",$response->status_line);
          print "\n::::::\n";
          $tagfileout = 0; #should die?
        }
        else 
        {
          $body = $response->decoded_content;
          # this rewrites to avoid some bad TeX spacing
          if ($type eq 'tex') {$body =~ s/(\w\,*)~(\w)/$1 $2/g; $body =~ s/,/, /g; } #allow normal spacing and line breaking
        }
      }
      else 
      {# new format author list; may be in UTF-8... or not
        $body = $response->content;
        my $contentType = $response->header('Content-Type');
        if ($response->header('Content-Encoding') ) 
        {
          my $contentCharset = $response->header('Content-Encoding');
          if ($contentCharset =~ /^utf-8|UTF-8|utf8|UTF8/) 
          {
            $body = Encode::decode('UTF-8',$body); # it's already in UTF-8, but this sets the flag.
            #$body = $newbody;
            ## diagnostic
            #$isUTF = utf8::is_utf8($body);
            #printf("Is UTF?: %0X\n",$isUTF);
          }
        }
      }
      # save the authorlist in append mode
      if (($type eq 'xml') || ($type eq 'revtex')) 
      {
        open(OUTF, ">>:encoding(UTF-8)", $tagfileout) || die("Can't open outputfile: $!");
      }
      else 
      {
        open(OUTF, ">>", $tagfileout) || die("Can't open outputfile: $!");
      }
      print OUTF $body;
      close(OUTF);
    }
  return $tagfileout;
}
&myCollabF
