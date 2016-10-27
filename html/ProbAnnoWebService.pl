#!/bin/env perl
$|=1;
use strict;
use warnings;

use HTTP::Request::Common;
use LWP::UserAgent;

MAIN: {
    my $ua=LWP::UserAgent->new;
    my $uri='http://probanno.systemsbiology.net/cgi-pub/ProbAnno';

    # using separate chromosome and coordinate fields
    # get variant information for four sample positions on chromosome 1  for hg38
    my $uploaded_file = "/local/local_webservices/probanno/ProbAnno-Standalone/genomes/test.fasta";
    my $template = "GramNegative";
    my $taxid = "1051631";

    my $req_args=[
                        template => $template,
			#uploaded_file => $uploaded_file,   # not working properly
			taxid => $taxid,
                 ];

    print "Retrieving information for taxonomy ID $taxid, template $template using param fields\n";
    #print "Retrieving information for file $uploaded_file, template $template\n";

    #sleep inserted for readability of print statement
    #do not use sleep statement in normal web service queries
    sleep 1;  


    my $res=$ua->request(POST $uri, $req_args);
    print $res->content;
 
    # using list field
    my $list = "template:GramNegative\ntaxid:224308";

    my $list_args = [
                        list => $list,
                    ];

    print "Retrieving information for taxonomy ID $taxid, template $template using list field\n";
    sleep 1;  
    my $list_res=$ua->request(POST $uri, $req_args);
    print $list_res->content;


}
