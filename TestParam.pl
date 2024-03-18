#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Verbose     = $v; undef $v;
my $Help        = $h; undef $h;
my $HelpXmlParam= $H; undef $H;
my $HelpXml     = $X; undef $X;
my $Format      = $F; undef $F;

use strict;

# This script uses the CheckParam.pl script 
my $CheckParamScript  = "share/Scripts/CheckParam.pl";
$CheckParamScript = "../../$CheckParamScript" unless -x $CheckParamScript;
die "Could not find $CheckParamScript" unless  -x $CheckParamScript;

# The -H and -X flags are transferred to CheckParam.pl
exec("$CheckParamScript -X") if $HelpXml;
exec("$CheckParamScript -H") if $HelpXmlParam;

if($Help){
    print "
Purpose:

    Check the correctness of the input parameter file using the XML description
    in PARAM.XML. The grid size returned by Config.pl -g, and the PRECISION 
    defined in Makefile.conf are also used during the check, if available.
    This simple script calls the general script in share/Scripts/CheckParam.pl.

Usage:

    TestParam.pl [-h] [-H] [-X] [-v] [PARAMFILE]

  -h            print help message and stop

  -H            print help about the XML tags used in PARAM.XML files and stop

  -X            print a short introduction to the XML language and stop

  -v            print verbose information

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'

Examples:

  Check the default parameter file run/PARAM.in:

      TestParam.pl

  Check another parameter file:

      TestParam.pl run/test.000/PARAM.expand

";
    exit 0;
}

my $XmlFile           = 'PARAM.XML';
my $ConfigScript      = 'Config.pl';
my $NameComp = 'SP';
my $Precision;       # passed to CheckParam.pl
my $SettingsPassed;  # passed to CheckParam.pl

if(-x $ConfigScript){
    my $Settings = `./Config.pl`;

    die "ERROR in TestParam.pl: $ConfigScript did not provide the settings\n" unless $Settings;

    foreach (split(/\n/,$Settings)){
	$Precision = $1 if /(single|double) precision/;
	# Extract : NAME=VALUE, NAME=VALUE
	if(s/.*://){
	    $SettingsPassed .= "$1:$2," while s/(\w+)\s*=\s*(\d+)//;
	}
    }
    chop $SettingsPassed;
}else{
    warn "WARNING could not execute $ConfigScript\n";
}

my @command = (
    $CheckParamScript,
    "-S",
    "-F=$Format",
    "-c=$NameComp",
    "-v=$Verbose",
    "-p=$Precision",
    "-s=$SettingsPassed",
    "-x=$XmlFile",
    @ARGV);

print "@command\n" if $Verbose;
system(@command);
