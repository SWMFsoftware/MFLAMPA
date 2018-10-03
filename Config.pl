#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

use strict;

our $Component = "SP";
our $Code = "MFLAMPA";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments= @ARGV;


my $config     = "share/Scripts/Config.pl";
# get util and share
my $GITCLONE = "git clone"; my $GITDIR = "herot:/GIT/FRAMEWORK/";

if (-f $config or -f "../../$config"){
}else{
    `$GITCLONE $GITDIR/share.git; $GITCLONE $GITDIR/util.git`;
}



if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# Variables inherited from share/Scripts/Config.pl
our %Remaining; # Unprocessed arguments
our $ERROR;
our $WARNING;
our $Help;
our $Verbose;
our $Show;
our $ShowGridSize;
our $NewGridSize;
our $NewGhostCell;


&print_help if $Help;

my $Src = 'src';

# Grid size variables
my $NameSizeFile = "$Src/ModSize.f90";
my $GridSize;
my $nP;

# Read previous grid size
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};
    if(/^-g$/)                {$ShowGridSize=1;                next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


# Set new grid size
&set_grid_size if ($NewGridSize and $NewGridSize ne $GridSize);

# Show current settings
my $Settings = &current_settings; print $Settings if $Show;

# Show grid
print "$GridSize\n" if $ShowGridSize;

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/;
	$nP   =$1 if /\bnParticleMax\s*=[^0-9]*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read nParticleMax from $NameSizeFile\n" 
	unless length($nP);                         

    $GridSize = "$nP";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

    if($GridSize =~ /^[1-9]\d*$/){
	$nP = $GridSize;
    }elsif($GridSize){
	die "$ERROR -g=$GridSize must be only one integer\n";
    }
    # Check the grid size (to be set)
    die "$ERROR nParticleMax=$nP must be positive\n" if $nP<=0;

    print "Writing new grid size $GridSize into ".
	"$NameSizeFile ...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nParticleMax\s*=[^0-9]*)(\d+)/$1$nP/i;
	print;
    }
}

#############################################################################

sub print_help{

    print "
Additional options for MFLAMPA/Config.pl:

-g=nP,nLon,nLat 
                Set grid size. 
                nP is maximum number of particles per field line,
                nLon, nLat are the grid size at the origin surface.
\n";
    exit 0;
}

#############################################################################

sub current_settings{

    $Settings .= 
	"Number of particles per line   : nParticleMax=$nP\n";
}

