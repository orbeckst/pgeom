#!/usr/bin/perl -w
# $Id: porecharger.pl,v 1.2 2005/11/15 16:09:27 oliver Exp $

use strict;
use Messages; 

use Getopt::Long;
use Data::Dumper;

my $Version_ID = q($Id: porecharger.pl,v 1.2 2005/11/15 16:09:27 oliver Exp $ );

# imported from Messages
# DEBUG level says how chatty we are:
#  <0  dead silent
#   0  no diagnostic messages, only output
#   1  a few additional important informations
#   3  verbose
#   5  extra verbose
#  10  start with real debugging stuff
#  15  many intermediate results
#  16  print input files line-by-line
#  20  everything   
$DEBUG = 1;

my $DataFn;
my $OutputFn = "chl_Q.itp";
my ($opt_verbose, $opt_version, $opt_help);

my $CHARGE = 0.380; # backbone carboxyl/amide
my $ListFn = "/dev/stdin";
my %name = ("plus" => "DC", "minus" => "DO");
my (@groups,@itp);   
# groups: arrays of charge groups ([+,-,+,-],[+,-],...)
# @itp: whole itp file


# how to indicate a comment line in the data file (regex)
my $commentchar = '[;#]';


sub print_usage {
    my $Prog = Progname;
    print <<EOF;
usage: $Prog [OPTIONS] FILE

Read chl.itp file FILE and replace selected atoms with positively and
negatively charged pseudo atoms (list is read from the charge list
file or stdin if non is given).

The OPTIONS can be abbreviated; defaults are given in [].	 
--help              help
--verbose           (alias for DEBUGLEVEL=3)
--version           print version id
--debug=DEBUGLEVEL  0 (almost) quiet, <0 silent,
                    1..3 verbose, 
                    >5 really debugging stuff (20 max) [$DEBUG]

--output=FILE       modified itp [$OutputFn]
--chargelist=FILE   atom numbers to be changed (format below)
--charge=FLOAT      absolute value on the charged atoms [$CHARGE]
                    (default is the partial charge of an carboxyl/amide 
                    in the ffgmx forcefield)

chargelist format
=================
# comment (also initial ';' allowed)
# each line is a separate chargegroup and must be neutral
# charges alternate
#    +         -         +         -
atom1_pos atom2_neg
atom3_pos atom4_neg atom5_pos atom6_neg

EOF
    exit;
}




######################################
#
# MAIN
#
######################################


# options
&GetOptions( "output=s" => \$OutputFn,
             "charge=f" => \$CHARGE,
             "chargelist=s" => \$ListFn,
	     "debug=i" => \$DEBUG,
	     "verbose" => \$opt_verbose,
	     "version" => \$opt_version,
             "help" => \$opt_help) or &print_usage;
if ($opt_help) {&print_usage};
if ($opt_version) { msg (0, "%s\n", $Version_ID); exit 0 };
if ($opt_verbose) { $DEBUG = 3 };


if ($ARGV[0]) {
    $DataFn = $ARGV[0];
    -e $DataFn or 
        die "The  file \`".$DataFn."' to read from does not exist, stopped ";
};

# read parameters
open LIST, "< $ListFn" or die "Need atom numbers from \`".$ListFn."', stopped";
GROUP: while (<LIST>) {
    chomp;
    next GROUP if /^ *$commentchar/ || /^ *$/;
    push @groups, [ split ];
}
close(LIST);

msg(15, "Dump of the parsed chargelist file `$ListFn':\n%s\n",
    Dumper(\@groups));

# read whole itp 
#[ atoms ]
#;   nr     type    resnr    res atom     cgnr charge         mass
#     1      CH4        1    MTH CH4         1  0.000

@itp = <>;

#print Dumper(\@itp);
my ($g,$a1,$i,$q);
my ($nr, $atype, $resnr, $res, $atom, $cgnr, $charge, $mass);
my ($newatom, $newcgnr);
my ($ngrp, $qTot)  = (0,0);

foreach $g (@groups) {
    # determine number of charge group, cgnr (the way that chl.itp is
    # made, each atom has its own cgnr. We simply take the first one.  
    $a1 = $g->[0];    # first atom in the group
    $ngrp++;
    foreach (@itp) {
        next if (/moleculetype/ .. /atoms/);  # skip fwd to [ atoms ] section
        next if (/^ *;/ || /^$/);
        if  (/^ *$a1 +[a-zA-Z]+/) {
            msg (15, "Finding cgnr for atom $a1: ".$_) ;
            ($nr, $atype, $resnr, $res, $atom, $cgnr, $charge, $mass) = split;
            msg (14, "atom[$nr] type=$atype cgnr=$cgnr\n");
            $newcgnr = $cgnr;
            last;
        };
    };
    # found cgnr, but no sanity checks...
    # now do the replacements
    for ($i=0; $i < @{$g}; $i++) {
        $q = ($i % 2 == 0) ? $CHARGE : -$CHARGE;  # alternating + - + - ...
        $newatom = (($q > 0) ? $name{"plus"} : $name{"minus"}) . $ngrp;
        foreach (@itp) {
            # now search for nr and replace charge -> q, cgnr, atom
            next if (/^ *;/);
            next unless (/^\[ atoms \]/ .. /^ *$/); 
            next if (/^\[ atoms \]/ || /^ *$/);
            ($nr, $atype, $resnr, $res, $atom, $cgnr, $charge, $mass) = split;
            # msg(15, "replace: a1=$a1 i=$i: $nr $resnr $res $atom $cgnr\n");
            if ($nr == $g->[$i]) {
                $atom = $newatom;
                $cgnr = $newcgnr;
                $charge = $q;
                $qTot += $q;
                msg(15, "replace: a1=$a1 i=$i: $nr $resnr $res $atom $cgnr\n");
            }
            # modify array in place
            # (printf string from pgeom::mol_data.c itp_w_atom() )
            $_ = sprintf "%6d %8s %8d %6s %-6s %6d %6.3f\n", 
                     $nr, $atype, $resnr, $res, $atom, $cgnr, $charge;
        }
    }
};

if ($qTot != 0) {warning("Total charge is NOT ZERO (qTot=$qTot)!\n")};

# now output the final itp
open OUTPUT, "> $OutputFn" or die "Error opening `$OutputFn', stopped";
print OUTPUT "; charges q=+/-$CHARGE added to original itp $DataFn\n\n";
foreach (@itp) {
    print OUTPUT $_;
}
close(OUTPUT);

exit(0);	
