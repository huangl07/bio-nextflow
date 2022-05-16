#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin qw ($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($fIn,$fOut);
GetOptions(
    "help|?"=>\&USAGE,
    "input:s"=>\$fIn,
    "output:s"=>\$fOut,
) or &USAGE;
&USAGE unless ($fIn and $fOut);
my %type=(
    'lmxll'=>'A.1',
    'nnxnp'=>'B3.7',
    'hkxhk'=>'D2.15',
    'abxcd'=>'D1.10',
    'efxeg'=>'A.2',
    'ccxab'=>'D2.14',
    'abxcc'=>'D1.9'
);
my %alle=(
    'ee'=>'a',
    'ef'=>'ab',
    'eg'=>'ac',
    'fg'=>'bc',
    'lm'=>'ab',
    'll'=>'a',
    'nn'=>'a',
    'np'=>'ab',
    'aa'=>'a',
    'bb'=>'b',
    '--'=>'-',
    'hk'=>'ab',
    'hh'=>'h',
    'kk'=>'k',
    'kh'=>'ab',
    'pn'=>'ab',
    'ml'=>'ab'
);
open IN,$fIn;
open Out,">$fOut";
my $indi;
my $out;
my @marker;
while(<IN>){
    chomp;
    next if ($_ eq "" || /^$/ ||/^#/);
        my ($markerID,$type,@geno)=split(/\s+/,$_);
        $type =~ s/<|>//g;
        shift @geno if($geno[0] =~ /{/);
        $indi=scalar @geno;
        my $ntype=$type{$type};
        my @out;
        for(my $i=0;$i<@geno;$i++){
            if(!exists $alle{$geno[$i]}){
                print $geno[$i],"\n";
                die;
            }
            push @out,$alle{$geno[$i]};
        }
        push @marker,join("\t",$markerID,$ntype,@out);
}
my $nloc=scalar @marker;
print Out join(" ",$indi,$nloc,0),"\n";
print Out join("\n",@marker);
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -input	<file>	input file name
  -output	<file>	input keys of file name

  -h         Help

USAGE
        print $usage;
        exit;
}
