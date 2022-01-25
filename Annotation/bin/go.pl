#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,

			) or &USAGE;
&USAGE unless ($fIn and $fOut );
open In,$fIn;
open Out,">$fOut";
while (<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my @a=split;
    for(my $i=1;$i<@a;$i++){
        $a[$i]=~s/\s+//g;
        $a[$i]=~s/;//g;
        print Out $a[0],"\t",$a[$i],"\n";
    }
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	this Script will split fasta file into n files

	eg:
	perl $Script -i demo.fa -d Nr -o ./ -k demo

Usage:
  Options:
  -i	<file>	input fa file name
  -o	<dir>	output dir 
 
  -h         Help

USAGE
        print $usage;
        exit;
}
