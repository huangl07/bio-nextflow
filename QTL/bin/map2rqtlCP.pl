#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($loc,$map,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"l:s"=>\$loc,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($loc and $fOut);

open In,$loc;
open Out,">$fOut";
my @out;
my $nind;
my $nloc;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	my ($markerID,$type,$phase,@indi)=split(/\s+/,$_);
    $nloc++;
    $nind=scalar @indi;
    push @out,$_;
}
close In;
print Out "popt=CP\nnind=$nind\nnloc=$nloc\nname=pop\n";
print Out join("\n",@out);
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"l:s"=>\$loc,
	"m:s"=>\$map,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
