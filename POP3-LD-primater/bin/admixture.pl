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
&USAGE unless ($fIn and $fOut);
my @log=glob("$fIn/*.log");
my %cv;

foreach my $log( @log){
	my $k=(split(/\./,basename($log)))[1];
	open In,$log;
	while(<In>){
		chomp;
		next if($_ eq ""||/^$/);
		if(/CV error/){
			my $cv=(split(/\s+/,$_))[-1];
			$cv{$k}=$cv;
		}
	}
	close In;
}

my @order= sort {$cv{$a} <=> $cv{$b}} keys %cv;
my $best=$order[0];
my $num1=$best-1;
my $num2=$best+1;
`cp $fIn/pop.$best.result $fOut/pop.best$best.result`;
`cp $fIn/pop.$num1.result $fOut/pop.best-1.result` if(-f "$fIn/pop.$num1.result");
`cp $fIn/pop.$num2.result $fOut/pop.best+1.result` if(-f "$fIn/pop.$num2.result");
open Out,">$fOut/CV.list";
foreach my $k(sort {$a<=>$b} keys %cv){
	print Out $k,"\t",$cv{$k},"\n";
}
close Out;
`Rscript ${Bin}/admixture.R --best $fOut/pop.best$best.result --best1 $fOut/pop.best-1.result --best2 $fOut/pop.best+1.result --outfile $fOut/structure --bestk $best`;
`Rscript ${Bin}/cverror.R --infile CV.list --outfile $fOut/CV`;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
