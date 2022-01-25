#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$output,$gff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"output:s"=>\$output
			) or &USAGE;
&USAGE unless ($ref and $output);

open In,$ref;
$/=">";
my %seq;
my $level="scaffold";
my %length;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	if ($id =~ /Chr/||$id =~ /chr/) {
		$level="chromosome";
	}
	my $line=join("",@seq);
	my $gc=$line;
	$gc=~s/A|T//g;
	$gc=~s/a|t//g;
	$seq{full}+=length($line);
	$seq{gc}+=length($gc);
	$seq{num}++;
	$length{$id}=length($line);
}
close In;
my $sum=0;
my $N50=0;
foreach my $l (sort {$b <=> $a} values %length) {
	$sum+=$l;
	if ($sum / $seq{full} > 0.5) {
		$N50=$l;
		last;
	}
}
open Out,">$output";
print Out "AssembelLevel\t","$level\n";
print Out "SeqNum\t",$seq{num},"\n";
print Out "TotalLength\t",$seq{full},"\n";
print Out "N50\t",$N50,"\n";
print Out "GC%\t",int($seq{gc}/$seq{full}*10000)/100,"%\n";
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
Usage:
  Options:
  -ref	<file>	input reference file name
  -gff	<file>	input gff file 
  -output	<file>	input seq stat name
  -h         Help

USAGE
        print $usage;
        exit;
}
