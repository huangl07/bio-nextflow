#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output
			) or &USAGE;
&USAGE unless ($input and $output);
open In,"bgzip -dc $input|";
my @indi;
my %type;
my $totalsv=0;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=split(/\t/,$_);
		$totalsv++;
		my $SVtype;
		if ($info =~ /SVTYPE=([^;]*);/) {
			$SVtype=$1;
			for (my $i=0;$i<@samples;$i++) {
				my $gt=split(/\:/,$samples[0]);
				$stat{$indi[$i]}{$SVtype}++ if($gt ne "0/0" && $gt ne "./." );
			}
		}
	}
}
close In;
open Out,">$output";
print Out "#totalsv $totalsv\n";
print Out "#sampleID\tINS\tDEL\tDUP\tBND\n";
foreach my $sample (keys %stat) {
	print Out join("\t",$sample,"\t",$stat{$sample}{"INS"},"\t",$stat{$sample}{"DEL"},"\t",$stat{$sample}{"DUP"},"\t",$stat{$sample}{"BND"}),"\n";
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
Usage:
  Options:
  -input	<file>	input reference file name
  -output	<file>	input gff file name
  -h         Help

USAGE
        print $usage;
        exit;
}
