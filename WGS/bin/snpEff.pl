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
open In,$input;
my $totalsnp=0;
my $type;
my %table;
my @head;
my %stat;
my $vtype="snp";
while (<In>) {
	s/ ,/,/g;
	s/, /,/g;
	chomp;
	next if ($_ eq ""||/^$/);
	if (/Number_of_variants_processed/) {
		$totalsnp=(split(/,/,$_))[-1];
	}elsif(/\# Count by effects/){
		$type="effects";
	}elsif(/\# Count by genomic region/){
		$type="region";
	}elsif(/\# Ts\/Tv : All variants/){
		$type="TsTv";
	}elsif(/# Hom\/Het table/){
		$type="Home";
	}elsif(/^#/){
		$type="Other";
	}
	if (/InDel lengths/) {
		$vtype="indel";
	}
	if ($type eq "effects" || $type eq "region") {
		my @info=split(/\,/,$_);
		push @{$table{$type}},join("\t",@info);
	}elsif ($type eq "TsTv") {
		if (/Sample/) {
			(undef,@head)=split(/,/,$_);
		}else{
			my ($types,@value)=split(/,/,$_);
			for (my $i=0;$i<@head;$i++) {
				$stat{$head[$i]}{$types}=$value[$i];
			}
		}
		
	}elsif ($type eq "Home") {
		if (/Sample_names/) {
			(undef,@head)=split(/,/,$_);
		}else{
			my ($types,@value)=split(/,/,$_);
			for (my $i=0;$i<@head;$i++) {
				$stat{$head[$i]}{$types}=$value[$i];
			}
		}
	}
}
open Out,">$output.effects.xls";
print Out join("\n",@{$table{"effects"}}),"\n";
close Out;
open Out,">$output.region.xls";
print Out join("\n",@{$table{"region"}}),"\n";
close Out;
open Out,">$output.sample.xls";
if ($vtype eq "snp") {
	print Out "#total $totalsnp\n";
	print Out join("\t","#sampleID","Total","Ts/Tv","Het","Hom"),"\n";
	foreach my $sample (sort keys %stat) {
		next if ($sample =~ /Total/);
		my $line = join("\t",$sample,$stat{$sample}{"Het"}+$stat{$sample}{"Hom"},$stat{$sample}{"Ts/Tv"},$stat{$sample}{"Het"},$stat{$sample}{"Hom"});
		print Out $line,"\n";
	}
	close Out;
}else{
	print Out "#total $totalsnp\n";
	print Out join("\t","#sampleID","Total","Het","Hom"),"\n";
	foreach my $sample (sort keys %stat) {
		next if ($sample =~ /Total/);
		my $line = join("\t",$sample,$stat{$sample}{"Het"}+$stat{$sample}{"Hom"},$stat{$sample}{"Het"},$stat{$sample}{"Hom"});
		print Out $line,"\n";
	}
	close Out;
}

close In;


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
