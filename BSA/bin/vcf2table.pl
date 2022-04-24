#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$P1,$P2,$B1,$B2,$Pldep,$Bldep,$Bhdep,$Phdep,$popt,$Vtype,$group);
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fIn,
	"out:s"=>\$fOut,
	"group:s"=>\$group,
	"popt:s"=>\$popt,
	"vtype:s"=>\$Vtype
				) or &USAGE;
&USAGE unless ($fIn and $fOut );
$Pldep||=10;
$Phdep||=1000;
$Bldep||=10;
$Bhdep||=1000;
$popt||="F2";
$Vtype||="ALL";
open In,$group;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sampleID,$group,$low,$high)=split(/\s+/,$_);
	if($group eq "P1"){
		$P1=$sampleID;
		$Pldep=$low;
		$Phdep=$high;
	}elsif($group eq "P2"){
		$P2=$sampleID;
		$Pldep=$low;
		$Phdep=$high;
	}elsif($group eq "B1"){
		$B1=$sampleID;
		$Bldep=$low;
		$Bhdep=$high;
	}elsif($group eq "B2"){
		$B2=$sampleID;
		$Bldep=$low;
		$Bhdep=$high;
	}
}
close In;
if($B1 eq "-" && $B2 eq "-"){
	die "error group file die $B1 $B2";
}elsif($B2 eq "-"){
	die "$B2";
}

if ($fIn =~ /\.gz$/){
	open In,"gzip -dc $fIn|";
}else{
	open In,$fIn;
}
open Out,">$fOut";
my @Sample;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\t/,$_);
		push @Sample,@info;
		my @out;
		push @out,join("\t","$P1.GT","$P1.AD","$P1.DP") if ($P1 ne "-");
		push @out,join("\t","$P2.GT","$P2.AD","$P2.DP") if ($P2 ne "-");
		push @out,join("\t","$B1.GT","$B1.AD","$B1.DP") if ($B1 ne "-");
		push @out,join("\t","$B2.GT","$B2.AD","$B2.DP") if ($B2 ne "-");
		print Out join("\t","CHROM\tPOS\tRef\tAlt\tVtype",@out,"Ann"),"\n";
	}else{
		my (%samples,$vtype,@ann,$baseinfo,$allenum);
		$vtype=readvcf($_,\%samples,\@Sample,\@ann,\$baseinfo,\$allenum);
		if ($vtype ne $Vtype && $Vtype ne "ALL"){
			$stat{vtype}++;
			next;
		}
		if ($B1 ne "-" && ($samples{$B1}{GT} eq "NN" || $samples{$B1}{DP} <= $Bldep|| $samples{$B1}{DP} >= $Bhdep)){
			$stat{B1missingOrdepth}++;
			next;
		}
		if ($samples{$B2}{GT} eq "NN" || $samples{$B2}{DP} <= $Bldep || $samples{$B2}{DP} >= $Bhdep){
			$stat{B2missingORdepth}++;
			next;
		}
		if ($P1 ne "-" && ($samples{$P1}{GT} eq "NN"||$samples{$P1}{DP} <= $Pldep || $samples{$P1}{DP} >= $Phdep)){
			$stat{P1missingORdepth}++;
			next;
		}
		if ($P2 ne "-" && ($samples{$P2}{GT} eq "NN"||$samples{$P2}{DP} <= $Pldep || $samples{$P2}{DP} >= $Phdep)){
			$stat{P2missingORdepth}++;
			next;
		}
		if ($allenum > 2){
			$stat{allenumMoreThan2}++;
			next;
		};
		if ($popt ne "F1" && defined $P1 && defined $P2 && $samples{$P1}{GT} eq $samples{$P2}{GT}){
			$stat{PMsameGenotype}++;
			next;
		};
		my @g1=split("/",$samples{$P1}{GT});
		my @g2=split("/",$samples{$P2}{GT}) if (defined $P2);
		if ($popt ne "F1" && ((defined $P1 && $g1[0] ne $g1[1]) || (defined $P2 && $g2[0] ne $g2[1]))){
			$stat{popErrF2buthete}++;
			next;
		};
		if ($popt eq "F1" && ((defined $P2 && $g1[0] eq $g1[1]) ||(defined $P2 && $g2[0] eq $g2[1]))){
			$stat{popErrF1buthomo}++;
			next;
		}

		my @out;
		push @out,join("\t",$samples{$P1}{GT},$samples{$P1}{AD},$samples{$P1}{DP}) if (defined $P1);
		push @out,join("\t",$samples{$P2}{GT},$samples{$P2}{AD},$samples{$P2}{DP}) if (defined $P2);
		push @out,join("\t",$samples{$B1}{GT},$samples{$B1}{AD},$samples{$B1}{DP}) if (defined $B1);
		push @out,join("\t",$samples{$B2}{GT},$samples{$B2}{AD},$samples{$B2}{DP}) if (defined $B2);
		print Out join("\t",$baseinfo,$vtype,@out,join(";",@ann)),"\n";
	}
}
close In;
close Out;
open Out,">$fOut.stat";
foreach my $key (sort keys %stat){
	print Out $key,"\t",$stat{$key},"\n";
}
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub readvcf{
	my ($line,$sample,$Sample,$anninfo,$baseinfo,$allenum)=@_;
	my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\s+/,$line);
	$$baseinfo=join("\t",$chrom,$pos,$ref,$alt);
	my @alles=split(",",join(",",$ref,$alt));
	$$allenum=scalar @alles;
	if($info=~/ANN=([^\;]*)/g){
		my @ann=split(/\,/,$1);
		for (my $i=0;$i<@ann;$i++) {
			my @str=split(/\|/,$ann[$i]);
			$str[0]||="--";
			$str[1]||="--";
			$str[2]||="--";
			$str[3]||="--";
			$str[4]||="--";
			$str[10]||="--";
			my $ann=join("|",$str[0],$str[1],$str[2],$str[3],$str[4],$str[10]);
			push @{$anninfo},$ann;
		}
	}
	my %len;
	for (my $i=0;$i<@alles;$i++) {
		$len{length($alles[$i])}=1;
	}
	my $type="SNP";
	if (scalar keys %len > 1) {
		$type="INDEL";
	}
	my @format=split(/\:/,$format);
	for (my $i=0;$i<@info;$i++) {
		my @infos=split(/\:/,$info[$i]);
		for (my $j=0;$j<@infos;$j++) {
			if ($format[$j] eq "GT") {
				if ($infos[$j] =~ /\./) {
					$$sample{$$Sample[$i]}{$format[$j]}="NN";
					$$sample{$$Sample[$i]}{DP}=0;
					$$sample{$$Sample[$i]}{AD}=0;
				}else{
					my @gt=split(/\//,$infos[$j]);
					$$sample{$$Sample[$i]}{$format[$j]}=join("/",sort($alles[$gt[0]],$alles[$gt[1]]));
				}
			}
			if ($format[$j] eq "AD") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
			if ($format[$j] eq "DP") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
		}
	}
	return $type;
}
sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-vcf	<file>	input file 
	-out	<file>	output file
	--group	<file>	input group file
		mbid must be given

		wpid    P1	phdep   pldep
		wbid    P2	mbip    bhdep  

	-popt	<str>	population type default F2
	-vtype	<str>	varaint type default ALL
USAGE
	print $usage;
	exit;
}
