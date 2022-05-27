#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$pop,$queue,$shell,$num);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fIn,
	"gro:s"=>\$pop,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);

mkdir $fOut if (!-d $fOut);
$fOut=ABSOLUTE_DIR($fOut);
my %pop;
open In,$pop;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($id,$popid)=split(/\s+/,$_);
	$pop{$id}=$popid;
}
close In;
if($fIn !~ /.gz$/){
    open In,$fIn;
}else{
    open In,"bgzip -dc $fIn|";
}
my %seq;
my @Indi;
my %BASE=(
	"AA"=>"A","GG"=>"G","CC"=>"C","TT"=>"T",
	"AT"=>"W","AG"=>"R","AC"=>"M",
	"CG"=>"S","CT"=>"Y",
	"GT"=>"K"
);
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);
	}else{
		my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@indi)=split(/\t/,$_);
		my @ale=split(/\,/,join(",",$REF,$ALT));
		for (my $i=0;$i<@indi;$i++) {
			my $geno=(split(/\:/,$indi[$i]))[0];
			my ($g1,$g2)=split(/\/|\|/,$geno);
			if ($g1 eq ".") {
					if (!exists $seq{$Indi[$i]}{$CHROM}) {
						$seq{$Indi[$i]}{$CHROM}="N";
					}else{
						$seq{$Indi[$i]}{$CHROM}.="N";
					}
			}else{
				if (!exists $BASE{join("",sort($ale[$g1],$ale[$g2]))}){
					if (!exists $seq{$Indi[$i]}{$CHROM}) {
						$seq{$Indi[$i]}{$CHROM}="N";
					}else{
						$seq{$Indi[$i]}{$CHROM}.="N";
					}
					next;
				};
				if (!exists $seq{$Indi[$i]}{$CHROM}) {
					$seq{$Indi[$i]}{$CHROM}=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
				}else{
					$seq{$Indi[$i]}{$CHROM}.=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
				}
			}
		}
	}
}
close In;
my $head=0;
my %poppsmc;
open List,">$fOut/psmc.list";
foreach my $id (sort keys %seq) {
	open FA,">$fOut/$id.fasta";
	foreach my $chr (sort keys %{$seq{$id}}) {
		print FA ">$chr\n$seq{$id}{$chr}\n";
	}
	close FA;
    print List $pop{$id},"\t",$id,"\t","$fOut/$id.fasta","\n";
}
close List;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
	"vcf:s"=>\$fIn,
	"gro:s"=>\$pop,
	"out:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
