#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$num);
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
if($fIn !~ /.gz$/){
    open In,$fIn;
}else{
    open In,"bgzip -dc $fIn|";
}
open Out,">$fOut";
my @indi;
my $head;
my %filehand;
my $n=0;
my %nchro;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/ ||/^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	#next if ($Filter ne "PASS"  && $Filter ne "FILTER");
	if (/^#/) {
		push @indi,@geno;
		$head="rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t".join("\t",@indi);
        print Out $head,"\n";
	}else{
		my @ale=split(/\,/,join(",",$ref,$alt));
        $nchro{$chr}=1;
		my %ale;
		my @out;
		for (my $i=0;$i<@geno;$i++) {
			my ($g1,$g2)=split(/\/|\|/,(split(/\:/,$geno[$i]))[0]);
			if ($g1 eq "." || $g2 eq "."){
				push @out,"NN";
				next;
			}else{
				push @out,join("",sort($ale[$g1],$ale[$g2]));
			}
			$ale{$ale[$g1]}=1;
			$ale{$ale[$g2]}=1;
		}
		next if (scalar keys %ale != 2);
        $n++;
		my $snpid="snp$n";
		my $alleles=join("/",sort keys %ale);
		my $chrom=scalar keys %nchro;
		my $pos=$pos;
		my $strand="+";
		my $assembly="Majorbio";
		my $center="Majorbio";
		my $protLSID="Majorbio";
		my $assayLSID="Majorbio";
		my $panel="Majorbio";
		my $QCcode="Majorbio";
		print Out join("\t",$snpid,$alleles,$chrom,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode,@out),"\n";
	}
}
close In;
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
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
