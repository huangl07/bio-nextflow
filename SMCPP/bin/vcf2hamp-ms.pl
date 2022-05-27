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
if($fIn !~ /.gz$/){
    open In,$fIn;
}else{
    open In,"bgzip -dc $fIn|";
}
open Out,">$fOut";
my @indi;
my $nchr=0;
my %nchr;
my $nsnp=0;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/ ||/^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	#next if ($Filter ne "PASS"  && $Filter ne "FILTER");
	if (/^#/) {
		push @indi,@geno;
		print Out "#CHROM\tPOS\t".join("\t",@indi),"\n";
	}else{
		my @ale=split(/\,/,join(",",$ref,$alt));
		if (!exists $nchr{$chr}){
			$nchr++;
			$nchr{$chr}=$nchr;
		}
		my $chro=$nchr{$chr};
		my $nsnp++;
		my %ale;
		my @out;
		for (my $i=0;$i<@geno;$i++) {
			my ($g1,$g2)=split(/\/|\|/,(split(/\:/,$geno[$i]))[0]);
			if ($g1 eq "." || $g2 eq "."){
				push @out,"N";
				next;
			}else{
				push @out,$ale[$g2];
			}
			$ale{$ale[$g1]}=1;
			$ale{$ale[$g2]}=1;
		}
		next if (scalar keys %ale >=3);
		my $pos=$pos;
		print Out  join("\t","chr_".$chro,$pos,@out),"\n";
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	vcf thanslate to msmc format
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
