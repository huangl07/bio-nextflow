#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$vcf,$group);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "vcf:s"=>\$vcf,
    "g:s"=>\$group,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($vcf and $group and $fOut);
open In,$group;
my %fgroup;
open SH,">$fOut/poplddecay.sh";
open List,">$fOut/pop.list";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($id,$gid)=split(/\s+/,$_);
    if(!exists $fgroup{$gid}){
        open $fgroup{$gid},">$fOut/$gid.list";
        print SH "PopLDdecay --InVCF $vcf --OutStat $fOut/$gid --MAF 0.05 --Miss 0.3 -SubPop $fOut/$gid.list\n";
        print List "$gid\t$fOut/$gid.stat.gz\n";
    }
    print {$fgroup{$gid}} $id,"\n";

}
close In;
close SH;
close List;
`sh $fOut/poplddecay.sh`;
`Rscript $Bin/ld-decay.R --list $fOut/pop.list --outfile $fOut/pop`;
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
    "vcf:s"=>\$vcf,
    "g:s"=>\$group,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
