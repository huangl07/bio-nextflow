#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$chr,$dsh,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"group:s"=>\$group,
    "out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($group);
open In,$group;
my %group;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($sampleID,$groupID)=split(/\s+/,$_);
    push @{$group{$groupID}},$sampleID;
}
close In;
open Out,">$fOut";
foreach my $id(sort keys %group){
    print Out join("\t",$id,join(",",@{$group{$id}})),"\n";
}
close Out;
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	run smc
	eg:
	perl $Script -group -chr -out -dsh -h

Usage:
  Options:
  -group    <file>	input group list files#tab format eg:sample1    group;
  -out	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
