#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$gff,$dOut,$group,$condition);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
    "group:s"=>\$group,
	"condition:s"=>\$condition,
	"output:s"=>\$output,
    "dOut:s"=>\$dOut
			) or &USAGE;
&USAGE unless ($input and $output);
$dOut||="./";
if(!-d $dOut){mkdir $dOut;}
$dOut=ABSOLUTE_DIR($dOut);
$input=ABSOLUTE_DIR($input);
open In,$group;
my %gid;
while (<In>){
    chomp;
    next if ($_ eq ""||/^$/ ||/^#/);
    my ($id,$gid)=split(/\s+/,$_);
    my @bam=glob("$input/$id*bam");
    die $id if (scalar @bam !=1);
    push @{$gid{$gid}},$bam[0];
}
close In;
my %gfile;
foreach my $gid(sort keys %gid){
    open Out,">$dOut/$gid.bam.list";
    print Out join(",",@{$gid{$gid}}),"\n";
    close Out;
    $gfile{$gid}="$dOut/$gid.bam.list"
}
open Out,">$output";
if(!$condition){
    my @gid=sort keys %gid;
    for(my $i=0;$i<@gid;$i++){
        for(my $j=$i+1;$j<@gid;$j++){
            print Out "-b1  $gfile{$gid[$i]} -b2  $gfile{$gid[$j]} -od $gid[$i]vs$gid[$j]\n";
        }
    }
}else{
    open In,$condition;
    while(<In>){
        chomp;
        next if ($_ eq ""||/^$/ ||/^#/);
        my ($g1,$g2)=split;
        print $g1,"\n" if (!exists $gid{$g1});
        print Out "-b1 $gfile{$g1} -b2 $gfile{$g2} -od $g1"."vs"."$g2\n";
    }
    close In;

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
  -ref	<file>	input reference file name
  -gff	<file>	input gff file 
  -output	<file>	input seq stat name
  -h         Help

USAGE
        print $usage;
        exit;
}
