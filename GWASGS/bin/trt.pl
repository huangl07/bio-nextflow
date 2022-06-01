#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt);
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
mkdir $fOut if(!-d $fOut);
$fOut=ABSOLUTE_DIR($fOut);
open In,$fIn;
my %trt;
my @trtname;
my %filehand;
my $split="\t";
if($fIn=~/csv$/){
    $split=",";
}
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/genotype/){
        (undef,@trtname)=split($split,$_);
        foreach my $trt(@trtname){
            open $filehand{$trt},">$fOut/$trt.txt";
            print {$filehand{$trt}} "sampleID\t$trt\n";
        }
    }else{
        my ($id,@trt)=split($split,$_);
        for(my $i=0;$i<@trt;$i++){
            print {$filehand{$trtname[$i]}} $id,"\t",$trt[$i],"\n";
        }
    }
}
close In;
open Out,">$fOut/trt.list";
foreach my $trt(keys %filehand){
    close $filehand{$trt};
    print Out "$trt\t","$fOut/$trt.txt\n";
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
