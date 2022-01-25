#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$gc);
use Data::Dumper;
use JSON;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
my @stat=glob("*.stat");
my %stat;
open Out,">$fOut";
print Out "#sampleID\tReadsNumber\tAllRatio(%)\tUniqRatio(%)\n";
foreach my $stat(@stat){
    my $id=(split(/\./,basename($stat)))[0];
    open In,$stat;
    while (<In>){
        s/^\s+//g;
        chomp;
        next if ($_ eq ""||/^$/);
        if(/reads; of these:/){
            $stat{$id}{all}=(split(/\s+/,$_))[0];
        }
        if(/\(([^%]*)%\) aligned concordantly exactly 1 time/){
             $stat{$id}{one}=(split(/\s+/,$_))[0];
        }
        if(/\(([^%]*)%\) aligned discordantly exactly 1 time/){
             $stat{$id}{one}+=(split(/\s+/,$_))[0];
        }
        if(/\(([^%]*)%\) aligned exactly 1 time/){
             $stat{$id}{one}+=(split(/\s+/,$_))[0]/2;
        }
        if(/overall alignment rate/){
            s/\%//g;
             $stat{$id}{ratio}=(split(/\s+/,$_))[0];
        }

    }
    close In;
    print Out $id,"\t",$stat{$id}{all},"\t",$stat{$id}{ratio},"\t",sprintf("%.4f",$stat{$id}{one}/$stat{$id}{all})*100,"\n";
}
close Out;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub average{
	my ($a,$b,$c,$d)=@_;
	my $sum=0;
	my $n=0;
	for (my $i=0;$i<@{$a};$i++) {
		$sum+=$$a[$i];
		$sum+=$$b[$i];
		$n++;
	}
	for (my $i=0;$i < @{$c} ;$i++) {
		$sum+=$$c[$i];
		$sum+=$$d[$i];
		$n++;
	}
	return($sum/$n);
}
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
  -i	<file>	input dict name
  -o	<file>	out qc base file name
  -h         Help

USAGE
        print $usage;
        exit;
}
