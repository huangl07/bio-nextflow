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
mkdir $fOut if(!-d $fOut);
$fOut=ABSOLUTE_DIR($fOut);
open In,$fIn;
open List,">$fOut/group.list";
my %filehand;
while(<In>){
    chomp;
    next if($_ eq ""||/^$/);
    my ($id,$gid)=split(/\s+/,$_);
    if(!exists $filehand{$gid}){
        print List "$gid\t$fOut/$gid.list\n";
        open $filehand{$gid},">$fOut/$gid.list";
    }
    print {$filehand{$gid}} $id,"\n";
}
close In;
close List;

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
  -i	<file>	input dict name
  -o	<file>	out qc base file name
  -h         Help

USAGE
        print $usage;
        exit;
}
