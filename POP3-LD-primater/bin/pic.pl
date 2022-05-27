#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fgroup);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
    "g:s"=>\$fgroup,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fgroup;
my %group;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($sampleID,$groupID)=split;
    $group{$sampleID}=$groupID;
}
close In;
open In,$fIn;
my %PIC;
my @Indi;
my $head=0;
open Out,">$fOut\.out";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/ ||/^##/);
    if(/^#/){
        (undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\s+/,$_);
    }else{
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\s+/,$_);
        my %stat;
        my %sum;
        for(my $i=0;$i<@indi;$i++){
            my $gid=$group{$Indi[$i]};
            my $gt=(split(":",$indi[$i]))[0];
            next if ($gt eq "./." || $gt eq ".|.");
            my @gt=split(/\/|\|/,$gt);
            if(!exists $group{$Indi[$i]}){
                print $Indi[$i];die;
            }
            $stat{$gid}{$gt[0]}++;
            $stat{$gid}{$gt[1]}++;
            $sum{$gid}+=2;
        }
        my @head;
        push @head,sort keys %stat;
        my @out;
        foreach my $g(sort keys %stat){
            my @alle=keys %{$stat{$g}};
            my $sum=$sum{$g};
            my $minus=0;
            my @p;
            for(my $i=0;$i<@alle;$i++){
                my $p=$stat{$g}{$alle[$i]}/$sum;
                push @p,$p;
            }
            my $minus1=0;
            my $minus2=0;
            for(my $i=0;$i<@p;$i++){
                $minus1+=$p[$i]*$p[$i];
                for(my $j=$i+1;$j<@p;$j++){
                    $minus2+=2*$p[$i]*$p[$i]*$p[$j]*$p[$j];
                }
            }
            my $PIC = 1- $minus1-$minus2;
            $PIC{$g}{value}+=$PIC;
            $PIC{$g}{num}+=1;
            push @out,$PIC;
        }
        if($head == 0){
            print Out "#chr\tpos\t",join("\t",@head),"\n";
            $head =1;
        }
    }
}
close In;
close Out;
open Out,">$fOut\.stat";
foreach my $g(sort keys %PIC){
    print Out $g,"\t",$PIC{$g}{value}/$PIC{$g}{num},"\n";
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
