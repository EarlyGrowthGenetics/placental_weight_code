#!/usr/bin/perl
use warnings;
use strict;

sub classify_pair{
  my $f=$_[0];
  my $m=$_[1];
  my $F=$_[2];
  my $thresh=$_[3];
  my $classification="Unclassified";
  
  # which are > thresh
  my $fthresh=$$F[$f+2]<$thresh ? 1 : 0;
  my $mthresh=$$F[$m+2]<$thresh ? 1 : 0;
  
  # only 1 passes threshold
  if($fthresh==1 && $mthresh==0){
    $classification=(($$F[$f]-1.96*$$F[$f+1])>($$F[$m]+1.96*$$F[$m+1]) || (($$F[$f]+1.96*$$F[$f+1])<$$F[$m]-1.96*$$F[$m+1])) ? "Fetal" : "Unclassified";
  }
  if($fthresh==0 && $mthresh==1){
    $classification=(($$F[$f]-1.96*$$F[$f+1])>($$F[$m]+1.96*$$F[$m+1]) || (($$F[$f]+1.96*$$F[$f+1])<$$F[$m]-1.96*$$F[$m+1])) ? "Maternal" : "Unclassified";
  }

  # both pass threshold
  if($fthresh==1 && $mthresh==1){
    $classification="Fetal & Maternal";
  }
 
  return $classification;
}

open(my $in,"snps") or die $!;
my $line=<$in>;
my $f=15;
my $m=18;
my $p=21;
my $fe=24;
my $me=27;
my $pe=30;
my $fp=33;
my $mp=36;
my @F=split(' ',$line);
print join("\t",$F[$f],$F[$m],$F[$p],$F[$fe],$F[$me],$F[$pe],$F[$fp],$F[$mp])."\n";
while(<$in>){
  chomp;
  my @F=split(' ');
  my $pair=classify_pair($fp,$mp,\@F,0.001315789);
  print $F[3]."\t".$pair."\n";
}
close($in);
