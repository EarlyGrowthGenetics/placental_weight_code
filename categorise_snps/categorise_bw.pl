#!/usr/bin/perl
use warnings;
use strict;
# script to categorise BW and PW SNPs

sub classify_pair{
  my $pwf=$_[0];
  my $pwm=$_[1];
  my $bwf=$_[2];
  my $bwm=$_[3];
  my $F=$_[4];
  my $thresh=$_[5];
  my $classification="PW";
  
  # convert missing SNP to 0 
  if($$F[$bwf] eq "NA"){
    $$F[$bwf]=0;
    $$F[$bwf+1]=0;
    $$F[$bwf+2]=1;
  }
  if($$F[$bwm] eq "NA"){
    $$F[$bwm]=0;
    $$F[$bwm+1]=0;
    $$F[$bwm+2]=1;
  }

  # are either BW P > thresh
  my $fthresh=$$F[$bwf+2]<$thresh ? 1 : 0;
  my $mthresh=$$F[$bwm+2]<$thresh ? 1 : 0;
  $classification=($fthresh==1 || $mthresh==1) ? "PW & BW" : "PW";
  
  # if it's both, which beta is bigger
  my $PW=$$F[$pwf+2]<$$F[$pwm+2] ? $pwf : $pwm;
  my $BW=$$F[$bwf+2]<$$F[$bwm+2] ? $bwf : $bwm;
  my $BW1=$$F[$bwf+2]<$$F[$bwm+2] ? $bwm : $bwf;
  #print $$F[$pwf+2]."\t".$$F[$pwm+2]."\t".$$F[$bwf+2]."\t".$$F[$bwm+2]."\n";
  #if(($fthresh==1 || $mthresh==1) && $beta>0){
  my $bwlo=$$F[$BW]-1.96*$$F[$BW+1];
  my $bwhi=$$F[$BW]+1.96*$$F[$BW+1];
  my $bw1lo=$$F[$BW1]-1.96*$$F[$BW1+1];
  my $bw1hi=$$F[$BW1]+1.96*$$F[$BW1+1];
  my $pwlo=$$F[$PW]-1.96*$$F[$PW+1];
  my $pwhi=$$F[$PW]+1.96*$$F[$PW+1];
  $classification=($bwhi > $pwlo && $bwlo < $pwhi) || ($bwlo < $pwhi && $bwhi > $pwlo) ? "unclassified" : "single";	# if the strongest associated beta overlaps the PW beta
  $classification=$classification eq "unclassified" && $$F[$BW+2]<$thresh ? "PW & BW overlapping" : $classification;	# if they overlap, is the strongest less than the threshold
  $classification=$classification eq "unclassified" && !($mthresh || $fthresh) ? "PW overlapping" : $classification;
  if($classification eq "single"){
    $classification=$$F[$BW]>$$F[$PW] ? "BW stronger" : "PW stronger";
    $classification=$$F[$BW]<0 && $$F[$BW+2]<$thresh ? "PW & BW opposite directions" : $classification;
  }
  return $classification;
}

open(my $in,"pw_bw_snps1") or die $!;
open(my $out,">","snp_classifications");
my $line=<$in>;
chomp $line;
my $pwf=6;
my $pwm=9;
my $pwp=12;
my $pwf_wlm=15;
my $pwm_wlm=18;
my $pwp_wlm=21;
my $pwf_wlm_poe=24;
my $pwm_wlm_poe=27;
my $pwp_wlm_poe=30;
my $pwf_wlm_pair=33;
my $pwm_wlm_pair=46;
my $bwf=39;
my $bwm=42;
my $bwf_sem=45;
my $bwm_sem=48;
my @F=split(' ',$line);
print join("\t",$F[$pwf],$F[$pwm],$F[$pwp],$F[$pwf_wlm],$F[$pwm_wlm],$F[$pwp_wlm],$F[$pwf_wlm_poe],$F[$pwm_wlm_poe],$F[$pwp_wlm_poe],$F[$pwf_wlm_pair],$F[$pwm_wlm_pair],$F[$bwf],$F[$bwm],$F[$bwf_sem],$F[$bwm_sem])."\n";
while(<$in>){
  chomp;
  my @F=split(' ');
  my $beta_wlm=classify_pair($pwf_wlm,$pwm_wlm,$bwf_sem,$bwm_sem,\@F,0.001315789);
  print join("\t",$F[3],$beta_wlm."\n");
  print $out $_."\t".join("\t",$beta_wlm."\n");
}
close($in);
