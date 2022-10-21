#!/usr/bin/perl
# script to extract and align PW and BW SNPs from respective datasets
use warnings;
use strict;
# function to get column numbers for each person
sub get_colnum{
  my $array=$_[0];
  my $beta=$_[1];
  my $se=$_[2];
  my $p=$_[3];
  my $bi=$_[4];
  my $si=$_[5];
  my $pi=$_[6];
  for(my $i=0;$i<@$array;$i++){
    if($$array[$i] eq $beta){
      $bi=$i;
    }
    if($$array[$i] eq $se){
      $si=$i;
    }
    if($$array[$i] eq $p){
      $pi=$i;
    }
  }
}

sub extract_bw{
  my $file=$_[0];
  my $snp=$_[1];
  my $ea=$_[2];
  my $oa=$_[3];
  my $b=$_[4];
  my $s=$_[5];
  my $p=$_[6];
  my $rsid=$_[7];
  my $effect=$_[8];
  my $other=$_[9];
  my %found=();
  open(my $in,"zcat $file | ");
  while(<$in>){
    chomp;
    my @F=split(' ');
    if(exists $$rsid{$F[$snp]}){
      if(lc($F[$ea]) eq lc($$effect{$F[$snp]}) && lc($F[$oa]) eq lc($$other{$F[$snp]})){
        $$rsid{$F[$snp]}.="\t".$F[$b]."\t".$F[$s]."\t".$F[$p];
        $found{$F[$snp]}=1;
      }elsif(lc($F[$oa]) eq lc($$effect{$F[$snp]}) && lc($F[$ea]) eq lc($$other{$F[$snp]})){
        $$rsid{$F[$snp]}.="\t".-$F[$b]."\t".$F[$s]."\t".$F[$p];
        $found{$F[$snp]}=1;
      }
    }
  }
  close($in);
  for my $key (keys %$rsid){
    if(!exists $found{$key}){
      $$rsid{$key}.="\tNA\tNA\tNA";
    }
  }
}

open(my $in,"snps_to_get") or die $!;	# read in list of SNPs to extract
my %rsid=();
my %snp=();
my %ea=();
my %oa=();
while(<$in>){
  chomp;
  my @F=split(' ');
  $rsid{$F[2]}=1;
  $snp{$F[0].":".$F[1]}=1;
  $ea{$F[2]}=$F[3];
  $oa{$F[2]}=$F[4];
}
close($in);
# first read in PW data
open($in,"zcat ../frozen_meta_analysis_files/rsid_added/wlm_sex_gest_plus_p.gz | ");
my $head=<$in>;
chomp $head;
my @H=split(' ',$head);
# find column numbers of betas, ses and ps
my $fb=4;
my $fse=5;
my $fp=5;
my $mb=9;
my $mse=10;
my $mp=11;
my $pb=14;
my $pse=15;
my $pp=16;
my $fwb=18;
my $fwse=19;
my $fwp=36;
my $mwb=20;
my $mwse=21;
my $mwp=37;
my $pwb=22;
my $pwse=23;
my $pwp=38;
my $fwpb=24;
my $fwpse=25;
my $fwpp=39;
my $mwpb=26;
my $mwpse=27;
my $mwpp=40;
my $pwpb=28;
my $pwpse=29;
my $pwpp=41;
my $fwpab=30;
my $fwpase=31;
my $fwpap=42;
my $mwpab=32;
my $mwpase=33;
my $mwpap=43;
while(<$in>){
  chomp;
  my @F=split(' ');
  if(exists $rsid{$F[0]}){
    my $a=2;
    if(lc($F[1]) eq lc($ea{$F[0]}) && lc($F[2]) eq lc($oa{$F[0]})){
      $a=1;
    }elsif(lc($F[2]) eq lc($ea{$F[0]}) && lc($F[1]) eq lc($oa{$F[0]})){
      $a=0;
    }
    if($a==1){
      $rsid{$F[0]}=join("\t",$F[$fb],$F[$fse],$F[$fp],$F[$mb],$F[$mse],$F[$mp],$F[$pb],$F[$pse],$F[$pp],$F[$fwb],$F[$fwse],$F[$fwp],$F[$mwb],$F[$mwse],$F[$mwp],,$F[$pwb],$F[$pwse],$F[$pwp],$F[$fwpb],$F[$fwpse],$F[$fwpp],$F[$mwpb],$F[$mwpse],$F[$mwpp],,$F[$pwpb],$F[$pwpse],$F[$pwpp],$F[$fwpab],$F[$fwpase],$F[$fwpap],$F[$mwpab],$F[$mwpase],$F[$mwpap]);
    }elsif($a==0){
      $rsid{$F[0]}=join("\t",-$F[$fb],$F[$fse],$F[$fp],-$F[$mb],$F[$mse],$F[$mp],-$F[$pb],$F[$pse],$F[$pp],-$F[$fwb],$F[$fwse],$F[$fwp],-$F[$mwb],$F[$mwse],$F[$mwp],-$F[$pwb],$F[$pwse],$F[$pwp],-$F[$fwpb],$F[$fwpse],$F[$fwpp],-$F[$mwpb],$F[$mwpse],$F[$mwpp],-$F[$pwpb],$F[$pwpse],$F[$pwpp],-$F[$fwpab],$F[$fwpase],$F[$fwpap],-$F[$mwpab],$F[$mwpase],$F[$mwpap]);
    }
  }
}
close($in);
extract_bw("~/EGG_Summary_Data/BW5/Fetal_BW_European_meta.NG2019.txt.gz",11,3,4,6,7,8,\%rsid,\%ea,\%oa);
extract_bw("~/EGG_Summary_Data/BW5/Maternal_BW_European_meta.NG2019.txt.gz",0,3,4,6,7,8,\%rsid,\%ea,\%oa);
extract_bw("~/EGG_Summary_Data/BW5/Fetal_Effect_European_meta_NG2019.txt.gz",1,4,5,7,8,9,\%rsid,\%ea,\%oa);
extract_bw("~/EGG_Summary_Data/BW5/Maternal_Effect_European_meta_NG2019.txt.gz",1,4,5,7,8,9,\%rsid,\%ea,\%oa);
open($in,"pw_bw_snps") or die $!;
open(my $out,">","pw_bw_snps1");
my $line=<$in>;
print $out "Person\tChr\tPos\trsid\tEA\tOA\tbeta_fetal\tse_fetal\tp_fetal";
print $out "\tbeta_maternal\tse_maternal\tp_maternal";
print $out "\tbeta_paternal\tse_paternal\tp_paternal";
print $out "\tbeta_fetal_wlm\tse_fetal_wlm\tp_fetal_wlm";
print $out "\tbeta_maternal_wlm\tse_maternal_wlm\tp_maternal_wlm";
print $out "\tbeta_paternal_wlm\tse_paternal_wlm\tp_paternal_wlm";
print $out "\tbeta_fetal_wlm_poe\tse_fetal_wlm_poe\tp_fetal_wlm_poe";
print $out "\tbeta_maternal_wlm_poe\tse_maternal_wlm_poe\tp_maternal_wlm_poe";
print $out "\tbeta_paternal_wlm_poe\tse_paternal_wlm_poe\tp_paternal_wlm_poe";
print $out "\tbeta_fetal_wlm_pair\tse_fetal_wlm_pair\tp_fetal_wlm_pair";
print $out "\tbeta_maternal_wlm_pair\tse_maternal_wlm_pair\tp_maternal_wlm_pair";
print $out "\tbeta_fetal_bw\tse_fetal_bw\tp_fetal_bw";
print $out "\tbeta_maternal_bw\tse_maternal_bw\tp_maternal_bw";
print $out "\tbeta_fetal_bw_sem\tse_fetal_bw_sem\tp_fetal_bw_sem";
print $out "\tbeta_maternal_bw_sem\tse_maternal_bw_sem\tp_maternal_bw_sem\n";
while(<$in>){
  chomp;
  my @F=split(' ');
  print $out join("\t",$F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$rsid{$F[3]}."\n");
}
close($out);
close($in);
