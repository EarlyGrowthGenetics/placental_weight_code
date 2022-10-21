#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $poe=0;
my $pair=0;
my $cm=0.2653;
my $cp=0.197;
my $mp=0.0082;
GetOptions(
  "poe!"  => \$poe,
  "pair!" => \$pair,
  "mp=f"  => \$mp,
  "cm=f"  => \$cm,
  "cp=f"  => \$cp,
);
if($mp && $cp && $cm){
my $dir="../frozen_results_files/";
open(my $in,"zcat ".$dir."pw_fetal_sex_gest1.tbl.gz | ") or die $!;
my %fetal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $fetal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[8]."\t".$F[9]."\t".$F[15];
}
close($in);
open($in,"zcat ".$dir."pw_maternal_sex_gest1.tbl.gz | ") or die $!;
my %maternal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $maternal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[8]."\t".$F[9]."\t".$F[15];
}
close($in);
open($in,"zcat ".$dir."pw_paternal_sex_gest1.tbl.gz | ") or die $!;
my %paternal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $paternal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[8]."\t".$F[9]."\t".$F[15];
}
close($in);
open(my $out," | gzip -c > wlm_sex_gest.gz");
print $out "MarkerName\tEffect_allele\tOther_allele\teaf_fetal\tbeta_fetal\tse_fetal\tp_fetal\tn_fetal\teaf_maternal\tbeta_maternal\tse_maternal\tp_maternal\tn_maternal\teaf_paternal\tbeta_paternal\tse_paternal\tp_paternal\tn_paternal\twlm_beta_fetal\twlm_se_fetal\twlm_beta_maternal\twlm_se_maternal\twlm_beta_paternal\twlm_se_paternal";
if($poe){
  print $out "\twlm_beta_fetal_poe\twlm_se_fetal_poe\twlm_beta_maternal_poe\twlm_se_maternal_poe\twlm_beta_poe\twlm_se_poe";
}
if($pair){
  print $out "\twlm_beta_fetal_pair\twlm_se_fetal_pair\twlm_beta_maternal_pair\twlm_se_maternal_pair";
}
print $out "\n";
foreach my $key (keys %fetal){
  if(exists $maternal{$key} && exists $paternal{$key}){
    my @F=split(' ',$fetal{$key});
    my @M=split(' ',$maternal{$key});
    my @P=split(' ',$paternal{$key});
    print $out $key."\t".$fetal{$key};
    if($F[0] eq $M[0] && $F[1] eq $M[1]){
      print $out "\t".join("\t",$M[2],$M[3],$M[4],$M[5],$M[6]);
    }elsif($F[0] eq $M[1] && $F[1] eq $M[0]){
      $M[3]=-$M[3];
      print $out "\t".join("\t",$M[2],$M[3],$M[4],$M[5],$M[6]);
    }else{
      print "maternal\t".$key."\n";
    }
    if($F[0] eq $P[0] && $F[1] eq $P[1]){
      print $out "\t".join("\t",$P[2],$P[3],$P[4],$P[5],$P[6]);
    }elsif($F[0] eq $P[1] && $F[1] eq $P[0]){
      $P[3]=-$P[3];
      print $out "\t".join("\t",$P[2],$P[3],$P[4],$P[5],$P[6]);
    }else{
      print "paternal\t".$key."\n";
    }
    my $beta=2*$F[3]-$M[3]-$P[3];
    my $se=sqrt(4*$F[4]**2+$M[4]**2+$P[4]**2+2*$mp*sqrt(($M[4]**2)*($P[4]**2))-4*$cm*sqrt(($M[4]**2)*($F[4]**2))-4*$cp*sqrt(($P[4]**2)*($F[4]**2)));
    print $out "\t".$beta."\t".$se;
    $beta=(3*$M[3]-2*$F[3]+$P[3])/2;
    $se=sqrt((9*$M[4]**2)/4+$F[4]**2+($P[4]**2)/4-$cp*sqrt(($F[4]**2)*($P[4]**2))-3*$cm*sqrt(($M[4]**2)*($F[4]**2))+3*$mp*sqrt(($P[4]**2)*($M[4]**2))/2);
    print $out "\t".$beta."\t".$se;
    $beta=(3*$P[3]-2*$F[3]+$M[3])/2;
    $se=sqrt((9*$P[4]**2)/4+$F[4]**2+($M[4]**2)/4-$cm*sqrt(($F[4]**2)*($M[4]**2))-3*$cp*sqrt(($P[4]**2)*($F[4]**2))+3*$mp*sqrt(($P[4]**2)*($M[4]**2))/2);
    print $out "\t".$beta."\t".$se;
    if($poe){
      $beta=4*$F[3]-2*$M[3]-4*$P[3];
      $se=sqrt(16*$F[4]**2+4*$M[4]**2+16*$P[4]**2+16*$mp*sqrt(($M[4]**2)*($P[4]**2))-16*$cm*sqrt(($M[4]**2)*($F[4]**2))-32*$cp*sqrt(($P[4]**2)*($F[4]**2)));
      print $out "\t".$beta."\t".$se;
      $beta=(2*$M[3]-2*$F[3]+2*$P[3]);
      $se=sqrt((4*$M[4]**2)+4*$F[4]**2+(4*$P[4]**2)-8*$cp*sqrt(($F[4]**2)*($P[4]**2))-8*$cm*sqrt(($M[4]**2)*($F[4]**2))+8*$mp*sqrt(($P[4]**2)*($M[4]**2))/2);
      print $out "\t".$beta."\t".$se;
      $beta=(6*$P[3]-4*$F[3]+2*$M[3]);
      $se=sqrt((36*$P[4]**2)+16*$F[4]**2+4*($M[4]**2)-16*$cm*sqrt(($F[4]**2)*($M[4]**2))-48*$cp*sqrt(($P[4]**2)*($F[4]**2))+24*$mp*sqrt(($P[4]**2)*($M[4]**2))/2);
      print $out "\t".$beta."\t".$se;
    }
    if($pair){
      $beta=(4*$F[3]-2*$M[3])/3;
      $se=sqrt(16*$F[4]**2+4*$M[4]**2-16*$cm*sqrt(($M[4]**2)*($F[4]**2)))/3;
      print $out "\t".$beta."\t".$se;
      $beta=(4*$M[3]-2*$F[3])/3;
      $se=sqrt((16*$M[4]**2)+4*$F[4]**2-16*$cm*sqrt(($M[4]**2)*($F[4]**2)))/3;
      print $out "\t".$beta."\t".$se;
    }
    print $out "\n";
  }
}
close($in);
}else{
  print "\n\tUsage:\n\t--cm\toverlap of child/maternal\n\t--cp\toverlap of child/paternal\n\t--mp\toverlap of maternal/paternal\n\t--poe\toptionan\n\t--pair\toptional\n\n";
}
