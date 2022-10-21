#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $locus=undef;
my $analysis=undef;
GetOptions(
  "locus=s"    => \$locus,
  "analysis=s" => \$analysis,
);
# Set filename and suffix (whether we're after sex or Sex and GA adjusted analysis)
my $file=undef;
my $suff=undef;
if($analysis eq "sex"){
  $file="sex_loci_pairs";
  $suff="locus";
}elsif($analysis eq "gest"){
  $file="sex_gest_loci_pairs";
  $suff="gest_locus";
}
# Open list of SNPs and get the line corresponding to the SNP we're interested in
open(my $in,$file) or die $!;
my $line=undef;
for(my $i=0;$i<$locus;$i++){
  $line=<$in>;
}
close($in);
chomp $line;
my @F=split(' ',$line);
# Open person1 summary stats file and save summary stats to a hash
open($in,"PW_".$suff."_person1_".$locus) or die $!;
<$in>;
my %hash=();
while(<$in>){
  chomp;
  my @G=split(' ',$_,3);
  $hash{$G[1]}=$G[2]
}
close($in);
# Open person2 summary stats file and output joined stats
open($in,"PW_".$suff."_person2_".$locus) or die $!;
<$in>;
open(my $out,">",$suff."_".$locus."_joined") or die $!;
print $out "Marker\tpos\teaf\tp1_beta\tp1_se\tp1_p\tp1_n\tp2_beta\tp2_se\tp2_p\tp2_n\tp2_eaf\n";
while(<$in>){
  chomp;
  my @G=split(' ');
  if(exists $hash{$G[1]}){
    my @H=split(' ',$hash{$G[1]});
    if(uc($G[7]) eq uc($H[5]) && uc($G[8]) eq uc($H[6])){
      print $out $G[0].":".$G[1]."\t".$G[1]."\t".$H[0]."\t".$H[1]."\t".$H[2]."\t".$H[3]."\t".$H[4]."\t".$G[3]."\t".$G[4]."\t".$G[5]."\t".$G[6]."\t".$G[2]."\n";
    }elsif(uc($G[7]) eq uc($H[6]) && uc($G[8]) eq uc($H[5])){
      my $eaf=1-$H[0];
      my $beta=0-$H[1];
      print $out $G[0].":".$G[1]."\t".$G[1]."\t".$eaf."\t".$beta."\t".$H[2]."\t".$H[3]."\t".$H[4]."\t".$G[3]."\t".$G[4]."\t".$G[5]."\t".$G[6]."\t".$G[2]."\n";
    }
  }
}
close($out);
close($in);
