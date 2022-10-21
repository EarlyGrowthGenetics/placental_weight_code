#!/usr/bin/perl
# script to extract SNPs from meta-analysis results for wlm adjusted analysis
use warnings;
use strict;
use Getopt::Long;
my $infile=undef;
my $outfile=undef;
my $snpfile=undef;
my $trait=undef;
my $beta=undef;
my $se=undef;
GetOptions(
  "snps=s"     => \$snpfile,
  "infile=s"   => \$infile,
  "outfile=s"  => \$outfile,
  "trait=s"    => \$trait,
  "beta=s"     => \$beta,
  "se=s"       => \$se,
);
if($infile && $outfile && $snpfile && $trait && $beta && $se){
  my %snps=();
  open(my $in,$snpfile) or die $!;
  while(<$in>){
    chomp;
    my @F=split(' ');
    if(!exists $snps{$F[0].":".uc($F[1]).":".uc($F[2])} && !exists $snps{$F[0].":".uc($F[2]).":".uc($F[1])}){
      $snps{$F[0].":".uc($F[1]).":".uc($F[2])}=1;
    }
  }
  close($in);
  open($in,"zcat ".$infile." | ");
  open(my $out,">>",$outfile);
  <$in>;
  while(<$in>){
    chomp;
    my @F=split(' ');
    if(exists $snps{$F[0].":".uc($F[1]).":".uc($F[2])}){
      print $out $trait."\t".$F[0]."\t".uc($F[1])."\t".$F[$beta]."\t".$F[$se]."\tNA\n";
      delete $snps{$F[0].":".uc($F[1]).":".uc($F[2])};
    }elsif(exists $snps{$F[0].":".uc($F[2]).":".uc($F[1])}){
      print $out $trait."\t".$F[0]."\t".uc($F[1])."\t".$F[$beta]."\t".$F[$se]."\tNA\n";
      delete $snps{$F[0].":".uc($F[2]).":".uc($F[1])};
    }
  }
  close($out);
  close($in);
  foreach my $key (keys %snps){
    print $key."\n";
  }
}else{
  print "\n\tUsage:\n\t--infile\tinput summary stats filename\n\t--outfile\toutput filename (trait_name rsid effect_allele beta se p\n\t--snps\tSNPs filename (SNP_name a1 a2)\n\t--trait\ttrait name\n\t--beta\tcolumn containing beta\n\t--se\tcolumn containing SE\n\n";
}
