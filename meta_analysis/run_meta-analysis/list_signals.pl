#!/usr/bin/perl
# script to filter GWAS signals
# Robin Beaumont 2020
# r.beaumont@exeter.ac.uk
use warnings;
use strict;
use Getopt::Long;
my $files=undef;
my $thresh=5e-8;
my $n_thresh=5000;
my $study_thresh=2;
GetOptions(
  "files=s"  => \$files,
);

sub find_chr{
  my $chr=$_[0];
  my $array=$_[1];
  my $sorted=$_[2];
  my %kept=();
  foreach my $line (@$array){
    my @F=split(' ',$line,3);
    if($F[0] == $chr){
     $kept{$F[1]}=$F[1]."\t".$F[2];
    }
  }
  foreach my $key (sort {$a <=> $b} keys %kept){
    push @$sorted, $kept{$key};
  }
}

sub list_top_hits {
  my $hits=$_[0];
  my $chr=$_[1];
  my $count=$_[2];
  my $out=$_[3];
  $count++;
  my @F=split(' ',shift(@$hits));
  my $prev=$F[0];
  my $sig=join("\t",$chr,$F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[11],$F[12],$F[13],$F[14],$F[15]);
  my @a=(uc($F[1]),uc($F[2]));
  my @a1=sort @a;
  print $out join("\t",$chr.":".$F[0].":".$a1[0].":".$a1[1],$F[9])."\n";
  my $sigp=$F[9];
  my $nsnp=1;
  foreach my $val (@$hits){
    @F=split(' ',$val);
    @a=(uc($F[1]),uc($F[2]));
    @a1=sort @a;
    print $out join("\t",$chr.":".$F[0].":".$a1[0].":".$a1[1],$F[9])."\n";
    $nsnp++;
    if(($F[0]-$prev)>500000){
        print $sig."\t".$nsnp."\n";
        $sigp=$F[9];
        $sig=join("\t",$chr,$F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[11],$F[12],$F[13],$F[14],$F[15]);
        $count++;
        $nsnp=1;
    }else{
        # same signal
        if($F[9]<$sigp){
            $sig=join("\t",$chr,$F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[11],$F[12],$F[13],$F[14],$F[15]);
            $sigp=$F[9];
        }
    }
    $prev=$F[0];
  }
  print $sig."\t".$nsnp."\n";
  return $count;
}

if($files){
  my @F=split(",",$files);
  foreach my $file (@F){
    open(my $in,"zcat $file | ") or die $!;
    open(my $out,">",$file.".gwas_sig_ordered");
    open(my $outf," | gzip -c > ".$file.".filtered.gz");
#    print $out "CHR\tPOS\tEA\tOA\tFreq1\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tN\n";
    print $out "Markername\tP\n";
    my @GWAS=();
    my $line=<$in>;
    print $outf $line;
    while(<$in>){
      chomp;
      my @G=split(' ');
      my @H=split('',$G[10]);
      my $count=0;
      for(my $i=0;$i<@H;$i++){
        if($H[$i] ne "?"){
          $count++;
        }
      }
      if($G[15]>$n_thresh && $count>=$study_thresh){
        print $outf $_."\n";
      }
      if($G[9]<$thresh && $G[15]>$n_thresh && $count>=$study_thresh){
        my @H=split(":",$G[0]);
        my @Q=split(' ',$_,2);
        if($H[0] eq "X"){
          $H[0]=23;
        }
        push @GWAS, $H[0]."\t".$H[1]."\t".$Q[1];
      }
    }
    close($in);
    my $count=0;
    print "CHR\tPOS\tEA\tOA\tFreq1\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tN\tNSNP\n";
    for(my $i=1;$i<24;$i++){
      my @hits=();
      find_chr($i,\@GWAS,\@hits);
      if(@hits>0){
        $count=list_top_hits(\@hits,$i,$count,$out);
      }
    }
    print $count."\n";
    close($outf);
    close($out);
  }
}else{
  print "\n\t--files\tcomma separated list of files to read and generate gwas SNP lists\n\n";
}
