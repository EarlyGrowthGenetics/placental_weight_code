#!/usr/bin/perl
use warnings;
use strict;
sub format_results{
	my $infile=$_[0];
	my $outfile=$_[1];
	open(my $in,"zcat $infile | ") or die $!;
	open(my $out," | gzip -c > $outfile");
	my $line=<$in>;
	print $out $line;
	while(<$in>){
		chomp;
		my @F=split(' ',$_,2);
		my @G=split(":",$F[0]);
		print $out $G[0]."\t".$F[1]."\n";
	}
	close($out);
	close($in);
}
format_results("/slade/local/GWAS_datasets/Preeclampsia_2020/EGAF00004682218/fetal_all_chrALL_STERR_EU.1tbl.gz","fetal_preeclampsia_filtered.gz");
format_results("/slade/local/GWAS_datasets/Preeclampsia_2020/EGAF00004682220/mat_all_chrALL_STERR_EU.1tbl.gz","maternal_preeclampsia_filtered.gz");
