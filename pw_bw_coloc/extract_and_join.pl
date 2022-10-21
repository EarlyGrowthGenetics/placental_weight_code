#!/usr/bin/perl
use warnings;
use strict;
my $pfetal="../../frozen_results_files/pw_fetal_sex_gest_1.gz";
my $pmaternal="../../frozen_results_files/pw_maternal_sex_gest_1.gz";
my $ppaternal="../../frozen_results_files/pw_paternal_sex_gest.gz";
my $dfetal="~/decode_BW_gwas/Birthweight2021.gz";
my $dmaternal="~/decode_BW_gwas/Birthweight_offspring_mothers2021.gz";
my @b37=();	# b37 position
my @b38=();	# b38 position
my @pwperson=();	# person to get pw snps for
my @dperson=();	# person to get decode snps for
open(my $in1,"decode_loci_positions_b37") or die $!;
while(<$in1>){
	chomp;
	my @L=split(' ');
	push @b37, $L[0].":".$L[1];
	push @pwperson, $L[4];	# keep a record of which GWAS we need to pull this out from
	if($L[3] eq "fetal" || $L[3] eq "maternal"){	# decide whether we want fetal or maternal from decode for SNPs coming up in both based on PW person
		push @dperson, $L[3];
	}else{
		push @dperson, lc($L[4]);
	}
}
close($in1);
open($in1,"decode_loci_positions_b38") or die $!;
while(<$in1>){
	chomp;
	my @L=split(' ');
	push @b38, $L[0].":".$L[1];
}
close($in1);
# loop through the SNPs and pull out the summary stats
for(my $i=0;$i<@b37;$i++){
	my @L=split(":",$b37[$i]);
	my @M=split(":",$b38[$i]);
	print $L[0]."\t".$L[1]."\t".$M[0]."\t".$M[1]."\n";
	my $in=undef;
	if($pwperson[$i] eq "Fetal"){
		open($in,"zcat ".$pfetal." | ") or die $!;
	}elsif($pwperson[$i] eq "Maternal"){
		open($in,"zcat ".$pmaternal." | ") or die $!;
	}elsif($pwperson[$i] eq "Paternal"){
		open($in,"zcat ".$ppaternal." | ") or die $!;
	}else{
		die "Error: PW person not recognised\n";
	}
	<$in>;
	my %alleles=();
	my %egg=();
	my @egg_position=();
	while(<$in>){
		chomp;
		my @F=split(' ');
		if($F[16] eq $L[0] && $F[17]>$L[1]-500000 && $F[17]<$L[1]+500000){
			my @a1=sort($F[1], $F[2]);
			$alleles{$F[17]}=$a1[0].":".$a1[1];
			$egg{$F[17]}=join("\t",$F[2],$F[3],$F[7],$F[8],$F[9],$F[15]);
		}
	}
	close($in);
	my $loc=$i+1;
	open($in,"locus_".$loc) or die $!;	# read in b37 positions
	while(<$in>){
		chomp;
		my @F=split(' ');
		push @egg_position, $F[1];
	}
	close($in);
	open($in,"locus_".$loc."_lifted") or die $!;	# read in the b38 positions for the SNPs within the region
	my %lifted=();
	while(<$in>){
		chomp;
		my @F=split(' ');
		$lifted{$F[1]}=shift @egg_position;
	}
	close($in);
	print "deCODE\n";
	if($dperson[$i] eq "fetal"){
		open($in,"zcat ".$dfetal." | ") or die $!;
	}elsif($dperson[$i] eq "maternal"){
		open($in,"zcat ".$dmaternal." | ") or die $!;
	}else{
		die "Error: deCODE person not recognised\n";
	}
	open(my $out,">","locus_".$loc."_joined");
	<$in>;
	print $out "Marker\tpos\teaf\tpw_beta\tpw_se\tpw_p\tpw_n\tbw_beta\tbw_p\tbw_n\tbw_eaf\n";
	my $idx=0;
	while(<$in>){
		chomp;
		my @F=split(' ');
		if($F[0] eq "chr".$L[0] && exists $lifted{$F[1]}){
			if(exists $egg{$lifted{$F[1]}}){
				my @a1=sort(lc($F[3]),lc($F[4]));
				if($alleles{$lifted{$F[1]}} eq $a1[0].":".$a1[1]){
					my @G=split(' ',$egg{$lifted{$F[1]}},2);
					print $out join("\t",$F[0].":".$lifted{$F[1]},$lifted{$F[1]},$G[1]);
					$F[5]=$F[5]/100;
					if($G[0] eq lc($F[4])){
						print $out "\t".join("\t",$F[8],$F[9],"321223",$F[5])."\n";
					}else{
						print $out "\t".join("\t",-$F[8],$F[9],"321223",1-$F[5])."\n";
					}
				}else{
					print "Allele error for ".$F[1].". Expected ".$alleles{$lifted{$F[1]}}." got ".$a1[0].":".$a1[1]."\n";
				}
			}
			$idx++;
		}
	}
	close($out);
	close($in);
}
