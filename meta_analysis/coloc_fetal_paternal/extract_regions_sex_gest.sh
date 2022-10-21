#!/bin/bash

file="sex_gest_loci_pairs"
lines=`wc -l $file | awk '{print $1}'`
for i in $(seq 1 ${lines})
do
  line=`head -n $i $file | tail -n 1`
  IFS=$' ' read -r -a locus <<< $line
# person1
  if [[ ${locus[3]} == "Fetal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_fetal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person1_${i}
  elif [[ ${locus[3]} == "Maternal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_maternal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person1_${i}
  elif [[ ${locus[3]} == "Paternal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_paternal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person1_${i}
  else
    echo "Error: "${locus[3]}
  fi
# person2
  if [[ ${locus[7]} == "Fetal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_fetal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person2_${i}
  elif [[ ${locus[7]} == "Maternal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_maternal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person2_${i}
  elif [[ ${locus[7]} == "Paternal" ]]
  then
    zcat ../../frozen_meta_analysis_files/pw_paternal_sex1.tbl.gz | awk -v chr=${locus[0]} -v pos=${locus[1]} 'NR==1{print "CHR\tPOS\teaf\tbeta\tse\tp\tn\tea\toa";next} NR>1{split($1,F,":");if(F[1]==chr && F[2]>(pos-1000000) && F[2]<(pos+1000000)){print F[1]"\t"F[2]"\t"$4"\t"$8"\t"$9"\t"$10"\t"$16"\t"$2"\t"$3}}' > PW_gest_locus_person2_${i}
  else
    echo "Error: "${locus[7]}
  fi
done
