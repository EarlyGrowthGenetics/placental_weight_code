
#!/bin/bash

#the top two lines create a chr:pos column for the summary stats to be filtered on
#zcat /child_gest/LDSC_Files/Maternal//meta_results/pw_fetal_sex_gest.gz|awk '{print $17, $18}' OFS=":"  > new_column.txt

#zcat /child_gest/LDSC_Files/Maternal//meta_results/pw_fetal_sex_gest.gz| paste -d '\t' - new_column.txt > new_meta.txt
#gzip new_meta.txt

zcat /child_gest/LDSC_Files/Maternal/Clean/EGG_HRC_BW6.PW.mother.sex_gest.RaineStudy.european.CW.20200616.txt.gz|head -n 1 > /child_gest/LDSC_Files/Maternal/top_hits_files/Raine.txt

zgrep -f /child_gest/LDSC_Files/maternal_top.csv /child_gest/LDSC_Files/Maternal/Clean/EGG_HRC_BW6.PW.mother.sex_gest.RaineStudy.european.CW.20200616.txt.gz >> /child_gest/LDSC_Files/Maternal/top_hits_files/Raine.txt

