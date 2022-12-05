paste -d ' ' ../configs/na_inclafr.samples  ../configs/na_inclafr.pops ../res_inclafr/admix/K3_s1/sub0.q > names_and_Q.txt
awk '$3>0.95 {print $1"\tAFR"}' names_and_Q.txt > AFR.map
awk '$5>0.95 {print $1"\tEUR"}' names_and_Q.txt > EUR.map
awk '$4>0.95 {print $1"\tNA"}' names_and_Q.txt > NA.map
cat  AFR.map EUR.map NA.map > REF.map
~/software/bcftools/bcftools-1.16/bin/bcftools view --threads 20 -S <(cut -f 1 REF.map) -o REF.vcf.gz -O z chr1.vcf.gz
diff <(bcftools query -l REF.vcf.gz ) <(cut -f 1 REF.map )
~/software/bcftools/bcftools-1.16/bin/bcftools view --threads 20 -s HG02006,HG02089 -o ADM.vcf.gz -O z chr1.vcf.gz
/home/krishang/software/rfmix/rfmix/rfmix --reanalyze-reference \
	-e 10 -f ADM.vcf.gz -r REF.vcf.gz -m REF.map \
	-g hapmap-phase2-genetic-map.tsv -o res --n-threads=10 --chromosome=chr1
