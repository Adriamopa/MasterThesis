#!/bin/bash

wd="$(pwd)"
pops="AFR EUR EAS"
invs=$(cat "$wd/../VCFs/30X/Invs2Imp")

> $wd/Results_Stats/MAF_Freq_AFR.txt
> $wd/Results_Stats/MAF_Freq_EUR.txt
> $wd/Results_Stats/MAF_Freq_EAS.txt

for inv in $invs
do
	for pop in $pops
	do
		state=$(bcftools view -H "$wd/../VCFs/Inversion_Regions/Unphased/$inv/$pop.vcf.gz" | grep $inv | cut -f 10- | sed 's_/_|_g' | tr "|" "\t" | tr "\t" "\n" | sort | uniq -c)
		num=$(echo "$state" | wc -l)
		if [ $num -gt 1 ]
		then
			freqs=($(echo "$state" | awk '{print $1}'))
			inver=${freqs[1]}
			no_inv=${freqs[0]}
			total_freq=$(LC_NUMERIC=C printf %.4f $(echo "scale=10; $inver / $((inver + no_inv))" | bc -l))		
			if [ $(echo "$total_freq > 0.5" |bc -l) -eq 1 ]
			then
				total_freq=$(echo "1 - $total_freq" |bc -l)
			fi
			echo -e "$inv: $total_freq\n"
			echo -e "$inv\t$total_freq" >> $wd/Results_Stats/MAF_Freq_$pop.txt
		elif [ $num -eq 1 ]
		then
			freq=($(echo "$state" | awk '{print $2}'))
			echo -e "$inv: $freq\n"
			echo -e "$inv\t$freq" >> $wd/Results_Stats/MAF_Freq_$pop.txt
		else
			echo -e "Error in Inversion $inv and Population $Pop!!\n"
		fi
	done
done
