#!/bin/bash

#SBATCH --job-name=cpdx-paleo
#SBATCH --output=%x-%j.out
#SBATCH --account=wasser
#SBATCH --partition=wasser
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-10:00:00
#SBATCH --mem=10G
#SBATCH --chdir=/gscratch/stf/kimh11/analyses/cpdx-paleomix
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kimh11@uw.edu
#SBATCH --export=all


#Run Paleomix pipeline using each makefile
paleomix bam run makefile-museum.yaml --jre-option=-Xmx20g
paleomix bam run makefile-frozen.yaml --jre-option=-Xmx20g
paleomix bam run makefile-dnazoo.yaml --jre-option=-Xmx20g
paleomix bam run makefile-hu.yaml --jre-option=-Xmx20g
paleomix bam run makefile-tri.yaml --jre-option=-Xmx20g

#Run Picard Downsample in parallel
cat downsample.sh | parallel -j 20

#Run FgBio clipbam in parallel
cat clipbam.sh | parallel -j 2

#Index the new bam files
for file in *x.clipped.bam;
do samtools index $file;
done;


#Run angsd on each species separately

REF="/gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta"
TODO="-doMajorMinor 1 -doMaf 3 -doCounts 1 -doGeno 8 -doBcf 1 -doPost 2 -doGlf 2"


#crassicaudata
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 3 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b crassicaudata.txt -GL 1 $FILTERS $TODO -ref $REF -out ./crass/crassicaudata-all.{/}"

#javanica-17
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 2 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ./SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b javanica.txt -GL 1 $FILTERS $TODO -ref $REF -out ./jav/javanica-all.{/}"

#gig
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 3 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ./SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b gigantea.txt -GL 1 $FILTERS $TODO -ref $REF -out ./gig/gigantea-all.{/}"

#pentadactyla
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 16 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ./SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b pentadactyla.txt -GL 1 $FILTERS $TODO -ref $REF -out ./penta/pentadactyla-all.{/}"

#tri
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 8 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ./SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b tricuspis.txt -GL 1 $FILTERS $TODO -ref $REF -out ./tri/tricuspis-all.{/}"

#tetra
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 5 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b tetradactyla.txt  -GL 1 $FILTERS $TODO -ref $REF -out ./tetra/tetradactyla-all.{/}"

#temm
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 4 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.05"
ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b temminckii.txt -GL 1 $FILTERS $TODO -ref $REF -out ./temm/temminckii-all.{/}"


#Concatenate all the bcf files together
bcftools concat *.bcf --naive -Ob -o javanica-all.bcf
bcftools concat *.bcf --naive -Ob -o pentadactyla-all.bcf
bcftools concat *.bcf --naive -Ob -o crassicaudata-all.all.bcf
bcftools concat *.bcf --naive -Ob -o tricuspis-all.bcf
bcftools concat *.bcf --naive -Ob -o tetradactyla-all.bcf
bcftools concat *.bcf --naive -Ob -o gigantea-all.bcf
bcftools concat *.bcf --naive -Ob -o temminckii-all.bcf

#Convert all bcf files to vcfs
for file in *.bcf;
do echo "bcftools view  -Oz -o ${file%.bcf}.vcf.gz $file";
done

#Remove all transitions
#This uses a script by Ruoyun Hui that was modifed (https://github.com/ryhui/imputation-pipeline/blob/master/pmd_filter.py)
for file in *l.vcf.gz;
do python pmd_filter_edit.py $file | bgzip -c > ./pmd-edited/${file%.vcf.gz}-edit-noTs.vcf.gz;
done


#Filter each file using hard filters
for file in *noTs.vcf.gz;
do vcftools --gzvcf $file --maf 0.03 --max-missing 0.95 --minQ 20 --min-meanDP 3 --max-meanDP 25 --minDP 3 --maxDP 25 --recode --recode-INFO-all --out ${file%.vcf.gz}-maf0.03-maxmiss0.95-minQ20-minmeanDP3-maxmeanDP25-minDP3-maxDP25;
done

#Calculate He
for file in *noTs-maf0.03-maxmiss0.95-minQ20-minmeanDP3-maxmeanDP25-minDP3-maxDP25.recode.vcf;
do vcftools --vcf $file --maf 0.05 --het --out ${file%recode.vcf.gz};
done

#Calculate pi
for file in *noTs-maf0.03-maxmiss0.95-minQ20-minmeanDP3-maxmeanDP25-minDP3-maxDP25.recode.vcf;
do vcftools --vcf $file --window-pi 10000 --out ${file%recode.vcf.gz};
done


#Convert vcf to plink files
for file in *.recode.vcf;
do plink --vcf $file --recode --out ${file%.recode.vcf} --allow-extra-chr;
done

#Calculate roh
for file in *ped;
do plink --file ${file%.ped} --homozyg --allow-extra-chr --homozyg-snp 20 --homozyg-kb 200 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 20 --homozyg-window-het 5 --homozyg-window-missing 5 --homozyg-window-threshold 0.01  --out test-param3-${file%.ped};
done




