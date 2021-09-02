#Run freebayes

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/crass-georef.txt"

/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > crass-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/gig-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > gig-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/jav-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > jav-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/pent-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > pent-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/temm-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > temm-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/tetra-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > tetra-freebayes-geo-ref.vcf

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/tri-georef.txt"
/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/freebayes-parallel <(/gscratch/wasser/kimh11/conda3/envs/freebayesenv/bin/fasta_generate_regions.py /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta 50000) 25 \ --gvcf --limit-coverage 200 --use-reference-allele --fasta-reference /gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --bam-list $BAM > tri-freebayes-geo-ref.vcf

#remove indels
for file in *.vcf;
do vcftools --vcf $file --remove-indels --recode --recode-INFO-all --out ${file%.vcf}_SNPs_only;
done

for file in *recode.vcf; do bgzip $file; done

for file in *.vcf.gz; do bcftools index $file; done

#Run angsd

REF="/gscratch/wasser/kimh11/reference-genomes/ManJav2.0/GCF_014570535.1_YNU_ManJav_2.0_genomic.fasta"
#crass
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 3 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.2"

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/crass-georef.txt"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b $BAM -GL 1 $FILTERS $TODO -ref $REF -out ./crass/crassicaudata-georef.{/}"


#gig
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 3 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.2"

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/gig-georef.txt"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b $BAM -GL 1 $FILTERS $TODO -ref $REF -out ./gig/gigantea-georef.{/}"


#javanica
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 11 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.08"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b /gscratch/wasser/kimh11/reference-genomes/sample-list/jav-georef.txt -GL 1 $FILTERS $TODO -ref $REF -out ./javanica/javanica-georef.{/}"


#pent
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 8 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.11"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b /gscratch/wasser/kimh11/reference-genomes/sample-list/pent-georef.txt -GL 1 $FILTERS $TODO -ref $REF -out ./penta/pentadactyla-georef.{/}"


#temm
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 4 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.19"

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/temm-georef.txt"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b $BAM -GL 1 $FILTERS $TODO -ref $REF -out ./temm/temmickii-georef.{/}"


#tetra
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 5 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.16"

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/tetra-georef.txt"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b $BAM -GL 1 $FILTERS $TODO -ref $REF -out ./tetra/tetradactyla-georef.{/}"


#tri
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 7 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -baq 1 -snp_pval 1e-6 -minMaf 0.1"

BAM="/gscratch/wasser/kimh11/reference-genomes/sample-list/tri-georef.txt"

ls ../SPLIT_RF/* | parallel -j $SLURM_CPUS_PER_TASK "angsd -rf {} -b $BAM -GL 1 $FILTERS $TODO -ref $REF -out ./tri/tridactyla-georef.{/}"

#Run intersect

JAVFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/jav_SNPs_only.recode.vcf.gz
JAVANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/javanica/javanica-georef-angsd.vcf.gz

bcftools isec $JAVFREEBAYES $JAVANGSD -p jav --threads $SLURM_CPUS_PER_TASK --output-type z


PENTFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/pent_SNPs_only.recode.vcf.gz
PENTANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/penta/pentadactyla-georef-angsd.vcf.gz

bcftools isec $PENTFREEBAYES $PENTANGSD -p pent --threads $SLURM_CPUS_PER_TASK --output-type z

TRIFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/tri_SNPs_only.recode.vcf.gz
TRIANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/tri/tridactyla-georef-angsd.vcf.gz

bcftools isec $TRIFREEBAYES $TRIANGSD -p tri --threads $SLURM_CPUS_PER_TASK --output-type z

TEMMFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/temm_SNPs_only.recode.vcf.gz
TEMMANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/temm/temminckii-georef-angsd.vcf.gz

bcftools isec $TEMMFREEBAYES $TEMMANGSD -p temm --threads $SLURM_CPUS_PER_TASK --output-type z

GIGFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/gig_SNPs_only.recode.vcf.gz
GIGANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/gig/gigantea-georef-angsd.vcf.gz

bcftools isec $GIGFREEBAYES $GIGANGSD -p gig --threads $SLURM_CPUS_PER_TASK --output-type z

TETRAFREEBAYES=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-freebayes/tetra_SNPs_only.recode.vcf.gz
TETRAANGSD=/mmfs1/gscratch/wasser/kimh11/analysis/final-bams/species/species-georef-angsd/tetra/tetradactyla-georef-angsd.vcf.gz

bcftools isec $TETRAFREEBAYES $TETRAANGSD -p tetra --threads $SLURM_CPUS_PER_TASK --output-type z

#Filter intersect files
for file in */0002.vcf.gz;
do realpath $file; vcftools --gzvcf $file --maf 0.03 --max-missing 0.95 --minQ 20 --min-meanDP 5 --max-meanDP 25 --minDP 5 --maxDP 25 --recode --recode-INFO-all --out ${file}-filtered};
done


#Find common SNPs
bcftools isec crass/0002.vcf.gz jav/0002.vcf.gz pent/0002.vcf.gz temm/0002.vcf.gz tri/0002.vcf.gz tetra/0002.vcf.gz gig/0002.vcf.gz -p all-in3 --nfiles +2 -c all --threads $SLURM_CPUS_PER_TASK --output-type z

#Thin SNPs
vcftools --vcf combined-allin2-bcf-sorted --thin 20000 --out combined-allin2-bcf-sorted-thin10k

#Run tres and SCAT

