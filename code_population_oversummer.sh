module load apps/samtools/1.16.1-gnu485
samtools view -@ 20 -ubhSt CYR34_Chr_18.fasta.fai add_10.sam > add_10.bam
##CleanSam
java -Xmx128G -jar picard-2.27.5.jar CleanSam INPUT=add_10.bam OUTPUT=add_10.clean.bam VALIDATION_STRINGENCY=SILENT
##SortSam
java -Xmx128G -jar picard-2.27.5.jar SortSam INPUT=add_10.clean.bam OUTPUT=add_10.clean.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
##MarkDuplicates
java -Xmx128G -jar picard-2.27.5.jar MarkDuplicates INPUT=add_10.clean.sorted.bam OUTPUT=add_10.final.bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE METRICS_FILE=./DUP_METRICS.104E.OUT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
##samtools index
samtools index add_10.final.bam

#GATK 1 round 
gatk --java-options "-Xmx128g" HaplotypeCaller --native-pair-hmm-threads 48 -R /public/home/liyx/yuxiang.li/001_database/01_Pst_genome/CYR34/CYR34_Chr_18.fasta -I add_10.final.bam -O add_10.snps.indels.vcf
vcftools --vcf add_10.snps.indels.vcf --min-alleles 2 --max-alleles 2 --minQ 1000 --recode --recode-INFO-all --out add_10.snps.indels.gold

#GATK 2 round
gatk --java-options "-Xmx128g" IndexFeatureFile -I add_10.snps.indels.gold.recode.vcf
gatk --java-options "-Xmx128g" BaseRecalibrator -R /public/home/liyx/yuxiang.li/001_database/01_Pst_genome/CYR34/CYR34_Chr_18.fasta -I /public/home/liyx/yuxiang.li/002_primer_design/bam_file/add_10.final.bam --known-sites add_10.snps.indels.gold.recode.vcf -O add_10.recal_data.table
gatk --java-options "-Xmx128g" ApplyBQSR -R /public/home/liyx/yuxiang.li/001_database/01_Pst_genome/CYR34/CYR34_Chr_18.fasta -I /public/home/liyx/yuxiang.li/002_primer_design/bam_file/add_10.final.bam --bqsr-recal-file add_10.recal_data.table -O add_10.recal.bam
gatk --java-options "-Xmx128g" HaplotypeCaller --native-pair-hmm-threads 64 -R /public/home/liyx/yuxiang.li/001_database/01_Pst_genome/CYR34/CYR34_Chr_18.fasta --emit-ref-confidence GVCF -I add_10.recal.bam -O add_10.g.vcf

#CombineGVCF
gatk CombineGVCFs --java-options "-Xmx256g" -O 28_combine.g.vcf.gz -R CYR34_Chr_18.fasta -V 03_07.g.vcf -V 104E.g.vcf -V 11_140.g.vcf -V 134E.g.vcf -V 15_34.g.vcf -V 93_210.g.vcf -V DK_0911.g.vcf -V GS08.g.vcf -V J0085F.g.vcf -V j02-022.g.vcf -V PST-78.g.vcf -V QH06.g.vcf -V Race_31.g.vcf -V Race_Yr9.g.vcf -V SC01.g.vcf -V TB03.g.vcf -V add_1.g.vcf -V add_2.g.vcf -V add_3.g.vcf -V add_4.g.vcf -V add_5.g.vcf -V add_6.g.vcf -V add_7.g.vcf -V add_8.g.vcf -V add_9.g.vcf -V add_10.g.vcf -V add_11.g.vcf -V add_12.g.vcf

#Genotype
gatk GenotypeGVCFs --java-options "-Xmx512g" -O 28_genotype.snps.indels.vcf -R /public/home/liyx/yuxiang.li/001_database/01_Pst_genome/CYR34/CYR34_Chr_18.fasta -V 28_combine.g.vcf.gz

#Filtering
1. vcftools --vcf 28_genotype.snps.indels.vcf --recode --recode-INFO-all --remove-indels --min-meanDP 20 --max-meanDP 1000 --min-alleles 2 --max-alleles 2 --out 1

2. /public/home/liyx/yuxiang.li/005_plink/plink --allow-extra-chr --recode --chr-set 18 no-xy no-mt --vcf 1.recode.vcf --out 2.vcf
awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4}' 2.vcf.map > 2.vcf.pos.map
mv 2.vcf.ped 2.vcf.pos.ped
/public/home/liyx/yuxiang.li/005_plink/plink --file 2.vcf.pos --indep-pairwise 50 1 0.2 --out 3 --allow-extra-chr --chr-set 18 no-xy no-mt  
sed 's/_/\t/g' 3.prune.in > 3_selected.in
bgzip 1.recode.vcf
tabix 1.recode.vcf.gz
bcftools view -R 3_selected.in  1.recode.vcf.gz > 4.vcf

3. gatk VariantFiltration -V 4.vcf -O 5.vcf --genotype-filter-expression "isHet == 1" --genotype-filter-name "Het_site"
cat 5.vcf| grep -v 'Het_site' > 6.vcf

4. vcftools --vcf 6.vcf --recode --recode-INFO-all  --mac 3 --maf 0.1 --minQ 100 --max-missing 0.9 --out 7

5. vcftools --vcf 7.recode.vcf  --recode --recode-INFO-all --thin 1000  --out 9