
echo "mapping ${1}..."

/master/xli/software/bwa3/bwa-0.7.15/bwa mem /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta /master/xli/raw.data/BSA5_NF54xNHP4026_lifecycle/pre.blood/19047-10/${1}.R1.fastq.gz /master/xli/raw.data/BSA5_NF54xNHP4026_lifecycle/pre.blood/19047-10/${1}.R2.fastq.gz -t 12 -M -R "@RG\tID:${1}\tLB:${1}\tPL:ILLUMINA\tPM:HISEQ\tSM:${1}/" > SAM/${1}.sam

java -jar /master/xli/software/picard/picard.jar SortSam \
     INPUT=SAM/${1}.sam \
     OUTPUT=sorted.bam/${1}.sorted.bam \
     SORT_ORDER=coordinate

java -jar /master/xli/software/picard/picard.jar MarkDuplicates \
     INPUT=sorted.bam/${1}.sorted.bam \
     OUTPUT=dedup.sorted.bam/${1}.dedup.sorted.bam \
     METRICS_FILE=metrics/${1}.metrics.txt

cd dedup.sorted.bam

java -jar /master/xli/software/picard/picard.jar BuildBamIndex \
     INPUT=${1}.dedup.sorted.bam

cd ..

echo "BQSR $file..."	 
	java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar\
		 -T BaseRecalibrator\
		 -R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
		 -I dedup.sorted.bam/${1}.dedup.sorted.bam\
		 -knownSites /master/xli/Index/Known_sites/3d7_hb3.combined.final.karo.sort.vcf\
		 -knownSites /master/xli/Index/Known_sites/7g8_gb4.combined.final.karo.sort.vcf\
		 -knownSites /master/xli/Index/Known_sites/hb3_dd2.combined.final.karo.sort.vcf\
		 -o BQSR/${1}.recal.table
	java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar\
		 -T PrintReads\
		 -R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
		 -I dedup.sorted.bam/${1}.dedup.sorted.bam\
		 -BQSR BQSR/${1}.recal.table\
		 -o recal.bam/${1}.recal.bam 
		 
echo "Variant calling ${1}..."
	java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar\
         -T HaplotypeCaller\
         -R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
         -I recal.bam/${1}.recal.bam\
         --emitRefConfidence BP_RESOLUTION\
         -o g.vcf/${1}.g.vcf\
         --useNewAFCalculator \
         --sample_ploidy 1 \
         --dontUseSoftClippedBases \
         --num_cpu_threads_per_data_thread 16
        echo "${1} Variant calling done!"
