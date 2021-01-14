#combine g.vcf
java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
--variant AB-BSA-222-NHP4026.g.vcf \
--variant AB_BSA_221.g.vcf \
--variant AB_BSA_222.g.vcf \
--variant AB_BSA_223.g.vcf \
--variant AB_BSA_225.g.vcf \
--variant AB_BSA_229.g.vcf \
--variant AB_BSA_230.g.vcf \
--variant AB_BSA_231.g.vcf \
--variant AB_BSA_233.g.vcf \
--variant AB_BSA_237.g.vcf \
--variant AB_BSA_238.g.vcf \
--variant AB_BSA_239.g.vcf \
--variant AB_BSA_244.g.vcf \
--variant AB_BSA_245.g.vcf \
--variant AB_BSA_246.g.vcf \
--variant AB_BSA_247.g.vcf \
--variant AB_BSA_248.g.vcf \
--variant AB_BSA_249.g.vcf \
--variant AB_BSA_251.g.vcf \
--variant AB_BSA_252.g.vcf \
--variant AB_BSA_253.g.vcf \
--variant AB_BSA_254.g.vcf \
--variant AB_BSA_255.g.vcf \
--variant AB_BSA_256.g.vcf \
--variant AB_BSA_257.g.vcf \
--variant AB_BSA_258.g.vcf \
--variant AB_BSA_259.g.vcf \
--variant AB_BSA_260.g.vcf \
--variant AB_BSA_261.g.vcf \
--variant AB_BSA_262.g.vcf \
--variant AB_BSA_263.g.vcf \
--variant AB_BSA_264.g.vcf \
--variant AB_BSA_265.g.vcf \
--variant AB_BSA_266.g.vcf \
--variant AB_BSA_267.g.vcf \
-o BSA5.NF54-NHP4026.g.vcf

# genotype
java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
-V BSA5.NF54-NHP4026.g.vcf \
--useNewAFCalculator \
--sample_ploidy 1 \
-nt 10 \
-o BSA5.NF54-NHP4026.vcf



## SelectVariants
java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
-V BSA5.NF54-NHP4026.vcf \
-nt 10 \
-selectType SNP -o BSA5.NF54-NHP4026.SNP.vcf

# find SNPs called by both parents
java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/NF54ref/NF54.SNP.fasta \
   -V BSA5.NF54-NHP4026.SNP.vcf \
   --concordance /data/infectious/malaria_XUE/cross/BSA2_lifecycle_drug/NF54_NHP4026/P.genotype/NF54_NHP4026.VQSR.SNP.filter.clean.vcf \
   --restrictAllelesTo BIALLELIC -nt 10 \
   -o BSA5.NF54-NHP4026.SNP.filter.vcf

# vcf to table  

java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /data/infectious/malaria_XUE/cross/Parents/mergeAllparents/MAL31.REF/Mal31.SNP.fasta \
-V BSA5.NF54-NHP4026.SNP.filter.vcf \
-L /master/xli/Index/Known_sites/Core_Genome2.intervals \
-F CHROM -F POS -F REF -F ALT \
-GF AD -GF DP -GF GQ -GF PL \
-o BSA5.NF54-NHP4026.SNP.filter_coreGenome.table




