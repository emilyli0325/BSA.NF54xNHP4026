### '0/1' to non-call
# HetFilter

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -V Parent.12Samples.SNP.di.ann.vcf \
   --sample_name AB-BSA-220-NF54/ --sample_name AB-BSA-222-NHP4026/ --sample_name MKK2835.2G/ --sample_name NHP1337.12C/ \
   -o NF54.NHP4026.MKK2835.NHP1337.SNP.vcf 

   
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
--variant NF54.NHP4026.MKK2835.NHP1337.SNP.vcf \
--genotypeFilterExpression "isHet == 1" \
--genotypeFilterName "HetFilter" \
--setFilteredGtToNocall \
--out NF54.NHP4026.MKK2835.NHP1337.SNP.rmHet.vcf



java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V NF54.NHP4026.MKK2835.NHP1337.SNP.rmHet.vcf \
-F CHROM -F POS -F REF -F ALT \
-GF GT -GF AD -GF DP -GF GQ -GF PL  \
-o NF54.NHP4026.MKK2835.NHP1337.SNP.rmHet.table 



########################################################################################## 
##########################################################################################  
hmmIBD

setwd("D:/Dropbox (TX Biomed)/Emily/P01/5.5.BSA5_NF54xNHP4026_all/hmmIBD")

## data clean
Input <- read.delim("NF54.NHP4026.MKK2835.NHP1337.SNP.rmHet.table", sep="\t",header=T)

Input[Input == "./."]<-NA
Input<- Input[complete.cases(Input), ]
> dim(Input)
[1] 54746    24

GT <- Input[,seq(5,24,5)]

> dim(GT)
[1] 54746     4

GT.A <- function(X){unlist(strsplit(as.character(X),"\\/"))[1]}

GT.Aallele <- matrix(ncol=4,nrow=54746)
for(j in 1:4){
GT.Aallele[,j] <- sapply(GT[,j], GT.A, simplify="array")
}

GT.Aallele[GT.Aallele== "."]<-NA
colnames(GT.Aallele)<- colnames(GT)

GT.recode<-NULL

REF <- as.character(Input$REF)
ALT <- as.character(Input$ALT)

for(i in c(1:4)){
    X<-as.character(GT.Aallele[,i])
    Z<-NULL; 
	Z[which(X==REF & is.na(X)==F)]<- 0; 
	Z[which(X==ALT & is.na(X)==F)]<- 1; 
	Z[which(is.na(X)==T)]<- "-1";
    GT.recode <- cbind(GT.recode,Z)
}

GT.recode <- cbind(Input[,1:2],GT.recode)

colnames(GT.recode) <- c("chrom", "pos",gsub("\\.GT|\\.\\.GT","",colnames(GT)[1:4]))

GT.recode.cleanChr <- GT.recode
GT.recode.cleanChr$chrom <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",GT.recode$chrom)

GT.recode.cleanChr<- GT.recode.cleanChr[complete.cases(GT.recode.cleanChr), ]

write.table(GT.recode.cleanChr,"HmmIBD.Input.txt",row.names=F,col.names=T,sep="\t", quote=FALSE)



### hmmIBD: find crossovers
/master/xli/software/hmmIBD/hmmIBD-master/hmmIBD -i HmmIBD.Input.txt -o hmmIBD



