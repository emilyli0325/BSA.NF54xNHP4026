setwd("C:/Users/xli/Documents/Emily/P01/5.1.BSA5_NF54xNHP4026_lifecycle/3.results")

library(magrittr)
library(dplyr)
library("QTLseqr")
library("ggpubr")

format_genomic <- function(...) {
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            #prefix <- c("","Kb","Mb")
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...)
                #  ,prefix[i]
            )
      }
}

BSA <- read.table("BSA5.NF54-NHP4026.SNP.filter_coreGenome.table", sep = '\t',header = TRUE)
> dim(BSA)
[1] 15114   144

AD <- BSA[,seq(5,144,4)]
DP <- BSA[,seq(6,144,4)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)
[1] 15114    35

refFre.AD <- matrix(ncol=35,nrow=15114)
for(j in 1:35){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}
refFre.AD[DP<20]<-NA
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(BSA[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub("..AD", "", colnames(refFre.AD))


######################################################
######### remove outliers ############################
######################################################
outliersMAD <- function(data, MADCutOff = 2, replace = NA, values = FALSE, bConstant = 1.4826, digits = 2) {
    absMADAway <- abs(   (data - median(data, na.rm = T))  /  mad(data, constant = bConstant, na.rm = T)  )
   
    data[absMADAway > MADCutOff] <- replace
    
    if (values == TRUE) { 
        return(round(absMADAway, digits)) #if values == TRUE, return number of mads for each value
    } else {
        return(round(data, digits)) #otherwise, return values with outliers replaced
    }
}

outlierByMAD <- function (x, k){
	n <- length(x)
     y <- x
     
 for (i in (k + 1):(n - k)) {
 	
    data <- x[(i - k):(i + k)]
    y[i] <- outliersMAD(data)[k+1]}
    
   return(y)
                               }

AB.BSA.222.NHP4026.filter <- outlierByMAD(refFre.AD$AB.BSA.222.NHP4026,50)
AB_BSA_221.filter <- outlierByMAD(refFre.AD$AB_BSA_221,50)
AB_BSA_222.filter <- outlierByMAD(refFre.AD$AB_BSA_222,50)
AB_BSA_223.filter <- outlierByMAD(refFre.AD$AB_BSA_223,50)
AB_BSA_225.filter <- outlierByMAD(refFre.AD$AB_BSA_225,50)
AB_BSA_229.filter <- outlierByMAD(refFre.AD$AB_BSA_229,50)
AB_BSA_230.filter <- outlierByMAD(refFre.AD$AB_BSA_230,50)
AB_BSA_231.filter <- outlierByMAD(refFre.AD$AB_BSA_231,50)
AB_BSA_233.filter <- outlierByMAD(refFre.AD$AB_BSA_233,50)
AB_BSA_237.filter <- outlierByMAD(refFre.AD$AB_BSA_237,50)
AB_BSA_238.filter <- outlierByMAD(refFre.AD$AB_BSA_238,50)
AB_BSA_239.filter <- outlierByMAD(refFre.AD$AB_BSA_239,50)
AB_BSA_244.filter <- outlierByMAD(refFre.AD$AB_BSA_244,50)
AB_BSA_245.filter <- outlierByMAD(refFre.AD$AB_BSA_245,50)
AB_BSA_246.filter <- outlierByMAD(refFre.AD$AB_BSA_246,50)
AB_BSA_247.filter <- outlierByMAD(refFre.AD$AB_BSA_247,50)
AB_BSA_248.filter <- outlierByMAD(refFre.AD$AB_BSA_248,50)
AB_BSA_249.filter <- outlierByMAD(refFre.AD$AB_BSA_249,50)
AB_BSA_251.filter <- outlierByMAD(refFre.AD$AB_BSA_251,50)
AB_BSA_252.filter <- outlierByMAD(refFre.AD$AB_BSA_252,50)
AB_BSA_253.filter <- outlierByMAD(refFre.AD$AB_BSA_253,50)
AB_BSA_254.filter <- outlierByMAD(refFre.AD$AB_BSA_254,50)
AB_BSA_255.filter <- outlierByMAD(refFre.AD$AB_BSA_255,50)
AB_BSA_256.filter <- outlierByMAD(refFre.AD$AB_BSA_256,50)
AB_BSA_257.filter <- outlierByMAD(refFre.AD$AB_BSA_257,50)
AB_BSA_258.filter <- outlierByMAD(refFre.AD$AB_BSA_258,50)
AB_BSA_259.filter <- outlierByMAD(refFre.AD$AB_BSA_259,50)
AB_BSA_260.filter <- outlierByMAD(refFre.AD$AB_BSA_260,50)
AB_BSA_261.filter <- outlierByMAD(refFre.AD$AB_BSA_261,50)
AB_BSA_262.filter <- outlierByMAD(refFre.AD$AB_BSA_262,50)
AB_BSA_263.filter <- outlierByMAD(refFre.AD$AB_BSA_263,50)
AB_BSA_264.filter <- outlierByMAD(refFre.AD$AB_BSA_264,50)
AB_BSA_265.filter <- outlierByMAD(refFre.AD$AB_BSA_265,50)
AB_BSA_266.filter <- outlierByMAD(refFre.AD$AB_BSA_266,50)
AB_BSA_267.filter <- outlierByMAD(refFre.AD$AB_BSA_267,50)


refFre.AD.LC <- cbind(refFre.AD[,1:4], AB.BSA.222.NHP4026.filter, AB_BSA_221.filter,AB_BSA_222.filter,AB_BSA_223.filter,AB_BSA_225.filter,AB_BSA_229.filter,AB_BSA_230.filter,AB_BSA_231.filter,AB_BSA_233.filter,AB_BSA_237.filter,AB_BSA_238.filter,AB_BSA_239.filter,AB_BSA_244.filter,AB_BSA_245.filter,AB_BSA_246.filter,AB_BSA_247.filter,AB_BSA_248.filter,AB_BSA_249.filter,AB_BSA_251.filter,AB_BSA_252.filter,AB_BSA_253.filter,AB_BSA_254.filter,AB_BSA_255.filter,AB_BSA_256.filter,AB_BSA_257.filter,AB_BSA_258.filter,AB_BSA_259.filter,AB_BSA_260.filter,AB_BSA_261.filter,AB_BSA_262.filter,AB_BSA_263.filter,AB_BSA_264.filter,AB_BSA_265.filter,AB_BSA_266.filter,AB_BSA_267.filter)

> dim(refFre.AD.LC)
[1] 15114    39


################################################################  
###### compare whole genome allele frequency  ##################  
################################################################  
library(reshape2)
library(ggplot2)
library(ggridges)
refFre.LC <- setNames(melt(refFre.AD.LC[,5:39]), c('BSAs', 'Allele frequency of NF54'))

library(doBy)
colnames(refFre.LC)[2]<- c("AlleleFrequency")
refFre.LC.filter <- refFre.LC[rowSums(is.na(refFre.LC)) == 0,]
sum <- summaryBy(AlleleFrequency ~ BSAs, data = refFre.LC.filter, FUN = list(mean, median))


write.csv(sum, file = "Whole.Genome.NF54-AF.sum_coreGenome.csv",row.names=FALSE)



ggplot(refFre.LC.filter[refFre.LC.filter$BSAs=="AB_BSA_221.filter",], aes(x=AlleleFrequency))+ geom_density(color="black", fill="lightblue") + theme_classic()+ theme(axis.text.x = element_text(size=10, angle=90),legend.position = c(0.8, 0.2)) + ylab("Density")+ xlab("NF54 allele frequency") + xlim(0,1)

ggplot(refFre.AD.LC, aes(x=AB_BSA_221.tricube))+ geom_density(color="black", fill="lightblue") + theme_classic()+ theme(axis.text.x = element_text(size=10, angle=90),legend.position = c(0.8, 0.2)) + ylab("Density")+ xlab("NF54 allele frequency") + xlim(0,1)


################################################################
######################## Allele frequency ######################
################################################################

#######tricubesmooth#######		
tricubeStat<-function(POS,Stat,windowSize=1e5,...)
{
if(windowSize<=0)
stop("Apositivesmoothingwindowisrequired")
stats::predict(locfit::locfit(Stat~locfit::lp(POS,h=windowSize,deg=0),...),POS)
}


SAnalysis.1<-function(SNPset,windowSize=1e5,...)
{SNPset<-SNPset%>%
dplyr::group_by(CHROM)%>%

dplyr::mutate(AB_BSA_267.tricube =tricubeStat(POS=POS,Stat=AB_BSA_267.filter,windowSize,...))


return(as.data.frame(SNPset))};refFre.AD.LC<-SAnalysis.1(refFre.AD.LC)	



> dim(refFre.AD.LC)
[1] 15114    74

write.csv(refFre.AD.LC, file = "BSA5.AF.filter.tricube_coregenome.csv",row.names=FALSE)

refFre.AD.LC <- read.csv(file = "BSA5.AF.filter.tricube_coregenome.csv",header = TRUE)


refFre.AD.LC$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",refFre.AD.LC$CHROM)
refFre.AD.LC$CHROM <- as.numeric(as.character(refFre.AD.LC$CHROM))
refFre.AD.LC <- refFre.AD.LC[!is.na(refFre.AD.LC$AB.BSA.222.NHP4026.filter),]
> dim(refFre.AD.LC)
[1] 13450    74

################################################################
# all

p0 <- ggplot(data=refFre.AD.LC)+scale_x_continuous(breaks=seq(from=0,to=max(refFre.AD.LC$POS),by=10^(floor(log10(max(refFre.AD.LC$POS))))),labels=format_genomic())+ylim(0,1)+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()

p0 + geom_point(aes_string(x = "POS", y = "AB_BSA_248.filter"),color = "orange",size=0.5) + geom_line(aes_string(x = "POS", y = "AB_BSA_221.tricube"),color = "black",size=1)

p221 <- p0 + geom_point(aes_string(x = "POS", y = "AB_BSA_221.filter"),color = "orange",size=0.5) + geom_line(aes_string(x = "POS", y = "AB_BSA_221.tricube"),color = "black",size=1)


pdf('cage2.pdf', width=10, height=18)
ggarrange(p221,p229,p237,p244,p248,p252,p253,p254,p255,
          labels = c("Oocyst d4", "Oocyst d10", "spz", "Liver","in vivo", "serum 1","serum 2","Albumax 1","Albumax 2"),
          ncol = 1, nrow = 9)
dev.off()

pdf('cage3.pdf', width=10, height=18)
ggarrange(p222, p230, p238, p245, p249, p256, p257, p258, p259,
          labels = c("Oocyst d4", "Oocyst d10", "spz", "Liver","in vivo", "serum 1","serum 2","Albumax 1","Albumax 2"),
          ncol = 1, nrow = 9)
dev.off()

pdf('cage4.pdf', width=10, height=18)
ggarrange(p223, p231, p239, p246, p260, p261, p262, p263,
          labels = c("Oocyst d4", "Oocyst d10", "spz", "Liver", "serum 1","serum 2","Albumax 1","Albumax 2"),
          ncol = 1, nrow = 8)
dev.off()

pdf('cage6.pdf', width=10, height=18)
ggarrange(p225, p233, p247, p251, p264, p265, p266, p267,
          labels = c("Oocyst d4", "Oocyst d10", "Liver","in vivo", "serum 1","serum 2","Albumax 1","Albumax 2"),
          ncol = 1, nrow = 8)
dev.off()

################################################################

p.Oo10_invivo_invitro4 <- p0 + ylab("NF54 allele frequency") + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_229.tricube"),color = "black",size=1) + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_230.tricube"),color = "black",size=1) + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_233.tricube"),color = "black",size=1) + 		  
		  geom_line(aes_string(x = "POS", y = "AB_BSA_248.tricube"),color = "orange",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_249.tricube"),color = "orange",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_251.tricube"),color = "orange",size=1) + 	
		  geom_line(aes_string(x = "POS", y = "AB_BSA_260.tricube"),color = "green",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_256.tricube"),color = "green",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_264.tricube"),color = "green",size=1) 
		  		  
pdf('p.Oo10_invivo_invitro4_3.pdf', width=10, height=3)
p.Oo10_invivo_invitro4
dev.off()
		  
################################################################		  
#Cross NF54xNHP4026
p.serum_Albumax <- p0 + ylab("NF54 allele frequency") + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_256.tricube"),color = "black",size=1) + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_260.tricube"),color = "black",size=1) + 
		  geom_line(aes_string(x = "POS", y = "AB_BSA_264.tricube"),color = "black",size=1) + 		  
		  geom_line(aes_string(x = "POS", y = "AB_BSA_258.tricube"),color = "orange",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_262.tricube"),color = "orange",size=1) +
		  geom_line(aes_string(x = "POS", y = "AB_BSA_266.tricube"),color = "orange",size=1)
		  
pdf('p.serum_Albumax.pdf', width=10, height=3)
p.serum_Albumax
dev.off()



################################################################

dfcage2.invivo <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.filter_coreGenome2.table",highBulk = "AB_BSA_229",lowBulk = "AB_BSA_248")

df_filt <-filterSNPs(SNPset = dfcage2.invivo,minTotalDepth = 80,maxTotalDepth = 1000,minSampleDepth = 30,minGQ = 99)

df_filt <- runQTLseqAnalysis(SNPset = df_filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(200, 200),replications = 10000,intervals = c(95, 99))

df_filt <- runGprimeAnalysis(df_filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)

plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()


plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()

plotQTLStats(SNPset = df_filt, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+ylim(-0.4, 0.4)



################################################################

# cage2
p2 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_221.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_229.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_248.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_244.tricube"),color = "orange",size=1)

p2.2 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_252.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_253.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_254.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_255.tricube"),color = "orange",size=1)

# cage3
p3 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_222.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_230.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_249.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_245.tricube"),color = "orange",size=1)
p3.2 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_256.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_257.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_258.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_259.tricube"),color = "orange",size=1)

# cage4
p4 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_223.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_231.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_246.tricube"),color = "orange",size=1)
p4.2 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_260.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_261.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_262.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_263.tricube"),color = "orange",size=1)
# cage6
p6 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_225.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_233.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_251.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_247.tricube"),color = "orange",size=1)
p6.2 <- p0 +geom_line(aes_string(x = "POS", y = "AB_BSA_264.tricube"),color = "black",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_265.tricube"),color = "red",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_266.tricube"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "AB_BSA_267.tricube"),color = "orange",size=1)
