setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(warn=1)
options(stringsAsFactors = FALSE)

rm(list=ls())

library("stringr")
library("xlsx")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggrepel")
library("ggbeeswarm")
library("pwr")
library("boot")
library("ROCR")

# ----- Global control variables
xlOut <- FALSE # Toggle report output
graphOut <- FALSE # Toggle graph output
stripOut <- FALSE # Toggle stripchart output
volcanoOut <- FALSE # Toggle volcano plot output
logitOut <- FALSE # Toggle logistic regression plot output
aucOut <- FALSE # Toggle auc plot output

thrFC <- 1 # Threshold log2 fold-change 
thrAlpha <- 0.05 # Alpha threshold for comparative statistics 
thrLogit <- 0.5 # Probability threshold for the logistic regression 
optPower <- 0.8 # Optimal statistical power

stockfmolug <- 40

filepaths <- list(
			pepInfo="C:\\Users\\lem\\Dropbox\\LTRI\\2 Tools\\peptide_info.csv",
			pepconc="peptconc.csv",
			skyline="Sample_Skyline_Results.csv",
			xloutput="R_report_onetail.xlsx"
			)

# ----- Functions
itranslate <- function(x, y, newvar, by){
	x[newvar] <- y[match(x[[by]], y[[by]]), newvar]
	return(x)
}

firstcap <- function(x) {
	substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	return(x)
}

spreadby <- function(new_dt, old_dt, id, spreader, dat){
	id <- enquo(id)
	spreader <- enquo(spreader)
	dat <- enquo(dat)
	n_id <- quo_name(id)
	n_spreader <- quo_name(spreader)
	n_dat <- quo_name(dat)

	new_dt <- old_dt %>% 
		select(!!id, !!spreader, !!dat) %>% 
		mutate(!!n_spreader:=paste0(!!spreader, firstcap(n_dat))) %>% 
		spread_(key=n_spreader, value=n_dat) %>% 
		left_join(new_dt, ., by=n_id)
	return(new_dt)
}


# Equation for hyperbolic cut-off threshold in volcano plots
cthresh <- function(x){
	c  <-  0.5
	c/(abs(x-thrFC))+(-log10(thrAlpha))
} 

# ----- Data wrangling

# Load peptide information database
pepInfo <- read.csv(filepaths$pepInfo)
pepconc <- read.csv(filepaths$pepconc)

# Load Skyline data
ds <- read.csv(filepaths$skyline, na.strings="#N/A")

# Rename columns for consistency
ds <- ds %>% 
	rename(
		peptide = Peptide.Sequence,
		peptideMod = Peptide.Modified.Sequence,
		protein = Protein.Name,
		subjectID = SubjectID,
		bioreplicate = BioReplicate,
		condition = Condition,
		sampleType = Sample.Type,
		replicateName = Replicate.Name,
		isotope = Isotope.Label.Type,
		precursorMz = Precursor.Mz,
		precursorCharge = Precursor.Charge,
		productMz = Product.Mz,
		productCharge = Product.Charge,
		fragmentIon = Fragment.Ion,
		retentionTime = Retention.Time,
		area = Area,
		background = Background,
		ratioLightToHeavy = RatioLightToHeavy,
		fmoluL = Analyte.Concentration,
		ratioToStandard = Ratio.To.Standard,
		rdotp = DotProductLightToHeavy,
		peakRank = Peak.Rank)

# Populate protein name from pepInfo database
ds <- itranslate(ds, pepInfo, "protein", "peptide")
ds <- itranslate(ds, pepconc, "fmolug", "peptide")

# Create columns for IS concentration, trial, Precursor and Transition names
ds <- ds %>% 
	mutate(trial = as.numeric(str_match(.$replicateName, "rep([[:digit:]]+)")[,2])) %>% 
	mutate(precursorName = paste0(.$peptideMod, "_", .$precursorCharge, "+")) %>% 
	mutate(transitionName = paste0(.$precursorName, "-->", .$fragmentIon, "_", .$productCharge, "+")) %>% 
	mutate(transitionID = paste0(.$replicateName, "_", .$transitionName)) %>% 
	mutate(expID=paste(subjectID, condition)) %>% 
	mutate(expTrID=paste(.$expID, transitionName)) %>% 
	mutate(fmoluL = fmolug*stockfmolug)


#Reorder columns for convenience
ds <- ds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "precursorMz", "precursorCharge", "productMz", "productCharge", "fragmentIon", "isotope", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "retentionTime", "area", "background", "ratioLightToHeavy", "fmolug", "fmoluL", "ratioToStandard", "rdotp", "peakRank")]

# Create lists of unique values of concentrations, transitions and precursors
UfmoluL <- unique(ds$fmoluL)
Ufmolug <- unique(ds$fmolug)
Utransitions <- unique(ds$transitionName)
Uprecursors <- unique(ds$precursorName)
Upeptides <- unique(ds$peptide)
Uproteins <- unique(ds$protein)
Uisotope <- unique(ds$isotope)
UsubjectID <- unique(ds$subjectID)

# Create a tidy version of ds
# Spread columns by isotope label type
tds <- ds %>% 
	select(-c(isotope, area, background, peakRank, retentionTime, precursorMz, productMz))

tds <- spreadby(tds, ds, transitionID, isotope, area)
tds <- spreadby(tds, ds, transitionID, isotope, background)
tds <- spreadby(tds, ds, transitionID, isotope, retentionTime)
tds <- spreadby(tds, ds, transitionID, isotope, peakRank)
tds <- spreadby(tds, ds, transitionID, isotope, precursorMz)
tds <- spreadby(tds, ds, transitionID, isotope, productMz)


# Remove duplicates created by spreading
tds <- tds %>% 
	distinct(transitionID, .keep_all=TRUE)

# Reorder columns for convenience
tds <- tds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "heavyPrecursorMz", "lightPrecursorMz", "heavyProductMz", "lightProductMz", "precursorCharge", "productCharge", "fragmentIon", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "heavyRetentionTime", "lightRetentionTime", "fmolug", "fmoluL", "heavyArea", "lightArea", "heavyBackground", "lightBackground", "heavyPeakRank", "lightPeakRank", "rdotp", "ratioLightToHeavy", "ratioToStandard")]

#------------- Retention Time Analysis
# Calculate mean and standard deviation of Retention time for each precursor and compare differences between light and heavy
rtDS <- tds %>% 
	select(precursorName:ratioToStandard) %>% 
	group_by(precursorName) %>% 
	summarise(
		meanHeavyRetentionTime=mean(heavyRetentionTime, na.rm=TRUE),
		cvHeavyRetentionTime=sd(heavyRetentionTime, na.rm=TRUE)/mean(heavyRetentionTime, na.rm=TRUE),
		meanLightRetentionTime=mean(lightRetentionTime, na.rm=TRUE),
		cvLightRetentionTime=sd(lightRetentionTime, na.rm=TRUE)/mean(lightRetentionTime, na.rm=TRUE)) %>% 
	mutate(diffLightHeavy=abs(.$meanHeavyRetentionTime-.$meanLightRetentionTime)) 


rtDS <- itranslate(rtDS, ds, "protein", "precursorName")
rtDS <- itranslate(rtDS, ds, "peptide", "precursorName")

# Sort precursors by heavy retention time
rtDS <- rtDS %>% arrange(meanHeavyRetentionTime)

# Sort for convenience
rtDS <- rtDS %>% select(peptide, everything())
rtDS <- rtDS %>% select(protein, everything())

#------------- Replicate summary
repDS <- tds %>% 
	group_by(expTrID) %>% 
	summarise(
		meanHeavyRetentionTime = mean(heavyRetentionTime),
		cvHeavyRetentionTime = sd(heavyRetentionTime)/mean(heavyRetentionTime),
		meanLightRetentionTime = mean(lightRetentionTime),
		cvLightRetentionTime = sd(lightRetentionTime)/mean(lightRetentionTime),
		meanHeavyArea = mean(heavyArea),
		cvHeavyArea = sd(heavyArea)/mean(heavyArea),
		meanLightArea = mean(lightArea),
		cvLightArea = sd(lightArea)/mean(lightArea),
		meanHeavyBackground = mean(heavyBackground),
		cvHeavyBackground = sd(heavyBackground)/mean(heavyBackground),
		meanLightBackground = mean(lightBackground),
		cvLightBackground = sd(lightBackground)/mean(lightBackground),
		meanHeavyPeakRank = mean(heavyPeakRank),
		cvHeavyPeakRank = sd(heavyPeakRank)/mean(heavyPeakRank),
		meanLightPeakRank = mean(lightPeakRank),
		cvLightPeakRank = sd(lightPeakRank)/mean(lightPeakRank),
		meanRdotp = mean(rdotp),
		cvRdotp = sd(rdotp)/mean(rdotp),
		meanRatioLightToHeavy = mean(ratioLightToHeavy),
		cvRatioLightToHeavy = sd(ratioLightToHeavy)/mean(ratioLightToHeavy))

# Sort by heavy retention time
repDS <- repDS %>% arrange(meanHeavyRetentionTime)

# Add corresponding protein and peptide names for convenience
repDS <- itranslate(repDS, ds, "fmoluL", "expTrID")
repDS <- itranslate(repDS, ds, "protein", "expTrID")
repDS <- itranslate(repDS, ds, "peptide", "expTrID")
repDS <- itranslate(repDS, ds, "condition", "expTrID")
repDS <- itranslate(repDS, ds, "subjectID", "expTrID")

# Sort for convenience
repDS <- repDS %>% select(peptide, everything())
repDS <- repDS %>% select(protein, everything())
repDS <- repDS %>% select(condition, everything())
repDS <- repDS %>% select(subjectID, everything())

# ----- Peptide concentrations summary

peptideDS <- repDS %>% 
	group_by(subjectID, condition, peptide) %>% 
	summarise(meanLTH=mean(meanRatioLightToHeavy),
		meanHeavyArea=mean(meanHeavyArea),
		meanLightArea=mean(meanLightArea),
		fmoluL=mean(fmoluL)) %>% 
	ungroup()

peptideDS <- peptideDS %>% 	
	mutate(meanLTH=ifelse(.$meanLTH==0, .$meanLightArea/.$meanHeavyArea, .$meanLTH)) %>% 
	select(-c(meanLightArea, meanHeavyArea)) %>% 
	mutate(concentration_fmul=meanLTH*fmoluL) %>% 
	mutate(molw=pepInfo$mw[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(concentration_ugml=concentration_fmul*molw*1e-6) %>% 
	mutate(pepID=paste0(subjectID, peptide)) %>% 
	mutate(protein=pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(conclass=ifelse(condition=="POST",1,0))

# Sort for convenience
peptideDS <- peptideDS[c("pepID", "protein", "stageSpecificity", "peptide", "subjectID", "condition",  "meanLTH", "fmoluL", "concentration_fmul", "molw", "concentration_ugml" ,"conclass")]


# ----- Base data wrangling ends


# ----- Ranking of patient samples by protein concentrations
rankDS <- peptideDS %>%
	filter(condition == "PREvas") %>% 
	select(protein, subjectID, concentration_fmul) %>% 
	group_by(protein) %>% 
	mutate(rank=rank(-concentration_fmul)) %>% 
	ungroup() %>%
	select(-c(concentration_fmul)) %>%  
	spread(protein, rank) %>% 
	mutate(rankSum=rowSums(.[2:21])) %>% 
	arrange(rankSum)

# # ----- Summary of peptide concentrations by subjectID (spread version of the peptideDS)
subjectDS <- peptideDS %>% 
	select(pepID, condition, concentration_fmul) %>% 
	mutate(condition=paste0("concentration_fmul", condition)) %>% 
	spread(condition, concentration_fmul) %>%
	select(pepID, concentration_fmulPRE, concentration_fmulPOST) %>% 
	left_join(peptideDS, ., by="pepID") %>% 
	distinct(pepID, .keep_all=TRUE) %>% 
	select(-c(pepID, concentration_fmul, condition, meanLTH, fmoluL)) %>% 
	mutate(fc=concentration_fmulPRE/concentration_fmulPOST) %>% 
	mutate(protein= pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity= pepInfo$stageSpecificity[match(.$peptide, pepInfo$peptide)])

# Sort for convenience
subjectDS <- subjectDS[c("subjectID", "protein", "peptide", "stageSpecificity", "concentration_fmulPRE", "concentration_fmulPOST", "fc")]

# ----- Summary of biomarker performance
biomarkerDS <- subjectDS %>% 
	group_by(peptide) %>% 
	summarise(
		meanPRE=mean(concentration_fmulPRE),
		meanPOST=mean(concentration_fmulPOST),
		N=min(length(concentration_fmulPRE), length(concentration_fmulPOST)),
		sdPRE=sd(concentration_fmulPRE),
		sdPOST=sd(concentration_fmulPOST),
		sdPooled=sqrt((sd(concentration_fmulPRE)^2+sd(concentration_fmulPOST)^2)/2),
		sePRE=sd(concentration_fmulPRE)/length(concentration_fmulPRE),
		sePOST=sd(concentration_fmulPOST)/length(concentration_fmulPOST),
		p.ttest=t.test(concentration_fmulPRE, concentration_fmulPOST, alternative="greater", paired=TRUE)[["p.value"]],
		t.ttest=t.test(concentration_fmulPRE, concentration_fmulPOST, alternative="greater", paired=TRUE)[["statistic"]],
		cohenD=t.test(concentration_fmulPRE, concentration_fmulPOST,paired=TRUE)[["statistic"]]/sqrt(min(length(concentration_fmulPRE), length(concentration_fmulPOST))),
		p.utest=wilcox.test(concentration_fmulPRE, concentration_fmulPOST,paired=TRUE)[["p.value"]],
		meanFC=mean(fc)
		) %>% 
	mutate(protein=pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(.$protein, pepInfo$protein)]) %>% 
	mutate(log2FC=log2(meanFC)) %>% 
	mutate(p.bh=p.adjust(p.ttest, method = "BH")) %>% 
	mutate(p.utBH=p.adjust(p.utest, method = "BH")) %>% 
	mutate(p.tPower=pwr.t.test(power=NULL, n=.$N, d=.$cohenD, sig.level=thrAlpha, type="paired")[["power"]]) %>% 
	mutate(p.optTPower=sapply(.$cohenD, function(d) round(pwr.t.test(power=0.8, n=NULL, d=d, sig.level=0.05, type="paired")[["n"]],1))) %>% 
	
	mutate(thrT.test=abs(.$log2FC) > thrFC & -log10(.$p.ttest) > cthresh(.$log2FC)) %>% 
	mutate(thrBH=abs(.$log2FC) > thrFC & -log10(.$p.bh) > cthresh(.$log2FC)) %>%
	mutate(thrutBH=abs(.$log2FC) > thrFC & -log10(.$p.utBH) > cthresh(.$log2FC)) 
	

biomarkerDS <- biomarkerDS[c("peptide", "protein", "stageSpecificity", "N", "meanFC", "log2FC", "meanPRE", "meanPOST", "sdPRE", "sdPOST", "sdPooled", "sePRE", "sePOST", "t.ttest", "p.ttest", "p.bh", "cohenD", "p.tPower", "p.optTPower", "thrT.test", "thrBH", "p.utest", "p.utBH", "thrutBH")]

# ----- AUROC analysis

set.seed(1354)
aucDS <- data.frame()

for (target in Upeptides){
	protTarget <- pepInfo$protein[match(target, pepInfo$peptide)]
	#Split into training and test set
	fullset <- peptideDS %>% 
		filter(peptide==target) 

	train <- fullset %>% 
		sample_frac(0.3, replace=FALSE)

	test <- fullset %>% 
		setdiff(train)

	#Train logistic regression classifier
	model <- glm(formula=conclass~concentration_ugml, data=fullset, family=binomial)
	modelInt <- coefficients(model)[[1]]
	modelSlope <- coefficients(model)[[2]]
	modelConcAtThr <- (logit(thrLogit)-modelInt)/modelSlope

	# inv.logit(modelInt + modelSlope * x)
	# anov <- anova(model, test="Chisq")

	#Test the classifier on the test set
	pred <- predict(model, fullset, type="response")
	predT <- ifelse(pred>thrLogit, 1, 0)

	#Confusion Matrix
	confmat <- table(predT, fullset$conclass)
	accuracy <- round(sum(diag(confmat))/sum(confmat),2)
	sensitivity <- round(confmat[2,2]/sum(confmat[,2]),2)
	specificity <- round(confmat[1,1]/sum(confmat[,1]),2)

	# Logistic regression graph
	if (graphOut & logitOut){
		p <- fullset %>% 
			ggplot(aes(x=concentration_ugml, y=conclass))+
			geom_jitter(size=2, width = 0, height = 0.02)+
			labs(title=paste(protTarget, " (", substr(target, 1, 5), "...)", sep=""), subtitle=sprintf("Accuracy: %s, Sensitivity: %s, Specificity: %s", accuracy, sensitivity, specificity))+
			theme(text = element_text(size=15), legend.position="none")+
			xlab(expression("Concentration ("*mu*"g/mL)")) + ylab("Probability of POST classification")+
			stat_function(fun=function(x) inv.logit(modelInt+modelSlope*x), geom="line", colour="deepskyblue4", size=1)+
			geom_vline(xintercept=modelConcAtThr, linetype="longdash")+
			scale_x_log10()

		ggsave(paste("g3_", protTarget, "_classifier.jpg", sep=""), plot=p, device="jpg")
	}

	#ROCR analysis
	pr <- prediction(pred, fullset$conclass)
	prf <- performance(pr, measure = "tpr", x.measure = "fpr")
	accprf <- performance(pr, measure = "acc")

	# Make plot data
	plotdat <- data.frame(fp=prf@x.values[[1]],tp=prf@y.values[[1]],cut=prf@alpha.values[[1]], acc=accprf@y.values[[1]])
	 
	auc <- performance(pr, measure = "auc")
	auc <- round(auc@y.values[[1]],2)

	
	if (graphOut){
		# AUC plots
		if (aucOut){
			p <- plotdat %>% 
				ggplot(aes(x=fp, y=tp, colour=cut))+
				geom_line()+ geom_point()+
				scale_colour_gradient2(low = "green3", mid = "orange3",
					high = "red4", midpoint = 0.5, name="Threshold")+
				labs(title=paste(protTarget, " (", substr(target, 1, 5), "...)", sep=""), subtitle=paste("AUC: ", auc, sep=""))+
				theme(text = element_text(size=15), legend.title = element_text(size=10))+
				xlab("False Positive Rate") + ylab("True Positive Rate")+
				stat_function(fun=function(x) x, geom="line", colour="black", linetype="dotted")+
				geom_text_repel(aes(label=round(acc,2)))

			ggsave(paste("g4_", protTarget, "_AUC.jpg", sep=""), plot=p, device="jpg")
		}	
	}

	drow <- data.frame(peptide=target, concentration_ugmlAtThr=modelConcAtThr,
			accuracy=accuracy, sensitivity=sensitivity , specificity=specificity, auc=auc)
	aucDS <- rbind(aucDS, drow)
}


# Merge aucDS with biomarkerDS
biomarkerDS <- biomarkerDS %>% 
	left_join(., aucDS, by="peptide")
rm("drow")


bioSummary <- biomarkerDS %>% 
	select(peptide, protein, stageSpecificity, meanFC, p.ttest, p.bh, cohenD, p.tPower, p.optTPower, concentration_ugmlAtThr, accuracy, sensitivity, specificity, auc) %>% 
	rename(
		T.test.p = p.ttest,
		p.adjusted.BH = p.bh,
		effect.Size = cohenD,
		power = p.tPower,
		N.at.power0.8 = p.optTPower,
		concentration_ugml.cutoff = concentration_ugmlAtThr
		) %>% 
	mutate(score=(((meanFC/max(meanFC)+power+accuracy+auc)/4))) %>% 
	arrange(desc(score))
set.seed(42)
# ----- Graphs
if (graphOut & volcanoOut){
	# Volcano plot
	p <- biomarkerDS %>% 
		ggplot(aes(x=log2FC, y=-log10(p.bh), colour=thrBH)) +
		geom_point(size=3) +
		theme_bw()+
		scale_color_manual(values=c("FALSE" = "grey44", "TRUE" = "black"))+
		scale_y_continuous(breaks = round(seq(1, 4, by = 1),1))+
		scale_x_continuous(breaks = round(seq(0, 9, by = 1),1))+
		coord_cartesian(ylim=c(1.4, 4), xlim=c(1, 9))+
		stat_function(fun=cthresh, geom="line", colour="gray34", xlim=c(1,10), linetype="longdash")+
		# geom_vline(xintercept=thrFC, linetype="dotted")+
		# geom_hline(yintercept=-log10(thrAlpha), linetype="dotted")+
		geom_text_repel(aes(label=.$protein), size=4, segment.color="black", segment.alpha=0, force=5, point.padding = unit(0.1, "lines"), max.iter=3e2, nudge_x = ifelse(biomarkerDS$thrBH == FALSE, -1, 0.1))+
		theme(aspect.ratio=1, text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"), legend.position="none")+
		xlab("log2 Mean Fold Change") + ylab("-log10 Adjusted p-value")
		
	ggsave("g2_volcano_bh.jpg", plot=p, device="jpg", width=5, height=5*1)

	p <- biomarkerDS %>% 
		ggplot(aes(x=log2FC, y=-log10(p.utBH), colour=thrutBH)) +
		geom_point(size=3) +
		scale_color_manual(values=c("FALSE" = "sienna4", "TRUE" = "green4"))+
		scale_y_continuous(breaks = round(seq(1, 6, by = 1),1))+
		scale_x_continuous(breaks = round(seq(0, 9, by = 1),1))+
		coord_cartesian(ylim=c(1, 6), xlim=c(1, 9))+
		geom_vline(xintercept=thrFC, linetype="dotted")+
		geom_hline(yintercept=-log10(thrAlpha), linetype="dotted")+
		geom_text_repel(aes(label=.$protein), size=4, segment.alpha=0.5, force=5)+
		theme(text = element_text(size=15), legend.position="none")+
		xlab("log2 Mean Fold Change") + ylab("-log10 B-H adj. U test p-value (alpha=0.05)")+
		stat_function(fun=cthresh, geom="line", colour="gray34", xlim=c(1,9), linetype="longdash")

	ggsave("g2_volcano_utest.jpg", plot=p, device="jpg", width=7, height=7)
}

peptideDS <- peptideDS %>% 
	mutate(concentration_ugmlAtThr=biomarkerDS$concentration_ugmlAtThr[match(.$protein, biomarkerDS$protein)])

# Stripchars for each peptide (with and without SubjectID labels)
if (stripOut){
	p <- peptideDS %>% 
		ggplot(aes(x=condition, y=concentration_ugml)) +
		geom_hline(aes(yintercept=concentration_ugmlAtThr), linetype="longdash")+
		theme_bw()+
		geom_beeswarm(size=0.5, cex=1) +
		# geom_jitter(position=position_jitter(0.05)) +
		# labs(title=sprintf("%s (%s)", protTarget, stagel), subtitle=sprintf("%s...", substr(target,1,5)))+
		scale_x_discrete(limits=c("PRE", "POST"))+
		xlab("Condition") + ylab(expression("Concentration ("*mu*"g/mL)"))+
		theme(aspect.ratio=0.8, text = element_text(size=8), axis.text = element_text(size=6, color="black"), axis.text.x = element_text(angle=90), axis.line=element_line(color="black"), panel.border=element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
		stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, color="red", size=0.2)+
		scale_y_log10()+
		facet_wrap(~ protein, scales="free_y", nrow=4)
			
	ggsave(paste("g1_facet_beeswarm",".jpg", sep=""), plot=p, device="jpg", width=7, height=7*0.8)
}


#------------- Output
if(xlOut){
		write.xlsx(as.data.frame(rtDS), filepaths$xloutput, sheetName="rtDS", row.names=FALSE)
	for (sheet in c("rankDS", "bioSummary")){
		write.xlsx(as.data.frame(get(sheet)), filepaths$xloutput, sheetName=sheet, row.names=FALSE, append=TRUE)
	}
}else{
	View(bioSummary)
}