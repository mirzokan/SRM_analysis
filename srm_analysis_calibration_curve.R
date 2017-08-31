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
xlOut <- TRUE # Toggle report output
graphOut <- FALSE # Toggle graph output
stripOut <- FALSE # Toggle stripchart output
volcanoOut <- FALSE # Toggle volcano plot output
logitOut <- FALSE # Toggle logistic regression plot output
aucOut <- FALSE # Toggle auc plot output

thrFC <- 1 # Threshold log2 fold-change 
thrAlpha <- 0.05 # Alpha threshold for comparative statistics 
thrLogit <- 0.5 # Probability threshold for the logistic regression 
optPower <- 0.8 # Optimal statistical power

stockuguL <- 40

filepaths <- list(pepInfo="C:\\Users\\lem\\Dropbox\\LTRI\\2 Tools\\peptide_info.csv",
			skyline="Mirzo_Skyline_Results.csv",
			xloutput="R_report_stcurve.xlsx",
			pepconc="peptconc.csv"
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
		peakRank = Peak.Rank
		)

# Populate protein name from pepInfo database
ds$protein <- pepInfo$protein[match(ds$peptide, pepInfo$peptide)]
ds$fmolug <- pepconc$fmolug[match(ds$peptide, pepconc$peptide)]

# Create columns for IS concentration, trial, Precursor and Transition names
ds <- ds %>% 
	mutate(trial = as.numeric(str_match(.$replicateName, "rep([[:digit:]]+)")[,2])) %>% 
	mutate(curvePoint = str_match(.$replicateName, "_P([[:digit:]]+_*[[:digit:]]*)-")[,2]) %>% 
	mutate(curvePoint = as.numeric(gsub("_", ".", curvePoint))) %>% 
	mutate(fmolug = fmolug) %>% 
	mutate(fmoluL = fmolug*stockuguL) %>% 
	mutate(precursorName = paste0(.$peptideMod, "_", .$precursorCharge, "+")) %>% 
	mutate(transitionName = paste0(.$precursorName, "-->", .$fragmentIon, "_", .$productCharge, "+")) %>% 
	mutate(transitionID = paste0(.$replicateName, "_", .$transitionName)) %>% 
	mutate(expID=paste(subjectID, condition)) %>% 
	mutate(expTrID=paste(curvePoint, transitionName))


#Reorder columns for convenience
ds <- ds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "precursorMz", "precursorCharge", "productMz", "productCharge", "fragmentIon", "isotope", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "retentionTime", "area", "background", "ratioLightToHeavy", "curvePoint", "fmoluL", "fmolug", "ratioToStandard", "rdotp", "peakRank")]

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
tds <- tds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "heavyPrecursorMz", "lightPrecursorMz", "heavyProductMz", "lightProductMz", "precursorCharge", "productCharge", "fragmentIon", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "heavyRetentionTime", "lightRetentionTime", "curvePoint", "fmoluL", "fmolug", "heavyArea", "lightArea", "heavyBackground", "lightBackground", "heavyPeakRank", "lightPeakRank", "rdotp", "ratioLightToHeavy", "ratioToStandard")]

# Take care of ratioLightToHeavy = 0
tds <- tds %>% 
	mutate(ratioLightToHeavy=ifelse(ratioLightToHeavy==0, lightArea/heavyArea, ratioLightToHeavy))

# Create lists of unique values of concentrations, transitions and precursors
Ufmolug <- unique(ds$fmolug)
UfmoluL <- unique(ds$fmoluL)
Utransitions <- unique(ds$transitionName)
Uprecursors <- unique(ds$precursorName)
Upeptides <- unique(ds$peptide)
Uproteins <- unique(ds$protein)
Uisotope <- unique(ds$isotope)
UsubjectID <- unique(ds$subjectID)
UcurvePoint <- unique(ds$curvePoint)

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
	mutate(diffLightHeavy=abs(meanHeavyRetentionTime-meanLightRetentionTime))

# Sort precursors by heavy retention time
rtDS <- rtDS %>% arrange(meanHeavyRetentionTime)

# Add corresponding protein and peptide names for convenience
rtDS$protein <- ds$protein[match(rtDS$precursorName, ds$precursorName)]
rtDS$peptide <- ds$peptide[match(rtDS$precursorName, ds$precursorName)]

# Sort for convenience
rtDS <- rtDS %>% select(peptide, everything())
rtDS <- rtDS %>% select(protein, everything())

#------------- Replicate summary
# repDS <- tds %>% 
# 	group_by(expTrID) %>% 
# 	summarise(
# 		meanHeavyRetentionTime = mean(heavyRetentionTime),
# 		cvHeavyRetentionTime = sd(heavyRetentionTime)/mean(heavyRetentionTime),
# 		meanLightRetentionTime = mean(lightRetentionTime),
# 		cvLightRetentionTime = sd(lightRetentionTime)/mean(lightRetentionTime),
# 		meanHeavyArea = mean(heavyArea),
# 		cvHeavyArea = sd(heavyArea)/mean(heavyArea),
# 		meanLightArea = mean(lightArea),
# 		cvLightArea = sd(lightArea)/mean(lightArea),
# 		meanHeavyBackground = mean(heavyBackground),
# 		cvHeavyBackground = sd(heavyBackground)/mean(heavyBackground),
# 		meanLightBackground = mean(lightBackground),
# 		cvLightBackground = sd(lightBackground)/mean(lightBackground),
# 		meanHeavyPeakRank = mean(heavyPeakRank),
# 		cvHeavyPeakRank = sd(heavyPeakRank)/mean(heavyPeakRank),
# 		meanLightPeakRank = mean(lightPeakRank),
# 		cvLightPeakRank = sd(lightPeakRank)/mean(lightPeakRank),
# 		meanRdotp = mean(rdotp),
# 		cvRdotp = sd(rdotp)/mean(rdotp),
# 		meanRatioLightToHeavy = ifelse(mean(ratioLightToHeavy)==0, mean(lightArea/heavyArea), mean(ratioLightToHeavy)),
# 		cvRatioLightToHeavy = ifelse(mean(ratioLightToHeavy)==0, sd(lightArea/heavyArea)/mean(lightArea/heavyArea), sd(ratioLightToHeavy)/mean(ratioLightToHeavy)),
# 		meanRatioHeavyToLight = ifelse(mean(ratioLightToHeavy)==0, mean(heavyArea/lightArea), mean(1/ratioLightToHeavy)),
# 		cvRatioHeavyToLight = ifelse(mean(ratioLightToHeavy)==0, sd(heavyArea/lightArea)/mean(heavyArea/lightArea), sd(1/ratioLightToHeavy)/mean(1/ratioLightToHeavy)))


# # Sort by heavy retention time
# repDS <- repDS %>% arrange(meanHeavyRetentionTime)

# # Add corresponding protein and peptide names for convenience
# repDS$fmoluL <- ds$fmoluL[match(repDS$expTrID, ds$expTrID)]
# repDS$fmolug <- ds$fmolug[match(repDS$expTrID, ds$expTrID)]
# repDS$protein <- ds$protein[match(repDS$expTrID, ds$expTrID)]
# repDS$peptide <- ds$peptide[match(repDS$expTrID, ds$expTrID)]
# repDS$curvePoint <- ds$curvePoint[match(repDS$expTrID, ds$expTrID)]

# # Sort for convenience
# repDS <- repDS %>% select(peptide, everything())
# repDS <- repDS %>% select(protein, everything())
# repDS <- repDS %>% select(curvePoint, everything())

repDS <- tds %>% 
	group_by(curvePoint, peptide, trial) %>% 
	summarise(
		fmoluL = mean(fmoluL, na.rm=TRUE),
		rdotp = mean(rdotp, na.rm=TRUE),
		LTH = mean(ratioLightToHeavy, na.rm=TRUE)
		) %>% 
	mutate(HTL=1/LTH) %>% 
	mutate(theorHfmul=fmoluL*curvePoint) %>% 
	mutate(hConcentration_fmul=fmoluL*HTL) %>% 
	mutate(molw=pepInfo$mw[match(peptide, pepInfo$peptide)]) %>% 
	mutate(theorHugml=theorHfmul*molw*1e-6) %>% 
	mutate(hConcentration_ugml=hConcentration_fmul*molw*1e-6) %>% 
	mutate(protein=pepInfo$protein[match(peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(peptide, pepInfo$peptide)])


# if (graphOut){
# 	for (i in Upeptides){
# 		lProtein <- pepInfo$protein[match(i, pepInfo$peptide)]
# 		lStagel <- pepInfo$stageSpecificity[match(i, pepInfo$peptide)]
# 		p <- repDS %>% 
# 			filter(peptide==i) %>% 
# 			ggplot(aes(x=theorHugml, y=HTL))+
# 			geom_point(size=2)+
# 			scale_y_log10()+
# 			scale_x_log10()+
# 			labs(title=lProtein, subtitle=sprintf("%s (%s...)", lStagel, substr(i, 1,5)))+
# 			geom_smooth(method='lm')+
# 			geom_hline(yintercept=1, linetype="longdash")

# 		ggsave(paste("g1_", lProtein, "_stcurve.jpg", sep=""), plot=p, device="jpg")
# 	}
# }
	
stcurveDS <- repDS %>% 
	group_by(curvePoint, peptide) %>% 
	summarise(
		meanHTL=mean(HTL, na.rm=TRUE),
		sdHTL=sd(HTL, na.rm=TRUE),
		cvHTL=sd(HTL, na.rm=TRUE)/mean(HTL, na.rm=TRUE),
		theorHugml=mean(theorHugml, na.rm=TRUE)
		) %>% 
	mutate(protein=pepInfo$protein[match(peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(protein, pepInfo$protein)]) %>% 
	mutate(sdthr=(cvHTL<0.2))

stcurveDS <- stcurveDS %>% 
	group_by(peptide) %>% 
	summarise(maxFcv=max(ifelse(!sdthr,curvePoint,0))) %>% 
	left_join(stcurveDS, ., by="peptide") %>% 
	mutate(linlim=pepconc$linlim[match(peptide, pepconc$peptide)]) %>% 
	mutate(linthr=curvePoint>=linlim) %>% 
	mutate(sdthr=ifelse(curvePoint<=maxFcv, FALSE, sdthr)) %>% 
	mutate(cthr=linthr&sdthr) %>% 
	select(-c(maxFcv, linlim))


if (graphOut){
	for (i in Upeptides){
		lProtein <- pepInfo$protein[match(i, pepInfo$peptide)]
		lStagel <- pepInfo$stageSpecificity[match(i, pepInfo$peptide)]
		p <- stcurveDS %>% 
			filter(peptide==i) %>% 
			ggplot(aes(x=theorHugml, y=meanHTL, colour=cthr))+
			scale_color_manual(values=c("FALSE" = "sienna4", "TRUE" = "green4"))+
			geom_point(size=2)+
			geom_errorbar(aes(ymin=meanHTL-sdHTL, ymax=meanHTL+sdHTL), width=.2)+
			scale_y_log10(breaks=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000), labels = scales::comma)+
			scale_x_log10()+
			labs(title=lProtein, subtitle=sprintf("%s (%s...)", lStagel, substr(i, 1,5)))+
			geom_smooth(method="lm", se=FALSE)+
			# geom_text_repel(aes(label=curvePoint))+
			xlab(expression("Concentration of heavy IS ("*mu*"g/mL)")) + ylab("Heavy-to-light ratio")+
			theme(text = element_text(size=15), legend.position="none")+
			geom_hline(yintercept=1, linetype="longdash")

		ggsave(paste("g2_", lProtein, "_stcurve1p.jpg", sep=""), plot=p, device="jpg")
	}
}


loqDS <- stcurveDS %>% 
	filter(cthr==TRUE) %>% 
	group_by(peptide) %>% 
	summarise(
		loqPoint=min(curvePoint, na.rm=TRUE),
		loq_ugmL=min(theorHugml)
		) %>% 
	mutate(protein=pepInfo$protein[match(peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(peptide, pepInfo$peptide)])

loqDS <- loqDS %>% select(stageSpecificity, everything())
loqDS <- loqDS %>% select(peptide, everything())
loqDS <- loqDS %>% select(protein, everything())

# Linear regression

# for (i in Upeptides){
# lProtein <- pepInfo$protein[match(i, pepInfo$peptide)]
# lStagel <- pepInfo$stageSpecificity[match(i, pepInfo$peptide)]

# 	curvesub <- stcurveDS %>% 
# 		filter(peptide==i & sdthr)

# 	linest <- lm(meanHTL~theorHugml, data=curvesub)

# 	# print(i)
# 	# print(summary(linest)$adj.r.squared)


# 	if (graphOut){
# 		p <- stcurveDS %>% 
# 			filter(peptide==i) %>% 
# 			ggplot(aes(x=theorHugml, y=meanHTL, colour=sdthr))+
# 			scale_color_manual(values=c("FALSE" = "sienna4", "TRUE" = "green4"))+
# 			geom_point(size=2)+
# 			geom_errorbar(aes(ymin=meanHTL-sdHTL, ymax=meanHTL+sdHTL), width=.2)+
# 			scale_y_log10()+
# 			scale_x_log10()+
# 			labs(title=lProtein, subtitle=sprintf("%s (%s...)", lStagel, substr(i, 1,5)))+
# 			geom_smooth(method="lm")+
# 			geom_hline(yintercept=1, linetype="longdash")

# 		ggsave(paste("g3_", lProtein, "_stcurve1p.jpg", sep=""), plot=p, device="jpg")
		
# 	}

# }




# # ----- Opimal fmolug
# fmolugDS <- tds %>% 
# 	filter(curvePoint==1) %>% 
# 	group_by(peptide) %>% 
# 	summarise(
# 		fmolug=mean(fmolug, na.rm=TRUE),
# 		fmoluL=mean(fmoluL, na.rm=TRUE),
# 		meanLTH = mean(ratioLightToHeavy),
# 		cvLTH = sd(ratioLightToHeavy)/mean(ratioLightToHeavy)
# 		)  %>% 
# 	mutate(optFmolug=round((fmolug*meanLTH)*3, 2)) %>% 
# 	mutate(protein=tds$protein[match(.$peptide, tds$peptide)]) %>% 
# 	select(protein, everything())








# # ----- Ranking of patient samples by protein concentrations
# rankDS <- peptideDS %>%
# 	filter(condition == "PRE") %>% 
# 	select(protein, subjectID, concentration_fmul) %>% 
# 	group_by(protein) %>% 
# 	mutate(rank=rank(-concentration_fmul)) %>% 
# 	ungroup() %>%
# 	select(-c(concentration_fmul)) %>%  
# 	spread(protein, rank) %>% 
# 	mutate(rankSum=rowSums(.[2:21])) %>% 
# 	arrange(rankSum)

# # ----- Summary of peptide concentrations by subjectID (spread version of the peptideDS)
# subjectDS <- peptideDS %>% 
# 	select(pepID, condition, concentration_fmul) %>% 
# 	mutate(condition=paste0("concentration_fmul", condition)) %>% 
# 	spread(condition, concentration_fmul) %>%
# 	select(pepID, concentration_fmulPRE, concentration_fmulPOST) %>% 
# 	left_join(peptideDS, ., by="pepID") %>% 
# 	distinct(pepID, .keep_all=TRUE) %>% 
# 	select(-c(pepID, concentration_fmul, condition, meanLTH, fmoluL)) %>% 
# 	mutate(fc=concentration_fmulPRE/concentration_fmulPOST) %>% 
# 	mutate(protein= pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
# 	mutate(stageSpecificity= pepInfo$stageSpecificity[match(.$peptide, pepInfo$peptide)])

# # Sort for convenience
# subjectDS <- subjectDS[c("subjectID", "protein", "peptide", "stageSpecificity", "concentration_fmulPRE", "concentration_fmulPOST", "fc")]

# # ----- Summary of biomarker performance
# biomarkerDS <- subjectDS %>% 
# 	group_by(peptide) %>% 
# 	summarise(
# 		meanPRE=mean(concentration_fmulPRE),
# 		meanPOST=mean(concentration_fmulPOST),
# 		N=min(length(concentration_fmulPRE), length(concentration_fmulPOST)),
# 		sdPRE=sd(concentration_fmulPRE),
# 		sdPOST=sd(concentration_fmulPOST),
# 		sdPooled=sqrt((sd(concentration_fmulPRE)^2+sd(concentration_fmulPOST)^2)/2),
# 		sePRE=sd(concentration_fmulPRE)/length(concentration_fmulPRE),
# 		sePOST=sd(concentration_fmulPOST)/length(concentration_fmulPOST),
# 		p.ttest=t.test(concentration_fmulPRE, concentration_fmulPOST, alternative="greater", paired=TRUE)[["p.value"]],
# 		t.ttest=t.test(concentration_fmulPRE, concentration_fmulPOST, alternative="greater", paired=TRUE)[["statistic"]],
# 		cohenD=t.test(concentration_fmulPRE, concentration_fmulPOST,paired=TRUE)[["statistic"]]/sqrt(min(length(concentration_fmulPRE), length(concentration_fmulPOST))),
# 		p.utest=wilcox.test(concentration_fmulPRE, concentration_fmulPOST,paired=TRUE)[["p.value"]],
# 		meanFC=mean(fc)
# 		) %>% 
# 	mutate(protein=pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
# 	mutate(stageSpecificity=pepInfo$stageSpecificity[match(.$protein, pepInfo$protein)]) %>% 
# 	mutate(log2FC=log2(meanFC)) %>% 
# 	mutate(p.bh=p.adjust(p.ttest, method = "BH")) %>% 
# 	mutate(p.utBH=p.adjust(p.utest, method = "BH")) %>% 
# 	mutate(p.tPower=pwr.t.test(power=NULL, n=.$N, d=.$cohenD, sig.level=thrAlpha, type="paired")[["power"]]) %>% 
# 	mutate(p.optTPower=sapply(.$cohenD, function(d) round(pwr.t.test(power=0.8, n=NULL, d=d, sig.level=0.05, type="paired")[["n"]],1))) %>% 
	
# 	mutate(thrT.test=abs(.$log2FC) > thrFC & -log10(.$p.ttest) > cthresh(.$log2FC)) %>% 
# 	mutate(thrBH=abs(.$log2FC) > thrFC & -log10(.$p.bh) > cthresh(.$log2FC)) %>%
# 	mutate(thrutBH=abs(.$log2FC) > thrFC & -log10(.$p.utBH) > cthresh(.$log2FC)) 
	

# biomarkerDS <- biomarkerDS[c("peptide", "protein", "stageSpecificity", "N", "meanFC", "log2FC", "meanPRE", "meanPOST", "sdPRE", "sdPOST", "sdPooled", "sePRE", "sePOST", "t.ttest", "p.ttest", "p.bh", "cohenD", "p.tPower", "p.optTPower", "thrT.test", "thrBH", "p.utest", "p.utBH", "thrutBH")]

# # ----- AUROC analysis

# set.seed(1354)
# aucDS <- data.frame()

# for (target in Upeptides){
# 	protTarget <- pepInfo$protein[match(target, pepInfo$peptide)]
# 	#Split into training and test set
# 	fullset <- peptideDS %>% 
# 		filter(peptide==target) 

# 	train <- fullset %>% 
# 		sample_frac(0.3, replace=FALSE)

# 	test <- fullset %>% 
# 		setdiff(train)

# 	#Train logistic regression classifier
# 	model <- glm(formula=conclass~concentration_ugml, data=fullset, family=binomial)
# 	modelInt <- coefficients(model)[[1]]
# 	modelSlope <- coefficients(model)[[2]]
# 	modelConcAtThr <- (logit(thrLogit)-modelInt)/modelSlope

# 	# inv.logit(modelInt + modelSlope * x)
# 	# anov <- anova(model, test="Chisq")

# 	#Test the classifier on the test set
# 	pred <- predict(model, fullset, type="response")
# 	predT <- ifelse(pred>thrLogit, 1, 0)

# 	#Confusion Matrix
# 	confmat <- table(predT, fullset$conclass)
# 	accuracy <- round(sum(diag(confmat))/sum(confmat),2)
# 	sensitivity <- round(confmat[2,2]/sum(confmat[,2]),2)
# 	specificity <- round(confmat[1,1]/sum(confmat[,1]),2)

# 	# Logistic regression graph
# 	if (graphOut & logitOut){
# 		p <- fullset %>% 
# 			ggplot(aes(x=concentration_ugml, y=conclass))+
# 			geom_jitter(size=2, width = 0, height = 0.02)+
# 			labs(title=paste(protTarget, " (", substr(target, 1, 5), "...)", sep=""), subtitle=sprintf("Accuracy: %s, Sensitivity: %s, Specificity: %s", accuracy, sensitivity, specificity))+
# 			theme(text = element_text(size=15), legend.position="none")+
# 			xlab(expression("Concentration ("*mu*"g/mL)")) + ylab("Probability of POST classification")+
# 			stat_function(fun=function(x) inv.logit(modelInt+modelSlope*x), geom="line", colour="deepskyblue4", size=1)+
# 			geom_vline(xintercept=modelConcAtThr, linetype="longdash")+
# 			scale_x_log10()

# 		ggsave(paste("g3_", protTarget, "_classifier.jpg", sep=""), plot=p, device="jpg")
# 	}

# 	#ROCR analysis
# 	pr <- prediction(pred, fullset$conclass)
# 	prf <- performance(pr, measure = "tpr", x.measure = "fpr")
# 	accprf <- performance(pr, measure = "acc")

# 	# Make plot data
# 	plotdat <- data.frame(fp=prf@x.values[[1]],tp=prf@y.values[[1]],cut=prf@alpha.values[[1]], acc=accprf@y.values[[1]])
	 
# 	auc <- performance(pr, measure = "auc")
# 	auc <- round(auc@y.values[[1]],2)

	
# 	if (graphOut){
# 		# AUC plots
# 		if (aucOut){
# 			p <- plotdat %>% 
# 				ggplot(aes(x=fp, y=tp, colour=cut))+
# 				geom_line()+ geom_point()+
# 				scale_colour_gradient2(low = "green3", mid = "orange3",
# 					high = "red4", midpoint = 0.5, name="Threshold")+
# 				labs(title=paste(protTarget, " (", substr(target, 1, 5), "...)", sep=""), subtitle=paste("AUC: ", auc, sep=""))+
# 				theme(text = element_text(size=15), legend.title = element_text(size=10))+
# 				xlab("False Positive Rate") + ylab("True Positive Rate")+
# 				stat_function(fun=function(x) x, geom="line", colour="black", linetype="dotted")+
# 				geom_text_repel(aes(label=round(acc,2)))

# 			ggsave(paste("g4_", protTarget, "_AUC.jpg", sep=""), plot=p, device="jpg")
# 		}

# 		# Stripchars for each peptide (with and without SubjectID labels)
# 		if (stripOut){
# 			for (j in c(TRUE, FALSE)){
# 					stagel <- pepInfo$stageSpecificity[match(target, pepInfo$peptide)]
# 					p <- peptideDS %>% 
# 						filter(peptide==target) %>%
# 						ggplot(aes(x=condition, y=concentration_ugml)) +
# 						theme_bw()+
# 						geom_beeswarm(size=2, cex=2) +
# 						# geom_jitter(position=position_jitter(0.05)) +
# 						geom_hline(yintercept=modelConcAtThr, linetype="longdash")+
# 						labs(title=sprintf("%s (%s)", protTarget, stagel), subtitle=sprintf("%s...", substr(target,1,5)))+
# 						scale_x_discrete(limits=c("PRE", "POST"))+
# 						xlab("Condition") + ylab(expression("Concentration ("*mu*"g/mL)"))+
# 						theme(text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"),)+
# 						stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.2)+
# 						scale_y_log10()
							
# 					if (j){
# 						p <- p + geom_text_repel(aes(label=subjectID), size=3, segment.alpha=0.2)
# 						filepref <- "g_tl_"
# 					}else{
# 						filepref <- "g_nl_"
# 					}
					
# 					ggsave(paste(filepref, protTarget,"_beeswarm",".jpg", sep=""), plot=p, device="jpg", width=7, height=7)
# 			}
# 		}
# 	}

# 	drow <- data.frame(peptide=target, concentration_ugmlAtThr=modelConcAtThr,
# 			accuracy=accuracy, sensitivity=sensitivity , specificity=specificity, auc=auc)
# 	aucDS <- rbind(aucDS, drow)
# }

# # Merge aucDS with biomarkerDS
# biomarkerDS <- biomarkerDS %>% 
# 	left_join(., aucDS, by="peptide")
# rm("drow")


# bioSummary <- biomarkerDS %>% 
# 	select(peptide, protein, stageSpecificity, meanFC, p.ttest, p.bh, cohenD, p.tPower, p.optTPower, concentration_ugmlAtThr, accuracy, sensitivity, specificity, auc) %>% 
# 	rename(
# 		T.test.p = p.ttest,
# 		p.adjusted.BH = p.bh,
# 		effect.Size = cohenD,
# 		power = p.tPower,
# 		N.at.power0.8 = p.optTPower,
# 		concentration_ugml.cutoff = concentration_ugmlAtThr
# 		) %>% 
# 	mutate(score=(((meanFC/max(meanFC)+power+accuracy+auc)/4))) %>% 
# 	arrange(desc(score))
# set.seed(42)
# # ----- Graphs
# if (graphOut & volcanoOut){
# 	# Volcano plot
# 	p <- biomarkerDS %>% 
# 		ggplot(aes(x=log2FC, y=-log10(p.bh), colour=thrBH)) +
# 		geom_point(size=3) +
# 		theme_bw()+
# 		scale_color_manual(values=c("FALSE" = "grey44", "TRUE" = "black"))+
# 		scale_y_continuous(breaks = round(seq(1, 4, by = 1),1))+
# 		scale_x_continuous(breaks = round(seq(0, 9, by = 1),1))+
# 		coord_cartesian(ylim=c(1.4, 4), xlim=c(1, 9))+
# 		stat_function(fun=cthresh, geom="line", colour="gray34", xlim=c(1,10), linetype="longdash")+
# 		# geom_vline(xintercept=thrFC, linetype="dotted")+
# 		# geom_hline(yintercept=-log10(thrAlpha), linetype="dotted")+
# 		geom_text_repel(aes(label=.$protein), size=4, segment.color="black", segment.alpha=0, force=5, point.padding = unit(0.1, "lines"), max.iter=3e2, nudge_x = ifelse(biomarkerDS$thrBH == FALSE, -1, 0.1))+
# 		theme(aspect.ratio=0.6, text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"), legend.position="none")+
# 		xlab("log2 Mean Fold Change") + ylab("-log10 Adjusted p-value")
		
# 	ggsave("g2_volcano_bh.jpg", plot=p, device="jpg", width=7, height=4.2)

# 	# p <- biomarkerDS %>% 
# 	# 	ggplot(aes(x=log2FC, y=-log10(p.utBH), colour=thrutBH)) +
# 	# 	geom_point(size=3) +
# 	# 	scale_color_manual(values=c("FALSE" = "sienna4", "TRUE" = "green4"))+
# 	# 	scale_y_continuous(breaks = round(seq(1, 6, by = 1),1))+
# 	# 	scale_x_continuous(breaks = round(seq(0, 9, by = 1),1))+
# 	# 	coord_cartesian(ylim=c(1, 6), xlim=c(1, 9))+
# 	# 	geom_vline(xintercept=thrFC, linetype="dotted")+
# 	# 	geom_hline(yintercept=-log10(thrAlpha), linetype="dotted")+
# 	# 	geom_text_repel(aes(label=.$protein), size=4, segment.alpha=0.5, force=5)+
# 	# 	theme(text = element_text(size=15), legend.position="none")+
# 	# 	xlab("log2 Mean Fold Change") + ylab("-log10 B-H adj. U test p-value (alpha=0.05)")+
# 	# 	stat_function(fun=cthresh, geom="line", colour="gray34", xlim=c(1,9), linetype="longdash")

# 	# ggsave("g2_volcano_utest.jpg", plot=p, device="jpg", width=7, height=7)
# }

#------------- Output
if(xlOut){
		write.xlsx(as.data.frame(rtDS), filepaths$xloutput, sheetName="rtDS", row.names=FALSE)
	for (sheet in c("loqDS", "stcurveDS")){
		write.xlsx(as.data.frame(get(sheet)), filepaths$xloutput, sheetName=sheet, row.names=FALSE, append=TRUE)
	}
}else{
	# View(bioSummary)
}