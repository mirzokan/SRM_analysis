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
library("agricolae")

# ----- Global control variables
xlOut <- TRUE # Toggle report output
graphOut <- TRUE # Toggle graph output
stripOut <- TRUE # Toggle stripchart output
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
			skyline="Mirzo_Skyline_Results.csv",
			xloutput="R_report.xlsx"
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
ds <- itranslate(ds, pepconc, "loq_ugml", "peptide")

# Create columns for IS concentration, trial, Precursor and Transition names
ds <- ds %>% 
	mutate(trial = as.numeric(str_match(.$replicateName, "rep([[:digit:]]+)")[,2])) %>% 
	mutate(precursorName = paste0(.$peptideMod, "_", .$precursorCharge, "+")) %>% 
	mutate(transitionName = paste0(.$precursorName, "-->", .$fragmentIon, "_", .$productCharge, "+")) %>% 
	mutate(transitionID = paste0(.$replicateName, "_", .$transitionName)) %>% 
	mutate(expID=paste(subjectID, condition)) %>% 
	mutate(expTrID=paste(.$expID, transitionName)) %>% 
	mutate(fmoluL = fmolug*stockfmolug)

# Create lists of unique values of concentrations, transitions and precursors
UfmoluL <- unique(ds$fmoluL)
Ufmolug <- unique(ds$fmolug)
Utransitions <- unique(ds$transitionName)
Uprecursors <- unique(ds$precursorName)
Upeptides <- unique(ds$peptide)
Uproteins <- unique(ds$protein)
Uisotope <- unique(ds$isotope)
UsubjectID <- unique(ds$subjectID)

#Reorder columns for convenience
ds <- ds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "precursorMz", "precursorCharge", "productMz", "productCharge", "fragmentIon", "isotope", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "retentionTime", "area", "background", "ratioLightToHeavy", "fmolug", "fmoluL", "loq_ugml", "ratioToStandard", "rdotp", "peakRank")]


# Remove Outliers
ds <- ds %>% 
	filter(!(subjectID %in% c("PCASV-013", "PCASV-032", "PCASV-020")))


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
tds <- tds[c("replicateName", "subjectID", "bioreplicate", "condition", "sampleType", "trial", "protein", "peptide", "peptideMod", "heavyPrecursorMz", "lightPrecursorMz", "heavyProductMz", "lightProductMz", "precursorCharge", "productCharge", "fragmentIon", "precursorName", "transitionName", "transitionID", "expID", "expTrID", "heavyRetentionTime", "lightRetentionTime", "fmolug", "fmoluL", "loq_ugml","heavyArea", "lightArea", "heavyBackground", "lightBackground", "heavyPeakRank", "lightPeakRank", "rdotp", "ratioLightToHeavy", "ratioToStandard")]

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
repDS <- itranslate(repDS, ds, "loq_ugml", "expTrID")
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
		fmoluL=mean(fmoluL),
		loq_ugml=mean(loq_ugml)) %>% 
	ungroup()

peptideDS <- peptideDS %>% 	
	mutate(meanLTH=ifelse(.$meanLTH==0, .$meanLightArea/.$meanHeavyArea, .$meanLTH)) %>% 
	select(-c(meanLightArea, meanHeavyArea)) %>% 
	mutate(concentration_fmul=meanLTH*fmoluL) %>% 
	mutate(molw=pepInfo$mw[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(concentration_ugml=concentration_fmul*molw*1e-6) %>% 
	mutate(pepID=paste0(subjectID, peptide)) %>% 
	mutate(protein=pepInfo$protein[match(.$peptide, pepInfo$peptide)]) %>% 
	mutate(stageSpecificity=pepInfo$stageSpecificity[match(.$peptide, pepInfo$peptide)])

# Sort for convenience
peptideDS <- peptideDS[c("pepID", "protein", "stageSpecificity", "peptide", "subjectID", "condition",  "meanLTH", "fmoluL", "concentration_fmul", "molw", "concentration_ugml", "loq_ugml")]



# ----- Ranking of patient samples by protein concentrations
rankDS <- peptideDS %>%
	select(condition, protein, subjectID, concentration_fmul) %>% 
	group_by(condition, protein) %>% 
	mutate(rank=rank(-concentration_fmul)) %>% 
	ungroup() %>%
	select(-c(concentration_fmul)) %>%  
	spread(protein, rank) %>%  
	mutate(rankSum=rowSums(.[3:22])) %>% 
	select(condition, subjectID, rankSum) %>% 
	arrange(condition, rankSum)


# ----- Base data wrangling ends


# Faceted stripchars for each peptide (with and without SubjectID labels)
if (stripOut & graphOut){
	for (target in Upeptides){
		for (j in c(TRUE, FALSE)){
			protTarget <- pepInfo$protein[match(target, pepInfo$peptide)]
			stagel <- pepInfo$stageSpecificity[match(target, pepInfo$peptide)]
			p <- peptideDS %>%
				filter(peptide==target) %>%  
				ggplot(aes(x=condition, y=concentration_ugml)) +
				geom_hline(aes(yintercept=loq_ugml), linetype="longdash")+
				# theme_bw()+
				geom_beeswarm(size=2, cex=2) +
				# geom_jitter(position=position_jitter(0.05)) +
				labs(title=sprintf("%s (%s)", protTarget, stagel), subtitle=sprintf("%s...", substr(target,1,5)))+
				scale_x_discrete(limits=c("PREvas", "POSTvas", "OA", "HSG", "MA", "SCO"))+
				xlab("Condition") + ylab(expression("Concentration ("*mu*"g/mL)"))+
				theme(text = element_text(size=15))+
				# theme(text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"))+
				scale_y_log10()+
				stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.2, color="red")
				# facet_wrap(~ protein, scales="free_y", nrow=4)
		

			if (j){
				p <- p + geom_text_repel(aes(label=subjectID), size=3, segment.alpha=0.2)
				filepref <- "g_tl_"
			}else{
				p <- p + geom_boxplot(aes(group=condition), alpha=0.5)
				filepref <- "g_nl_"
			}
			
			ggsave(paste(filepref, protTarget,"_beeswarm",".jpg", sep=""), plot=p, device="jpg", width=7, height=7)
		}	
	}
}

# ----- ANOVA

anovaDS <- peptideDS %>% 
	select(peptide, subjectID, condition, concentration_ugml) %>% 
	filter(condition %in% c("HSG", "MA", "OA", "POSTvas", "SCO"))

anovares <- list()
tukres <- list()
duncanres <- list()

for (target in Upeptides){
	protTarget <- pepInfo$protein[match(target, pepInfo$peptide)]
	subDS <- anovaDS %>% filter(peptide==target)
	
	# Anova
	aovout <- aov(concentration_ugml~condition, data=subDS)
	anovares[[protTarget]]  <- summary(aovout)[[1]][["Pr(>F)"]][[1]]

	# Tukey
	tukres[[protTarget]] <- TukeyHSD(aovout)

	# Duncan
	duncanres[[protTarget]] <- duncan.test(aovout, "condition", group=FALSE)
}


res <- anovares[anovares<0.05]
repep <- attr(res, "names")

tukp <- list()
duncanp <- list()
for (pep in repep){
	# Tukey
	rowds <- as.data.frame(tukres[repep][[pep]][["condition"]][,4])
	tukp[[pep]] <- data.frame(peptide=pep ,combinations=row.names(rowds), p=rowds[,1])
	
	tukpds <- do.call("rbind", tukp)
	rownames(tukpds) <- c()

	tukpds <- tukpds %>% 
		filter(p<0.05)

	# Duncan
	rowds <- duncanres[repep][[pep]][["comparison"]]
	duncanp[[pep]] <- data.frame(peptide=pep ,combinations=row.names(rowds), p=rowds$pvalue)

	duncanpds <- do.call("rbind", duncanp)
	rownames(duncanpds) <- c()

	duncanpds <- duncanpds %>% 
		filter(p<0.05)

}

# # Faceted stripchars for each peptide (with and without SubjectID labels)
# if (stripOut & graphOut){
# 	p <- peptideDS %>% 
# 		ggplot(aes(x=condition, y=concentration_ugml)) +
# 		# geom_hline(aes(yintercept=concentration_ugmlAtThr), linetype="longdash")+
# 		theme_bw()+
# 		geom_beeswarm(size=0.5, cex=1) +
# 		# geom_jitter(position=position_jitter(0.05)) +
# 		# labs(title=sprintf("%s (%s)", protTarget, stagel), subtitle=sprintf("%s...", substr(target,1,5)))+
# 		# scale_x_discrete(limits=c("PRE", "POST"))+
# 		xlab("Condition") + ylab(expression("Concentration ("*mu*"g/mL)"))+
# 		theme(aspect.ratio=0.8, text = element_text(size=8), axis.text = element_text(size=6, color="black"), axis.text.x = element_text(angle=90), axis.line=element_line(color="black"), panel.border=element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
# 		# stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, color="red", size=0.2)+
# 		scale_y_log10()+
# 		facet_wrap(~ protein, scales="free_y", nrow=4)
			
# 	ggsave(paste("g1_facet_beeswarm",".jpg", sep=""), plot=p, device="jpg", width=7, height=7*0.8)
# }




#------------- Output
if(xlOut){
		write.xlsx(as.data.frame(rtDS), filepaths$xloutput, sheetName="rtDS", row.names=FALSE)
	for (sheet in c("rankDS", "tukpds", "duncanpds")){
		write.xlsx(as.data.frame(get(sheet)), filepaths$xloutput, sheetName=sheet, row.names=FALSE, append=TRUE)
	}
}