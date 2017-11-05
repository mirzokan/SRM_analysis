setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(warn=1)
options(stringsAsFactors = FALSE)

rm(list=ls())

library("stringr")
library("xlsx")
library("magrittr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggrepel")
library("ggbeeswarm")
library("rlang")
library("pwr")
library("boot")
library("ROCR")
library("agricolae")
library("dunn.test")

# ----- Global control variables
xl_out <- FALSE
graph_out <- FALSE
lod_out <- TRUE
strip_out <- TRUE
cormat_out <- TRUE
volcano_out <- FALSE
logit_out <- FALSE
auc_out <- FALSE
remove_botched <- TRUE

thrFC <- 1 # Threshold log2 fold-change 
thrAlpha <- 0.05 # Alpha threshold for comparative statistics 
thrLogit <- 0.5 # Probability threshold for the logistic regression 
optPower <- 0.8 # Optimal statistical power

stock_fmol_ug <- 40

lib_paths <- list(
			lib_pep_info="C:\\Users\\lem\\Dropbox\\LTRI\\2 Tools\\peptide_info.csv",
			lib_stage_info="C:\\Users\\lem\\Dropbox\\LTRI\\2 Tools\\stage_info.csv"
			)

# ----- Functions
source("C:\\Users\\lem\\Documents\\Git\\SRM_analysis\\srm_functions.R")


# ----- Load libraries
lib_pep_info <- read.csv(lib_paths$lib_pep_info)
lib_stage_info <- read.csv(lib_paths$lib_stage_info)



# ------ Wranglign one Set
set1 <- "patient_samples"

patient_samples_files <- list(
			skyline="Sample_Skyline_Results.csv",
			lib_pep_conc="peptconc.csv"
			)

patient_samples_botched_list <- c("PCASV-013", "PCASV-032", "PCASV-020")

df1 <- load_skyline(patient_samples_files, set=set1, remove_botched=TRUE, botched_list=patient_samples_botched_list)

df1_tidy <- pivot(df1, id="measurement_transitionID", key="isotope", values=c("area", "background", "max_height", "transition_rank", "retention_time", "precursor_mz", "product_mz"))

list_peptides <- unique(df1_tidy$peptide)

# ----- Base Wrangling Ends





#------------- Replicate averaging
replicates <- replicate_average(df1_tidy)

df1_replicate_means <- replicates[["means"]]
df1_replicate_cv <- replicates[["cv"]]

# ----- Peptide concentrations summary
df1_concentration <- peptide_concentrations(df1_replicate_means, df1_tidy)

#------------- Retention Time Analysis
report_retention_time <- rt_analysis(df1_tidy)


# ----- Lod
if (graph_out & lod_out){
	for (target in list_peptides){
		p <- lod_graph(df1_concentration, target,save=TRUE)
	}
}


# ------ Stripchars
if (strip_out & graph_out){
	for (target in list_peptides){
		
		strip_chart(df1_concentration, target, "condition", categories=c("PREvas", "POSTvas", "OA", "NOA"), save=TRUE, boxplot=TRUE) 
		# strip_chart(df1_concentration, target, "condition", categories=c("PREvas", "POSTvas", "OA", "NOA"), save=TRUE, point_label=TRUE) 
		strip_chart(df1_concentration, target, "histology", save=TRUE, boxplot=TRUE) 
		strip_chart(df1_concentration, target, "tese", categories=c("sperm", "spermatids", "no sperm"), save=TRUE, boxplot=TRUE) 
	}
	strip_chart(df1_concentration, "facet", "condition", categories=c("PREvas", "POSTvas", "OA", "NOA"), save=TRUE) 
}

# ----- Protein concentration correlation
if (graph_out & cormat_out){
	cormat <- protein_cross_correlation(df1_concentration, save=TRUE)	
}

# ----- Variance Analysis
noa_anova <- df1_concentration %>% 
	select(peptide, subjectID, histology, concentration_ug_ml) %>% 
	filter(!histology=="unknown") %>%
	anova_posthoc("histology") 


noa_kruskal <- df1_concentration %>% 
	select(peptide, subjectID, histology, concentration_ug_ml) %>% 
	filter(!histology=="unknown") %>%
	kruskal_dunn("histology") 


condition_kruskal <- df1_concentration %>% 
	select(peptide, subjectID, condition, concentration_ug_ml) %>% 
	kruskal_dunn("condition") 	

# ----- Epididymal

epi_anova <- df1_concentration %>% 
	filter(stage_specificity=="Epididymus") %>% 
	mutate(condition=ifelse(condition %in% c("POSTvas"), "OA", condition)) %>% 
	anova_posthoc("condition")

epi_kruskal <- df1_concentration %>% 
	filter(stage_specificity=="Epididymus") %>% 
	mutate(condition=ifelse(condition %in% c("POSTvas"), "OA", condition)) %>% 
	kruskal_dunn("condition") 	

# ----- Ranking of patient samples by protein concentrations
report_sample_ranks <- df1_concentration %>%
	select(protein, sampleID, concentration_fm_ul) %>% 
	group_by(protein) %>% 
	mutate(rank=rank(-concentration_fm_ul)) %>% 
	ungroup() %>%
	select(-concentration_fm_ul) %>%  
	spread(protein, rank) %>%
	mutate(rankSum=rowSums(.[2:21])) %>% 
	mutate(rank=dense_rank(rankSum)) %>%  
	arrange(rankSum)


# ----- Epi AUROC analysis

# epiDS %<>% 
# 	filter(!condition=="PREvas") %>% 
# 	mutate(conclass=ifelse(condition=="OA", 1L, 0L))

# set.seed(1354)
# epiAucDS <- data.frame()

# for (protTarget in epirepep){
# 	target <- lib_pep_info$peptide[match(protTarget, lib_pep_info$protein)]

# 	#Split into training and test set
# 	fullset <- epiDS %>% 
# 		filter(peptide==target) 

# 	train <- fullset %>% 
# 		sample_frac(0.3, replace=FALSE)

# 	test <- fullset %>% 
# 		setdiff(train)

# 	#Train logistic regression classifier
# 	model <- glm(formula=conclass~concentration_ug_ml, data=fullset, family=binomial)
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
# 	if (graph_out & logit_out){
# 		p <- fullset %>% 
# 			ggplot(aes(x=concentration_ug_ml, y=conclass))+
# 			geom_jitter(size=2, width = 0, height = 0.02)+
# 			labs(title=paste0(protTarget, " (", substr(target, 1, 5), "...)"), subtitle=sprintf("Accuracy: %s, Sensitivity: %s, Specificity: %s", accuracy, sensitivity, specificity))+
# 			theme(text = element_text(size=15), legend.position="none")+
# 			xlab(expression("Concentration ("*mu*"g/mL)")) + ylab("Probability of OA classification")+
# 			stat_function(fun=function(x) inv.logit(modelInt+modelSlope*x), geom="line", colour="deepskyblue4", size=1)+
# 			geom_vline(xintercept=modelConcAtThr, linetype="longdash")+
# 			scale_x_log10()

# 		ggsave(paste0("g4_", protTarget, "_epi_classifier.jpg"), plot=p, device="jpg")
# 	}

# 	#ROCR analysis
# 	pr <- prediction(pred, fullset$conclass)
# 	prf <- performance(pr, measure = "tpr", x.measure = "fpr")
# 	accprf <- performance(pr, measure = "acc")

# 	# Make plot data
# 	plotdat <- data.frame(fp=prf@x.values[[1]],tp=prf@y.values[[1]],cut=prf@alpha.values[[1]], acc=accprf@y.values[[1]])
	 
# 	auc <- performance(pr, measure = "auc")
# 	auc <- round(auc@y.values[[1]],2)

	
# 	if (graph_out){
# 		# AUC plots
# 		if (auc_out){
# 			p <- plotdat %>% 
# 				ggplot(aes(x=fp, y=tp, colour=cut))+
# 				geom_line()+ geom_point()+
# 				scale_colour_gradient2(low = "green3", mid = "orange3",
# 					high = "red4", midpoint = 0.5, name="Threshold")+
# 				labs(title=paste0(protTarget, " (", substr(target, 1, 5), "...)"), subtitle=paste0("AUC: ", auc))+
# 				theme(text = element_text(size=15), legend.title = element_text(size=10))+
# 				xlab("False Positive Rate") + ylab("True Positive Rate")+
# 				stat_function(fun=function(x) x, geom="line", colour="black", linetype="dotted")+
# 				geom_text_repel(aes(label=round(acc,2)))

# 			ggsave(paste0("g4_", protTarget, "_AUC.jpg"), plot=p, device="jpg")
# 		}	
# 	}

# 	drow <- data.frame(peptide=target, concentration_ug_mlAtThr=modelConcAtThr,
# 			accuracy=accuracy, sensitivity=sensitivity , specificity=specificity, auc=auc)
# 	epiAucDS <- rbind(epiAucDS, drow)


# 	# if (graph_out & strip_out){
# 	if (FALSE){
# 		for (target in epipeps){
# 			for (j in c(TRUE, FALSE)){
# 				protTarget <- lib_pep_info$protein[match(target, lib_pep_info$peptide)]
# 				stagel <- lib_pep_info$stage_specificity[match(target, lib_pep_info$peptide)]

# 				subD <- epiDS %>%
# 					filter(peptide==target)

# 				thresh <- ifelse(target %in% epiAucDS$peptide, epiAucDS$concentration_ug_mlAtThr[epiAucDS$peptide==target], 0)

# 				p <- subD %>% 
# 					ggplot(aes(x=condition, y=concentration_ug_ml)) +
# 					geom_hline(aes(yintercept=lod), linetype="longdash")+
# 					geom_hline(aes(yintercept=thresh), color="red", linetype="longdash")+
# 					# theme_bw()+
# 					geom_beeswarm(size=2, cex=2) +
# 					# geom_jitter(position=position_jitter(0.05)) +
# 					labs(title=sprintf("%s (%s)", protTarget, stagel), subtitle=sprintf("%s...", substr(target,1,5)))+
# 					scale_x_discrete(limits=c("OA", "NOA"))+
# 					xlab("Condition") + ylab(expression("Concentration ("*mu*"g/mL)"))+
# 					theme(text = element_text(size=15))+
# 					# theme(text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"))+
# 					scale_y_log10()+
# 					stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.2, color="red")
# 					# facet_wrap(~ protein, scales="free_y", nrow=4)
			

# 				if (j){
# 					p <- p + geom_text_repel(data=subD[subD$concentration_ug_ml>subD$lod,], aes(label=run), size=3, segment.alpha=0.2)
# 					filepref <- "g3_tl_"
# 				}else{
# 					# p <- p + geom_boxplot(aes(group=condition), alpha=0.5)
# 					filepref <- "g3_nl_"
# 				}
				
# 				ggsave(paste0(filepref, protTarget,"_epi_beeswarm",".jpg"), plot=p, device="jpg", width=7, height=7)
# 			}	
# 		}
# 	}
# }






# # ------------- Output
# if(xl_out){
# 		write.xlsx(as.data.frame(report_retention_time), lib_paths$xloutput, sheetName="report_retention_time", row.names=FALSE)
# 	for (sheet in c("rankDS", "tukpds", "duncanpds")){
# 		write.xlsx(as.data.frame(get(sheet)), lib_paths$xloutput, sheetName=sheet, row.names=FALSE, append=TRUE)
# 	}
# }
