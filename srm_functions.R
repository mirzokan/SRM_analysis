

# Equation for hyperbolic cut-off threshold in volcano plots
cthresh <- function(x){
	c  <-  0.5
	c/(abs(x-thrFC))+(-log10(thrAlpha))
} 

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

firstcap <- function(x) {
	substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	return(x)
}

pop_columns <- function(df, columns){
	columns_reversed <- rev(columns)
	for (column in columns_reversed){
		column_s <- sym(column)
		df %<>% select(!!column_s, everything()) 
	}
	return(df)
}

pivot <- function(df, id, key, values, rename=TRUE){
	id_sym <- sym(id)
	key_sym <- sym(key)
	key_name <- quo_name(key)

	remove_columns <- syms(c(key, values))
	new_df <- df %>% select(-c(!!!remove_columns))

	for (value in values){
		value_sym <- sym(value)

		intermediate_df <- df %>%
			select(!!id_sym, !!key_sym, !!value_sym)
		
		if(rename){
			intermediate_df %<>% mutate(!!key_name:=paste0((!!key_sym), "_", value))
		}
		
		new_df <- intermediate_df %>% 
			spread(key=(!!key_sym), value=(!!value_sym)) %>% 
			left_join(new_df, ., by=id)
	}

	new_df %<>% distinct(!!id_sym, .keep_all=TRUE)

	return(new_df)
}

#  ----- Wrangling

load_skyline <- function(file_paths, remove_botched=FALSE, botched_list=c()){

	df <- read.csv(file_paths$skyline, na.strings="#N/A")
	lib_pep_conc <- read.csv(file_paths$lib_pep_conc)

	# ----- Data wrangling
	df %<>% rename(
		replicate_name = Replicate.Name,
		fmol_ul = Analyte.Concentration,
		peptide = Peptide.Sequence,
		peptide_mods = Peptide.Modified.Sequence,
		sample_type = Sample.Type,
		isotope = Isotope.Label.Type,
		precursor_mz = Precursor.Mz,
		precursor_charge = Precursor.Charge,
		product_mz = Product.Mz,
		product_charge = Product.Charge,
		fragment_ion = Fragment.Ion,
		retention_time = Retention.Time,
		area = Area,
		max_height = Max.Height,
		background = Background,
		lth = RatioLightToHeavy,
		rdotp = DotProductLightToHeavy,
		transition_rank = Peak.Rank)

	df <- lib_pep_info %>% 
		select(peptide, protein) %>% 
		left_join(df, ., by="peptide")

	df <- lib_pep_conc %>% 
		select(peptide, fmol_ug, loq_ug_ml, lod) %>% 
		left_join(df, ., by="peptide")

	df %<>%  
		mutate(run = as.numeric(str_match(.$replicate_name, "^([[:digit:]]+)_")[,2])) %>% 
		mutate(trial = as.numeric(str_match(.$replicate_name, "rep([[:digit:]]+)")[,2])) %>% 
		mutate(precursor_name = paste0(.$peptide_mod, "_", .$precursor_charge, "+")) %>% 
		mutate(transition_name = paste0(.$precursor_name, "-->", .$fragment_ion, "_", .$product_charge, "+")) %>% 
		mutate(sampleID = paste(subjectID, condition, timepoint, sep="_")) %>% 
		mutate(experimentID = paste(.$sampleID, curvepoint, sep="_")) %>% 
		mutate(experiment_transitionID = paste0(experimentID, "_", .$transition_name)) %>% 
		mutate(experiment_peptideID = paste0(experimentID, "_", .$peptide)) %>% 
		mutate(transitionID = paste0(run, "_", .$transition_name)) %>% 
		mutate(fmol_ul = fmol_ug*stock_fmol_ug)

	# Remove botched samples
	if (remove_botched){
		df %<>%  
			filter(!(subjectID %in% botched_list))
	}

	df %<>%  
		mutate(area=ifelse(is.na(area), 0, area)) %>% 
		mutate(background=ifelse(background==0, 1, background)) %>% 
		mutate(background=ifelse(is.na(background), 1, background)) %>% 
		mutate(lth=ifelse(is.na(lth), 0, lth)) %>% 
		mutate(rdotp=ifelse(is.na(rdotp), 0, rdotp)) %>% 
		mutate(max_height=ifelse(is.na(max_height), 0, max_height)) %>% 
		mutate(transition_rank=ifelse(is.na(transition_rank), 0, transition_rank)) %>% 
		mutate(retention_time=ifelse(is.na(retention_time), 0, retention_time))

	return(df)
}


replicate_average <- function(df_tidy){
	df_replicate_means <- df_tidy %>% 
		group_by(experiment_transitionID) %>% 
		summarise(
			heavy_retention_time = mean(heavy_retention_time),
			light_retention_time = mean(light_retention_time),
			heavy_area = mean(heavy_area),
			light_area = mean(light_area),
			heavy_background = mean(heavy_background),
			light_background = mean(light_background),
			heavy_max_height = mean(heavy_max_height),
			light_max_height = mean(light_max_height),
			rdotp = mean(rdotp),
			lth = mean(lth))

	df_replicate_means <- df_tidy %>% 
		select(experiment_transitionID, experiment_peptideID, experimentID, heavy_transition_rank, peptide) %>% 
		distinct(experiment_transitionID, .keep_all=TRUE) %>% 
		left_join(df_replicate_means, ., by="experiment_transitionID") 

	df_replicate_cv <- df_tidy %>% 
		group_by(experiment_transitionID) %>% 
		summarise(
			heavy_retention_time = percent(sd(heavy_retention_time)/mean(heavy_retention_time)),
			light_retention_time = percent(sd(light_retention_time)/mean(light_retention_time)),
			heavy_area = percent(sd(heavy_area)/mean(heavy_area)),
			light_area = percent(sd(light_area)/mean(light_area)),
			heavy_background = percent(sd(heavy_background)/mean(heavy_background)),
			light_background = percent(sd(light_background)/mean(light_background)),
			heavy_max_height = percent(sd(heavy_max_height)/mean(heavy_max_height)),
			light_max_height = percent(sd(light_max_height)/mean(light_max_height)),
			rdotp = percent(sd(rdotp)/mean(rdotp)),
			lth = percent(sd(lth)/mean(lth)))

		output <- list(
					means=df_replicate_means,
					cv=df_replicate_cv
			)
		return(output)
}


peptide_concentrations <- function(df_replicate_means){
	df_concentration <- df_replicate_means %>% 
		filter(heavy_transition_rank==1) %>% 
		group_by(experiment_peptideID) %>% 
		summarise(
			lth = mean(lth),
			heavy_area = mean(heavy_area),
			light_area = mean(light_area),
			heavy_background = mean(heavy_background),
			light_background = mean(light_background),
			heavy_max_height = mean(heavy_max_height),
			light_max_height = mean(light_max_height),
			rdotp = mean(rdotp),
			lh_rt_diff = abs(mean(light_retention_time) - mean(heavy_retention_time))) %>% 
		ungroup()


	df_concentration <- df %>% 
		select(run, condition, protein, histology, tese, texlevel, subjectID, sampleID, peptide, experiment_peptideID, fmol_ul, loq_ug_ml, lod) %>% 
		distinct(experiment_peptideID, .keep_all=TRUE) %>% 
		left_join(df_concentration, ., by="experiment_peptideID") 


		df_concentration %<>%  
		mutate(heavy_snr=heavy_area/heavy_background) %>% 
		mutate(light_snr=light_area/light_background) %>% 
		select(-c(light_area, heavy_area)) %>% 
		mutate(concentration_fm_ul=lth*fmol_ul) %>% 
		mutate(molecular_weigth=lib_pep_info$mw[match(.$peptide, lib_pep_info$peptide)]) %>% 
		mutate(concentration_ug_ml=concentration_fm_ul*molecular_weigth*1e-6) %>% 
		mutate(protein=lib_pep_info$protein[match(.$peptide, lib_pep_info$peptide)]) %>% 
		mutate(stage_specificity=lib_pep_info$stage_specificity[match(.$peptide, lib_pep_info$peptide)]) %>%
		mutate(lh_rt_diffRel=ifelse(is.na(lh_rt_diff), 1, lh_rt_diff/max(.$lh_rt_diff))) %>% 
		pop_columns(c("run", "peptide", "protein", "stage_specificity", "subjectID", "condition", "histology", "tese", "texlevel", "rdotp", "concentration_fm_ul", "concentration_ug_ml")) %>% 
		arrange(peptide, -concentration_ug_ml)
	return(df_concentration)
}



#  ----- Analysis

lod_graph <- function(df_concentration, target, save=FALSE){
	
	target_protein <- lib_pep_info$protein[match(target, lib_pep_info$peptide)]
	stage_label <- lib_pep_info$stage_specificity[match(target, lib_pep_info$peptide)]

	subD <- df_concentration %>% 
		filter(peptide==target)

	p <- subD %>% 
		ggplot(aes(x=concentration_ug_ml, y=rdotp, color=lh_rt_diff))+
		geom_point()+
		labs(title=sprintf("%s (%s)", target_protein, stage_label), subtitle=sprintf("%s...", substr(target,1,5)))+
		geom_vline(aes(xintercept=lod), linetype="longdash")+
		scale_x_log10()+ scale_y_log10(breaks=c(0.2, 0.3, 0.5, 0.8, 0.9, 1), limits = c(0.2, 1.3))+
		xlab(expression("Concentration ("*mu*"g/mL)"))+
		ylab("dot Product")+
		scale_color_gradientn(colours=c("blue", "red", "grey"), values=c(0, 0.1, 0.3), limits=c(0,1), name="RT Difference\n")

		if(save){
			ggsave(paste0("g_lod_", target_protein, ".jpg"), plot=p, device="jpg", width=7, height=7)
		}

		return(p)
}

strip_chart <- function(df_concentration, target, x, categories=c(), bw=FALSE, point_label=FALSE, boxplot=FALSE, save=FALSE){
	
	x_sym <- sym(x)
	
	subD <- df_concentration %>% 
		filter(!((!!x_sym)=="unknown"))

	if(!target=="facet"){
	
		target_protein <- lib_pep_info$protein[match(target, lib_pep_info$peptide)]
		stage_label <- lib_pep_info$stage_specificity[match(target, lib_pep_info$peptide)]

		subD %<>% 
			filter(peptide==target)
	}

	p <- subD %>% 
		ggplot(aes_string(x=x, y="concentration_ug_ml")) +
		geom_hline(aes(yintercept=lod), linetype="longdash")+
		geom_beeswarm(size=2, cex=2, groupOnX=TRUE) +
		xlab(firstcap(x)) + ylab(expression("Concentration ("*mu*"g/mL)"))+
		theme(text = element_text(size=15))+
		theme(axis.text.x = element_text(angle = 90, hjust = 1))+
		scale_y_log10()+
		stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.2, color="red")
		

		if (target=="facet"){
			p <- p + facet_wrap(~ protein, scales="free_y", nrow=4)
			target_protein <- "facet_"
			g_width <- 14
			g_height <- 7
		} else{
			p <- p + labs(title=sprintf("%s (%s)", target_protein, stage_label), subtitle=sprintf("%s...", substr(target,1,5)))
			g_width <- 7
			g_height <- 7
		}

		if (!length(categories)==0){
			p <- p + scale_x_discrete(limits=categories)
		}

		if (bw){
			p <- p + theme_bw()+
				theme(text = element_text(size=15), axis.text = element_text(size=13, color="black"), axis.line=element_line(color="black"), panel.border=element_rect(color="black"))
		}

		if (point_label){
			p <- p + geom_text_repel(data=subD, aes(label=run), size=3, segment.alpha=0.2)
			filepref <- "label_"
		} else {
		filepref <- "nolabel_"
		}

		if (boxplot){
			p <- p + geom_boxplot(aes_string(group=x), alpha=0.5)
		}
		
		if (save){
			ggsave(paste0("g_stripchart_", x, "_", filepref, target_protein,".jpg"), plot=p, device="jpg", width=g_width, height=g_height)
		}

		return(p)
}

mlod <- function(df_concentration, pepnum, lodp=40, zoom=c(0.2,1.3), lab=0, xzoom=NULL){
	target <- list_peptides[[pepnum]]

	subD <- df_concentration %>% 
		filter(peptide==target)

	nlod <- subD$concentration_ug_ml[subD$run==lodp]

	p <- lod_graph(target)
	
	p <- p + geom_vline(xintercept=nlod, color="red", linetype="longdash")+
	coord_cartesian(xlim = xzoom, ylim = zoom)

	if (lab==1){
		# label all
		p <- p + geom_text_repel(data=subD, aes(label=run))
	} else if (lab==2){
		# label lod and new lod
		p <- p + geom_text_repel(data=subD[(round(subD$concentration_ug_ml,4)==round(subD$lod,4))|(subD$concentration_ug_ml==nlod),], aes(label=run))
	} else if (lab==3){
		# label lod, new lod, rtdiff less than 5%
		p <- p + geom_text_repel(data=subD[(round(subD$concentration_ug_ml,4)==round(subD$lod,4))|(subD$concentration_ug_ml==nlod)|(subD$lh_rt_diff<0.05),], aes(label=run))
	}

	print(p)
	print(nlod)
}

protein_cross_correlation <- function(df_concentration, save=FALSE){
	cormat  <-  df_concentration %>% 
		select(run, protein, concentration_ug_ml) %>% 
		pivot(., id="run", key="protein", values="concentration_ug_ml", rename=FALSE) %>% 
		select(-run) %>%
		round(2) %>%
		cor() %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("protein") %>% 
		gather(combination, correlation, -protein) %>% 
		arrange(protein, desc(combination))

	p <- ggplot(data = cormat, aes(x=protein, y=combination, fill=correlation)) + 
	  geom_tile()+
	  geom_text(aes(label=round(correlation, 2)), size=2) +
	  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
	     midpoint = 0, limit = c(-1,1))+
	  theme(axis.text.x = element_text(angle = 90, hjust = 1))

	if(save){
		ggsave("g_correlation_matrix.jpg", plot=p, device="jpg", width=7, height=7)
	} else {
		return(p)
	}
}


anova_posthoc <- function(df, group_by){

	gb_sym <- sym(group_by)
	peptide_list <- unique(df$peptide)

	subDf <- df %>% 
		select(peptide, subjectID, !!gb_sym, concentration_ug_ml) 

	duncan_model <- list()
	tukey_model <- list()

	anova_pvalues <- list()
	
	duncan_results <- data.frame()
	tukey_results <- data.frame()

	for (target in peptide_list){
		target_protein <- lib_pep_info$protein[match(target, lib_pep_info$peptide)]
		pepSubDf <- subDf %>% filter(peptide==target)
		
		# Anova
		anova_model <- aov(as.formula(paste0("concentration_ug_ml~", group_by)), data=pepSubDf)
		
		# Extract anova p-values
		anova_pvalues[[target_protein]]  <- summary(anova_model)[[1]][["Pr(>F)"]][[1]]

		# run Duncan
		duncan_model[[target_protein]] <- duncan.test(anova_model, group_by, group=FALSE)
		
		# run Tukey
		tukey_model[[target_protein]] <- TukeyHSD(anova_model)
	}

	significant_anova <- anova_pvalues[anova_pvalues<thrAlpha]
	significant_peptides <- attr(significant_anova, "names")

	
	for (target in significant_peptides){
		# Duncan
		duncan_results <- duncan_model[significant_peptides][[target]][["comparison"]] %>% 
			as_tibble %>% 
			rename(p=pvalue) %>% 
			select(p) %>% 
			tibble::rownames_to_column("combinations") %>%
			mutate(peptide=target) %>% 
			pop_columns("peptide") %>% 
			bind_rows(duncan_results,.)

		# Tukey
		tukey_results <- tukey_model[significant_peptides][[target]][[group_by]][,4] %>%
			as_tibble() %>%
			rename(p=value) %>% 
			tibble::rownames_to_column("combinations") %>%
			mutate(peptide=target) %>% 
			pop_columns("peptide") %>% 
			bind_rows(tukey_results,.)
	}

	output <- list(anova=significant_anova, duncan=duncan_results, tukey=tukey_results)
	return(output)
}

kruskal_dunn <- function(df, group_by){

	gb_sym <- sym(group_by)
	peptide_list <- unique(df$peptide)

	df[group_by] <- as.factor(df[[group_by]])

	subDf <- df %>% 
		select(peptide, subjectID, !!gb_sym, concentration_ug_ml) 

	dunn_model <- list()

	kruskal_pvalues <- list()
	
	dunn_results <- data.frame()

	for (target in peptide_list){
		target_protein <- lib_pep_info$protein[match(target, lib_pep_info$peptide)]
		pepSubDf <- subDf %>% filter(peptide==target)
		
		# Kruskal-Wallis
		kruskal_model <- kruskal.test(as.formula(paste0("concentration_ug_ml~", group_by)), data=pepSubDf)
		# Extract kruskal p-values
		kruskal_pvalues[[target_protein]]  <- kruskal_model[["p.value"]]


		# Dunn
		invisible(capture.output(dunn_model[[target_protein]] <- dunn.test(pepSubDf$concentration_ug_ml, g=pepSubDf[[group_by]], kw=FALSE, method="bh", table=FALSE, list=FALSE)))
	}

	significant_kruskal <- kruskal_pvalues[kruskal_pvalues<thrAlpha]
	significant_peptides <- attr(significant_kruskal, "names")

	
	for (target in significant_peptides){

		# Dunn
		dunn_results <- data.frame(peptide=target, combinations=dunn_model[[target]][["comparisons"]], p=dunn_model[[target]][["P.adjusted"]]) %>% 
		bind_rows(dunn_results,.)
	}

	out <- list(kruskal=significant_kruskal, dunn=dunn_results)
	return(out)
}




# Diagnostic plots

normal_hist <- function(df, var, binwidth=1){
	p <- ggplot(df, aes_string(var)) +
		geom_histogram(aes(y=..density..), binwidth=binwidth)+
		stat_function(fun=dnorm, args=list(mean = mean(df[[var]], na.rm=TRUE), sd = sd(df[[var]], na.rm=TRUE)))

	return(p)
}

qq_plot <- function(df, var, colour=NULL){
	p <- ggplot(df, aes_string(sample=var, colour=colour))+
		stat_qq()

	return(p)
}