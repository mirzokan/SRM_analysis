# Description

R script that helps to analyse SRM experiments for biomarker discovery. The script takes an output report from Skyline, wrangles the data to a more comprehensible form and calculates some vital statistics. Current features:

* Produces a volcano plot to find biomarker candidates with significant fold-change between treatments
* Uses logistical regression to find the best threshold value for each biomarker candidate
* Calculates biomarker performance parameters, such as accuracy, sensitivity, specificity, and area under receiver operating characteristic curve (AUROC)