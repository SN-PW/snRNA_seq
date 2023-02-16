Reanalysis of expressions of Asct2 (Slc1a5), Smct1 (Slc5a8) and Smct2 (Slc5a12) in renal proximal tubules at the single cell level. 
Expressions of Asct2, Smct1, and Smct2 in early Ischemia Reperfusion Injury (4 h and 12 h IRI) were compared to the controls.

The dataset of single nucleus RNA sequencing (snRNA-seq) was originally published by Kirita et al. PNAS 117(27)15874 – 15883 (https://doi.org/10.1073/pnas.2005477117). 

The corresponding data available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139107 include:    
“GSE139107_MouseIRI_4hours.dge.txt”   
“GSE139107_MouseIRI_12hours.dge.txt”   
“GSE139107_MouseIRI_control.dge.txt”    
“GSE139107_MouseIRI.metadata.txt”

To create Fig. S7, “20230209_PTS_separation.ipynb” is executed using the following files to make “SpeedTest” files.
-Extract_control_GSE139107_MouseIRI.metadata.txt
-Extract_4hours_GSE139107_MouseIRI.metadata.txt
-Extract_12hours_GSE139107_MouseIRI.metadata.txt
Subsequently, exported “SpeedTest” files should be used in the R script file named “20230210_SMCTs_Integrated_Dotplot.R”.
