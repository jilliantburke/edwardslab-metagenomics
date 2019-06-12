# edwardslab-metagenomics
This respository contains the codes used to create a predictive model of the source location (freshwater, marine water, human, animal, etc.) of metagenomic data. These metagenomes can be either RNAseq or shotgun DNAseq, but not amplicon/16s rRNA seq. 

This README document is intended to provide a follow-along style protocol for re-creating this model. Relevant codes are included in this repository. For more detailed explinations of the algorithm and reasoning behind the use of this particular method, see the pdf file in this repo called "thesis_final_version.docx." This document also contains the expected output of the graphs from all of included R scripts. Email me if you have questions or comments: jillian.t.burke@gmail.com


Summary of steps to creating the model

Preprocessing of training set: 
0. Process the training set to include a fairly balanced number of samples from each source category

Step 1 - Generate input features: 
1A. Run the PARTIE tool
1B. Estimates of the human DNA percentage
1C. Run the FOCUS tool & calculate the percentages of bacterial reads in the PARTIE coming from each source category
1D. Parse the outputs of steps 1A-1C to include only the metagenomes from step 0. (those in the training set) and combine into a single table. 

Step 2 - Create Random Forest Model:
2A. Scale the data from step 1D
2B. Generate a random forest classifier using the scaled input data

Step 3 - Summarize model (visalization):
3A. Obtain the error rates for each level of the model (sources) and create a plot of these 

Step 4 - optional PCA analysis: 
4A. Scale the data from step 1D and perform PCA 
4B. Isolate single levels of the model (sources) and plot these  


Detailed Steps to creating the model

Preprocessing of training set: 
0. Process the training set to include a fairly balanced number of samples from each source category
To do this: 

Step 1 - Generate input features: 
1A. Run the PARTIE tool:
A spreadsheet containing the PARTIE output for many metagenomes from the SRA database is located here:
https://raw.githubusercontent.com/linsalrob/partie/master/SRA_PARTIE_DATA.txt

This spreadsheet was created Rob Edwards, ask him for more detail, or see the PARTIE documentation manual.

1B. Estimates of the human DNA percentage
1C. Run the FOCUS tool & calculate the percentages of bacterial reads in the PARTIE coming from each source category
1D. Parse the outputs of steps 1A-1C to include only the metagenomes from step 0. (those in the training set) and combine into a single table. 

Step 2 - Create Random Forest Model and visualize results :
2A. Scale the data from step 1D
2B. Generate a random forest classifier using the scaled input data
2C. Obtain the error rates for each level of the model (sources) and create a plot of these 

These steps are performed using the script "random_forest_unbalanced.R";the input file here is the output of step 1D, the output is a picture file (png) of the bargraph of the within group error rates. You will need to change the output and input file paths manually in the code.

Step 3 - optional PCA analysis: 
3A. Scale the data from step 1D and perform PCA 
3B. Isolate single levels of the model (sources) and plot these  
These steps are performed by running the code: "pca_by_environment.R";  the input file here is the output of step 1D, the output is a picture file (png) of the PCA dotplot. You will need to change the output and input file paths as well as the name of the source environment that you wish to plot manually in the code. 



