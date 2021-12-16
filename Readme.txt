A new classification framework for single-cell RNA-seq data by combining the statistical method and deep learning

Code

Abstract

Genetic cell classification is extremely important no matter in medical or biological fields. All of the data processing and analysis for this paper were done in R and Python. As for Python, these packages are required tensorflow-gpu==1.15.0 keras==2.3.1 numpy==1.19.1 scanpy==1.3.7 scikit-learn==0.23.1. The version of Python should be greater than or equal to 3.6. As for R, edgeR, limma, PoiClaClu, reticulate, sSeq, stringr are required. 

Description

The files contain all codes of studies and analysis on datasets. 
Specifically, we make an explain for each document in the following.

###############  DeepZINB ######################

(1). Main.R---an r code for tutorial

(2). ZINBL.R ---an r code for some basic functions including the data processing 

(3). ZINB.py --- a python code for training the neural network 

(4). ZINBAutoEncoder.py--- a python code for architecture of  neural network 


###############  function ######################
(1). ZINBScore
       Description
	Combining with the trained parameters, it is used to calculate the discriminant score and obtain the labels of test set.

       Usage 
	ZINBScore(X_train_, y_train, X_test_, y_test, batch_size, py_path)

       Arguments 
	X_train_ 	   Matrix or data.The input training sample.

	y_train 	   A vector shows the class of the input sample

	X_test_ 	   Matrix or data.The input test sample.

	y_test 	   A vector shows the class of the output sample

	batch_size   Number of training batches of the neural network.	

	py_path 	    The path of python.exe.

      Value
	The prediction of X_test and y_test.

(2).  select2
       Description
	This function is used to select the top genes for the binary dataset.  For datasets with three or four classes, you can choose function select3 or select4.

        Usage 
	select2(dat, gene_no_list)

        Arguments 
	dat 	       Matrix or data. The  binary input dataset.

	gene_no_list  The number of top genes in the ranking.

        Value

	Matrix with the top genes in the ranking. 
   
      
############### datasets ######################

(1). realdata_xtr 
       Description
	This real data set is used as a training dataset, including 595 samples and 23654 genes, with three classes.(GSE113069) 

       Usage 
	A dataset contains 595 samples and 23654 genes.

       Source
	Cembrowski MS, Wang L, Lemire AL, Copeland M et al. The subiculum is a patchwork of discrete subregions. Elife 2018 Oct 30;7.

(2). realdata_xte
       Description
	This real data set is used as a test dataset, including 595 samples and 23654 genes, with three classes.(GSE113069) 

       Usage 
	A dataset contains 595 samples and 23654 genes.

       Source
	Cembrowski MS, Wang L, Lemire AL, Copeland M et al. The subiculum is a patchwork of discrete subregions. Elife 2018 Oct 30;7.

(3). realdata_ytr 
       Description
	Category label of the training real dataset.

       Source
	Cembrowski MS, Wang L, Lemire AL, Copeland M et al. The subiculum is a patchwork of discrete subregions. Elife 2018 Oct 30;7.

(4). realdata_yte
       Description
	Category label of the test real dataset 

       Source
	Cembrowski MS, Wang L, Lemire AL, Copeland M et al. The subiculum is a patchwork of discrete subregions. Elife 2018 Oct 30;7.




Optional Information

The following R packages are necessary to successfully run the codes:

(1) reticulate (https://cran.r-project.org/web/packages/reticulate)
(2) PoiClaClu (https://cran.rstudio.com/web/packages/PoiClaClu)
(3) sSeq (http://www.bioconductor.org/packages/release/bioc/html/sSeq.html)
(4) limma (http://www.bioconductor.org/packages/release/bioc/html/limma.html)
(5) edgeR(http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
(6) stringr(https://cran.r-project.org/web/packages/stringr)

The following python packages are necessary to successfully run the codes:
(7)tensorflow-gpu==1.15.0 
(8)keras==2.3.1 
(9)numpy==1.19.1
(10)scanpy==1.3.7
(11)scikit-learn==0.23.1

Instructions for Use

The code provided can be used to reproduce the results of DeepZINB.

First, we need to install and load the following R and Python packages.

Second, to reproduce the results, put the files provided in the 
working directory and execute the following commands in R:source("ZINBL.R")

Third, 

run: Main.R






