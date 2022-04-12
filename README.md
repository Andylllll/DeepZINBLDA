# `LazAE`

> a deep learning framework with statistical method to classify RNA-seq data.

## Description

Genetic cell classification is extremely important no matter in medical or biological fields. We design an antoencoder network to learn the parameters of the model, and use the zero-inflated negative binomial distribution as the loss function for better feature selection. Then we put the trained parameters into the zero-inflated negative binomial linear discriminant model, which based on the Bayesian model. 


## Requirements

It is required to install the following dependencies in order to be able to run the code of scDLC

- [Anaconda3](https://www.anaconda.com/products/individual)  
- [R>=3.6.0](https://cran.r-project.org/)  
- [python 3](https://www.python.org/downloads/)  
  [sklearn](https://pypi.org/project/sklearn/0.0/)
  [numpy 1.19.1](https://pypi.org/project/numpy/1.19.1/)
  [tensorflow-gpu 1.15.0](https://pypi.org/project/tensorflow-gpu/1.15.0/)
  [keras 2.3.1](https://pypi.org/project/keras/2.3.1/)
  [scanpy 1.3.7](https://pypi.org/project/scanpy/1.3.7/)
 - [R>=4.1.0](https://www.r-project.org/)  
  [reticulate](https://cran.r-project.org/web/packages/reticulate)
  [PoiClaClu](https://cran.rstudio.com/web/packages/PoiClaClu)
  [sSeq](http://www.bioconductor.org/packages/release/bioc/html/sSeq.html)
  [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
  [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
  [stringr](https://cran.r-project.org/web/packages/stringr)
  

## Data

Simulation data and real data are contained in the data folder. You can see the detailed description in the html of this R package.


## Usage

ZINBScore(X_train_, y_train, X_test_, y_test, batch_size, py_path)

## Arguments
-X_train :Matrix or data.The input training sample.

-X_test :Matrix or data.The input test sample.

-y_train :A vector shows the class of the input sample

-y_test :A vector shows the class of the output sample

-batch_size 	:Number of training batches of the neural network.

-py_path :The path of python.exe.
## Value
The prediction of X_test and y_test.
