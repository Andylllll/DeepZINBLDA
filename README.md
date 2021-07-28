# `DeepZINBLDA`

> a deep learning framework with statistical method to classify RNA-seq data.

## Description

Genetic cell classification is extremely important no matter in medical or biological fields. We design an antoencoder network to learn the parameters of the model, and use the zero-inflated negative binomial distribution as the loss function for better feature selection. Then we put the trained parameters into the zero-inflated negative binomial linear discriminant model, which based on the Bayesian model. 


## Requirements

It is required to install the following dependencies in order to be able to run the code of scDLC

- [Anaconda3](https://www.anaconda.com/products/individual)  
- [R>=3.6.0](https://cran.r-project.org/)  
- [python 3](https://www.python.org/downloads/)  
  [sklearn](https://pypi.org/project/sklearn/0.0/)，[numpy 1.19.1](https://pypi.org/project/numpy/1.19.1/)，[tensorflow-gpu 1.15.0](https://pypi.org/project/tensorflow-gpu/1.15.0/),[keras 2.3.1](https://pypi.org/project/keras/2.3.1/)，[scanpy 1.3.7](https://pypi.org/project/scanpy/1.3.7/)
  
  

## Data

Simulation data and real data are contained in the data folder. You can see the detailed description in the html of this R package.


## Usage

DeepZINB(X_train, X_test, y_train, y_test, batch_size, path_py)

## Arguments
-X_train :Matrix or data.The input training sample.
-X_test :Matrix or data.The input test sample.
-y_train :A vector shows the class of the input sample
-y_test :A vector shows the class of the output sample
-batch_size 	:Number of training batches of the neural network.
-path_py :The path of python.exe.


