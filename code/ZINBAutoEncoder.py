import pandas as pd
from time import time
import tensorflow as tf
import os
import numpy as np
from keras.models import Model
import keras.backend as K
from keras.engine.base_layer import InputSpec
from keras.engine.topology import Layer
from keras.layers import Dense, Input, GaussianNoise, Layer, Activation
from keras.models import Model
from keras.optimizers import SGD, Adam
from keras.utils.vis_utils import plot_model
from keras.callbacks import EarlyStopping
from layers import ConstantDispersionLayer, SliceLayer,ColWiseMultLayer
from sklearn.model_selection import train_test_split
from loss import poisson_loss, NB, ZINB
from preprocess import read_dataset, normalize
import h5py
import scanpy as sc
from sklearn.model_selection import StratifiedShuffleSplit
from keras.utils.np_utils import to_categorical
from numpy.random import seed
seed(2211)
from tensorflow import set_random_seed
set_random_seed(2211)

MeanAct = lambda x: tf.clip_by_value(K.exp(x), 1e-5, 1e6)
DispAct = lambda x: tf.clip_by_value(tf.nn.softplus(x), 1e-4, 1e4)

def autoencoder(dims, init='glorot_uniform', act='relu'):
    n_stacks = len(dims) - 1

    sf_layer = Input(shape=(1,), name='size_factors')
    x = Input(shape=(dims[0],), name='counts')
    h = x
  
    for i in range(n_stacks-1):
        h = Dense(dims[i + 1], kernel_initializer=init, name='encoder_%d' % i)(h)  
        h = Activation(act)(h)

    h = Dense(dims[-1], kernel_initializer=init, name='encoder_hidden')(h)  # dims[-1]=32


    for i in range(n_stacks-1, 0, -1):
        h = Dense(dims[i], activation=act, kernel_initializer=init, name='decoder_%d' % i)(h) #对称生成decoder层


    pi = Dense(dims[0], activation='sigmoid', kernel_initializer=init, name='pi')(h)
    disp = Dense(dims[0], activation=DispAct, kernel_initializer=init, name='dispersion')(h)
    mean = Dense(dims[0], activation=MeanAct, kernel_initializer=init, name='mean')(h)

    output = ColWiseMultLayer(name='output')([mean, sf_layer])#[batch_size,基因数]
    output = SliceLayer(0, name='slice')([mean, disp, pi])

    return Model(inputs=[x, sf_layer], outputs=output)


class SCDeepCluster(object):
    def __init__(self,
                 dims,
                 noise_sd=0,
                 alpha=1.0,
                 ridge=0,
                 debug=False):

        super(SCDeepCluster, self).__init__()

        self.dims = dims 
        self.input_dim = dims[0] 
        self.n_stacks = len(self.dims) - 1 

        self.alpha = alpha
        self.act = 'relu'
        self.ridge = ridge
        self.debug = debug
        self.autoencoder = autoencoder(self.dims, act = self.act)
        
        ae_layers = [l for l in self.autoencoder.layers]
        hidden = self.autoencoder.input[0]
        for i in range(1, len(ae_layers)):
            if "noise" in ae_layers[i].name:
                next
            elif "dropout" in ae_layers[i].name:
                next
            else:
                hidden = ae_layers[i](hidden)
            if "encoder_hidden" in ae_layers[i].name:  # only get encoder layers
                break
        self.encoder = Model(inputs=self.autoencoder.input, outputs=hidden)

        pi = self.autoencoder.get_layer(name='pi').output
        disp = self.autoencoder.get_layer(name='dispersion').output
        mean = self.autoencoder.get_layer(name='mean').output
        zinb = ZINB(pi, theta=disp, ridge_lambda=self.ridge, debug=self.debug) 
        self.loss = zinb.loss
        self.parameter = Model(inputs=[self.autoencoder.input[0], self.autoencoder.input[1]],
                           outputs=[pi,disp,mean])
        self.model = Model(inputs=[self.autoencoder.input[0], self.autoencoder.input[1]],
                           outputs=self.autoencoder.output)

        self.pretrained = False
    
    def pretrain(self, x, y, batch_size=10, epochs=20, optimizer='adam', ae_file=None):
        print('...Pretraining autoencoder...')
        self.autoencoder.compile(loss=self.loss, optimizer=optimizer)
        es = EarlyStopping(monitor="loss", patience=50, verbose=0)
        self.autoencoder.fit(x=x, y=y, batch_size=batch_size, epochs=epochs, callbacks=[es],verbose=0)

        self.pretrained = True
        
    def load_weights(self, weights_path):  # load weights of scDeepCluster model
        self.model.load_weights(weights_path)

    def extract_feature(self, x):  # extract features from before clustering layer
        return self.encoder.predict(x)
    
    def extract_parameter(self, x): 
        return self.parameter.predict(x)

