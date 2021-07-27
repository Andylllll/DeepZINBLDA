import pandas as pd
from time import time
import tensorflow as tf
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
from ZINBAutoEncoder import SCDeepCluster
from sklearn.neighbors import KNeighborsClassifier
from loss import poisson_loss, NB, ZINB
from preprocess import read_dataset, normalize
import h5py
from normal import normalize_tr
from normal import normalize_te
import scanpy as sc
from sklearn.model_selection import StratifiedShuffleSplit
from keras.utils.np_utils import to_categorical
from numpy.random import seed
seed(2211)
from tensorflow import set_random_seed
set_random_seed(2211)
import sklearn.svm as svm
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
import warnings
warnings.filterwarnings('ignore')



def DeepZINBL(X_train_,y_train,X_test_,y_test,batch_size):
        
        pretrain_epochs=800
        optimizer1 = Adam(amsgrad=True)
        sf,train_raw,X_train,var_names,y_train=normalize_tr(X_train_,y_train)
        sf_test,test_raw,X_test,y_test=normalize_te(X_test_,var_names,y_test)
        batch_size=batch_size
        update_interval = int(X_train.shape[0]/batch_size)
        input_size=X_train.shape[1]
        gamma=1
        scDeepCluster = SCDeepCluster(dims=[input_size,64,32,16])
        ae_weight_file='ae_weights.h5'
        ae_weights='ae_weights.h5'
        save_dir='results/scDeepCluster'
        scDeepCluster.pretrain(x=[X_train, sf], y=train_raw, batch_size=batch_size, epochs=pretrain_epochs, optimizer=optimizer1, ae_file=ae_weight_file)
        pi,disp,mean=scDeepCluster.extract_parameter(x=[X_train, sf])
        pi_t,disp_t,mean_t=scDeepCluster.extract_parameter(x=[X_test, sf_test])
        parameter={'pi': pi,  'disp' :disp,    'pi_t': pi_t,   'disp_t':disp_t,    'X_train':train_raw,   'X_test':test_raw, 'y_train':y_train,   'y_test':y_test}
        K.clear_session()
        return   parameter
    





