import numpy as np
import configparser
import glob
import os
import scipy.io
from skimage.morphology import skeletonize

config = configparser.ConfigParser()
config.read("./NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
train_dir = '../mat/train/'
def train_prep(case_name):
    exp_preds_event = export_preds_dir + case_name
    npz_list = glob.glob(exp_preds_event + "/*npz")
    true_buf = []
    false_buf = []
    for ppi_file in npz_list:
        print(ppi_file)
        data = np.load(ppi_file)
        evalbox = data['evalbox']
        evalline = skeletonize(evalbox)
        # evalline = evalbox.astype(bool)
        true_set = data['inputNF'][evalline]
        true_set=true_set[np.min(np.isfinite(true_set),axis=1)]
        false_set = data['inputNF'][np.logical_not(data['evalbox'])].reshape(-1, data['inputNF'].shape[2])
        false_set=false_set[np.min(np.isfinite(false_set),axis=1)]
        if true_set.size>0:
            true_buf.append(true_set)
            num_pick = int(np.ceil(true_set.shape[0]*16/6))
            idx = np.random.choice(false_set.shape[0], num_pick, replace=False)
            false_buf.append(np.take(false_set, idx, axis=0))
            print(true_buf[-1].shape,false_buf[-1].shape)
    true_set = np.concatenate(true_buf,axis=0)
    false_set = np.concatenate(false_buf,axis=0)
    trainNF = np.concatenate((np.concatenate((true_set,np.ones((true_set.shape[0],1))),axis=1)
        ,np.concatenate((false_set,np.zeros((false_set.shape[0],1))),axis=1))
        ,axis=0)
    scipy.io.savemat(train_dir+'train_set_'+case_name+'.mat', {"trainNF":trainNF})
    # scipy.io.savemat(train_dir+'train_set_wide_'+case_name+'.mat', {"trainNF":trainNF})


if __name__ == '__main__':
    train_prep(config["Settings"]["case_name"])