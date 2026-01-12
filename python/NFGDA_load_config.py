import configparser
import math_kit as mk
import numpy as np
import os
from scipy.interpolate import LinearNDInterpolator

config = configparser.ConfigParser()
config.read("NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
os.makedirs(export_preds_dir,exist_ok=True)
evalbox_on = config.getboolean('Settings', 'evalbox_on')
export_forecast_dir = config["Settings"]["export_forecast_dir"]
V06_dir = config["Settings"]["V06_dir"]
radar_id = config["Settings"]["radar_id"]
fig_dir = config["Settings"]["fig_dir"]
label_on = config.getboolean('labels', 'label_on')
if label_on:
    label_loc = list(map(float,config.get("labels", "loc").split(",")))
    radar_loc = list(map(float,config.get("labels", "rloc").split(",")))
    sitex, sitey = mk.geopoints_to_relative_xy(radar_loc,label_loc)

Cx, Cy = np.meshgrid(np.arange(-100,100.5,0.5),np.arange(-100,100.5,0.5))
r = np.sqrt(Cx**2+Cy**2)
rmask = r>=100

# radar_id = 'KABX'
# V06_dir = '../V06/runtime/'+radar_id
os.makedirs(V06_dir,exist_ok=True)
nf_dir = V06_dir+'npz/'
os.makedirs(nf_dir,exist_ok=True)

PARROT_mask_on = False

thrREF = -5
thrdREF = 0.3
RegR = np.arange(0,400)/4
RegAZ = np.arange(0,360,0.5)*np.pi/180
RegPolarX = RegR[:,np.newaxis] * np.sin(RegAZ[np.newaxis,:])
RegPolarY = RegR[:,np.newaxis] * np.cos(RegAZ[np.newaxis,:])
interpolator = LinearNDInterpolator((RegPolarX.reshape(-1),RegPolarY.reshape(-1)), np.zeros(RegPolarX.shape).reshape(-1))

###### Beta Cell magic numbers ##########
cellthresh = 5
cbcellthrsh = 0.8
cellcsrthresh=0.5
crsize = 5
cellINT = crsize + 2
widecellINT =crsize+4
avgINT = 8
s2xnum = [10, 15]
s2ynum = [-3, 1]
s2xdel = s2xnum[1]-s2xnum[0]
s2ydel = s2ynum[1]-s2ynum[0]
s2g = s2ydel/s2xdel
s2gc = s2ynum[1]-s2g*s2xnum[1]
Celldp = np.load("Celldp.npy")
Celldpw = np.load("Celldpw.npy")
###### Beta Cell magic numbers ##########

######### FTC Beta Z, dZ displacements ###########
datacy = np.arange(-8,9).reshape(1,-1)
datacx = np.zeros((1,17))
datac = np.swapaxes(np.array([datacy,datacx]),0,2)

datasy = np.array([*np.arange(-7,0,2),0,*np.arange(1,8,2),*np.arange(-7,0,2),0,*np.arange(1,8,2)]).reshape(1,-1)
datasx = np.array([-4*np.ones((9)),4*np.ones((9))]).reshape(1,-1)
datas = np.swapaxes(np.array([datasy,datasx]),0,2)

class path_struct():
    def __init__(self):
        self.nf_dir = nf_dir
        self.V06_dir = V06_dir
        self.nf_preds_dir = export_preds_dir
        self.nf_forecast_dir = export_forecast_dir

path_config = path_struct()