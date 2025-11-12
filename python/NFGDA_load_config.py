import configparser
from tminlib import math_kit as mk
import numpy as np

config = configparser.ConfigParser()
config.read("NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
export_forecast_dir = config["Settings"]["export_forecast_dir"]
fig_dir = config["Settings"]["fig_dir"]
label_on = config.getboolean('labels', 'label_on')
if label_on:
    label_loc = list(map(float,config.get("labels", "loc").split(",")))
    radar_loc = list(map(float,config.get("labels", "rloc").split(",")))
    sitex, sitey = mk.geopoints_to_relative_xy(radar_loc,label_loc)

Cx, Cy = np.meshgrid(np.arange(-100,100.5,0.5),np.arange(-100,100.5,0.5))