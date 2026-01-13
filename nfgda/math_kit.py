"""
tminlib.math_kit
============

Module of math toolbox.

Last updated: 2025-11-07
Author: Min-Duan Tzeng
Radar Innovation Laboratory (RIL),
School of Electrical and Computer Engineering,
The University of Oklahoma,
Norman, Oklahoma
E-mail: tzengmin2@gmail.com ; tmin@ou.edu ; b02209039@ntu.edu.tw

Functions:
    fftconvolve()
    circularconvolve()
    angdiff()
    geopoints_to_relative_xy()
    window_blend()
"""
import numpy as np

def fftconvolve(ax,ay,axis=-1):
    nex=ax.shape[axis]+ay.shape[axis]-1
    return np.fft.ifft(np.fft.fft(ax,n=nex,axis=axis)
        *np.fft.fft(ay,n=nex,axis=axis),axis=axis)

def circularconvolve(ax,ay,axis=-1):
    return np.fft.ifft(np.fft.fft(ax,axis=axis)
        *np.fft.fft(ay,axis=axis),axis=axis)

def angdiff(ar,ab):
    return np.angle(np.exp(1j*ar/180*np.pi)/np.exp(1j*ab/180*np.pi))*180/np.pi

def geopoints_to_relative_xy(geo_1,geo_2):
    lat1, lon1 = geo_1
    lat2, lon2 = geo_2
    lat1=lat1*np.pi/180
    lat2=lat2*np.pi/180
    lon1=lon1*np.pi/180
    lon2=lon2*np.pi/180
    dlat = lat2-lat1
    dlon = lon2-lon1
    Re = 6371.0 * 1000.0
    a = np.sin(dlat/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    s = 2*Re*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    theta=np.pi/2-np.arctan2(np.sin(dlon)*np.cos(lat2),
        np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(dlon))
    x = s*np.cos(theta)
    y = s*np.sin(theta)
    return  x, y

def window_blend(wd,blend):
    return  (1.0 - blend) * wd + blend