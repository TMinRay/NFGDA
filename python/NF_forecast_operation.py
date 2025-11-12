import numpy as np
import configparser
import glob
import os

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import tminlib.colorlevel as cl
from tminlib import plot_helper as phelp
from tminlib import math_kit as mk
from datetime import datetime, timedelta
from skimage.morphology import skeletonize, disk, binary_dilation, remove_small_objects, dilation
from skimage.measure import label
import scipy.io
from scipy.signal import medfilt2d
from scipy.spatial.distance import pdist, squareform
# from scipy.spatial import cKDTree
from scipy.ndimage import convolve
import datetime
from NFGDA_load_config import *
from pathlib import Path

VM=cl.VarMap()
varname_table=VM.varname_table
varunit_table=VM.varunit_table

kernel = np.ones((3,3), dtype=int)
ele_t_const, mean_t_const = 29.04, 28.43
r = np.sqrt(Cx**2+Cy**2)
rmask = r>=100

except_text = ""  # start empty log

def log_print(*args, **kwargs):
    """Behaves like print, but also saves the message."""
    global except_text
    s = " ".join(str(a) for a in args)
    except_text += s + "\n"
    print(*args, **kwargs) 

def rotation_matrix_2d(theta):
    """
    Returns a 2×2 rotation matrix for rotating points counterclockwise by angle theta (in radians).
    """
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

def find_roation_coord(points):
    points_T = points.T  # shape = (N, 2)
    dist_matrix = squareform(pdist(points_T))
    i, j = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)
    max_dist = dist_matrix[i, j]
    delta = points[:,i] - points[:,j]
    angle = np.arctan2(delta[1], delta[0])
    return angle, points[:,i]

def rotation_polyfit(points,n,fitn=None):
    # print(points.shape)
    if points.shape[1]>2:
        ang, origin = find_roation_coord(points)
        rot_points = np.matmul(rotation_matrix_2d(-ang),points-origin[:,np.newaxis])
        try:
            coeffs = np.polyfit(rot_points[0,:], rot_points[1,:], n)
        except:
            log_print(rot_points[0,:], rot_points[1,:])
            return points
        if fitn is None:
            fx = np.arange(np.min(rot_points[0,:]),np.max(rot_points[0,:])+0.25,0.25)
        else:
            fx = np.linspace(np.min(rot_points[0,:]),np.max(rot_points[0,:]),fitn)
        fy = np.polyval(coeffs, fx)
        return np.matmul(rotation_matrix_2d(ang),np.array([fx,fy]))+origin[:,np.newaxis]
    else:
        return points
        # if fitn is None:
        #     return points
        # else:
        #     if points.shape[1]==1:
        #         points = np.repeat(points, 2, axis=1)
        #     t = np.linspace(0, 1, fitn)
        #     return (1 - t) * points[:, [0]] + t * points[:, [1]]


def points_to_binary_grid(xy_points, gshape, origin, dg):
    grid = np.zeros(gshape, dtype=bool)

    # Round or floor the points to nearest integers
    xy_int = ((xy_points-origin)//dg).astype(int)

    # Clip to ensure they fall within bounds
    xy_int[0, :] = np.clip(xy_int[0, :], 0, gshape[1] - 1)
    xy_int[1, :] = np.clip(xy_int[1, :], 0, gshape[0] - 1)

    grid[xy_int[1, :], xy_int[0, :]] = True # grid [y, x]
    return grid

GFG_NONE    = 0b00000
GFG_LABEL   = 0b00001
GFG_NFGDA   = 0b00010
GFG_PREDICT = 0b00100
class GFGroups:
    ogn = np.array([Cx[0,0],Cy[0,0]])[:,np.newaxis]
    shape = Cx.shape
    def __init__(self, arc_anchors, timestamp=None, datakind = GFG_NONE):
        self.arc_anchors = arc_anchors
        self.timestamp = timestamp
        self.cur_motions = np.full(arc_anchors.shape, np.nan)
        self.pre_motions = np.full(arc_anchors.shape, np.nan)
        self.next_gp = np.full(arc_anchors.shape[0], np.nan)
        self.datakind = datakind
    def anchors_to_arcs(self):
        self.arc_points = []
        if len(self.arc_anchors)!=0:
            for nf_anchor in self.arc_anchors:
                self.arc_points.append(rotation_polyfit(nf_anchor,2))
        return self.arc_points
    def anchors_to_arcs_map(self):
        if len(self.arc_anchors)!=0:
            self.anchors_to_arcs()
            fit_points = np.concatenate(self.arc_points, axis=1)
            self.arcs_map = points_to_binary_grid(fit_points,self.shape,self.ogn,0.5)
            return self.arcs_map
        else:
            return np.zeros(self.shape,dtype=bool)

    def save(self, file_path: str | Path):
        """Save the GFGroups data to an .npz file"""
        np.savez(
            file_path,
            arc_anchors=self.arc_anchors,
            timestamp=self.timestamp,
            cur_motions=self.cur_motions,
            pre_motions=self.pre_motions,
            next_gp=self.next_gp,
            datakind = self.datakind,
        )

    @classmethod
    def load(cls, file_path: str | Path) -> "GFGroups":
        """Load a GFGroups object from an .npz file"""
        file_path = Path(file_path)
        data = np.load(file_path, allow_pickle=True)
        obj = cls.__new__(cls)  # create instance without calling __init__
        # assign directly
        obj.arc_anchors = data["arc_anchors"]
        obj.timestamp = data["timestamp"]
        obj.cur_motions = data["cur_motions"]
        obj.pre_motions = data["pre_motions"]
        obj.next_gp = data["next_gp"]
        obj.datakind = data["datakind"]
        return obj

class DataGFG(GFGroups):
    n_anchors = 10
    def __init__(self, data, binary_mask, kind = GFG_NONE):
        arc_anchors = []
        neighbors = convolve(binary_mask.astype(int), kernel, mode='constant') - binary_mask
        branch = (binary_mask & (neighbors >= 3))
        # Temporarily remove branch points
        skel_wo_branches = binary_mask.copy()
        skel_wo_branches[branch] = 0
        self.groups = label(skel_wo_branches, connectivity=2)
        # self.groups = label(binary_mask, connectivity=2)
        reduce = []
        for im in range(1,np.max(self.groups)+1):
            mask = self.groups == im
            gf_points = np.array([Cx[mask],Cy[mask]])
            if gf_points.shape[1]>2:
                rot_points = rotation_polyfit(gf_points,2,self.n_anchors)
                arc_anchors.append(rot_points)
            else:
                reduce.append(im)
        for im in reduce[::-1]:
            self.groups[self.groups==im] = 0
            self.groups[self.groups>=im] -= 1
        super().__init__(np.array(arc_anchors),data['timestamp'],kind)

class Prediction_Connection:
    def __init__(self, endidx, flip, igps, egps):
        self.endidx = endidx
        self.flip = flip
        self.igp_anchor = np.array(igps.arc_anchors)
        self.egp_anchor = np.array(egps.arc_anchors)
        self.motions = np.full(self.igp_anchor.shape,np.nan)
        self.speeds = np.full(self.igp_anchor.shape[0],np.nan)
        self.dt = (egps.timestamp - igps.timestamp)/np.timedelta64(60, 's')
        self.make_motion()
        self.copy_motion(igps, egps)
    def make_motion(self):
        for ii, ie in enumerate(self.endidx):
            ep = self.egp_anchor[ie]
            if self.flip[ii]:
                ep = np.fliplr(ep)
            self.motions[ii] = ep-self.igp_anchor[ii]
        self.motions = self.motions/self.dt 
        if self.motions.ndim>1:
            self.speeds = np.sqrt(np.sum(self.motions**2,axis=1))
    def copy_motion(self, igps, egps):
        igps.next_gp = self.endidx
        igps.cur_motions = self.motions
        for ii, ie in enumerate(set(self.endidx)):
            egps.pre_motions[ie] = np.mean(self.motions[self.endidx==ie].reshape(-1,2,egps.n_anchors),axis=0)
            # if np.sum(self.endidx==ie)>1:
            #     print(egps.pre_motions[ie],self.motions[self.endidx==ie].shape,self.flip[self.endidx==ie])

class Prediction_Worker:
    def __init__(self, gps):
        self.gps = gps
        self.connects = {}
    def update_velocitys(self, k=None):
        if k is None:
            for ig in range(len(self.gps)-1):
                A = self.gps[ig]
                B = self.gps[ig+1]
                # endidx, flip, dist = self.GFG_motion(A.arc_anchors,B.arc_anchors)
                # self.connects[ig] = Prediction_Connection(endidx, flip, A, B, ig)
                self.GFG_motion(A,B,ig)
        else:
            A = self.gps[k]
            B = self.gps[k+1]
            # endidx, flip, dist = self.GFG_motion(A.arc_anchors,B.arc_anchors)
            # self.connects[k] = Prediction_Connection(endidx, flip, A, B)
            self.GFG_motion(A,B,k)
    def GFG_motion(self,sgfg,egfg,conn):
        A = sgfg.arc_anchors
        B = egfg.arc_anchors
        if (A.ndim != 3) or (B.ndim != 3):
            # print(A,B)
            self.connects[conn] = Prediction_Connection([], [], sgfg, egfg)
        else:
            A_exp = A[:, None, :, :]
            B_exp = B[None, :, :, :]

            # Compute normal distances (mean over last axis=2, then mean over coordinates)
            dist_norm = np.linalg.norm(A_exp - B_exp, axis=2)
            
            # Compute flipped distances
            B_flip = B[:, :, ::-1]  # flip along the sample axis (10)
            B_flip_exp = B_flip[None, :, :, :]
            dist_flip = np.linalg.norm(A_exp - B_flip_exp, axis=2)
            
            # Choose the smaller distance
            dist_final = np.max(dist_norm, axis=-1).copy()
            flip_arc = np.sum(dist_norm, axis=-1) > np.sum(dist_flip, axis=-1)
            dist_final[flip_arc]=np.max(dist_flip, axis=-1)[flip_arc]
            
            endpoint = np.argmin(dist_final,1)
            self.connects[conn] = Prediction_Connection(endpoint, flip_arc[np.arange(endpoint.size),endpoint], sgfg, egfg)
            # return endpoint, flip_arc[np.arange(endpoint.size),endpoint], dist_final[np.arange(endpoint.size),endpoint]
    def prediction(self,startf,dt,mode='element'):
        if mode == 'mean':
            anchors = self.connects[startf].igp_anchor+np.mean(self.connects[startf].motions,axis=2)[:,:,np.newaxis]*dt
            tstp = self.gps[startf].timestamp + dt*np.timedelta64(1, 'm')
        elif mode == 'element':
            anchors = self.connects[startf].igp_anchor+self.connects[startf].motions*dt
            tstp = self.gps[startf].timestamp + dt*np.timedelta64(1, 'm')
        return GFGroups(anchors,tstp)

def exp_weight(dt,t_const):
    return np.exp(-dt/t_const)

def nfgda_forecast(case_name):
    exp_preds_event = export_preds_dir + case_name
    savedir = os.path.join(export_forecast_dir[:-1]+'-operation/', case_name)
    os.makedirs(savedir,exist_ok=True)

    npz_list = glob.glob(exp_preds_event + "/*npz")

    evs = []
    gps = []
    data = []
    for ifn in npz_list:
        data.append(np.load(ifn))
        evs.append(DataGFG(data[-1],skeletonize(data[-1]['evalbox'])))
        # gps.append(DataGFG(data[-1],data[-1]['nfout']))
    predict_worker = Prediction_Worker(gps)
    eval_worker = Prediction_Worker(evs)
    worker = eval_worker
    for iframe in range(len(npz_list)-1):
        worker.update_velocitys(iframe)
    max_predict = 20
    for pframe in range(len(npz_list)-1):
        ppi_id = os.path.basename(npz_list[pframe])
        pgf = np.zeros(Cx.shape)
        ps = 0
        # fig, axs = plt.subplots(1, 2, figsize=(7/0.7, 2.5/0.7),dpi=250, gridspec_kw=dict(left=0.08, right=1-0.085, top=1-0.08, bottom=0.06, wspace=0.25, hspace=0.16))
        pdata = np.ma.masked_where(rmask,data[pframe]['inputNF'][:,:,1])
        # axs[0].pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
        # axs[0].contour(Cx,Cy,data[pframe]['evalbox'],colors='k')
        tvec=[]
        ele_map=[]
        mean_map=[]
        for iframe in range(pframe-max_predict,pframe):
            if iframe<0:
                continue
            print(f'[{iframe}] -> [{pframe}]')
            dt = (worker.gps[pframe].timestamp-worker.gps[iframe].timestamp)/np.timedelta64(60, 's')
            if worker.connects[iframe].igp_anchor.ndim>1:
                tvec.append(worker.gps[iframe].timestamp)
                # ele_w = exp_weight(dt,ele_t_const)
                # mean_w = exp_weight(dt,mean_t_const)
                # # ele_w = 1
                # # mean_w = 1
                # ps += ele_w + mean_w
                end = worker.prediction(iframe, dt)
                ele_map.append(end.anchors_to_arcs_map())

                # axs[0].plot(end.arc_anchors[:,0,:].T,end.arc_anchors[:,1,:].T,'.-',color=((0.5+0.5*(pframe-iframe)/max_predict),0,0),label='element',alpha=0.5)
                # pgf += ele_w*binary_dilation(end.anchors_to_arcs_map(), footprint=disk(3)).astype(float)

                end = worker.prediction(iframe, dt,mode='mean')
                mean_map.append(end.anchors_to_arcs_map())

        #         axs[0].plot(end.arc_anchors[:,0,:].T,end.arc_anchors[:,1,:].T,'.-',color=(0,0,(0.5+0.5*(pframe-iframe)/max_predict)),label='mean',alpha=0.5)
        #         pgf += mean_w*binary_dilation(end.anchors_to_arcs_map(), footprint=disk(3)).astype(float)
        # pgf=pgf/ps*1e2

        # pcm=axs[1].pcolormesh(Cx,Cy,pgf,vmin=0,vmax=80,cmap='jet')
        # plt.colorbar(pcm)
        # axs[1].contour(Cx,Cy,data[pframe]['evalbox'],colors='k')
        # handles, labels = axs[0].get_legend_handles_labels()
        # by_label = dict(zip(labels, handles))
        # axs[0].legend(by_label.values(), by_label.keys())
        # # axs[1].set_title(gps[iframe].timestamp)
        # valid_time = worker.gps[pframe].timestamp
        # fig.suptitle(valid_time.astype(datetime.datetime).item().strftime('%y/%m/%d %H:%M:%S'))
        # for ia in range(2):
        #     axs[ia].set_xlim(-100,100)
        #     axs[ia].set_ylim(-100,100)

        # fig.savefig(os.path.join(savedir, ppi_id[:-4]+'.png'))
        # plt.close(fig)
        data_dict = {"ele_map": ele_map, "mean_map": mean_map, "start_timestamps":tvec}
        np.savez(os.path.join(savedir, ppi_id[:-4]+'.npz'), **data_dict)
    print("\n=== Full Log ===")
    print(except_text)

if __name__ == '__main__':
    nfgda_forecast(config["Settings"]["case_name"])