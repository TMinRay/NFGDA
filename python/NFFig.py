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
from skimage.morphology import skeletonize, disk, binary_dilation, binary_erosion, remove_small_objects, dilation
from skimage.measure import label
from scipy.signal import medfilt2d
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import cKDTree
import scipy.io

VM=cl.VarMap()
varname_table=VM.varname_table
varunit_table=VM.varunit_table
config = configparser.ConfigParser()
config.read("./NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
fig_dir = config["Settings"]["fig_dir"]
label_on = config.getboolean('labels', 'label_on')

if label_on:
    label_loc = list(map(float,config.get("labels", "loc").split(",")))
    radar_loc = list(map(float,config.get("labels", "rloc").split(",")))
    sitex, sitey = mk.geopoints_to_relative_xy(radar_loc,label_loc)
# export_statistics_dir = config["Settings"]["export_statistics_dir"]
# os.makedirs(export_statistics_dir,exist_ok=True)
export_preds_fapos_dir = export_preds_dir[:-1]+'_pos/'
os.makedirs(export_preds_fapos_dir,exist_ok=True)

def eval_nf(nfloc,evalline,evalbox,gd):
    nfpredict = dilation(nfloc, disk(5))
    Mhits = np.logical_and(evalline,nfpredict)
    Mmiss = np.logical_and(evalline,~nfpredict)
    q_arc = nf_arc(gd['Cx'],gd['Cy'],nfloc)
    Mfa = np.logical_and(np.logical_not(evalbox),q_arc)
    sch,scm,scf,scp = np.sum(Mhits), np.sum(Mmiss), np.sum(Mfa), np.sum(q_arc)
    HR = 1e2*sch/(scm+sch)
    FAR = 1e2*scf/(scp)
    return (HR, FAR), (Mhits, Mmiss, Mfa, q_arc)

def log_stat(fn, stats_list):
    grouped = list(zip(*stats_list))
    # Stack each group
    stacked = [np.stack(arr_list) for arr_list in grouped]
    # Create a dictionary with key names
    save_dict = {
    'hits':stacked[0], 
    'miss':stacked[1], 
    'fa':stacked[2], 
    'q_arc':stacked[3]
    }
    # Save to .npz
    np.savez(fn, **save_dict)
    return save_dict

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

def rotation_polyfit(points,n):
    # print(points.shape)
    if points.shape[1]>2:
        ang, origin = find_roation_coord(points)
        rot_points = np.matmul(rotation_matrix_2d(-ang),points-origin[:,np.newaxis])
        coeffs = np.polyfit(rot_points[0,:], rot_points[1,:], n)
        fx = np.arange(np.min(rot_points[0,:]),np.max(rot_points[0,:])+0.25,0.25)
        fy = np.polyval(coeffs, fx)
        return np.matmul(rotation_matrix_2d(ang),np.array([fx,fy]))+origin[:,np.newaxis]
    else:
        return points

def points_to_binary_grid(xy_points, gshape, origin, dg):
    grid = np.zeros(gshape, dtype=bool)

    # Round or floor the points to nearest integers
    xy_int = ((xy_points-origin)//dg).astype(int)

    # Clip to ensure they fall within bounds
    xy_int[0, :] = np.clip(xy_int[0, :], 0, gshape[1] - 1)
    xy_int[1, :] = np.clip(xy_int[1, :], 0, gshape[0] - 1)

    grid[xy_int[1, :], xy_int[0, :]] = True # grid [y, x]
    return grid

def nf_arc(xx,yy,bw):
    groups = label(bw, connectivity=2)
    fit_chunks = []
    for im in range(1,np.max(groups)+1):
        mask = groups == im
        gf_points = np.array([xx[mask],yy[mask]])
        rot_points = rotation_polyfit(gf_points,2)
        fit_chunks.append(rot_points)
    ogn = np.array([xx[0,0],yy[0,0]])[:,np.newaxis]
    if len(fit_chunks)!=0:
        fit_points = np.concatenate(fit_chunks, axis=1)
        return points_to_binary_grid(fit_points,xx.shape,ogn,0.5)
    else:
        return np.zeros(xx.shape,dtype=bool)

def clean_indices(idx,shp,edg):
    dim0 = idx[:,0]
    dim1 = idx[:,1]
    inbox = (dim0>=edg) & (dim0< shp[0]-edg) & (dim1>=edg) & (dim1< shp[1]-edg)
    return idx[inbox,:]

def post_proc(inGST):
    hGST = medfilt2d(inGST.astype(float), kernel_size=3)
    binary_mask = post_moving_avg(hGST) >= 0.6  # Thresholding
    pskel_nfout = binary_dilation(binary_mask, disk(5))
    skel_nfout = skeletonize(pskel_nfout*inGST)
    skel_nfout2 = remove_small_objects(skel_nfout, min_size=10, connectivity=2)
    return skel_nfout2

# gfv = [4, 32]
class GFSpace:
    def __init__(self, gfv = [18, 72]):
        self.gfv = gfv
        self.data =[]
        self.nf_pair = []
        self.nf_loc = []
        self.tstamp = []
    
    def load_nf(self,fn):
        self.data.append(np.load(fn))
        self.tstamp.append(get_tstamp(fn))
        self.nf_loc.append( {'x':self.data[-1]['xi2'][self.data[-1]['nfout']],'y':self.data[-1]['yi2'][self.data[-1]['nfout']],
                             'idx':np.arange(self.data[-1]['xi2'].size).reshape(self.data[-1]['xi2'].shape)[self.data[-1]['nfout']]})
        if len(self.nf_loc)>1:
            self.nf_pair.append(self.connect_nf(self.nf_loc[-2],self.nf_loc[-1],(self.tstamp[-1]-self.tstamp[-2]).total_seconds()))
        else:
            self.shp = self.data[0]['xi2'].shape
            
    def connect_nf(self,t1,t2,dt):
        A = np.column_stack((t1['x'].reshape(-1), t1['y'].reshape(-1)))
        B = np.column_stack((t2['x'].reshape(-1), t2['y'].reshape(-1)))
        # coord [dim_sample, dim_xy]
        # coord [dim_sample, 2]
        tree = cKDTree(A)
        dists_for, indices_for = tree.query(B)
        tree = cKDTree(B)
        dists_bac, indices_bac = tree.query(A)
        pair_pool = np.vstack((np.column_stack((indices_for,np.arange(B.shape[0]))),np.column_stack((np.arange(A.shape[0]),indices_bac)))).astype(int)
        dists_pool = np.concatenate((dists_for,dists_bac),axis=0)
        mask = np.logical_and(dists_pool > self.gfv[0]*dt/3600, dists_pool < self.gfv[1]*dt/3600)
        return np.concatenate((pair_pool[mask,:],np.zeros(np.sum(mask),dtype=bool).reshape(-1,1)),axis=1)
        
        # return pair_pool[mask,:].astype(int)
    def clean_short_track(self):
        for ic in range(len(self.nf_pair)-1):
            keep_mask = np.logical_or(np.isin(self.nf_pair[ic][:,1], np.unique(self.nf_pair[ic+1][:,0])),self.nf_pair[ic][:,2]>0)
            self.nf_pair[ic] = self.nf_pair[ic][keep_mask,:]
            self.nf_pair[ic+1][:,2] = np.logical_or(self.nf_pair[ic+1][:,2], np.isin(self.nf_pair[ic+1][:, 0], self.nf_pair[ic][:, 1]))
    
    def clean_random_track_motion(self):
        self.cal_motion()
        for ic in range(len(self.nf_pair)-1):
            curdir = self.nf_pair[ic][:,-2]+1j*self.nf_pair[ic][:,-1]
            nextdir = np.zeros(curdir.size)
            for ip in range(curdir.size):
                mask = self.nf_pair[ic+1][:,0]==self.nf_pair[ic][ip,1]
                nextdir[ip] = np.mean(self.nf_pair[ic+1][mask,-2]+1j*self.nf_pair[ic+1][mask,-1],axis=0)
            dirdiff = get_dirdiff(curdir,nextdir)
            self.nf_pair[ic]=np.concatenate((self.nf_pair[ic],dirdiff.reshape(-1,1)),axis=1)
        
    def cal_motion(self):
        for ic in range(len(self.nf_pair)):
            t1 = self.nf_loc[ic]
            t2 = self.nf_loc[ic+1]
            A = np.column_stack((t1['x'].reshape(-1), t1['y'].reshape(-1)))
            B = np.column_stack((t2['x'].reshape(-1), t2['y'].reshape(-1)))
            start = A[self.nf_pair[ic][:,0]]
            end = B[self.nf_pair[ic][:,1]]
            motion = end-start
            self.nf_pair[ic]=np.concatenate((self.nf_pair[ic],motion),axis=1)
    
    def get_cln_nf(self,ic):
        buf = np.zeros(self.shp,dtype=bool).reshape(-1)
        buf[self.nf_loc[ic]['idx'][self.nf_pair[ic][:,0].astype(int)]]=True
        return buf.reshape(self.shp)
        # return remove_small_objects(buf.reshape(self.shp), min_size=5, connectivity=2)

def get_tstamp(ppi_file):
    ppi_id = os.path.basename(ppi_file)
    ppi_name = ppi_id[11:]
    date_part = ppi_name[4:12]   # 5:12 in MATLAB → 4:12 in Python
    time_part = ppi_name[13:19]  # 14:19 in MATLAB → 13:19 in Python
    tstamp_date = datetime.strptime(date_part, "%Y%m%d")
    tstamp_time = datetime.strptime(time_part, "%H%M%S").time()
    tstamp = datetime.combine(tstamp_date.date(), tstamp_time)
    return tstamp

def get_dirdiff(dir1,dir2):
    return np.rad2deg(np.angle(dir1*np.conj(dir2)))

def nffig_proc(case_name,plot_on=False):
    exp_preds_event = export_preds_dir + case_name
    savedir = os.path.join(fig_dir, case_name)
    os.makedirs(savedir,exist_ok=True)

    export_preds_fapos_event = export_preds_fapos_dir + case_name
    os.makedirs(export_preds_fapos_event,exist_ok=True)

    npz_list = glob.glob(exp_preds_event + "/*npz")
    wgfspace = GFSpace([18,72])
    for ppi_file in npz_list:
        wgfspace.load_nf(ppi_file)
    # wgfspace.clean_short_track()
    # wgfspace.clean_random_track_motion()
    Cx = wgfspace.data[0]['xi2']
    Cy = wgfspace.data[0]['yi2']
    r = np.sqrt(Cx**2+Cy**2)
    rmask = r>=100
    
    tvec = wgfspace.tstamp[:-1]
    hr_pre = []
    hr_pos = []
    fa_pre = []
    fa_pos = []
    stat_pre = []
    stat_pos = []
    for ic,data in enumerate(wgfspace.data[:-1]):
        if evalbox_on:
            evalbox = data['evalbox']
        else:
            evalbox = np.zeros(Cx.shape)
        evalline = skeletonize(evalbox)
        gcoord = {'Cx':data['xi2'],'Cy':data['yi2']}
        (hr,fa), eval_pre = eval_nf(data['nfout'],evalline,evalbox,gcoord)
        hr_pre.append(hr)
        fa_pre.append(fa)
        stat_pre.append(eval_pre)

        proc_nf = wgfspace.get_cln_nf(ic)
        matout = export_preds_fapos_event + '/' + npz_list[ic].split('/')[-1]
        data_dict = {"xi2":data['xi2'],"yi2":data['yi2'],"REF":data['REF'], \
                    "nfout": proc_nf,"inputNF":data['inputNF'],
                    "evalbox":evalbox}
                    # ,'outputGST':data['outputGST']}
        scipy.io.savemat(matout, data_dict)

        (hr,fa), eval_pos = eval_nf(proc_nf,evalline,evalbox,gcoord)
        hr_pos.append(hr)
        fa_pos.append(fa)
        stat_pos.append(eval_pos)

        ppi_file = npz_list[ic]
        print(ppi_file)
        ppi_id = os.path.basename(ppi_file)
        ppi_name = ppi_id[11:]  # MATLAB 12:end is Python 11: (0-based)
        date_part = ppi_name[4:12]   # 5:12 in MATLAB → 4:12 in Python
        time_part = ppi_name[13:19]
    
        radar_id = ppi_name[0:4]  # 1:4 in MATLAB → 0:4 in Python
        tstamp_date = datetime.strptime(date_part, "%Y%m%d")
        tstamp_time = datetime.strptime(time_part, "%H%M%S").time()
        tstamp = datetime.combine(tstamp_date.date(), tstamp_time)
        ppi_desc = f"{radar_id}, {tstamp.strftime('%m/%d/%Y, %H:%M:%S %Z')}"
    
    
        fig, axs = plt.subplots(1, 2, figsize=(7/0.7, 2.5/0.7),dpi=250, gridspec_kw=dict(left=0.08, right=1-0.085, top=1-0.08, bottom=0.06, wspace=0.25, hspace=0.16))
        REF = data['inputNF'][:,:,1]
        pdata = np.ma.masked_where(rmask,REF)

        def plot_nf(ax,nfout,hrt,fat,evalmasks):
            # evalbox = data['evalbox']
            nfloc = np.logical_and(~rmask,nfout)
            # nfpredict = dilation(nfloc, disk(5))
            # Mhits = np.logical_and(evalline,nfpredict)
            # Mmiss = np.logical_and(evalline,~nfpredict)
            ax.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)

            ax.plot(Cx[nfloc],Cy[nfloc],'k.',markersize=0.8)
            ax.plot(Cx[np.logical_and(nfloc,evalbox)],Cy[np.logical_and(nfloc,evalbox)],'r.',markersize=0.8)
            ax.plot(Cx[evalmasks[3]],Cy[evalmasks[3]],'.',color=(1,0.5,0),markersize=0.8)
            ax.plot(Cx[evalmasks[2]],Cy[evalmasks[2]],'.',color=(1,0,1),markersize=0.8)
            ax.plot(Cx[evalmasks[1]],Cy[evalmasks[1]],'b.',markersize=0.8)
            ax.plot(Cx[evalmasks[0]],Cy[evalmasks[0]],'.',color=(0,1,0),markersize=0.8)
            ax.contour(Cx,Cy,evalbox,[0.5], colors='y',linewidths=0.8)
            ax.text(0.025, 0.975,  f'PLD = {hrt:.2f}%\nPFD = {fat:.2f}%', 
                     transform=ax.transAxes, verticalalignment='top', fontsize=7)
            # ax.text(0.025, 0.975,  f'PLD = {hrt:.2f}%\nPFD = {fat:.2f}%', 
            #          transform=ax.transAxes, verticalalignment='top', fontsize=7)
            ax.set_title(ppi_desc)
            ax.axis('equal')
            phelp.add_cbar(ax.collections[0],fig,ax,unit_text = varunit_table['Zh'],size='3%')
        plot_nf(axs[0],data['nfout'],hr_pre[-1],fa_pre[-1],eval_pre)
        plot_nf(axs[1],proc_nf,hr_pos[-1],fa_pos[-1],eval_pos)
        if label_on:
            axs[0].plot(sitex/1e3, sitey/1e3, 'r*', markersize=8)
            axs[1].plot(sitex/1e3, sitey/1e3, 'r*', markersize=8)
        fig.savefig(os.path.join(savedir, ppi_id[:-4]+'.png'))
        plt.close(fig)

    summ_pre = log_stat(os.path.join(savedir, 'stat_pre.npz'),stat_pre)
    summ_pos = log_stat(os.path.join(savedir, 'stat_pos.npz'),stat_pos)
    PLD = 1e2*np.sum(summ_pre['hits'])/(np.sum(summ_pre['hits'])+np.sum(summ_pre['miss']))
    PFD = 1e2*np.sum(summ_pre['fa'])/(np.sum(summ_pre['q_arc']))
    PLDp = 1e2*np.sum(summ_pos['hits'])/(np.sum(summ_pos['hits'])+np.sum(summ_pos['miss']))
    PFDp = 1e2*np.sum(summ_pos['fa'])/(np.sum(summ_pos['q_arc']))

    figst, axs = plt.subplots(1, 1, figsize=(4/0.65, 3/0.65),dpi=250, gridspec_kw=dict(left=0.08, right=1-0.085, top=1-0.08, bottom=0.06, wspace=0.25, hspace=0.16))

    axs.plot(tvec,hr_pre,'b-',label=f'PLD : {PLD:.1f}%')
    axs.plot(tvec,np.array(fa_pre),'r-',label=f'PFD : {PFD:.1f}%')
    axs.plot(tvec,hr_pos, 'b--',label=f'PLD* : {PLDp:.1f}%')
    axs.plot(tvec,np.array(fa_pos), 'r--',label=f'PFD* : {PFDp:.1f}%')
    axs.set_ylim(0,100.5)
    axs.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    figst.autofmt_xdate()
    axs.set_title(case_name,loc='left')
    # plt.xlabel('Time (Frame)')
    axs.set_ylabel('Percentage (%)')
    axs.grid()
    # lines = plt.gca().get_lines()
    plt.legend(ncol=2,loc='lower right',
    bbox_to_anchor=(1, 1.02),  # (x=1 means right end of axes, y=just above)
    borderaxespad=0,
    frameon=True)
    # figst.tight_layout()
    figst.savefig(os.path.join(savedir, 'far_com.png'),bbox_inches='tight')
    plt.close(figst)
    return wgfspace

if __name__ == '__main__':
    nffig_proc(config["Settings"]["case_name"])
# case_name = 