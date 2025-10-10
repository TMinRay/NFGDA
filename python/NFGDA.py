import time
tic = time.time()
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.interpolate import LinearNDInterpolator
import scipy.io
from scipy.signal import medfilt2d
# from scipy.ndimage import gaussian_filter
from skimage.morphology import skeletonize, disk, binary_dilation, remove_small_objects
import matplotlib.pyplot as plt
import glob
import sys
import os
import configparser
# from tminlib.utility import *
from datetime import datetime

config = configparser.ConfigParser()
config.read("NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
PARROT_mask_on = True

thrREF = -5
thrdREF = 0.3
RegR = np.arange(0,400)/4
RegAZ = np.arange(0,360,0.5)*np.pi/180
RegPolarX = RegR[:,np.newaxis] * np.sin(RegAZ[np.newaxis,:])
RegPolarY = RegR[:,np.newaxis] * np.cos(RegAZ[np.newaxis,:])
interpolator = LinearNDInterpolator((RegPolarX.reshape(-1),RegPolarY.reshape(-1)), np.zeros(RegPolarX.shape).reshape(-1))
Cx, Cy = np.meshgrid(np.arange(-100,100.5,0.5),np.arange(-100,100.5,0.5))

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
def rot_displace(dp,origindeg):
    dpvector = np.swapaxes(dp,1,2)
    origindeg = origindeg*np.pi/180
    backprocess = np.array([[np.cos(origindeg),np.sin(origindeg)], \
                            [-np.sin(origindeg),np.cos(origindeg)]])
    rotcord = np.matmul(backprocess,dpvector)
    rotidx = np.round(rotcord)
    return np.swapaxes(rotidx,1,2)

def make_ftc_cscore(c_para):
    cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
    def f(cbox):
        # params is captured from outer scope
        llscore = np.zeros(cbox.shape)
        # llscore = np.full(cbox.shape,np.nan)
        llscore[cbox<=cnum1] = gaussmf(cbox[cbox<=cnum1], csig1, cnum1)*cfactor1+cintersec1
        llscore[np.logical_and(cbox>cnum1, cbox<=cnum2)] = cyfill
        llscore[cbox>cnum2] = gaussmf(cbox[cbox>cnum2], csig2, cnum2)*cfactor2+cintersec2
        return llscore
    return f

def make_ftc_sscore(s_para):
    snum1, snum2, ssig1, sfactor1, sintersec1, ssig2, sfactor2, sintersec2, syfill = s_para
    def f(sbox):
        # params is captured from outer scope
        ssscore = np.zeros(sbox.shape)
        # ssscore = np.full(sbox.shape,np.nan)
        ssscore[sbox<snum1] = syfill
        con1 = np.logical_and(sbox>=snum1, sbox<=snum2)
        con2 = sbox>snum2
        ssscore[con1] = gaussmf(sbox[con1], ssig1, snum1)*sfactor1 + sintersec1
        ssscore[con2] = gaussmf(sbox[con2], ssig2, snum2)*sfactor2 + sintersec2
        return ssscore
    return f

class FTC_PLAN:
    def __init__(self,displace,scorefun,scale):
        self.displace = displace
        self.scorefun = scorefun
        self.numINT = np.max(np.abs(displace))
        self.scale = scale
    def gather_pixels(self,ar,center):
        # idx [direction, displacement, center_pixel, yx]
        idx = (center[np.newaxis,np.newaxis,:,:] + self.displace).astype(int)
        return ar[idx[:,:,:,0],idx[:,:,:,1]]
    def get_score(self,ar,center):
        cbox = self.gather_pixels(ar,center)
        pixel_score = self.scorefun(cbox)
        # pixel_score[np.isnan(pixel_score)] = -3
        return np.nansum(pixel_score,axis=1)

def gen_beta(a2,a2_thr,ftcs):
    center_indices = np.argwhere(a2>a2_thr)
    c_indices = clean_indices(center_indices, a2.shape, ftcs[0].numINT)
    total_score = ftcs[0].get_score(a2,c_indices)
    score_scale = ftcs[0].scale*ftcs[0].displace.shape[1]
    for ftcplan in ftcs[1:]:
        total_score += ftcplan.get_score(a2,c_indices)
        score_scale += ftcplan.scale*ftcplan.displace.shape[1]
    total_score = np.max(total_score,axis=0)
    scoremt = np.zeros(a2.shape)
    scoremt[c_indices[:,0],c_indices[:,1]] = total_score/score_scale
    return scoremt
######### FTC Beta Z, dZ displacements ###########
datacy = np.arange(-8,9).reshape(1,-1)
datacx = np.zeros((1,17))
datac = np.swapaxes(np.array([datacy,datacx]),0,2)

datasy = np.array([*np.arange(-7,0,2),0,*np.arange(1,8,2),*np.arange(-7,0,2),0,*np.arange(1,8,2)]).reshape(1,-1)
datasx = np.array([-4*np.ones((9)),4*np.ones((9))]).reshape(1,-1)
datas = np.swapaxes(np.array([datasy,datasx]),0,2)

ftcc=[]
ftcs=[]

rotdegree = 180/9
for irot in np.arange(0,180,rotdegree):
    ftcc.append(rot_displace(datac,irot))
    ftcs.append(rot_displace(datas,irot))
ftcc = np.array(ftcc)
ftcs = np.array(ftcs)
######### FTC Beta Z, dZ displacements ###########
mvdiscx = np.zeros((17,17))
mvdiscy = np.zeros((17,17))
for ix in range(17):
    mvdiscx[ix,:]=np.ceil(np.arange(-8,9)*np.sin(np.pi/2/8*(ix)))
    mvdiscy[ix,:]=np.ceil(np.arange(-8,9)*np.cos(np.pi/2/8*(ix)))
nccx = mvdiscx.shape[0]
mvdisc = np.swapaxes(np.array([mvdiscy,mvdiscx]),0,2)[:,:,np.newaxis,:]

def gaussmf(x, sigma, c):
    return np.exp(-((x - c) ** 2) / (2 * sigma ** 2))

def clean_indices(idx,shp,edg):
    dim0 = idx[:,0]
    dim1 = idx[:,1]
    inbox = (dim0>=edg) & (dim0< shp[0]-edg) & (dim1>=edg) & (dim1< shp[1]-edg)
    return idx[inbox,:]

def probor(ar):
    buf = np.zeros(ar.shape[:-1])
    for iv in range(ar.shape[-1]):
        buf = buf + ar[...,iv] - buf * ar[...,iv]
    return buf

def post_moving_avg(a2):
    center_indices = np.argwhere(a2>0)

    c_indices = clean_indices(center_indices, a2.shape, avgINT)
    cidx = (c_indices[np.newaxis,np.newaxis,...] + mvdisc).astype(int)
    cbox = a2[cidx[...,0],cidx[...,1]]
    # cbr = np.sum(cbox>thrREF,0)/datacx.size
    # cbox = cbox[:,cbr>0.5]
    cbr = np.sum(cbox>0,axis=0)/nccx
    validcenter = np.max(cbr>0.1,axis=0)
    cbox = cbox[:,:,validcenter,...]
    mc = np.nanmean(cbox,axis=0)
    mc[np.logical_not(cbr>0.1)]=0
    s_indices = c_indices[validcenter,:]
    result = np.zeros(a2.shape)
    result[s_indices[:,0],s_indices[:,1]] = np.max(mc,axis=0)
    return result

class NFModule:
    def __init__(self,fismat):
        buf = scipy.io.loadmat(fismat)
        # [1,1,rule,vars]
        self.c = buf['incoef'][:,1,:][np.newaxis,np.newaxis,...]
        self.sig = buf['incoef'][:,0,:][np.newaxis,np.newaxis,...]
        self.outcoef = buf['outcoef'][np.newaxis,np.newaxis,...]
        self.rulew = buf['rulelogic'][:,0]
        self.rulecon = buf['rulelogic'][:,1]
    def eval_fis(self,pxls):
        # pxls[x,y,vars] -> [x,y,rule,vars]
        x = pxls[:,:,np.newaxis,:]
        irr = np.exp(-(x-self.c)**2/(2*self.sig**2))
        # w[x,y,rule]
        w = np.zeros(irr.shape[:-1])
        for ir in range(self.rulecon.size):
            if self.rulecon[ir]==1:
                w[...,ir] = np.prod(irr[...,ir,:],axis=-1)
            else:
                w[...,ir] = probor(irr[...,ir,:])
        sw = np.sum(w,axis=-1)
        pad_window = tuple([(0,0),] * (x.ndim-1) + [(0,1),])
        orr = np.pad(x, pad_width = pad_window, mode='constant', constant_values=1)
        # orr = [x,ones(size(x,1),1)];
        unImpSugRuleOut = np.sum(orr*self.outcoef,axis=-1)
        return np.sum( unImpSugRuleOut*w, axis=-1 )/sw

# # def findclosest(arr, target):
# #     return np.argmin(np.abs(arr-target))

# # class GFSpace:
# #     def __init__(self, x, y, nfout, radius):
# #         self.gf = nfout.copy()
# #         self.radius = radius
# #         self.x = x
# #         self.y = y
# #         self.groups = np.full(self.gf.shape, np.nan)
# #         self.n_groups = 0
# #         self.areas = []
# #         self.centroids = []
# #         for ix in range(self.gf.shape[0]):
# #             for iy in range(self.gf.shape[1]):
# #                 if self.gf[ix,iy]:
# #                     area, cum_pos = self.dfs(ix, iy)
# #                     self.areas.append(area)
# #                     self.centroids.append(cum_pos/area)
# #                     self.n_groups += 1

# #     def dfs(self, r, c):
# #         print(r, c)
# #         self.groups[r,c] = self.n_groups
# #         self.gf[r,c] = False
# #         curr_X = self.x[r,c]
# #         curr_Y = self.y[r,c]
# #         positions_sum = np.array([curr_X,curr_Y])
# #         area = 1
# #         c_low = findclosest(self.x[0,:], curr_X - self.radius/2)
# #         c_high = findclosest(self.x[0,:], curr_X + self.radius/2)
# #         r_low = findclosest(self.y[:,0], curr_Y - self.radius/2)
# #         r_high = findclosest(self.y[:,0], curr_Y + self.radius/2)
# #         for ix in range(r_low,r_high+1):
# #             for iy in range(c_low,c_high+1):
                
# #                 newX = self.x[ix,iy]
# #                 newY = self.y[ix,iy]
# #                 if ((newX - curr_X)**2 + (newY - curr_Y)**2 <= self.radius**2) \
# #                     and self.gf[ix,iy] and np.isnan(self.groups[ix,iy]):
# #                     new_area, new_positions_sum = self.dfs(ix, iy)
# #                     area = area + new_area
# #                     positions_sum = positions_sum + new_positions_sum
# #         return area, positions_sum

# class GFSpace:
#     def __init__(self, x, y, nfout, radius):
#         self.radius = radius
#         self.xp = x[nfout]
#         self.yp = y[nfout]
#         self.pos = np.array([self.xp, self.yp])
#         self.gp = np.zeros(self.xp.shape,dtype=int)
#         self.n_groups = 0
#         for p in range(self.yp.size):
#             if self.gp[p] == 0:
#                 self.n_groups += 1
#                 self.buf = [p]
#                 self.n_proc = 0
#                 while len(self.buf)>self.n_proc:
#                     self.dfs(self.buf[self.n_proc],self.n_groups)
#                     self.n_proc += 1

#         self.groups = np.zeros(x.shape)
#         self.groups[nfout] = self.gp
#         self.areas = np.zeros((1,self.n_groups))
#         self.centroids = np.zeros((2,self.n_groups))
#         self.gp -= 1
#         for ix, iy, ig in zip(self.xp, self.yp, self.gp):
#             self.areas[0,ig] = self.areas[0,ig] + 1
#             self.centroids[:,ig] = self.centroids[:,ig] + np.array([ix,iy])
#         self.centroids = self.centroids/self.areas

#     def dfs(self,p,color):
#         in_this = np.logical_and(((self.xp - self.xp[p])**2 + (self.yp - self.yp[p])**2) < self.radius**2 , self.gp==0)
#         self.gp[in_this] = color
#         self.buf.extend(np.argwhere(in_this))
#         # for bp in np.argwhere(in_this):
#         #     self.dfs(bp,color)

fuzzGST = NFModule('NF00ref_YHWANG_fis4python.mat')

def nfgda_unit_step(nexrad_0,nexrad_1,process_tag):
    ifn = nexrad_1
    print(ifn)
    exp_preds_event = export_preds_dir + process_tag
    os.makedirs(exp_preds_event,exist_ok=True)
    label_path = os.path.join('../V06/',process_tag,process_tag+'_labels')

    buf = np.load(nexrad_0)
    PARROT0 = buf['PARROT']
    if PARROT_mask_on:
        PARROT0[buf['mask']] = np.nan
    PARROT0 = np.asfortranarray(PARROT0)

    PARROT_buf = np.load(nexrad_1)
    PARROT = PARROT_buf['PARROT']
    if PARROT_mask_on:
        PARROT[PARROT_buf['mask']] = np.nan
    PARROT = np.asfortranarray(PARROT)
    
    # z1 = PARROT[:,:,0].copy()
    # z0 = PARROT0[:,:,0].copy()
    # diffz = z1-z0
    diffz = PARROT[:,:,0] - PARROT0[:,:,0]
    
    PARITP = np.zeros((*Cx.shape,PARROT.shape[-1]))
    
    for iv in [0,1,3,4,5]:
        if iv == 3:
            sdphi=np.zeros((*RegPolarX.shape,5))
            phi = PARROT[:,:,iv]
            phi[phi<0] = np.nan
            phi[phi>360] = np.nan
            sdphi[4:-2,:,:]=sliding_window_view(phi[2:,:], 5, axis=0)
            interpolator.values = np.nanstd(sdphi,axis = 2, ddof=1).reshape(-1,1)
        else:
            interpolator.values = PARROT[:,:,iv].reshape(-1,1)
        PARITP[:,:,iv] = interpolator(Cx, Cy)
    # scipy.io.savemat('../mat/pyPARROT.mat', {"PARITP": PARITP})

    V_window = sliding_window_view(PARITP[:,:,1], (3, 3))
    V_window = V_window.reshape((*V_window.shape[:2],-1))
    cbr = np.sum(~np.isnan(V_window),axis = 2)/9
    SD_buf = np.zeros(V_window.shape[:2])
    SD_buf[cbr>=0.3] = np.nanstd(V_window[cbr>=0.3].reshape(-1,9),axis = 1, ddof=1)
    stda = np.zeros(Cx.shape)
    stda[1:-1,1:-1] = SD_buf
    
    ########## FTC beta #############
    ############# Beta Cell ############
    a2 = PARITP[:,:,0]
    center_indices = np.argwhere(a2>cellthresh)
    c_indices = clean_indices(center_indices, a2.shape, cellINT)
    cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
    cbox = a2[cidx[:,:,0],cidx[:,:,1]]
    cbr = np.sum( cbox>cellthresh,0)/Celldp.shape[0]
    cbox = cbox[:,cbr>cbcellthrsh]
    c_indices = c_indices[cbr>cbcellthrsh,:]
    llscore = np.zeros(cbox.shape)
    llscore[cbox<=s2xnum[0]] = s2ynum[0]
    pp = np.logical_and(cbox>=s2xnum[0], cbox<s2xnum[1])
    llscore[pp] = s2g*cbox[pp]+s2gc
    llscore[cbox>=s2xnum[1]] = s2ynum[1]
    clscore = np.nansum(llscore,0)/Celldp.shape[0]
    totscore = np.zeros(Cx.shape)
    totscore[c_indices[:,0],c_indices[:,1]] = clscore
    CELLline = medfilt2d(totscore, kernel_size=11)
    
    a2 = CELLline
    center_indices = np.argwhere(a2>cellcsrthresh)
    c_indices = clean_indices(center_indices, a2.shape, widecellINT)
    cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
    cbox = a2[cidx[:,:,0],cidx[:,:,1]]>cellcsrthresh
    cbr = np.sum( cbox>cellcsrthresh,0)/Celldp.shape[0]
    center_indices = c_indices[cbr<1,:]
    cidx = (center_indices[np.newaxis,:,:] + Celldpw).astype(int)
    # a2[cidx[:,:,0],cidx[:,:,1]] = 1
    widecellz = a2>0.5
    ############# Beta Cell ############
    ############# Beta Z, dZ ############
    rotgz = PARITP[:,:,0]
    interpolator.values = diffz.reshape(-1,1)
    rotitp = interpolator(Cx, Cy)
    
    # # cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
    # z_cfun = make_ftc_cscore([15, 20, 3, 3, -1, 12, 4, -2, 3])
    # z_sfun = make_ftc_sscore([0, 5, 5, 2,-1, 5, 3,-3, 1])
    z_cfun = make_ftc_cscore([3, 8, 3, 3, -1, 12, 4, -2, 3])
    z_sfun = make_ftc_sscore([-8, -3, 5, 2,-1, 5, 3,-3, 1])
    zftcs = [FTC_PLAN(ftcc,z_cfun,3), FTC_PLAN(ftcs,z_sfun,1)]
    zbeta = gen_beta(rotgz,thrREF,zftcs)
    # dz_cfun = make_ftc_cscore([5,15,4,3,-2,9,4,-3,2])
    # dz_sfun = make_ftc_sscore([-10,5,5,2,-1,8,2,-3,1])
    dz_cfun = make_ftc_cscore([0,10,4,3,-2,9,4,-3,2])
    dz_sfun = make_ftc_sscore([-10,-5,5,2,-1,8,2,-3,1])
    dzftcs = [FTC_PLAN(ftcc,dz_cfun,2), FTC_PLAN(ftcs,dz_sfun,1)]
    dzbeta = gen_beta(rotitp,thrdREF,dzftcs)
    ############# Beta Z, dZ ############
    zbeta[zbeta<0]=0
    dzbeta[dzbeta<0]=0
    pbeta = (zbeta+dzbeta)/2
    pbeta[np.isnan(PARITP[:,:,0])] = np.nan
    beta = pbeta-widecellz
    beta[beta<0] = 0
    ########## FTC beta #############

    ########## NFGDA eval ###########
    inputNF = np.zeros((*Cx.shape,6))
    inputNF[:,:,0] = beta
    inputNF[:,:,1] = PARITP[:,:,0] # reflectivity
    inputNF[:,:,2] = PARITP[:,:,4] # cross_correlation_ratio
    inputNF[:,:,3] = PARITP[:,:,5] # differential_reflectivity
    inputNF[:,:,4] = stda
    inputNF[:,:,5] = PARITP[:,:,3]

    pnan = np.isnan(inputNF)
    pnansum = np.max(pnan,2)
    inputNF[pnansum,:] = np.nan

    outputGST = fuzzGST.eval_fis(inputNF)
    ########## NFGDA raw output ###########
    ########## post-processing  ###########
    # hh = outputGST>=0.24
    hh = outputGST>=0.6
    hGST = medfilt2d(hh.astype(float), kernel_size=3)
    # smoothedhGST = gaussian_filter(hGST, sigma=1, mode='nearest')
    # skel_nfout = skeletonize(smoothedhGST > 0.3)

    binary_mask = post_moving_avg(hGST) >= 0.6  # Thresholding
    pskel_nfout = binary_dilation(binary_mask, disk(5))
    skel_nfout = skeletonize(pskel_nfout*hh)
    skel_nfout2 = remove_small_objects(skel_nfout, min_size=10, connectivity=2)

    matout = os.path.join(exp_preds_event,'nf_pred'+os.path.basename(ifn)[5:-3]+'mat')
    data_dict = {
    # "xi2":Cx,"yi2":Cy,"REF":PARITP[:,:,0], \
                "nfout": skel_nfout2,"inputNF":inputNF,
                "timestamp":PARROT_buf["timestamp"]}
    if evalbox_on:
        mhandpick = os.path.join(label_path,ifn.split('/')[-1][9:-4]+'.mat')
        try:
            handpick = scipy.io.loadmat(mhandpick)
            evalbox = handpick['evalbox']
        except:
            print(f'Warning: No {mhandpick} filling zeros.')
            evalbox = np.zeros(Cx.shape)
        interpolator.values = diffz.reshape(-1,1)
        diffz = interpolator(Cx, Cy)
        data_dict.update({"evalbox":evalbox, \
            "diffz": diffz, \
            'outputGST':outputGST})

    scipy.io.savemat(matout, data_dict)
    np.savez(matout[:-3]+'npz', **data_dict)

def nfgda_proc(case_name):
    # exp_preds_event = export_preds_dir + case_name
    # os.makedirs(exp_preds_event,exist_ok=True)
    # label_path = os.path.join('../V06/',case_name,case_name+'_labels')

    v6m_path = os.path.join('../mat/','POLAR',case_name)
    v6m_list = glob.glob(v6m_path + "/polar*npz")
    # PARROT_buf = np.load(v6m_list[0])
    # PARROT0 = PARROT_buf['PARROT']
    # if PARROT_mask_on:
    #     PARROT0[PARROT_buf['mask']] = np.nan
    # PARROT0 = np.asfortranarray(PARROT0)

    for iv in range(min((config.getint('Settings', 'i_end'),len(v6m_list)))-1):
        nfgda_unit_step(v6m_list[iv],v6m_list[iv+1],case_name)
    # for ifn in v6m_list[1:config.getint('Settings', 'i_end')+1]:
        # # ifn = v6m_list[iv6m+1]
        # print(ifn)

        # ########### polar coordinate (PARROT) to cartesian coordinate (PARITP) ########
        # PARROT_buf = np.load(ifn)
        # PARROT = PARROT_buf['PARROT']
        # if PARROT_mask_on:
        #     PARROT[PARROT_buf['mask']] = np.nan
        # PARROT = np.asfortranarray(PARROT)
        
        # z1 = PARROT[:,:,0].copy()
        # z0 = PARROT0[:,:,0].copy()
        # diffz = z1-z0
        
        # PARITP = np.zeros((*Cx.shape,PARROT.shape[-1]))
        
        # for iv in [0,1,3,4,5]:
        #     if iv == 3:
        #         sdphi=np.zeros((*RegPolarX.shape,5))
        #         phi = PARROT[:,:,iv]
        #         phi[phi<0] = np.nan
        #         phi[phi>360] = np.nan
        #         sdphi[4:-2,:,:]=sliding_window_view(phi[2:,:], 5, axis=0)
        #         interpolator.values = np.nanstd(sdphi,axis = 2, ddof=1).reshape(-1,1)
        #     else:
        #         # interpolator.values = PARROT[1:,:,iv].reshape(-1,1)
        #         interpolator.values = PARROT[:,:,iv].reshape(-1,1)
        #     PARITP[:,:,iv] = interpolator(Cx, Cy)
        # # scipy.io.savemat('../mat/pyPARROT.mat', {"PARITP": PARITP})

        # V_window = sliding_window_view(PARITP[:,:,1], (3, 3))
        # V_window = V_window.reshape((*V_window.shape[:2],-1))
        # cbr = np.sum(~np.isnan(V_window),axis = 2)/9
        # SD_buf = np.zeros(V_window.shape[:2])
        # SD_buf[cbr>=0.3] = np.nanstd(V_window[cbr>=0.3].reshape(-1,9),axis = 1, ddof=1)
        # stda = np.zeros(Cx.shape)
        # stda[1:-1,1:-1] = SD_buf
        
        # ########## FTC beta #############
        # ############# Beta Cell ############
        # a2 = PARITP[:,:,0]
        # center_indices = np.argwhere(a2>cellthresh)
        # c_indices = clean_indices(center_indices, a2.shape, cellINT)
        # cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
        # cbox = a2[cidx[:,:,0],cidx[:,:,1]]
        # cbr = np.sum( cbox>cellthresh,0)/Celldp.shape[0]
        # cbox = cbox[:,cbr>cbcellthrsh]
        # c_indices = c_indices[cbr>cbcellthrsh,:]
        # llscore = np.zeros(cbox.shape)
        # llscore[cbox<=s2xnum[0]] = s2ynum[0]
        # pp = np.logical_and(cbox>=s2xnum[0], cbox<s2xnum[1])
        # llscore[pp] = s2g*cbox[pp]+s2gc
        # llscore[cbox>=s2xnum[1]] = s2ynum[1]
        # clscore = np.nansum(llscore,0)/Celldp.shape[0]
        # totscore = np.zeros(Cx.shape)
        # totscore[c_indices[:,0],c_indices[:,1]] = clscore
        # CELLline = medfilt2d(totscore, kernel_size=11)
        
        # a2 = CELLline
        # center_indices = np.argwhere(a2>cellcsrthresh)
        # c_indices = clean_indices(center_indices, a2.shape, widecellINT)
        # cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
        # cbox = a2[cidx[:,:,0],cidx[:,:,1]]>cellcsrthresh
        # cbr = np.sum( cbox>cellcsrthresh,0)/Celldp.shape[0]
        # center_indices = c_indices[cbr<1,:]
        # cidx = (center_indices[np.newaxis,:,:] + Celldpw).astype(int)
        # # a2[cidx[:,:,0],cidx[:,:,1]] = 1
        # widecellz = a2>0.5
        # ############# Beta Cell ############
        # ############# Beta Z, dZ ############
        # rotgz = PARITP[:,:,0]
        # interpolator.values = diffz.reshape(-1,1)
        # rotitp = interpolator(Cx, Cy)
        
        # # # cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
        # # z_cfun = make_ftc_cscore([15, 20, 3, 3, -1, 12, 4, -2, 3])
        # # z_sfun = make_ftc_sscore([0, 5, 5, 2,-1, 5, 3,-3, 1])
        # z_cfun = make_ftc_cscore([3, 8, 3, 3, -1, 12, 4, -2, 3])
        # z_sfun = make_ftc_sscore([-8, -3, 5, 2,-1, 5, 3,-3, 1])
        # zftcs = [FTC_PLAN(ftcc,z_cfun,3), FTC_PLAN(ftcs,z_sfun,1)]
        # zbeta = gen_beta(rotgz,thrREF,zftcs)
        # # dz_cfun = make_ftc_cscore([5,15,4,3,-2,9,4,-3,2])
        # # dz_sfun = make_ftc_sscore([-10,5,5,2,-1,8,2,-3,1])
        # dz_cfun = make_ftc_cscore([0,10,4,3,-2,9,4,-3,2])
        # dz_sfun = make_ftc_sscore([-10,-5,5,2,-1,8,2,-3,1])
        # dzftcs = [FTC_PLAN(ftcc,dz_cfun,2), FTC_PLAN(ftcs,dz_sfun,1)]
        # dzbeta = gen_beta(rotitp,thrdREF,dzftcs)
        # ############# Beta Z, dZ ############
        # zbeta[zbeta<0]=0
        # dzbeta[dzbeta<0]=0
        # pbeta = (zbeta+dzbeta)/2
        # pbeta[np.isnan(PARITP[:,:,0])] = np.nan
        # beta = pbeta-widecellz
        # beta[beta<0] = 0
        # ########## FTC beta #############

        # ########## NFGDA eval ###########
        # inputNF = np.zeros((*Cx.shape,6))
        # inputNF[:,:,0] = beta
        # inputNF[:,:,1] = PARITP[:,:,0] # reflectivity
        # inputNF[:,:,2] = PARITP[:,:,4] # cross_correlation_ratio
        # inputNF[:,:,3] = PARITP[:,:,5] # differential_reflectivity
        # inputNF[:,:,4] = stda
        # inputNF[:,:,5] = PARITP[:,:,3]

        # pnan = np.isnan(inputNF)
        # pnansum = np.max(pnan,2)
        # inputNF[pnansum,:] = np.nan

        # outputGST = fuzzGST.eval_fis(inputNF)
        # ########## NFGDA raw output ###########
        # ########## post-processing  ###########
        # # hh = outputGST>=0.24
        # hh = outputGST>=0.6
        # hGST = medfilt2d(hh.astype(float), kernel_size=3)
        # # smoothedhGST = gaussian_filter(hGST, sigma=1, mode='nearest')
        # # skel_nfout = skeletonize(smoothedhGST > 0.3)

        # binary_mask = post_moving_avg(hGST) >= 0.6  # Thresholding
        # pskel_nfout = binary_dilation(binary_mask, disk(5))
        # skel_nfout = skeletonize(pskel_nfout*hh)
        # skel_nfout2 = remove_small_objects(skel_nfout, min_size=10, connectivity=2)

        # # radar_id = ifn.split('_')[-3][:4]
        # # date_part = ifn.split('_')[-3][4:]
        # # time_part = ifn.split('_')[-2]
        # # tstamp_date = datetime.strptime(date_part, "%Y%m%d")
        # # tstamp_time = datetime.strptime(time_part, "%H%M%S").time()

        # matout = os.path.join(exp_preds_event,'nf_pred'+os.path.basename(ifn)[5:-3]+'mat')
        # data_dict = {"xi2":Cx,"yi2":Cy,"REF":PARITP[:,:,0], \
        #             "nfout": skel_nfout2,"inputNF":inputNF,
        #             # "timestamp":np.datetime64(datetime.strptime(date_part+time_part, "%Y%m%d%H%M%S"))}
        #             "timestamp":PARROT_buf["timestamp"]}
        # if evalbox_on:
        #     mhandpick = os.path.join(label_path,ifn.split('/')[-1][9:-4]+'.mat')
        #     try:
        #         handpick = scipy.io.loadmat(mhandpick)
        #         evalbox = handpick['evalbox']
        #     except:
        #         print(f'Warning: No {mhandpick} filling zeros.')
        #         evalbox = np.zeros(Cx.shape)
        #     interpolator.values = diffz.reshape(-1,1)
        #     diffz = interpolator(Cx, Cy)
        #     data_dict.update({"evalbox":evalbox, \
        #         "diffz": diffz, \
        #         'outputGST':outputGST})

        # scipy.io.savemat(matout, data_dict)
        # np.savez(matout[:-3]+'npz', **data_dict)

        # PARROT0 = PARROT
        # # nf_history.append( GFSpace(Cx, Cy, skel_nfout2, 4) )
        # # gfworker = GFSpace(Cx, Cy, skel_nfout2, 4)
        # # print(gfworker.n_groups)
        # # plt.pcolormesh(gfworker.groups)
        # # plt.show()
        # # exit()

    toc = time.time()  # End timer
    print(f"Elapsed time: {toc - tic:.6f} seconds")


if __name__ == '__main__':
    nfgda_proc(config["Settings"]["case_name"])
