import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.interpolate import LinearNDInterpolator
import scipy.io
from scipy.signal import medfilt2d
import time
import matplotlib.pyplot as plt

angint = 0.5
rotdegree = 180/9
rotnum = int(np.round(180/rotdegree))
rotbackrad = np.deg2rad(rotdegree)
thrREF = 5
thrdREF = 0.3
cellthresh = 5
cbcellthrsh = 0.8
cellcsrthresh=0.5
crsize = 5
cellINT = crsize + 2
widecellINT =crsize+4

s2xnum = [10, 15]
s2ynum = [-3, 1]

s2xdel = s2xnum[1]-s2xnum[0]
s2ydel = s2ynum[1]-s2ynum[0]

s2g = s2ydel/s2xdel
s2gc = s2ynum[1]-s2g*s2xnum[1]

datacx = scipy.io.loadmat('../IMG/tmpCELLdatax.mat')
datacy = scipy.io.loadmat('../IMG/tmpCELLdatay.mat')
Celldp = np.moveaxis(np.array([datacy['datacy'],datacx['datacx']]),[0,1,2],[2,0,1])
# Celldp = np.moveaxis(np.array([datacx['datacx'],datacy['datacy']]),[0,1,2],[2,0,1])
np.save("Celldp.npy", Celldp)
Celldp = np.load("Celldp.npy")
datacx = scipy.io.loadmat('../IMG/tmpCELLdatax2.mat')
datacy = scipy.io.loadmat('../IMG/tmpCELLdatay2.mat')
Celldpw = np.moveaxis(np.array([datacy['datacy2'],datacx['datacx2']]),[0,1,2],[2,0,1])
# Celldpw = np.moveaxis(np.array([datacx['datacx2'],datacy['datacy2']]),[0,1,2],[2,0,1])
np.save("Celldpw.npy", Celldpw)
Celldpw = np.load("Celldpw.npy")
# exit()

# # plt.plot(Celldp[:,:,0],Celldp[:,:,1],'o')
# # plt.show()

datacy = np.arange(-8,9).reshape(1,-1)
datacx = np.zeros((1,17))
datac = np.swapaxes(np.array([datacy,datacx]),0,2)
# datac = np.swapaxes(np.array([datacx,datacy]),0,2)
datasy = np.array([*np.arange(-7,0,2),0,*np.arange(1,8,2),*np.arange(-7,0,2),0,*np.arange(1,8,2)]).reshape(1,-1)
datasx = np.array([-4*np.ones((9)),4*np.ones((9))]).reshape(1,-1)
datas = np.swapaxes(np.array([datasy,datasx]),0,2)
# datas = np.swapaxes(np.array([datasx,datasy]),0,2)

def gaussmf(x, sigma, c):
    return np.exp(-((x - c) ** 2) / (2 * sigma ** 2))

def clean_indices(idx,shp,edg):
    dim0 = idx[:,0]
    dim1 = idx[:,1]
    inbox = (dim0>=edg) & (dim0< shp[0]-edg) & (dim1>=edg) & (dim1< shp[1]-edg)
    return idx[inbox,:]

def gen_tot_score(a2,c_para,s_para,thrREF,numINT,scorediv):
    cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
    center_indices = np.argwhere(a2>thrREF)
    c_indices = clean_indices(center_indices, a2.shape, numINT)
    # dim0 = center_indices[:,0]
    # dim1 = center_indices[:,1]
    # inbox = (dim0>=numINT) & (dim0< a2.shape[0]-numINT) & (dim1>=numINT) & (dim1< a2.shape[1]-numINT)
    # c_indices = center_indices[inbox,:]
    cidx = (c_indices[np.newaxis,:,:] + datac).astype(int)

    cbox = a2[cidx[:,:,0],cidx[:,:,1]]
    cbr = np.sum(cbox>thrREF,0)/datacx.size
    cbox = cbox[:,cbr>0.5]

    s_indices = c_indices[cbr>0.5,:]
    sidx = (s_indices[np.newaxis,:,:] + datas).astype(int)
    sbox = a2[sidx[:,:,0],sidx[:,:,1]]

    llscore = np.zeros(cbox.shape)
    llscore[cbox<=cnum1] = gaussmf(cbox[cbox<=cnum1], csig1, cnum1)*cfactor1+cintersec1
    llscore[np.logical_and(cbox>cnum1, cbox<=cnum2)] = cyfill;
    llscore[cbox>cnum2] = gaussmf(cbox[cbox>cnum2], csig2, cnum2)*cfactor2+cintersec2;
    clscore = np.nansum(llscore,0)

    snum1, snum2, ssig1, sfactor1, sintersec1, ssig2, sfactor2, sintersec2, syfill = s_para

    con1 = np.logical_and(sbox>=snum1, sbox<=snum2)
    con2 = sbox>snum2

    ssscore = np.zeros(sbox.shape)
    ssscore[sbox<snum1] = syfill;
    ssscore[con1] = gaussmf(sbox[con1], ssig1, snum1)*sfactor1 + sintersec1;
    ssscore[con2]=gaussmf(sbox[con2], ssig2, snum2)*sfactor2 + sintersec2;
    sdscore = np.nansum(ssscore,0)

    pretotscore = sdscore+clscore
    pretotscore = pretotscore/scorediv
    pretotscore[pretotscore<0] = 0
    scoremt = np.zeros(a2.shape)
    scoremt[s_indices[:,0],s_indices[:,1]]=pretotscore
    return scoremt

def rot_score_back(a2,origindeg):
    backprocess = np.array([[np.cos(origindeg),np.sin(origindeg)], \
                            [-np.sin(origindeg),np.cos(origindeg)]])
    center_indices = np.argwhere(a2>0)
    rotcord = center_indices/2 - 100
    xyv = np.zeros(rotcord.shape)
    xyv[:,0] = rotcord[:,1]
    xyv[:,1] = rotcord[:,0]
    # print(rotcord.shape)
    # oldcord = np.matmul(backprocess,rotcord[:,:,np.newaxis])
    oldcord = np.matmul(backprocess,xyv[:,:,np.newaxis])
    oldidx = np.squeeze( np.round((oldcord+100)*2), axis=-1 )

    mappxl = np.logical_and( oldidx[:,0]>=0, np.logical_and(oldidx[:,0]<a2.shape[0],\
     np.logical_and(oldidx[:,1]>=0,oldidx[:,1]<a2.shape[1])))

    center_indices = center_indices[mappxl,:].astype(int)
    oldidx = oldidx[mappxl,:].astype(int)

    buf = np.zeros(a2.shape);
    buf[oldidx]=a2[center_indices]
    return buf




RegR = np.arange(0,400)/4
RegAZ = np.arange(0,360,0.5)*np.pi/180
RegPolarX = RegR[:,np.newaxis] * np.sin(RegAZ[np.newaxis,:])
RegPolarY = RegR[:,np.newaxis] * np.cos(RegAZ[np.newaxis,:])
Cx, Cy = np.meshgrid(np.arange(-100,100.5,0.5),np.arange(-100,100.5,0.5))

mat_data0 = scipy.io.loadmat('../mat/POLAR/KABX20200705_21/polar_03_KABX20200705_212755_V06.mat')
mat_data = scipy.io.loadmat('../mat/POLAR/KABX20200705_21/polar_04_KABX20200705_213306_V06.mat')
z1 = mat_data['PARROT'][:,:,0]
z0 = mat_data0['PARROT'][:,:,0]
z1[np.isnan(z1)] = 0
z0[np.isnan(z0)] = 0
diffz = z1-z0

tic = time.time()  # Start timer
PARITP = np.zeros((*Cx.shape,mat_data['PARROT'].shape[-1]))
interpolator = LinearNDInterpolator((RegPolarX[1:,:].reshape(-1),RegPolarY[1:,:].reshape(-1)), mat_data['PARROT'][1:,:,0].reshape(-1))

# for iv in range(mat_data['PARROT'].shape[-1]):
for iv in [0,1,3,4,5]:
# for iv in [1]:
    if iv == 3:
        sdphi=np.zeros((*RegPolarX.shape,5))
        phi = mat_data['PARROT'][:,:,iv]
        phi[phi<0] = np.nan
        phi[phi>360] = np.nan
        # NR = phi.shape[0]
        # for displaceR in range(-2,3):
        #     sdphi[4:-2,:,displaceR+2] = phi[4+displaceR:NR-2+displaceR,:]
        sdphi[4:-2,:,:]=sliding_window_view(phi[2:,:], 5, axis=0)
        interpolator.values = np.nanstd(sdphi,axis = 2, ddof=1)[1:,:].reshape(-1,1)
    else:
        interpolator.values = mat_data['PARROT'][1:,:,iv].reshape(-1,1)
    PARITP[:,:,iv] = interpolator(Cx, Cy)
toc = time.time()  # End timer
print(f"Elapsed time: {toc - tic:.6f} seconds")
scipy.io.savemat('../mat/pyPARROT.mat', {"PARITP": PARITP})

tic = time.time()  # Start timer

V_window = sliding_window_view(PARITP[:,:,1], (3, 3))
V_window = V_window.reshape((*V_window.shape[:2],-1))
cbr = np.sum(~np.isnan(V_window),axis = 2)/9
SD_buf = np.zeros(V_window.shape[:2])
SD_buf[cbr>=0.3] = np.nanstd(V_window[cbr>=0.3].reshape(-1,9),axis = 1, ddof=1)
stda = np.zeros(Cx.shape)
stda[1:-1,1:-1] = SD_buf
scipy.io.savemat('../mat/pystda.mat', {"stda": stda})

oriz = mat_data['PARROT'][:,:,0]
orirot = diffz
zoriginscore = np.zeros((*Cx.shape,rotnum))
originscore = np.zeros((*Cx.shape,rotnum))
for irot in range(rotnum):
    # print('----------')
    indi = int(rotdegree*irot/angint)
    origindeg = rotbackrad*irot
    rotz = np.roll(oriz, shift=indi, axis=1)
    interpolator.values = rotz.reshape(-1,1)
    rotgz = interpolator(Cx, Cy)
    ztotscore = gen_tot_score(rotgz, \
        [15, 20, 3, 3, -1, 12, 4, -2, 3], \
        [0, 5, 5, 2,-1, 5, 3,-3, 1], \
        thrREF,10,(3*17+1*18))
    zoriginscore[:,:,irot] = rot_score_back(ztotscore,origindeg)

    roted = np.roll(orirot, shift=indi, axis=1)
    interpolator.values = roted.reshape(-1,1)
    rotitp = interpolator(Cx, Cy)
    delztotscore = gen_tot_score(rotitp, \
        [5,10,4,3,-2,9,4,-3,2], \
        [-10,5,5,2,-1,8,2,-3,1], \
        thrdREF, 8, (2*17+1*18))
    originscore[:,:,irot] = rot_score_back(delztotscore,origindeg)

toc = time.time()  # End timer
print(f"Elapsed time: {toc - tic:.6f} seconds")


linez = np.max(zoriginscore,2)
linedelz = np.max(originscore,2)

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
clscore = clscore/Celldp.shape[0]
totscore = np.zeros(Cx.shape)
totscore[c_indices[:,0],c_indices[:,1]] = clscore
CELLline = medfilt2d(totscore, kernel_size=11)


# CELLline=medfilt2(totscore,[11 11]);
# exit()
a2 = CELLline

center_indices = np.argwhere(a2>cellcsrthresh)
c_indices = clean_indices(center_indices, a2.shape, widecellINT)

cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
cbox = a2[cidx[:,:,0],cidx[:,:,1]]>cellcsrthresh
cbr = np.sum( cbox>cellcsrthresh,0)/Celldp.shape[0]
center_indices = center_indices[cbr<1,:]
cidx = (c_indices[np.newaxis,:,:] + Celldpw).astype(int)

# cbox = cbox[:,cbr<1]
# c_indices = c_indices[cbr>cbcellthrsh,:]

a2[cidx[:,:,0],cidx[:,:,1]] = 1
widecellz = a2>0.5


# # %%%%%%%%%%%%%%      ../IMG/exe_3_img_exe.m
pbeta = (linez+linedelz)/2
pbeta[np.isnan(PARITP[:,:,0])] = np.nan
beta = pbeta-widecellz
beta[beta<0] = 0
plt.pcolormesh(linez)
plt.figure()
plt.pcolormesh(linedelz)
plt.figure()
plt.pcolormesh(widecellz)
plt.figure()
plt.pcolormesh(beta)
# plt.pcolor(PARITP[:,:,0])
plt.show()