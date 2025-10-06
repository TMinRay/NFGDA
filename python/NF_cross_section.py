import numpy as np
import configparser
import glob
import os
import scipy.io
from skimage.morphology import skeletonize
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import tminlib.colorlevel as cl
from tminlib import plot_helper as phelp
from scipy.ndimage import gaussian_filter
plt.style.use('dark_background')

def orientation_from_binary(mask, sigma=2.0, spacing=(1.0, 1.0), min_weight=1e-3):
    """
    Estimate local elongation direction for a binary mask.

    Parameters
    ----------
    mask : (H, W) bool/0-1 array
        True where the object exists.
    sigma : float
        Gaussian std (in pixels) for local window size. Larger -> smoother, more global.
    spacing : (dy, dx)
        Physical spacing of pixels for correct Cartesian scaling.
    min_weight : float
        Minimum local weight (foreground mass) to accept orientation.

    Returns
    -------
    theta : (H, W) float
        Dominant orientation angle in radians in [-pi/2, pi/2].
    coherence : (H, W) float
        Anisotropy measure in [0, 1]: (λ1-λ2)/(λ1+λ2).
    valid : (H, W) bool
        True where orientation is considered valid (enough local foreground).
    """
    mask = mask.astype(float)
    H, W = mask.shape
    dy, dx = spacing

    # Coordinate grids in physical units
    y, x = np.mgrid[0:H, 0:W]
    X = x * dx
    Y = y * dy

    # Weighted local sums via Gaussian filtering
    S   = gaussian_filter(mask, sigma)                         # total weight
    Sx  = gaussian_filter(mask * X, sigma)
    Sy  = gaussian_filter(mask * Y, sigma)
    Sxx = gaussian_filter(mask * X * X, sigma)
    Syy = gaussian_filter(mask * Y * Y, sigma)
    Sxy = gaussian_filter(mask * X * Y, sigma)

    # Avoid division by zero
    valid = S > min_weight

    # Local means
    mux = np.zeros_like(S); muy = np.zeros_like(S)
    mux[valid] = Sx[valid] / S[valid]
    muy[valid] = Sy[valid] / S[valid]

    # Covariance components (second central moments)
    Cxx = np.zeros_like(S); Cyy = np.zeros_like(S); Cxy = np.zeros_like(S)
    Cxx[valid] = Sxx[valid] / S[valid] - mux[valid]**2
    Cyy[valid] = Syy[valid] / S[valid] - muy[valid]**2
    Cxy[valid] = Sxy[valid] / S[valid] - mux[valid]*muy[valid]

    # Orientation: 0.5 * atan2(2Cxy, Cxx - Cyy), range [-pi/2, pi/2]
    theta = np.zeros_like(S, dtype=float)
    theta[valid] = 0.5 * np.arctan2(2 * Cxy[valid], (Cxx[valid] - Cyy[valid]))

    # Eigenvalues for coherence measure
    # λ = (Cxx + Cyy)/2 ± sqrt(((Cxx - Cyy)/2)^2 + Cxy^2)
    trace = Cxx + Cyy
    diff  = Cxx - Cyy
    root  = np.sqrt( (diff*diff)/4.0 + Cxy*Cxy )
    lam1  = (trace / 2.0) + root   # major
    lam2  = (trace / 2.0) - root   # minor

    coherence = np.zeros_like(S)
    denom = lam1 + lam2
    nz = valid & (denom > 0)
    coherence[nz] = (lam1[nz] - lam2[nz]) / denom[nz]

    return theta, coherence, valid

def cross_section(theta, x0, y0, xgrid, ygrid, length=20):
    """
    Extract cross-section perpendicular to elongate direction at (x0, y0).
    
    Parameters
    ----------
    theta : float
        Elongation direction angle at this point (radians).
    x0, y0 : float
        Physical coordinates of center point.
    xgrid, ygrid : 2D arrays
        Coordinate grids (same shape as mask).
    length : float
        Length of cross-section in same units as x0,y0.
    """
    # perpendicular direction
    angle = theta + np.pi/2
    dx, dy = np.cos(angle), np.sin(angle)
    s = np.arange(-(length//2),(length+1)//2)
    # start and end coordinates in physical space
    x, y = x0 - dx[np.newaxis,:]*s[:,np.newaxis], y0 - dy[np.newaxis,:]*s[:,np.newaxis]
    # x2, y2 = x0 + dx*length/2, y0 + dy*length/2
    dx_grid = xgrid[0,1] - xgrid[0,0]
    dy_grid = ygrid[1,0] - ygrid[0,0]
    i1 = np.clip(np.round((y - ygrid[0,0]) / dy_grid), 0, xgrid.shape[0]-1).astype(int)
    j1 = np.clip(np.round((x - xgrid[0,0]) / dx_grid), 0, xgrid.shape[1]-1).astype(int)
    # return x, y
    return i1, j1

config = configparser.ConfigParser()
config.read("./NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
train_dir = '../mat/train/'

case_name = 'KAMA20240227_21'

# case_name = 'KABX20200704_02'
# case_name = 'KABX20200721_19'
exp_preds_event = export_preds_dir + case_name
npz_list = glob.glob(exp_preds_event + "/*npz")
true_buf = []
false_buf = []
leng=15
ds = np.arange(-(leng//2),(leng+1)//2)
ppi_file=npz_list[8]
print(ppi_file)
data = np.load(ppi_file)
evalbox = data['evalbox']
evalline = skeletonize(evalbox)
Cx = data['xi2']
Cy = data['yi2']
# evalline[Cx>81]=0
# evalline[Cx<39]=0
evalline[Cx>1]=0
evalline[Cx<-51]=0
r = np.sqrt(Cx**2+Cy**2)
rmask = r>=100
# fig, axs = plt.subplots(2, 1, figsize=(4, 6.5),dpi=250, gridspec_kw=dict(left=0.18, right=1-0.085, top=1-0.11, bottom=0.08, wspace=0.25, hspace=0.16))
figp,axp=plt.subplots(1,1,figsize=(4.5, 4.5),dpi=400,gridspec_kw=dict(left=0.135, right=1-0.1, top=1-0.08, bottom=0.1, wspace=0.25, hspace=0.16))

REF = data['inputNF'][:,:,1]
pdata = np.ma.masked_where(rmask,REF)
# axs[0].pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
# axs[0].contour(Cx,Cy,evalbox,[0.5],colors='#ffff00',linewidths=1.5)
axp.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
axp.contour(Cx,Cy,evalbox,[0.5],colors='#ffff00',linewidths=1.5)

theta, coherence, valid = orientation_from_binary(evalline, sigma=3)

x,y =  cross_section(theta[evalline], data['xi2'][evalline], data['yi2'][evalline],data['xi2'],data['yi2'], length=leng)
ppi_id = os.path.basename(ppi_file)
ppi_name = ppi_id[11:]  # MATLAB 12:end is Python 11: (0-based)
date_part = ppi_name[4:12]   # 5:12 in MATLAB → 4:12 in Python
time_part = ppi_name[13:19]

radar_id = ppi_name[0:4]  # 1:4 in MATLAB → 0:4 in Python
tstamp_date = datetime.strptime(date_part, "%Y%m%d")
tstamp_time = datetime.strptime(time_part, "%H%M%S").time()
tstamp = datetime.combine(tstamp_date.date(), tstamp_time)
ppi_desc = f"{radar_id}, {tstamp.strftime('%m/%d/%Y, %H:%M:%S %Z')}"

# z = data['inputNF'][x,y,1].ravel()
cx = np.tile(ds[:, None], (1, x.shape[1])).ravel()
# dZ = data['diffz'][x,y].ravel()

z = data['inputNF'][x,y,1]
r = data['inputNF'][x,y,2].ravel()
d = data['inputNF'][x,y,3].ravel()
# box0 = axs[0].get_position()
# box1 = axs[1].get_position()


# plt.plot(data['xi2'][evalline],data['yi2'][evalline],'b.')
# axs[0].set_title(r'Z')
# fig.suptitle(ppi_desc)
# phelp.add_cbar(axs[0].collections[0],fig,axs[0],pad="9%",unit_text = 'dBZ')
# axs[0].axis('equal')
# axs[1].set_xlim(-leng//2, (leng+1)//2)
# axs[1].set_ylim(-20, 40)
# axs[1].set_xlabel('cross GF pixel')
# axs[1].set_ylabel('dBZ')
# axs[1].set_position([box0.x0+0.02, box1.y0, box0.width-0.11, box1.height])
# box1 = axs[1].get_position()
# axs[0].set_xlim(-25, 25)
# axs[0].set_ylim(-50, 0)
# axs[1].xticks(np.arange(-8,12,4))

# fig.suptitle(ppi_desc)
phelp.add_cbar(axp.collections[0],figp,axp,pad="9%",unit_text = 'dBZ')
axp.axis('equal')
axp.set_xlim(40, 80)
axp.set_ylim(30, 70)

figc,ax=plt.subplots(1,1,figsize=(4.5, 4.5),dpi=400,gridspec_kw=dict(left=0.15, right=1-0.085, top=1-0.08, bottom=0.1, wspace=0.25, hspace=0.16))
plt.plot([-1,1,1,-1,-1],[10,10,30,30,10],color=(116/255,180/255,128/255),linewidth=2,zorder=10)
plt.plot([-5,-3,-3,-5,-5],[-20,-20,0,0,-20],color=(221/255,240/255,99/255),linewidth=2,zorder=10)
plt.plot([5,3,3,5,5],[-20,-20,0,0,-20],color=(221/255,240/255,99/255),linewidth=2,zorder=10)
plt.xlim(-15//2, (15+1)//2)
plt.ylim(-20, 40)
plt.xlabel('Cross GF Pixels')
plt.ylabel('dBZ')
plt.xticks(np.arange(-8,12,4))
# plt.savefig('AMS/FTC_zone.png')
figp.savefig(f'tracking_points/AMS/p-box.png')

# for ic in np.arange(0,x.shape[1],5):
#     axp.plot(data['xi2'][x[:,ic],y[:,ic]],data['yi2'][x[:,ic],y[:,ic]],'r-',linewidth=2)
#     ax.plot(ds,z[:,ic],zorder=1,linewidth=1)
#     figc.savefig(f'tracking_points/AMS/c-{ic:0>3d}.png')
#     figp.savefig(f'tracking_points/AMS/p-{ic:0>3d}.png')
axp.plot(data['xi2'][x[:,::5],y[:,::5]],data['yi2'][x[:,::5],y[:,::5]],'r-',linewidth=2)
ax.plot(ds,z[:,::5],zorder=1,linewidth=1)
figc.savefig(f'tracking_points/AMS/c.png')
figp.savefig(f'tracking_points/AMS/p.png')
# # plt.title('')
plt.close(figp)
plt.close(figc)
# fig, axs = plt.subplots(2, 1, figsize=(4, 6.5),dpi=250, gridspec_kw=dict(left=0.18, right=1-0.085, top=1-0.11, bottom=0.08, wspace=0.25, hspace=0.16))
# REF = data['inputNF'][:,:,1]
# pdata = np.ma.masked_where(rmask,REF)
# axs[0].pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)

# axs[0].contour(Cx,Cy,evalbox,[0.5],colors='#ffff00',linewidths=1.5)
# axs[0].set_title(r'Z')
# fig.suptitle(ppi_desc)
# phelp.add_cbar(axs[0].collections[0],fig,axs[0],pad="9%",unit_text = 'dBZ')
# axs[0].axis('equal')
# # axs[1].set_xlabel('cross GF pixel')
# # axs[1].set_ylabel('dBZ')
# # z = data['inputNF'][x,y,1]
# box0 = axs[0].get_position()
# axs[0].plot(data['xi2'][x[:,::10],y[:,::10]],data['yi2'][x[:,::10],y[:,::10]],'r-')

# z = z.ravel()
# axs[1].hist2d(cx[np.isfinite(z)], z[np.isfinite(z)], bins=(np.arange(-leng//2, (leng+1)//2, 1)+0.5, np.arange(-20, 40, 1)), cmap='magma')
# phelp.add_cbar(axs[1].collections[0],fig,axs[1],pad="9%")
# axs[1].set_xlim(-leng//2, (leng+1)//2)
# axs[1].set_ylim(-20, 40)
# axs[1].set_xlabel('cross GF pixel')
# axs[1].set_ylabel('dBZ')
# axs[1].xticks(np.arange(-8,12,4))
# axs[1].set_position([box1.x0, box1.y0, box1.width/0.91, box1.height])
# fig.savefig(f'tracking_points/AMS/hist.png')
# # ax[1].set_title('Z')

z = z.ravel()
fig,ax=plt.subplots(1,1,figsize=(4.719, 4.5),dpi=400,gridspec_kw=dict(left=0.15, right=1-0.085, top=1-0.08, bottom=0.1, wspace=0.25, hspace=0.16))
# ax.hist2d(cx[np.isfinite(z)], z[np.isfinite(z)], bins=(np.arange(-leng//2, (leng+1)//2, 1)+0.5, np.arange(-20, 40, 1)), cmap='magma')
# ax.set_ylabel('dBZ')
# ax.set_ylim(-20, 40)
# fn=f'tracking_points/AMS/hist.png'

# ax.hist2d(cx[np.isfinite(r)], r[np.isfinite(r)], bins=(np.arange(-leng//2, (leng+1)//2, 1)+0.5, np.arange(0, 1.1, 0.05)), cmap='magma')
# ax.set_ylim(0.1, 1.05)
# fn=f'tracking_points/AMS/hist-R.png'

ax.hist2d(cx[np.isfinite(d)], d[np.isfinite(d)], bins=(np.arange(-leng//2, (leng+1)//2, 1)+0.5, np.arange(-9, 9, 0.2)), cmap='magma')
ax.set_ylim(-8, 8)
ax.set_ylabel('dB')
fn=f'tracking_points/AMS/hist-D.png'

phelp.add_cbar(ax.collections[0],fig,ax)
ax.set_xlim(-leng//2, (leng+1)//2)
ax.set_xlabel('Cross GF Pixels')

ax.set_xticks(np.arange(-8,12,4))
# ax.set_position([box1.x0, box1.y0, box1.width/0.91, box1.height])
# ax.plot([-1,1,1,-1,-1],[10,10,30,30,10],color=(116/255,180/255,128/255),linewidth=2,zorder=10)
# ax.plot([-5,-3,-3,-5,-5],[-20,-20,0,0,-20],color=(221/255,240/255,99/255),linewidth=2,zorder=10)
# ax.plot([5,3,3,5,5],[-20,-20,0,0,-20],color=(221/255,240/255,99/255),linewidth=2,zorder=10)
fig.savefig(fn)