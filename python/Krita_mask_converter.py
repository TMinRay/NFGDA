from PIL import Image
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from tminlib.utility import *
import configparser
import os
from skimage.filters import median
from skimage.morphology import erosion, disk

config = configparser.ConfigParser()
config.read("./NFGDA.ini")
export_preds_dir = config["Settings"]["export_preds_dir"]
# evalbox_on = config.getboolean('Settings', 'evalbox_on')
fig_dir = config["Settings"]["label_dir"]
case_name = config["Settings"]["case_name"]
savedir = os.path.join(fig_dir, case_name)
maskdir = os.path.join(savedir, 'mask')
labeldir = savedir+f'/{case_name}_labels'
os.makedirs(labeldir,exist_ok=True)
# if config["labels"]["label_on"]:
#     label_loc = list(map(float,config.get("labels", "loc").split(",")))
#     radar_loc = list(map(float,config.get("labels", "rloc").split(",")))
#     sitex, sitey = mk.geopoints_to_relative_xy(radar_loc,label_loc)
# # export_statistics_dir = config["Settings"]["export_statistics_dir"]
# # os.makedirs(export_statistics_dir,exist_ok=True)
# export_preds_fapos_dir = export_preds_dir[:-1]+'_pos/'
# os.makedirs(export_preds_fapos_dir,exist_ok=True)

# def is_red(rgb, threshold=200):
#     """Check if a pixel is 'red enough' (R high, G/B low)."""
#     r, g, b = rgb
#     return r > threshold and g < threshold // 2 and b < threshold // 2

def load_mask_image_as_red_binary(path, red_thresh=200):
    """Load an image and return binary mask where red pixels are True."""
    im = Image.open(path).convert("RGB")
    im_np = np.array(im)
    
    # Vectorized red detection
    r, g, b = im_np[..., 0], im_np[..., 1], im_np[..., 2]
    mask = (r > red_thresh) & (g < red_thresh // 2) & (b < red_thresh // 2)
    return mask  # shape (H, W)

def get_mask_values_from_xy(xx, yy, mask_image,
                            left=256, right=1024,
                            bottom=53, top=821,
                            xlim=(-100.25, 100.25),
                            ylim=(-100.25, 100.25)):
    """
    Given arrays of x and y coordinates (data space),
    return mask values from the mask image drawn over a plot.

    Parameters:
    - xx, yy: arrays of same shape (data coordinates)
    - mask_image: 2D NumPy array from PIL.Image (e.g., grayscale mask)
    - left, right, top, bottom: bounding box of plot in image pixels
    - xlim, ylim: data coordinate limits
    """
    # Normalize data coords to [0, 1] range
    x_frac = (xx - xlim[0]) / (xlim[1] - xlim[0])
    y_frac = (yy - ylim[0]) / (ylim[1] - ylim[0])

    # Convert to pixel coordinates in the image
    x_px = left + x_frac * (right - left)
    y_px = 892- 53 - y_frac * (top - bottom)  # y axis is flipped in image

    # Convert to int pixel index
    x_px = np.round(x_px).astype(int)
    y_px = np.round(y_px).astype(int)

    # Image size
    h, w = mask_image.shape

    # Clip to valid indices
    x_px = np.clip(x_px, 0, w - 1)
    y_px = np.clip(y_px, 0, h - 1)

    # Sample mask values
    return mask_image[y_px, x_px]

Cx, Cy = np.meshgrid(np.arange(-100,100.5,0.5),np.arange(-100,100.5,0.5))

evalboxfn = ls(os.path.join(maskdir, '*.png'))
exp_preds_event = export_preds_dir + case_name
npz_list = ls(exp_preds_event + "/*npz")
for ppi_file, krita_mask in zip(npz_list,evalboxfn):
    print(ppi_file, krita_mask)
    mask = load_mask_image_as_red_binary(krita_mask)
# fn = "/mnt/k/OU/NFGDA/NFGDA/python/tracking_points/label/KABX20250712_20/mask/frame0000.png"
    # mask_vals = get_mask_values_from_xy(Cx, Cy, mask)
    mask_vals = erosion(median(get_mask_values_from_xy(Cx, Cy, mask), footprint=np.ones((4,4), dtype=bool)), disk(1))
    scipy.io.savemat(os.path.join(labeldir, '_'.join(os.path.basename(ppi_file).split('_')[-3:])[:-4]+'.mat'), {'evalbox': mask_vals})

        # fig.savefig(os.path.join(labeldir, ))
# Your xx and yy arrays from pcolormesh (same shape as data)
# For example:

# xx, yy = np.meshgrid(np.linspace(-100.25, 100.25, 100),
#                      np.linspace(-100.25, 100.25, 100))

# Get mask values
    
# print(np.max(mask_vals))
# plt.pcolormesh(mask_vals)
# plt.show()
# # Optional threshold to create boolean mask
# binary_mask = mask_vals > 128