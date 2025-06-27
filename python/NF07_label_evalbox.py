### Gust Front Evalbox Identification From NEXRAD Data ###
### Created 2024/10/10 1602Z ###
### J. W. Thiesing 2024 ###
### M. D. Tzeng 03/19/2025 ###

import datetime as dt
import numpy as np
from scipy.io import savemat

import matplotlib.pyplot as plt
# import matplotlib.path as path
import matplotlib as mpl
from mpl_point_clicker import clicker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cartopy.crs as ccrs

import sys
import os.path
import os, pyart

def loadnexrad(downdir='./Data/NEXRAD/', ii=0):
	radarobjs = []
	d = list(os.walk(downdir))
	# d = []
	# for c in os.walk(downdir):
	# 	if c[0].replace(downdir, '') == '': continue
	# 	if '.DS_Store' in c[0]: continue
	# 	if '_labels' in c[0]: continue
	# 	d.append(c)
	# d.sort()
	parent = d[ii][0]
	files = d[ii][2]
	files.sort()
	for file in files:
		if file[-3:]=='V06':
			path = os.path.join(parent, file)
			# path = +'/'+
			if '.DS_Store' in path or '.mat' in path: continue
			# if os.path.isfile(parent+'_labels/'+file+'.mat'): continue
			radarobjs.append(path)
	return radarobjs

def plotradarobj(rpath, stlat, stlon, savepath='./'):

	radarobj = pyart.io.read(rpath)
	plt.rcParams.update({'font.size': 14})

	savepathsplit = savepath.split('/')
	savepath = os.path.join('/'.join(savepathsplit[:-1]),savepathsplit[-2]+'_labels',savepathsplit[-1])
	print(savepath)
	# if not os.path.exists('/'.join(savepath.split('/')[:-1])): os.makedirs('/'.join(savepath.split('/')[:-1]))
	os.makedirs('/'.join(savepath.split('/')[:-1]),exist_ok=True)

	global ax0
	global ax1
	fig, ax = plt.subplots(1,2, gridspec_kw={'wspace': 0.2, 'hspace': 0.1, 'left': 0.05, 'right': 0.95, 'bottom': 0.05, 'top': 0.95})
	ax0, ax1 = ax

	start_time = dt.datetime.strptime(radarobj.time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
	end_time = start_time + dt.timedelta(seconds=radarobj.time['data'][-1])

	plot_nexrad(ax0, radarobj, 'reflectivity', find_params('reflectivity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('reflectivity')))
	plot_nexrad(ax1, radarobj, 'velocity', find_params('velocity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('velocity')))
	plt.tight_layout()
	plt.show()
	print("How many gust fronts would you like to identify in this sweep?:")
	casecount = int(input())
	if casecount == 0:
		XX, YY = np.meshgrid(np.arange(-100,100.5,0.5), np.arange(-100,100.5,0.5))
		evalbox = np.zeros((401,401))
		savemat(savepath+'.mat', {'evalbox':evalbox,'xi2':XX,'yi2':YY})
		return

	fig, ax = plt.subplots(1,2, gridspec_kw={'wspace': 0.2, 'hspace': 0.1, 'left': 0.05, 'right': 0.95, 'bottom': 0.05, 'top': 0.95})
	ax0, ax1 = ax
	plot_nexrad(ax0, radarobj, 'reflectivity', find_params('reflectivity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('reflectivity')))
	plot_nexrad(ax1, radarobj, 'velocity', find_params('velocity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('velocity')))

	markerlabels = []
	for ii in range(0, casecount):
		cstr = 'GF'+str(ii+1)
		markerlabels.append(cstr)
	markers = ["o", "X", "*", "s", "P"]
	colors = ['fuchsia']*casecount

	klicker = clicker(
		ax0,
		markerlabels,
		markers=markers[:casecount],
		colors=colors
	)

	plt.tight_layout()
	plt.show(block=False)
	manager = plt.get_current_fig_manager()
	backend = plt.get_backend()
	if backend == "TkAgg":
		manager.window.wm_geometry("1920x1080+200+150")  # Tkinter
	elif "Qt" in backend:
		manager.window.setGeometry(200, 150, 1920, 1080)
	plt.show()

	coordsall = klicker.get_positions()

	pathsall = []
	for key in coordsall.keys():
		arr = coordsall[key].tolist()
		#arr.append([0,0])
		arr = np.array(arr)
		codes = [np.uint8(1)]
		for ii in range(0, len(arr)-1): codes.append(np.uint8(2))
		#codes.append(np.uint8(79))
		#print(codes)

		pathsall.append(mpl.path.Path(arr, codes=codes))

	# Path plotting code, optional
	plotpaths_debug = False
	if plotpaths_debug:
		import matplotlib.patches as patches
		fig, ax = plt.subplots(1,2, gridspec_kw={'wspace': 0.2, 'hspace': 0.1, 'left': 0.05, 'right': 0.95, 'bottom': 0.05, 'top': 0.95})
		ax0, ax1 = ax
		plot_nexrad(ax0, radarobj, 'reflectivity', find_params('reflectivity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('reflectivity')))
		plot_nexrad(ax1, radarobj, 'velocity', find_params('velocity', 15, 100, stlat, stlon, 0, 0., start_time, end_time, ('velocity')))
		for path in pathsall:
			patch = patches.PathPatch(path, facecolor='fuchsia', lw=6)
			ax0.add_patch(patch)
		plt.show()

	XX, YY = np.meshgrid(np.arange(-100,100.5,0.5), np.arange(-100,100.5,0.5))
	XX, YY = XX*1e3, YY*1e3
	evalbox = np.zeros((401,401))
	allpoints = []
	for x in range(0,401):
			for y in range(0,401):
				allpoints.append((XX[x,y], YY[x,y]))
	for path in pathsall:
		inside = path.contains_points(allpoints)
		inside = np.reshape(inside,(401,401))
		for x in range(0,401):
			for y in range(0,401):
				if inside[x,y]: evalbox[x,y] = 1

	savemat(savepath+'.mat', {'evalbox':evalbox,'xi2':XX/1e3,'yi2':YY/1e3})

	return

def plot_nexrad(ax, radarfile, field, params):
	fields_possible = ("reflectivity", "velocity")
	if not field in fields_possible:
		sys.exit(f"Field requested not found in fields possible. Possible fields: {fields_possible}. Exiting.")

	hasvel = []
	hasdp = []
	for ii in np.arange(0,radarfile.nsweeps,1):
		if np.any(radarfile.fields['velocity']['data'][radarfile.get_slice(ii)]):
			hasvel.append(ii)
		else:
			hasdp.append(ii)
	if field == "velocity":
		desiredsweep = hasvel[np.where(radarfile.fixed_angle['data'][hasvel] == np.min(radarfile.fixed_angle['data'][hasvel]))[0][0]]
	elif field == "reflectivity":
		desiredsweep = hasdp[np.where(radarfile.fixed_angle['data'][hasdp] == np.min(radarfile.fixed_angle['data'][hasdp]))[0][0]]
	#print(radarfile.fixed_angle['data'][hasvel], np.min(radarfile.fixed_angle['data'][hasvel]), np.where(radarfile.fixed_angle['data'][hasvel] == np.min(radarfile.fixed_angle['data'][hasvel])))
	#print(desiredsweep)

	R = radarfile.range['data']
	AZ = radarfile.azimuth['data'][radarfile.get_slice(desiredsweep)]
	EL = radarfile.elevation['data'][radarfile.get_slice(desiredsweep)]
	#print(hasvel, desiredsweep)
	fields = {
		'reflectivity': radarfile.fields['reflectivity']['data'][radarfile.get_slice(desiredsweep)],
		'velocity': radarfile.fields['velocity']['data'][radarfile.get_slice(desiredsweep)]
	}

	#print(field)
	final_params = params.copy()
	#radarfile.info()
	#print(radarfile.latitude['data'][0], radarfile.longitude['data'][0], params['lat'], params['lon'])
	final_params['lat'] = radarfile.latitude['data'][0]
	final_params['lon'] = radarfile.longitude['data'][0]
	#print(final_params['lon'], radarfile.longitude['data'][0])
	#print(radarfile.latitude['data'][0])
	if field == "reflectivity":
		# Reflectivity
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = -10., 75.
		final_params['productstr'] = "Reflectivity"
		final_params['unitstr'] = "[dBZ]"
		final_params['cbtickinterval'] = 5.

	dop_fixed_el = radarfile.fixed_angle['data'][desiredsweep]

	#selectedsweep = np.where(np.min(np.abs(dop_fixed_el-elev_disp)) == np.abs(dop_fixed_el-elev_disp))[0]
	#print(dop_fixed_el[selectedsweep])
	final_params['elevation'] = dop_fixed_el

	# Single sweep PPI plot
	Rpad = np.pad(R, ((0,1)), mode='edge')
	AZpad = np.pad(AZ, ((0,1)), mode='wrap') # NOTE: OFF TO BUGFIX
	ax = product_axes(ax, fields[field], Rpad, AZpad, EL, final_params, padded=True) # Plot lidar data

	return ax

earth_radius = 6371.0 * 1e3
dndh = 2.724e-5/1e3
k = 1./(1.+(earth_radius*dndh))

# Convert the 0 to 360 to -180 to 180 
def north0_to_arctheta(theta):
	return np.where(theta > 270, 450-theta, 90-theta)

def calculate_ground_dist(Reach, eleveach):
	return k*earth_radius*np.arcsin((Reach*np.cos(np.deg2rad(eleveach)))/(k*earth_radius + make_heights(Reach, eleveach, mask=np.zeros(Reach.shape)))) # From Rauber textbook

# Create cartesian axes from polar ranges & azimuths
def dis_angle_to_2Dxy(R,theta):
	X=R[:,np.newaxis].dot(np.cos(theta[np.newaxis,:]*np.pi/180))
	Y=R[:,np.newaxis].dot(np.sin(theta[np.newaxis,:]*np.pi/180))
	return X,Y

def mesh_XY_fixed(Reach,theta,eleveach):
	#theta = np.pad(theta, ((0,1)), mode='wrap')
	thetaeach = np.tile(theta,(eleveach.shape[1],1)).T
	S = calculate_ground_dist(Reach, eleveach)
	#X=R[:,np.newaxis].dot(np.cos(theta[np.newaxis,:]*np.pi/180))
	#Y=R[:,np.newaxis].dot(np.sin(theta[np.newaxis,:]*np.pi/180))
	X = S*np.cos(thetaeach*np.pi/180)
	Y = S*np.sin(thetaeach*np.pi/180)

	return X, Y

def make_heights(Reach, elev, mask, instrumentheight=0.0):
	heights = np.sqrt(Reach**2. + (k*earth_radius)**2. + 2.0*Reach*k*earth_radius*np.sin(np.deg2rad(elev))) - k*earth_radius + instrumentheight

	# ***NOTE***: For some reason, numpy keeps rounding these float64 arrays (I checked) to the nearest .0 or .5. I have no idea why, I did not consent to it, but it happens, and it's fairly inconsequential for the final result. Fitting algorithm error is probably significantly larger than a few centimeters' difference in height that's getting binned anyway. Just thought you'd like to know!

	#print(Reach.dtype, elev.dtype)
	#print(type(Reach[0][0]), type())
	#print(Reach, elev)
	#r = 22425.
	#eleva = 0.99
	#height = np.sqrt(r**2. + (k*ae)**2. + 2.0*r*k*ae*np.sin(eleva*np.pi/180.0)) - k*ae + instrumentheight
	#print(heights)
	heights = np.ma.array(heights, mask=mask)
	return heights

def product_axes(ax, field, ranges, azimuths, elev, params, padded=False, no_az_arctheta=False):
	# Polar to cartesian, create meshed axes
	if not no_az_arctheta:
		AZ_arctheta = north0_to_arctheta(azimuths)
	else:
		AZ_arctheta = azimuths
	if not params['use_fixed_mesh']:
		XX, YY = dis_angle_to_2Dxy(ranges,AZ_arctheta)
	else:
		eleveach = np.tile(elev,(field.shape[1],1)).T
		Reach = np.tile(ranges,(eleveach.shape[0],1))
		eleveach = np.pad(eleveach, ((0,1), (0,1)), mode='edge')
		Reach = np.pad(Reach, ((0,1), (0,0)), mode='edge')
		XX, YY = mesh_XY_fixed(Reach,AZ_arctheta,eleveach)
	#XX = np.pad(XX, ((0,1), (0,1)), mode='edge')
	#YY = np.pad(YY, ((0,1), (0,1)), mode='edge')

	# Georeference axes relative to lidar
	# XX = longc(XX, params['lat']) + params['lon']
	# YY = latc(YY) + params['lat']
	#print('XXfmin: ' + str(XX.min()) + '\nXXfmax: ' + str(XX.max()) + '\nYYfmin: ' + str(YY.min()) + '\nYYfmax: ' + str(YY.max()))

	# TODO: FUTURE: ADJUST WITH MOVING CENTER FROM GPS DATA

	# Georeference limits and tick intervals
	# params['leftc'] = longc(params['left'], params['lat'])*1e3 + params['lon']
	# params['rightc'] = longc(params['right'], params['lat'])*1e3 + params['lon']
	# params['bottomc'] = latc(params['bottom'])*1e3 + params['lat']
	# params['topc'] = latc(params['top'])*1e3 + params['lat']
	params['leftc'] = params['left']*1e3
	params['rightc'] = params['right']*1e3
	params['bottomc'] = params['bottom']*1e3
	params['topc'] = params['top']*1e3
	params['plottickintervalcx'] = params['plottickinterval']*1e3
	params['plottickintervalcy'] = params['plottickinterval']*1e3
	#print("corrected bounds: " + str(params['leftc']) + ' ' + str(params['rightc']) + ' ' + str(params['bottomc']) + ' ' + str(params['topc']))

	if not padded:
		field = np.ma.array(np.pad(field, ((0,1), (0,1)), mode='edge'))

	start_dati = params['start_dati']
	end_dati = params['end_dati']

	# Plot labels
	ax.set_xlabel('Zonal  Distance [km]', fontsize=params['fontsize']+2)
	ax.set_ylabel('Meridional Distance [km]',fontsize=params['fontsize']+2)
	if start_dati.date() == end_dati.date():
		ax.set_title(params['obs_method']+' '+start_dati.strftime('%Y-%m-%d %H:%M:%S')+' - '+end_dati.strftime('%H:%M:%S')+' '+str(params['elevation'])+'${^{o}}$',fontsize=params['fontsize']+2)
	else:
		ax.set_title(params['obs_method']+' '+start_dati.strftime('%Y-%m-%d %H:%M:%S')+' - '+end_dati.strftime('%Y-%m-%d %H:%M:%S')+' '+str(params['elevation'])+'${^{o}}$',fontsize=params['fontsize']+2)
	
	# Plot data!
	#cmesh=ax.pcolormesh(np.transpose(XX),np.transpose(YY), field, cmap=params['cmap'], vmin=params['vmin'], vmax=params['vmax'], transform=ccrs.PlateCarree())
	cmesh=ax.pcolormesh(XX, YY, field, cmap=params['cmap'], vmin=params['vmin'], vmax=params['vmax'])
	#cmesh=ax.pcolormesh(XX, YY, field, cmap=params['cmap'], vmin=params['vmin'], vmax=params['vmax'], transform=ccrs.PlateCarree())

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=plt.Axes)

	# Add colorbar
	cbar=plt.colorbar(cmesh,cax=cax,pad=0.065,label=params['productstr']+' '+params['unitstr'])

	# Add ticks
	xtickloc = np.arange(params['leftc'],params['rightc']+params['plottickintervalcx']/2,params['plottickintervalcx'])
	ytickloc = np.arange(params['bottomc'],params['topc']+params['plottickintervalcy']/2,params['plottickintervalcy'])
	ax.set_xticks(xtickloc)
	ax.set_yticks(ytickloc)
	ax.set_xticklabels(np.arange(params['left'],params['right']+params['plottickinterval']/2,params['plottickinterval']))
	ax.set_yticklabels(np.arange(params['bottom'],params['top']+params['plottickinterval']/2,params['plottickinterval']))
	ax.tick_params(axis='x', labelsize=params['fontsize'])
	ax.tick_params(axis='y', labelsize=params['fontsize'])

	ax.set_xlim(params['leftc'], params['rightc'])
	ax.set_ylim(params['bottomc'], params['topc'])
	ax.grid(linestyle = '--', linewidth = 1.25, c='gray', alpha=0.5)

	return ax


def find_params(name, velext, plotext, lat, lon, alt, el_fixed, start_time, end_time, fields_possible):
	if not name in fields_possible:
		sys.exit(f"Field requested not found in fields possible. Possible fields: {fields_possible}. Exiting.")

	if lat == -999.0:
		lat = 0.0
	if lon == -999.0:
		lon = 0.0

	# TODO: Fix velocity colormap, ask Cocoa
	# Huge parameters list
	product_parameters_VR = {
		'vmin': -velext, # Minimum value of the colorbar
		'vmax': velext, # Maximum value of the colorbar
		'cbtickinterval': 5, # Scan data colorbar tick interval
		'fontsize': 14,
		'figsize': (11, 11),
		'cmap': mpl.colormaps['pyart_Carbone42'],#mpl.colormaps['cool'],
		'elevation': el_fixed, # Fixed elevation angle of sweep, used in plot title
		'productstr': "Radial Velocity", # Product name
		'unitstr': "[m ${s^{-1}}$]", # Product units
		'left': -plotext, # Left plot extent, relative to lidar location in km
		'right': plotext, # Right plot extent, relative to lidar location in km
		'bottom': -plotext, # Bottom plot extent, relative to lidar location in km
		'top': plotext, # Top plot extent, relative to lidar location in km
		'plottickinterval': 0.25*plotext, #0.25*plotext, # Plot x/y tick interval
		'lat': lat, # lidar central latitude
		'lon': lon, # lidar central longitude
		'alt': alt, # lidar altitude
		'start_dati': start_time, # Sweep starting datetime
		'end_dati': end_time, # Sweep ending datetime
		'height_bounds': (0., 10000.), # Maximum and minimum height used, m
		'height_bin_width': 25., # Width of VWP algorithm height bins in m
		'minimum_samples_per_bin': 60, # Minimum total gates for a height bin to not be removed from data (for QC)
		'height_bins_per_cb_tick': 15, # Number of height bins (groups of x size) per height colorbar tick
		'plotfit': False, # Boolean, to plot the fitting algorithm in a specific layer or not
		'plotfitbin': False, # Layer to plot fitting algorithm of
		'fixed_height': 0., # When doing single elevation vector finding, this will be the bottom of the layer analyzed. 0.0 corresponds to "minimum height found in data"
		'fixed_height_depth': 120., # When doing single elevation vector finding, this is the height of the layer analyzed.
		'obs_method': "PPI", # Observation method string
		'filter_by_cv': True, # Forgoing actual data filtering, this instead filters each layer by the covariance of the fitting function.
		'cv_threshold': 0.15, # Threshold for covariance filtering 0.005
		'filter_by_phase_shift': False, # Filter using 180 degree phase shift method
		'phase_shift_threshold': 10., # Phase shift filter threshold
		'smooth_spacing': 2.5, # Spacing of the smoothing/averaging in phase shift filter, 360/smooth_spacing must be even
		'hodograph_grid_inc': 5., # Grid interval of the hodograph plot
		'use_fixed_mesh': True, # Uses refraction and elevation incorporated ground distance rather than just R
		'zebra_filter': False, # Toggles whether data uses my zebra dilter
	}

	if product_parameters_VR['filter_by_phase_shift']:
		product_parameters_VR['cv_threshold'] = 0.0075 #0.0125
	
	if name == "velocity" or name == "VR":
		# Velocity
		final_params = product_parameters_VR

	if name == "intensity":
		# Intensity, Halo:
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = -5, 1.
		final_params['productstr'] = "Intensity"
		final_params['unitstr'] = "[dB]"
		final_params['cbtickinterval'] = 1
	
	if name == "backscatter":
		# Backscatter, Halo:
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = 0.0, 0.00007
		final_params['productstr'] = "Attenuated Backscatter"
		final_params['unitstr'] = "[${km^{-1}}{sr^{-1}}$]"
		final_params['cbtickinterval'] = 0.00001

	if name == "lppd":
		# PWR, metroweather
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = 0, .5
		final_params['productstr'] = "Linear Peak Power Density"
		final_params['unitstr'] = "[%]"
		final_params['cbtickinterval'] = 0.05

	if name == "snr":
		# Signal to noise ratio, metroweather
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = 1., 5.
		final_params['productstr'] = "Signal-To-Noise Ratio"
		final_params['unitstr'] = "[dB]"
		final_params['cbtickinterval'] = 1

	if name == "wth":
		# Spectrum width, metroweather
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = 0, 5
		final_params['productstr'] = "Spectrum Width"
		final_params['cbtickinterval'] = 1

	if name == "noiselevel":
		# Noise level, metroweather
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = 0, 0.25
		final_params['productstr'] = "Noise Level"
		final_params['unitstr'] = "[dB]"
		final_params['cbtickinterval'] = 0.1

	if name == "reflectivity":
		final_params = product_parameters_VR.copy()
		final_params['cmap'] = mpl.colormaps["pyart_ChaseSpectral"]
		final_params['vmin'], final_params['vmax'] = -20, 80
		final_params['productstr'] = "Reflecitivity"
		final_params['unitstr'] = "[dBZ]"
		final_params['cbtickinterval'] = 5

	return final_params

if __name__ == "__main__":
	import argparse, sys
	import configparser

	config = configparser.ConfigParser()
	config.read("./NFGDA.ini")
	case_name = config["Settings"]["case_name"]
	downdir = '../V06/'+case_name+'/'
	print(downdir)
	radarobjs = loadnexrad(downdir=downdir)
	#print(radarobjs)
	for radarobj in radarobjs[:config.getint('Settings', 'i_end')]:
		print(radarobj)
		plotradarobj(radarobj, 0,0, savepath=radarobj)
