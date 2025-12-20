import asyncio
from concurrent.futures import ProcessPoolExecutor
from transitions.extensions import AsyncMachine
import time
import nexradaws
import sys
import datetime
import os

import pyart
import numpy as np

radar_id = 'KABX'
V06_dir = '../V06/runtime/'+radar_id
os.makedirs(V06_dir,exist_ok=True)
nf_dir = V06_dir+'/npz'
os.makedirs(nf_dir,exist_ok=True)

aws_int = nexradaws.NexradAwsInterface()

def tprint(*args, **kwargs):
    print(f"[{datetime.datetime.now():%H:%M:%S}]", *args, **kwargs)

# Simulated heavy processing task
def heavy_serial_task(data):
    tprint(f"[Worker] Start processing {data}")
    time.sleep(3)  # CPU-bound simulation
    return f"Result from {data}"

async def counter_loop(interval=120):
    """
    Counts seconds up to `interval` in the terminal.
    Updates in-place like: (count/120s)
    """
    for count in range(1, interval + 1):
        # print in-place
        print(f"\r({count}/{interval}s)", end="")
        await asyncio.sleep(1)  # non-blocking sleep
    print()  # move to next line after finishing

class HostDaemon:
    states = ['idle', 'checking', 'downloading']

    def __init__(self):
        self.machine = AsyncMachine(model=self, states=self.states, initial='idle')
        self.machine.add_transition('check', 'idle', 'checking', after='check_update')
        self.machine.add_transition('download', 'checking', 'downloading', after='handle_update')
        self.machine.add_transition('reset', '*', 'idle')
        self.executor = ProcessPoolExecutor(max_workers=2)
        self.last_nexrad = datetime.datetime.now()-datetime.timedelta(minutes=30)
        self.V06_dir = V06_dir

    async def check_update(self):
        tprint(f"[Daemon] Checking for [{radar_id}] updates... last =",self.last_nexrad - datetime.timedelta(seconds=1))
        scans = aws_int.get_avail_scans_in_range(self.last_nexrad, datetime.datetime.now(), radar_id)
        if len(scans)>0:
            self.last_nexrad = scans[-1].scan_time + datetime.timedelta(seconds=1)
            print(f"Find {len(scans)} volumes")
            self.new_nex = scans
            # for nxd in scans:
            #     print('  ',nxd.filename,nxd.scan_time)
            await self.download()
        else:
            await self.reset()

    async def handle_update(self):
        tprint("[Daemon] Submitting job to worker...")
        loop = asyncio.get_running_loop()
        # Submit heavy job to another process
        # for vol in self.new_nex:
        #     result = loop.run_in_executor(self.executor, self.get_nexrad, vol)
        futures=[]
        for vol in self.new_nex:
            futures.append(loop.run_in_executor(self.executor, get_nexrad, self.V06_dir, vol))

        # Wait for completion; results will be [None, None, ...]
        await asyncio.gather(*futures)
        await self.reset()

    async def run(self):
        while True:
            await self.check()
            await counter_loop(120)

def get_nexrad(dest,buf):
    fn = buf.filename
    if fn[-4:]=='_MDM':
        tprint(f"[Downloader] MDM! Skip: {fn}")
        return
    if not os.path.exists(os.path.join(dest,fn)):
        aws_int.download(buf, dest)
        tprint(f"[Downloader] Got Volume: {fn}")
        convert_v06_to_nf_input(os.path.join(dest,fn),nf_dir)
    else:
        tprint(f"[Downloader] Already downloaded. Skip: {fn}")


def ReadSliceElevation(radar, slice_idx):
    """ Copied from https://github.com/PreciousJatau47/VAD_correction/blob/master/RadarHCAUtils.py
    :param radar:
    :param slice_idx:
    :return:
    """
    sweep_ind = radar.get_slice(slice_idx)
    radar_el = radar.elevation['data'][sweep_ind]
    return radar_el
def ReadRadarSliceUpdate(radar, slice_idx):
    """ Copied from https://github.com/PreciousJatau47/VAD_correction/blob/master/RadarHCAUtils.py
    :param radar:
    :param slice_idx:
    :return:
    """
    radar_range = radar.range['data'] / 1000  # in km
    sweep_ind = radar.get_slice(slice_idx)
    radar_az_deg = radar.azimuth['data'][sweep_ind]  # in degrees
    radar_el = radar.elevation['data'][sweep_ind]

    ref_shape = radar.fields["reflectivity"]['data'][sweep_ind].shape
    placeholder_matrix = np.full(ref_shape, np.nan, dtype=np.float64)
    placeholder_mask = np.full(ref_shape, False, dtype=bool)

    data_slice = []
    labels_slice = list(radar.fields.keys())
    labels_slice.sort()
    mask_slice = []
    var_mask_slice = []

    for radar_product in labels_slice:
        if np.sum(radar.fields[radar_product]['data'][sweep_ind].mask == False) > 0:
            data_slice.append(radar.fields[radar_product]['data'][sweep_ind])
            mask_slice.append(True)
            var_mask_slice.append(radar.fields[radar_product]['data'][sweep_ind].mask)
        else:
            data_slice.append(placeholder_matrix)
            mask_slice.append(False)
            var_mask_slice.append(placeholder_mask)

    return radar_range, radar_az_deg, radar_el, data_slice.copy(), mask_slice.copy(), labels_slice, var_mask_slice

END_GATE = 400
NUM_AZ = 720
var_2_parrot_idx = {'reflectivity': 0, 'velocity': 1, 'spectrum_width': 2, 'differential_phase': 3,
                'cross_correlation_ratio': 4, 'differential_reflectivity': 5}
def convert_v06_to_nf_input(v06_file, nf_dir):

    l2_file = os.path.basename(v06_file)
    if not (l2_file.endswith('_V06') or l2_file.startswith('._')):
        print("[Converter] Skip: ", l2_file)
        return

    print("[Converter] Processing ", l2_file)

    # Output path.
    output_folder = nf_dir
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder,exist_ok=True)

    nf_input_file = l2_file.split('.')[0]+'.npz'

    py_path = os.path.join(nf_dir, nf_input_file)

    # read l2 data
    radar_obj = pyart.io.read_nexrad_archive(v06_file)

    # TODO(pjatau) erase below.
    nsweeps = radar_obj.nsweeps
    vcp = radar_obj.metadata['vcp_pattern']
    print("[Converter] VCP: ", vcp)

    # VCP 212.
    # slices 0-2-4 contain only dual-pol. super res.
    # slices 1-3-5 contain vel products. super res.
    # slices >= 6 contain all products. normal res.

    # Initialize data cube
    # PARROT = np.ma.full((END_GATE, NUM_AZ, 6), np.nan, dtype=np.float64)
    PARROT = np.ma.array(np.ma.array(np.full((END_GATE, NUM_AZ, 6), np.nan, dtype=np.float64)), mask=np.full((END_GATE, NUM_AZ, 6), True))
    in_parrot = np.full(6,False)

    for slice_idx in range(nsweeps):
        radar_el = ReadSliceElevation( radar_obj, slice_idx)
        scan_el = np.nanmedian(radar_el)
        # if abs(scan_el-target_el)>0.3:
        #     continue
        radar_range, az_sweep_deg, radar_el, data_slice, mask_slice, labels_slice, data_mask_slice = ReadRadarSliceUpdate(
            radar_obj, slice_idx)
        print("Processing elevation {} degrees".format(np.nanmedian(radar_el)))

        i_zero_az = np.argmin(np.abs(az_sweep_deg))
        az_shift = -i_zero_az

        var_idx_slice = {labels_slice[i]: i for i in range(len(labels_slice))}

        for var in var_2_parrot_idx.keys():
            i_var = var_idx_slice[var]
            i_parrot = var_2_parrot_idx[var]

            if not mask_slice[i_var] or in_parrot[i_parrot]:
                continue
            in_parrot[i_parrot] = True

            print("Processing {}. parrot idx {}".format(var, i_parrot))

            curr_data = data_slice[i_var][:, :END_GATE]
            curr_mask = data_mask_slice[i_var][:, :END_GATE]
            # curr_data[curr_mask] = np.nan  # (720, 400)
            curr_data = np.roll(a=curr_data, shift=az_shift, axis=0)
            PARROT[:, :, i_parrot] = curr_data.T
            PARROT[:, :, i_parrot].mask = curr_mask.T
        if np.min(in_parrot):
            timestamp=np.datetime64(pyart.graph.common.generate_radar_time_sweep(radar_obj,slice_idx))
            print("slice idx {} timestamp".format(slice_idx),timestamp)
            break
        print()
    print()
    # scipy.io.savemat(output_path, {"PARROT": PARROT})
    np.savez(py_path, PARROT=PARROT.data,mask=PARROT.mask,timestamp=timestamp)


if __name__ == "__main__":
    daemon = HostDaemon()
    asyncio.run(daemon.run())
