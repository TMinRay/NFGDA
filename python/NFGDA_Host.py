import asyncio
from concurrent.futures import ProcessPoolExecutor
from transitions.extensions import AsyncMachine
import sys
import datetime
import os
import numpy as np
import NF_Lib
from NF_Lib import tprint
from NFGDA_load_config import *

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

class NFGDADaemon:
    states = ['idle', 'checking', 'detecting']
    def __init__(self,host):
        self.host = host
        self.machine = AsyncMachine(model=self, states=self.states, initial='idle')
        self.machine.add_transition('check', 'idle', 'checking', after='check_update')
        self.machine.add_transition('detect', 'checking', 'detecting', after='handle_detect')
        self.executor = ProcessPoolExecutor(max_workers=4)
        self.cur_nfgda_idx = 0
    async def check_update(self):
        tprint(f"[NFGDA] Checking readiness of nfgda_idx={self.cur_nfgda_idx}")
        if self.host.nfgda_ready[self.cur_nfgda_idx]:
            await self.detect()
        else:
            await self.reset()

    async def handle_detect(self):
        tprint("[NFGDA] Submitting job to worker...")
        if self.host.nfgda_ready[(self.cur_nfgda_idx-1)%self.host.nfgda_ready.size] != '':
        loop = asyncio.get_running_loop()
        futures=[]
        while self.host.nfgda_ready[self.cur_nfgda_idx]:
            futures.append(loop.run_in_executor(self.executor, \
                NF_Lib.nfgda_unit_step, self.host.live_nexrad[(self.cur_nfgda_idx-1)%self.host.nfgda_ready.size], self.host.live_nexrad[self.cur_nfgda_idx],""))
            self.host.nfgda_ready[self.cur_nfgda_idx] = False
            self.cur_nfgda_idx = (self.cur_nfgda_idx+1)%self.host.nfgda_ready.size
        await asyncio.gather(*futures)
        await self.reset()

    async def run(self):
        while True:
            await self.check()
            await counter_loop(120)

class HostDaemon:
    states = ['idle', 'checking', 'downloading']
    def __init__(self):
        self.machine = AsyncMachine(model=self, states=self.states, initial='idle')
        self.machine.add_transition('check', 'idle', 'checking', after='check_update')
        self.machine.add_transition('download', 'checking', 'downloading', after='handle_download')
        self.machine.add_transition('reset', '*', 'idle')
        self.executor = ProcessPoolExecutor(max_workers=2)
        self.last_nexrad = datetime.datetime.now()-datetime.timedelta(minutes=30)
        self.path_config = path_config
        self.live_nexrad = np.full((20),'', dtype=object)
        self.nfgda_ready = np.full((20),False, dtype=bool)
        self.cur_nex_idx = 0

    async def check_update(self):
        tprint(f"[Host] Checking for [{radar_id}] updates... last =",self.last_nexrad - datetime.timedelta(seconds=1))
        try:
            scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad, datetime.datetime.now(), radar_id)
            # scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad - datetime.timedelta(days=1), datetime.datetime.now()- datetime.timedelta(days=1), radar_id)
        except TypeError:
            tprint(
                "[Host] aws_int TypeError: failed to retrieve NEXRAD scans "
                "(possible radar outage, AWS latency, or invalid time window)."
            )
            await self.reset()
            return
        if len(scans)>0:
            self.last_nexrad = scans[-1].scan_time + datetime.timedelta(seconds=1)
            print(f"Find {len(scans)} volumes")
            self.new_nex = scans
            # for nxd in scans:
            #     print('  ',nxd.filename,nxd.scan_time)
            await self.download()
            return
        else:
            await self.reset()
            return

    async def handle_download(self):
        tprint("[Daemon] Submitting job to worker...")
        loop = asyncio.get_running_loop()
        # Submit heavy job to another process
        # for vol in self.new_nex:
        #     result = loop.run_in_executor(self.executor, self.get_nexrad, vol)
        futures=[]
        for vol in self.new_nex:
            futures.append(loop.run_in_executor(self.executor, NF_Lib.get_nexrad, self.path_config, vol))
        # Wait for completion; results will be [None, None, ...]
        await asyncio.gather(*futures)
        for vol in self.new_nex:
            self.cur_nex_idx = (self.cur_nex_idx + 1) % self.live_nexrad.size
            self.live_nexrad[self.cur_nex_idx] = vol.filename
            self.nfgda_ready[self.cur_nex_idx] = True
        await self.reset()

    async def run(self):
        while True:
            await self.check()
            await counter_loop(120)

if __name__ == "__main__":
    daemon_host = HostDaemon()
    daemon_nfgda = NFGDADaemon(daemon_host)
    asyncio.run(daemon_host.run())
    asyncio.run(daemon_nfgda.run())
