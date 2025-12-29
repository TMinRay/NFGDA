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


    async def run(self):
        while True:
            await self.check()
            await counter_loop(120)

class HostDaemon:
    states = ['idle', 'checking', 'downloading']
    def __init__(self,
                dl_workers=2,
                nfgda_workers=4,
                dl_qsize=20,
                nfgda_qsize=20):
        # self.machine = AsyncMachine(model=self, states=self.states, initial='idle')
        # self.machine.add_transition('check', 'idle', 'checking', after='check_update')
        # self.machine.add_transition('download', 'checking', 'downloading', after='handle_download')
        # self.machine.add_transition('reset', '*', 'idle')
        self.running = True
        self.poll_seconds = 120
        # self.executor = ProcessPoolExecutor(max_workers=2)
        self.last_nexrad = datetime.datetime.now()-datetime.timedelta(minutes=30)
        self.path_config = path_config

        self.nexrad_buf_size = 20
        self.live_nexrad = np.full((self.nexrad_buf_size),'', dtype=object)
        # self.nfgda_ready = np.full((20),False, dtype=bool)
        self.nfgda_ready = [asyncio.Event() for _ in range(self.nexrad_buf_size)]
        self.cur_nex_idx = 0

        # pipeline queues
        self.download_q = asyncio.Queue(maxsize=dl_qsize)
        self.nfgda_q    = asyncio.Queue(maxsize=nfgda_qsize)
        self.result_q   = asyncio.Queue()

        # executors
        self.dl_pool = ProcessPoolExecutor(max_workers=dl_workers)
        self.ng_pool = ProcessPoolExecutor(max_workers=nfgda_workers)

        # concurrency caps (usually match executor workers)
        self.dl_sem = asyncio.Semaphore(dl_workers)
        self.ng_sem = asyncio.Semaphore(nfgda_workers)

        # tasks list for shutdown
        self._tasks = []

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
            # await self.reset()
            return
        if len(scans)>0:
            self.last_nexrad = scans[-1].scan_time + datetime.timedelta(seconds=1)
            print(f"Find {len(scans)} volumes")
            self.new_nex = scans

            for vol in scans:
                if vol.filename[-4:]=='_MDM':
                    tprint(f"[Downloader] MDM! Skip: {vol.filename}")
                    continue
                self.live_nexrad[self.cur_nex_idx] = vol.filename
                await self.download_q.put((vol,self.cur_nex_idx))
                self.cur_nex_idx = (self.cur_nex_idx + 1) % self.live_nexrad.size
            # for nxd in scans:
            #     print('  ',nxd.filename,nxd.scan_time)
        #     await self.download()
        #     return
        # else:
        #     await self.reset()
        #     return

                    #     self.cur_nfgda_idx = 0
                    # async def check_update(self):
                    #     tprint(f"[NFGDA] Checking readiness of nfgda_idx={self.cur_nfgda_idx}")
                    #     if self.host.nfgda_ready[self.cur_nfgda_idx]:
                    #         await self.detect()
                    #     else:
                    #         await self.reset()

    # async def handle_download(self):
    #     tprint("[Daemon] Submitting job to worker...")
    #     loop = asyncio.get_running_loop()
    #     # Submit heavy job to another process
    #     # for vol in self.new_nex:
    #     #     result = loop.run_in_executor(self.executor, self.get_nexrad, vol)
    #     futures=[]
    #     for vol in self.new_nex:
    #         futures.append(loop.run_in_executor(self.executor, NF_Lib.get_nexrad, self.path_config, vol))
    #     # Wait for completion; results will be [None, None, ...]
    #     await asyncio.gather(*futures)
    #     for vol in self.new_nex:
    #         self.cur_nex_idx = (self.cur_nex_idx + 1) % self.live_nexrad.size
    #         self.live_nexrad[self.cur_nex_idx] = vol.filename
    #         self.nfgda_ready[self.cur_nex_idx].set()
    #     await self.reset()

    async def download_worker(self):
        loop = asyncio.get_running_loop()
        try:
            while True:
                vol,idx = await self.download_q.get()
                try:
                    async with self.dl_sem:
                        await loop.run_in_executor(
                            self.dl_pool,
                            NF_Lib.get_nexrad,
                            self.path_config,
                            vol
                        )
                    self.nfgda_ready[idx].set()
                    if self.live_nexrad[(idx - 1) % self.live_nexrad.size] != '':
                        await self.nfgda_q.put((idx))
                finally:
                    self.download_q.task_done()
        except asyncio.CancelledError:
            raise

    async def nfgda_worker(self):
        loop = asyncio.get_running_loop()
        try:
            while True:
                idx = await self.nfgda_q.get()
                pre_idx = (idx - 1) % self.live_nexrad.size
                await self.nfgda_ready[pre_idx].wait()
                tprint
                try:
                    async with self.ng_sem:
                        await loop.run_in_executor(
                            self.ng_pool,
                            NF_Lib.nfgda_unit_step,
                            self.live_nexrad[pre_idx],
                            self.live_nexrad[idx]
                        )
                    self.nfgda_ready[pre_idx].clear()
                finally:
                    self.nfgda_q.task_done()
        except asyncio.CancelledError:
            raise

    # async def handle_detect(self):
    #     tprint("[NFGDA] Submitting job to worker...")
    #     if self.host.nfgda_ready[(self.cur_nfgda_idx-1)%self.host.nfgda_ready.size] != '':
    #     loop = asyncio.get_running_loop()
    #     futures=[]
    #     while self.host.nfgda_ready[self.cur_nfgda_idx]:
    #         futures.append(loop.run_in_executor(self.executor, \
    #             NF_Lib.nfgda_unit_step, self.host.live_nexrad[(self.cur_nfgda_idx-1)%self.host.nfgda_ready.size], self.host.live_nexrad[self.cur_nfgda_idx],""))
    #         self.host.nfgda_ready[self.cur_nfgda_idx] = False
    #         self.cur_nfgda_idx = (self.cur_nfgda_idx+1)%self.host.nfgda_ready.size
    #     await asyncio.gather(*futures)
    #     await self.reset()

    async def run(self):

        self._tasks = [
            asyncio.create_task(self.download_worker(), name="download_worker"),
            asyncio.create_task(self.nfgda_worker(),   name="nfgda_worker"),
            # asyncio.create_task(self.status_worker(),  name="status_worker"),
        ]

        try:
            while self.running:
                await self.check_update()                 # assigns work
                await counter_loop(self.poll_seconds)
        finally:
            await self.shutdown()

        async def shutdown(self):
            self.running = False

            for t in self._tasks:
                t.cancel()
            await asyncio.gather(*self._tasks, return_exceptions=True)

            self.dl_pool.shutdown(wait=False)
            self.ng_pool.shutdown(wait=False)
        # try:
        #     while self.running:
        #         await self.main_loop()
        # finally:
        #     await self.shutdown()

        # while True:
        #     await self.check()
        #     await counter_loop(120)

if __name__ == "__main__":
    daemon_host = HostDaemon()
    asyncio.run(daemon_host.run())
