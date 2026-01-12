import asyncio
from concurrent.futures import ProcessPoolExecutor
from transitions.extensions import AsyncMachine
import sys
import datetime
import os
import numpy as np
import NF_Lib
from NF_Lib import tprint, C
from NFGDA_load_config import *
import traceback

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
    def __init__(self,
                dl_workers=4,
                nfgda_workers=4,
                df_workers=4,
                dl_qsize=20,
                nfgda_qsize=20,
                df_qsize = 20):
        self.running = True
        self.pull_seconds = 120
        self.last_nexrad = datetime.datetime.now()-datetime.timedelta(minutes=90)
        self.path_config = path_config

        self.nexrad_buf_size = 20
        self.live_nexrad = np.full((self.nexrad_buf_size),'', dtype=object)
        self.nfgda_ready = [asyncio.Event() for _ in range(self.nexrad_buf_size)]
        self.df_ready = [asyncio.Event() for _ in range(self.nexrad_buf_size)]
        self.cur_nex_idx = 0

        # pipeline queues
        self.download_q   = asyncio.Queue(maxsize=dl_qsize)
        self.nfgda_q      = asyncio.Queue(maxsize=nfgda_qsize)
        self.d_forecast_q = asyncio.Queue(maxsize=df_qsize)
        self.result_q     = asyncio.Queue()

        # executors
        self.dl_pool = ProcessPoolExecutor(max_workers=dl_workers)
        self.ng_pool = ProcessPoolExecutor(max_workers=nfgda_workers)
        self.df_pool = ProcessPoolExecutor(max_workers=df_workers)

        # concurrency caps (usually match executor workers)
        self.dl_sem = asyncio.Semaphore(dl_workers)
        self.ng_sem = asyncio.Semaphore(nfgda_workers)
        self.df_sem = asyncio.Semaphore(df_workers)

        # tasks list for shutdown
        self._tasks = []

    async def check_update(self):
        tprint(ht_tag+f"Checking for [{radar_id}] updates... latest nexrad =",self.last_nexrad - datetime.timedelta(seconds=1))
        try:
            scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad, datetime.datetime.now(), radar_id)
            # scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad - datetime.timedelta(days=1), datetime.datetime.now()- datetime.timedelta(days=1), radar_id)
        except TypeError:
            tprint(
                ht_tag +
                f"aws_int TypeError: failed to retrieve NEXRAD scans "
                f"(possible radar outage, AWS latency, or invalid time window).{C.RESET}"
            )
            return
        if len(scans)>0:
            self.last_nexrad = scans[-1].scan_time + datetime.timedelta(seconds=1)
            tprint(dl_tag+
                f"Find {len(scans)} volumes.")
            self.new_nex = scans

            for vol in scans:
                if vol.filename[-4:]=='_MDM':
                    tprint(dl_tag+
                        f"MDM! Skip: {vol.filename}")
                    continue
                self.live_nexrad[self.cur_nex_idx] = vol.filename
                await self.download_q.put((vol,self.cur_nex_idx))
                self.cur_nex_idx = (self.cur_nex_idx + 1) % self.live_nexrad.size

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
                    tprint(dl_tag+
                        f'nfgda_ready [{idx}] set')
                    if self.live_nexrad[(idx - 1) % self.live_nexrad.size] != '':
                        # tprint(f'self.live_nexrad[({idx} - 1)]:{self.live_nexrad[(idx - 1) % self.live_nexrad.size]} \
                        #     nfgda_q.put [{idx}]')
                        await self.nfgda_q.put((idx))
                except:
                    traceback.print_exc()
                    tprint(dl_tag+
                        f'{C.RED_B}Fatal Error.{C.RESET}')
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
                tprint(ng_tag+f'{self.live_nexrad[idx].strip()} [{idx}] wait nfgda_ready[{pre_idx}]')
                await self.nfgda_ready[pre_idx].wait()
                try:
                    async with self.ng_sem:
                        await loop.run_in_executor(
                            self.ng_pool,
                            NF_Lib.nfgda_unit_step,
                            self.live_nexrad[pre_idx],
                            self.live_nexrad[idx]
                        )
                    tprint(ng_tag+f'{self.live_nexrad[idx].strip()}[{idx}] nfgda_ready[{pre_idx}] clear; df_ready[{idx}] set')
                    self.nfgda_ready[pre_idx].clear()
                    self.df_ready[idx].set()
                    await self.d_forecast_q.put((idx))
                except:
                    traceback.print_exc()
                    tprint(ng_tag+f'{C.RED_B}Fatal Error.{C.RESET}')
                finally:
                    self.nfgda_q.task_done()
        except asyncio.CancelledError:
            raise

    async def d_forecast_worker(self):
        loop = asyncio.get_running_loop()
        try:
            while True:
                idx = await self.d_forecast_q.get()
                next_idx = (idx + 1) % self.live_nexrad.size
                tprint(df_tag+f'{self.live_nexrad[idx].strip()}[{idx}] wait df_ready[{next_idx}]')
                await self.df_ready[next_idx].wait()
                try:
                    async with self.ng_sem:
                        await loop.run_in_executor(
                            self.df_pool,
                            NF_Lib.nfgda_forecast,
                            self.live_nexrad[idx],
                            self.live_nexrad[next_idx]
                        )
                    tprint(df_tag+
                        f'df_ready[{idx}] clear.')
                    self.df_ready[idx].clear()

                except:
                    traceback.print_exc()
                    tprint(df_tag+
                        f'{C.RED_B}Fatal Error.{C.RESET}')
                finally:
                    self.d_forecast_q.task_done()
        except asyncio.CancelledError:
            raise

    async def run(self):
        self._tasks = [
            asyncio.create_task(self.download_worker(), name="download_worker"),
            asyncio.create_task(self.nfgda_worker(),   name="nfgda_worker"),
            asyncio.create_task(self.d_forecast_worker(),   name="d_forecast_worker"),
        ]

        try:
            while self.running:
                await self.check_update()
                await counter_loop(self.pull_seconds)
        finally:
            await self.shutdown()

    async def shutdown(self):
        self.running = False

        for t in self._tasks:
            t.cancel()
        await asyncio.gather(*self._tasks, return_exceptions=True)

        self.dl_pool.shutdown(wait=False)
        self.ng_pool.shutdown(wait=False)
        self.df_pool.shutdown(wait=False)

if __name__ == "__main__":
    daemon_host = HostDaemon()
    asyncio.run(daemon_host.run())
