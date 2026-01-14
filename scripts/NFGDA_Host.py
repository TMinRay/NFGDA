import asyncio
from concurrent.futures import ProcessPoolExecutor
import sys
import datetime
import os
import numpy as np
import nfgda.NF_Lib as NF_Lib
from nfgda.NFGDA_load_config import *
import traceback
from functools import partial

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
                cstart = None,
                cend = None,
                dl_workers=4,
                nfgda_workers=4,
                df_workers=8,
                sf_workers=8):
        self.running = True
        self.pull_seconds = 120
        if cstart is None:
            self.last_nexrad = datetime.datetime.now(datetime.timezone.utc)-datetime.timedelta(minutes=90)
            self.exit_time = datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(days=360)
            self.real_time_mode = True
        else:
            self.last_nexrad = cstart-datetime.timedelta(minutes=1)
            self.exit_time = cend
            self.real_time_mode = False

        self.path_config = path_config

        self.nexrad_buf_size = 30
        self.live_nexrad = np.full((self.nexrad_buf_size),'', dtype=object)
        self.nfgda_ready = [asyncio.Event() for _ in range(self.nexrad_buf_size)]
        self.df_ready = [asyncio.Event() for _ in range(self.nexrad_buf_size)]
        self.cur_nex_idx = 0
        self.live_forecasts = [() for _ in range(self.nexrad_buf_size)]
        # pipeline queues
        self.download_q   = asyncio.Queue(maxsize=self.nexrad_buf_size)
        self.nfgda_q      = asyncio.Queue(maxsize=self.nexrad_buf_size)
        self.d_forecast_q = asyncio.Queue(maxsize=self.nexrad_buf_size)
        self.s_forecast_q = asyncio.Queue(maxsize=self.nexrad_buf_size)

        # executors
        self.dl_pool = ProcessPoolExecutor(max_workers=dl_workers)
        self.ng_pool = ProcessPoolExecutor(max_workers=nfgda_workers)
        self.df_pool = ProcessPoolExecutor(max_workers=df_workers)
        self.sf_pool = ProcessPoolExecutor(max_workers=sf_workers)

        # concurrency caps (usually match executor workers)
        self.dl_sem = asyncio.Semaphore(dl_workers)
        self.ng_sem = asyncio.Semaphore(nfgda_workers)
        self.df_sem = asyncio.Semaphore(df_workers)
        self.sf_sem = asyncio.Semaphore(sf_workers)

        # tasks list for shutdown
        self._tasks = []

    async def check_update(self):
        tprint(ht_tag+f"Checking for [{radar_id}] updates... latest nexrad =",self.last_nexrad - datetime.timedelta(seconds=1))
        try:
            if self.real_time_mode:
                scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad, datetime.datetime.now(datetime.timezone.utc), radar_id)
            else:
                scans = NF_Lib.aws_int.get_avail_scans_in_range(self.last_nexrad, self.last_nexrad+datetime.timedelta(minutes=20), radar_id)
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
            tprint(ht_tag,self.last_nexrad, self.exit_time,self.last_nexrad > self.exit_time)
            if self.last_nexrad > self.exit_time:
                await self.delay_shutdown()

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
                if (next_idx == self.cur_nex_idx) and not(self.real_time_mode):
                    tprint(df_tag+f'Last Detection {self.live_nexrad[idx].strip()}[{idx}] close forecast worker.')
                    self.d_forecast_q.task_done()
                    continue
                await self.df_ready[next_idx].wait()
                try:
                    suppress = self.d_forecast_q.qsize()>3 and self.real_time_mode
                    async with self.ng_sem:
                        results = await loop.run_in_executor(
                            self.df_pool,
                            partial(
                                NF_Lib.nfgda_forecast,
                                self.live_nexrad[idx],
                                self.live_nexrad[next_idx],
                                suppress_fig = suppress,
                            )
                        )
                    tprint(df_tag+
                        f'df_ready[{idx}] clear.')
                    self.df_ready[idx].clear()
                    self.live_forecasts[idx] = results
                    if not suppress:
                        await self.s_forecast_q.put(idx)
                except:
                    traceback.print_exc()
                    tprint(df_tag+
                        f'{C.RED_B}Fatal Error.{C.RESET}')
                finally:
                    self.d_forecast_q.task_done()
        except asyncio.CancelledError:
            raise

    async def s_forecast_worker(self):
        loop = asyncio.get_running_loop()
        try:
            while True:
                idx = await self.s_forecast_q.get()
                # next_idx = (idx + 1) % self.live_nexrad.size
                # tprint(df_tag+f'{self.live_nexrad[idx].strip()}[{idx}] wait df_ready[{next_idx}]')
                # await self.df_ready[next_idx].wait()
                try:
                    async with self.ng_sem:
                        results = await loop.run_in_executor(
                            self.df_pool,
                            NF_Lib.nfgda_stochastic_summary,
                            self.live_forecasts,
                            self.live_nexrad[idx]
                        )
                except:
                    traceback.print_exc()
                    tprint(df_tag+
                        f'{C.RED_B}Fatal Error.{C.RESET}')
                finally:
                    self.s_forecast_q.task_done()
        except asyncio.CancelledError:
            raise

    async def run(self):
        self._tasks = [
            asyncio.create_task(self.download_worker(), name="download_worker"),
            asyncio.create_task(self.nfgda_worker(),   name="nfgda_worker"),
            asyncio.create_task(self.d_forecast_worker(),   name="d_forecast_worker"),
            asyncio.create_task(self.s_forecast_worker(),   name="s_forecast_worker"),
        ]

        try:
            while self.running:
                await self.check_update()
                if self.real_time_mode:
                    await counter_loop(self.pull_seconds)
                else:
                    await counter_loop(5)
                    await self.wait_until_qsize([
                        self.download_q,
                        self.nfgda_q,
                        self.d_forecast_q,
                        self.s_forecast_q],
                        self.nexrad_buf_size//4)
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

    async def delay_shutdown(self, timeout=3600):
        self.running = False
        try:
            await asyncio.wait_for(
                asyncio.gather(
                    self.download_q.join(),
                    self.nfgda_q.join(),
                    self.d_forecast_q.join(),
                ),
                timeout,
            )
        except asyncio.TimeoutError:
            tprint(ht_tag+
                "Shutdown timed out:\n"
                f"Queues not drained: download={self.download_q.qsize()} nfgda={self.nfgda_q.qsize()}"
                f"forecast={self.d_forecast_q.qsize()}")
        await self.shutdown()

    async def wait_until_qsize(self, qs, max_remaining, timeout=None):
        async def _wait():
            while max(q.qsize() for q in qs) > max_remaining:
                await asyncio.sleep(1)
        if timeout is None:
            await _wait()
        else:
            await asyncio.wait_for(_wait(), timeout)

if __name__ == "__main__":
    daemon_host = HostDaemon(custom_start_time,custom_end_time)
    asyncio.run(daemon_host.run())
