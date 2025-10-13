import asyncio
from concurrent.futures import ProcessPoolExecutor
from transitions.extensions import AsyncMachine
import time
import nexradaws
import sys
import datetime

radar_id = 'KABX'
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
    states = ['idle', 'checking', 'processing']

    def __init__(self):
        self.machine = AsyncMachine(model=self, states=self.states, initial='idle')
        self.machine.add_transition('check', 'idle', 'checking', after='check_update')
        self.machine.add_transition('process', 'checking', 'processing', after='handle_update')
        self.machine.add_transition('reset', '*', 'idle')
        self.executor = ProcessPoolExecutor(max_workers=2)
        self.last_nexrad = datetime.datetime.now()-datetime.timedelta(minutes=30)

    async def check_update(self):
        tprint(f"[Daemon] Checking for [{radar_id}] updates... last =",self.last_nexrad - datetime.timedelta(seconds=1))
        scans = aws_int.get_avail_scans_in_range(self.last_nexrad, datetime.datetime.now(), radar_id)
        if len(scans)>0:
            self.last_nexrad = scans[-1].scan_time + datetime.timedelta(seconds=1)
        print("Find : ")
        for nxd in scans:
            print('  ',nxd.filename,nxd.scan_time)
        await self.reset()
        # await asyncio.sleep(1)
        # update_found = True
        # if update_found:
        #     await self.process()
        # else:
        #     await self.reset()

    async def handle_update(self):
        tprint("[Daemon] Submitting job to worker...")
        loop = asyncio.get_running_loop()
        # Submit heavy job to another process
        result = await loop.run_in_executor(self.executor, heavy_serial_task, "data_block_1")
        tprint(f"[Daemon] Got result: {result}")
        await self.reset()

    async def run(self):
        while True:
            await self.check()
            await counter_loop(120)

if __name__ == "__main__":
    daemon = HostDaemon()
    asyncio.run(daemon.run())
