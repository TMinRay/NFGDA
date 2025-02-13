import nexradaws
import os
from datetime import datetime
import sys

aws_int = nexradaws.NexradAwsInterface()
radar_id = 'KMAX'
# year, mth, day = 
# start_hr, start_min = 
start = datetime(2023, 8, 15, 22, 0)
end = datetime(2023, 8, 16, 4, 30)

# start = year, mth, day, start_hr, start_min
# end = datetime(year, mth, day, end_hr, end_min)
# end = datetime(year, mth, day+1, end_hr, end_min)

dest_base = r"../../V06"
# dest_base = r"./V06"
dest_folder = '{}{}{:02}{:02}_{:02}'.format(radar_id, year, mth, day, start_hr)
dest = os.path.join(dest_base, dest_folder)
dest = os.path.normpath(dest)

if not os.path.isdir(dest):
    os.makedirs(dest)
else:
    sys.exit('{} already exists.'.format(dest_folder))

scans = aws_int.get_avail_scans_in_range(start, end, radar_id)
aws_int.download(scans, dest)