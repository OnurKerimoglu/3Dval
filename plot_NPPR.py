import netCDF4
import cartopy


dataroot = '/home/daniel/levante_work2/'

simdict = {'JT3': 'sns144-EH-ERA5-4g-JT3-CS', 'NEU': 'sns144-EH-ERA5-4g-NEU-CS', 'NEC': 'sns144-EH-ERA5-4g-NEC-CS',
           'GCU': 'sns144-EH-ERA5-4g-GCU-CS', 'GCC': 'sns144-EH-ERA5-4g-GCC-CS'}

def main(runlist):
    for ri,rn in enumerate(runlist):
        