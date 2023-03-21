import netCDF4
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import point
from shapely import wkt
import io
import os
import datetime
# get this here: pip -m install git+https://github.com/lkschn/WQ_tools
from WQ_tools.nwdmFunctions import wfsbuild, readUrl
from WQ_tools.dwaqFunctions import get_modkey, get_modTime
from WQ_tools.plotFunctions import plotTS_modelNWDM

from NWDM_funcs import wfs2csv

outpath = '/home/daniel/levante_work/IR/Harmonization/NWDM/datafiles/'

stations = []
yint=[2014,2017]
vars = ['Chl','DIN','DIP','salt']
# olf = 4.0
wfs2csv(outpath,stations, yint, vars)