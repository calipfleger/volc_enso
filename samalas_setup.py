"""Setup utilities and constants for Samalas eruption analysis.

Derived from the initial cells of ``V1_1x_Sam.ipynb``.  These
parameters centralize common configuration such as file paths,
plotting options, region definitions and helper routines.
"""

from datetime import datetime
from itertools import product

import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
import cftime
from cftime import DatetimeNoLeap
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from collections import OrderedDict, namedtuple
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from shapely.geometry import LinearRing, Polygon
from netCDF4 import Dataset as nc
from netCDF4 import num2date

# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

def get_current_datetime():
    """Return strings for today's date and current time."""
    now = datetime.now()
    return now.strftime("%D"), now.strftime("%H:%M:%S")

# -----------------------------------------------------------------------------
# File paths and naming conventions
# -----------------------------------------------------------------------------

FILEPATH = '/glade/derecho/scratch/calipfleger/ilme/'
BASEDIR = '/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/'
BASEDIR_OCN = '/glade/derecho/scratch/samantha/iLME_seasvolc/CESM-CAM5-LME/atm/proc/tseries/monthly/'
FILE_PREFIX = '.cam.h0.'
I_MODEL = '/b.ie12.B1850C5CN.f19_g16.LME.'
I_MODEL1 = '/b.ie12.B1850CN.f19_g16.'
MODEL = '/b.e11.BLMTRC5CN.f19_g16.'

# -----------------------------------------------------------------------------
# Plotting parameters
# -----------------------------------------------------------------------------

VARPLOT = '\n 1258 Samalas Eruption_1x 10 Members'
VARIABLE = 'Surface Temperature Relative to 30 Year Climatology'
UNIT = 'C'
PLT_YLABEL = 'Surface Temperature (C)'
VARI = 'TS'
MONTHS = ['-5','-4','-3','-2', '-1','0', '1', '2', '3', '4', '5']
COLORS = ['blue', 'black', 'red' ,'green']
COLORS_RIBBON = ['lightblue' ,'lightgrey','mistyrose', 'lightgreen']
VAR_NAME = ['January_1x',  'April_1x', 'July_1x', 'October_1x']
VAR_NAMES = ['DJF', 'MAM', 'JJA', 'SON']

# -----------------------------------------------------------------------------
# Region and geometry definitions
# -----------------------------------------------------------------------------

REGBOX = [-5, 5, 190, 240]  # Nino 3.4 region of interest
WLATS = [-5, -5, 5, 5]
WLONS = [-240, -190, -190, -240]
WRING = LinearRing(list(zip(WLONS, WLATS)))

# -----------------------------------------------------------------------------
# CESM time handling helpers
# -----------------------------------------------------------------------------

dates_sam = [DatetimeNoLeap(year, month, 1) for year, month in product(range(1250, 1265), range(1, 13))]
da_sam = xr.DataArray(np.arange(180), coords=[dates_sam], dims=['time'], name='time')

years_cntl = np.linspace(0,348, 360)
dates_cntl = [DatetimeNoLeap(years_cntl, month, 1) for years_cntl, month in product(range(1228, 1258), range(1, 13))]
da_cntl = xr.DataArray(np.arange(360), coords=[dates_cntl], dims=['time'], name='time')

# -----------------------------------------------------------------------------
# NetCDF utility
# -----------------------------------------------------------------------------

def save_ensemble(x, var_name, filename):
    """Concatenate a list of xarray objects and save as NetCDF.

    Parameters
    ----------
    x : list of xr.DataArray or xr.Dataset
        Individual ensemble member objects.
    var_name : str
        Variable name used in the output path.
    filename : str
        Suffix for the output file name (without extension).
    """
    xrda = xr.concat(x, dim='member')
    xrda.to_netcdf(f"{FILEPATH}{var_name}{filename}.nc")
    return xrda
