import numpy as np
from datetime import datetime, timedelta

"""
    This program is used to read input data.
    Edit this program to set input file directories.

    Naming convention for input binary data files
    * variables at pressure levels
    {input-dir}/yyyymm/{var}_{plev}.yyyymmddhh.bin     # var: slp, ta, ua, va, topo
                                                       # plev: pressure level (hPa) in 4 digits (e.g., 0850, 0500, 0250)

    * surface height
    {input-dir}/topo.bin
"""
#******************************************
# Edit here (input file directories)
#------------------------------------------
slpbasedir = "/mnt/nas02/data/JRA55_GANAL"
tabasedir  = "/mnt/nas02/data/JRA55_GANAL"
uabasedir  = "/mnt/nas02/data/JRA55_GANAL"
vabasedir  = "/mnt/nas02/data/JRA55_GANAL"
topobasedir= "/mnt/nas02/data/JRA55_GANAL"
landbasedir= "/mnt/nas02/data/JRA55_GANAL"   # land sea mask: Sea=0, Land >0
tsbasedir  = "/mnt/nas02/data/JRA55_GANAL"   # SST

lats = np.arange(-90, 90+0.001, 1.25)
lons = np.arange(0, 358.75+0.001, 1.25)
miss= -999.0  # missing value

#******************************************

ny  = len(lats)
nx  = len(lons)

dbasedir = {
    "slp":slpbasedir,
    "ta" :tabasedir,
    "ua" :uabasedir,
    "va" :vabasedir,
    "sst":tsbasedir,
    "topo":topobasedir,
    "land":landbasedir,
    }

def ret_lats(model=None):
    return lats

def ret_lons(model=None):
    return lons

def ret_ny(model=None):
    return ny

def ret_nx(model=None):
    return nx

def ret_miss(model=None):
    return miss

def Load_6hrPlev(model, var, DTime, plev):
    year,mon,day,hour = DTime.timetuple()[:4]
    srcdir  = dbasedir[var] + "/%04d%02d"%(year,mon)
    srcpath = srcdir + "/%s_%04d.%04d%02d%02d%02d.bin"%(var,plev,year,mon,day,hour)


    a = np.fromfile(srcpath, "float32").reshape(ny,nx)
    return np.fromfile(srcpath, "float32").reshape(ny,nx)

def Load_6hrSfc(model, var, DTime):
    year,mon,day,hour = DTime.timetuple()[:4]
    srcdir  = dbasedir[var] + "/%04d%02d"%(year,mon)
    srcpath = srcdir + "/%s.%04d%02d%02d%02d.bin"%(var,year,mon,day,hour)
    return np.fromfile(srcpath, "float32").reshape(ny,nx)

def Load_monSfc(model, var, Year, Mon):
    # Search file
    srcdir  = dbasedir[var] + "/%04d%02d"%(Year,Mon)
    srcpath = srcdir + "/%s.%04d%02d.bin"%(var, Year, Mon)
    print(srcpath)
    return np.fromfile(srcpath, "float32").reshape(ny,nx)

def Load_const(model, var):
    srcdir  = dbasedir[var] 
    srcpath = srcdir + "/%s.bin"%(var)
    return  np.fromfile(srcpath, "float32").reshape(ny,nx)


