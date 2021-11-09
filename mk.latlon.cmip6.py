import numpy as np
import netCDF4
from datetime import datetime, timedelta
from glob import glob
import os, sys
"""
    This program is used to read input data.
"""
#******************************************
# Edit here (input file directories)
#------------------------------------------
#lmodel  = ["MIROC6."]
#lmodel = ["MIROC6.piControl.r1i1p1f1"]
#lmodel = ["MRI-ESM2-0.piControl.r1i1p1f1"]
lmodel = ["MRI-ESM2-0.historical.r1i1p1f1"]
#lmodel = ["MPI-ESM1-2-HR.piControl.r1i1p1f1"]

for model in lmodel:
    vname  = "orog"
    srcdir = "/mnt/nas02/data/CMIP6/%s"%(model)
    ssearch = srcdir + "/%s_fx*.nc"%(vname)
    lsrcpath = glob(ssearch)
    srcpath  = lsrcpath[0]
    nc = netCDF4.Dataset(srcpath)
    
    a1lat = nc.variables["lat"][:].data
    a1lon = nc.variables["lon"][:].data
    
    latpath = srcdir + "/lat.npy"
    lonpath = srcdir + "/lon.npy"
    np.save(latpath, a1lat)
    np.save(lonpath, a1lon)
    
    print(a1lat)
    print(a1lon)
    print(len(a1lat), len(a1lon))
    print(latpath) 
    
