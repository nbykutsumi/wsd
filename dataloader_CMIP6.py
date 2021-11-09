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
slpbasedir = "/mnt/nas02/data/CMIP6"
tabasedir  = "/mnt/nas02/data/CMIP6"
uabasedir  = "/mnt/nas02/data/CMIP6"
vabasedir  = "/mnt/nas02/data/CMIP6"
tsbasedir  = "/mnt/nas02/data/CMIP6"
topobasedir= "/mnt/nas02/data/CMIP6"
landbasedir= "/mnt/nas02/data/CMIP6"
prbasedir  = "/mnt/nas02/data/CMIP6"   # not used for detection


#******************************************

dbasedir = {
    "slp":slpbasedir,
    "ta" :tabasedir,
    "ua" :uabasedir,
    "va" :vabasedir,
    "sst":tsbasedir,
    "topo":topobasedir,
    "land":landbasedir,
    "pr" :prbasedir,     # not used for detection
    }

dvar = {
    "slp":"psl",
    "ta" :"ta",
    "ua" :"ua",
    "va" :"va",
    "sst":"ts",
    "topo":"orog",
    "land":"sftlf",
    "pr" :"pr",    # not used for detection
    }

def ret_lats(model):
    return np.load(slpbasedir + "/%s/lat.npy"%(model))
    # MIROC6: -88.92773535 ~ 88.92773535, d=~1.4007664

def ret_lons(model):
    return np.load(slpbasedir + "/%s/lon.npy"%(model))
    # MIROC6: 0 ~ 358.59375, d=1.40625

def ret_ny(model):
    return len(ret_lats(model))

def ret_nx(model):
    return len(ret_lons(model))
    # MIROC6:   (128, 256)
def ret_miss(model):
    modelname = model.split(".")[0]
    if   modelname=="MIROC6":       miss_in= 9.969209968386869e+36
    elif modelname=="MRI-ESM2-0":   miss_in= 9.969209968386869e+36
    elif modelname=="MPI-ESM1-2-HR":miss_in= 9.969209968386869e+36

    return miss_in 

def Load_6hrPlev(model, var, DTime, plev):
    modelname, expr, ens = model.split(".")
    vname = dvar[var]
    iplev = [850, 500, 250].index(plev)
    # Search file
    srcdir  = dbasedir[var] + "/%s"%(model)
    ssearch = srcdir + "/%s_6hrPlev*.nc"%(vname)
    lsrcpath = glob(ssearch)
    for srcpath in lsrcpath:
        stime = os.path.basename(srcpath).split("_")[6].split(".")[0]
        stime0, stime1 = stime.split("-")

        dtime0 = datetime.strptime(stime0, "%Y%m%d%H%M")
        dtime1 = datetime.strptime(stime1, "%Y%m%d%H%M")
        if (dtime0<=DTime)&(DTime<=dtime1):
            break
    nc = netCDF4.Dataset(srcpath)

    #print(nc.variables)
    #print(srcpath)
    # Find time index
    basetime = {
        ("MIROC6","piControl"):         datetime(3200,1,1),

        ("MRI-ESM2-0","piControl"):     datetime(1850,1,1),
        ("MRI-ESM2-0","historical"):    datetime(1850,1,1),

        ("MPI-ESM1-2-HR","piControl"):  datetime(1850,1,1),
        }[modelname,expr]

    dtime0 = basetime + timedelta(days=float(nc.variables["time"][0]))
    idxtime = int((DTime - dtime0).total_seconds()/21600)    # 6-hour = 21600 sec
    #print(DTime, dtime0)
    #print(idxtime)
    return nc.variables[vname][idxtime, iplev]


def Load_6hrSfc(model, var, DTime):
    modelname, expr, ens = model.split(".")
    vname = dvar[var]

    # Search file
    srcdir  = dbasedir[var] + "/%s"%(model)
    ssearch = srcdir + "/%s_6hrPlev*.nc"%(vname)
    lsrcpath = np.sort(glob(ssearch))
    for srcpath in lsrcpath:
        stime = os.path.basename(srcpath).split("_")[6].split(".")[0]
        stime0, stime1 = stime.split("-")

        dtime0 = datetime.strptime(stime0, "%Y%m%d%H%M")
        dtime1 = datetime.strptime(stime1, "%Y%m%d%H%M")
        if (dtime0<=DTime)&(DTime<=dtime1):
            break
    nc = netCDF4.Dataset(srcpath)
    #print(nc.variables)
    #print(srcpath)

    # Find time index
    basetime = {
        ("MIROC6","piControl"):         datetime(3200,1,1),
        ("MRI-ESM2-0","piControl"):     datetime(1850,1,1),
        ("MRI-ESM2-0","historical"):    datetime(1850,1,1),
        }[modelname,expr]

    dtime0 = basetime + timedelta(days=float(nc.variables["time"][0]))
    idxtime = int((DTime - dtime0).total_seconds()/21600)    # 6-hour = 21600 sec

    return nc.variables[vname][idxtime]
    #return nc.variables[vname].shape

def Load_monSfc(model, var, Year, Mon):
    modelname, expr, ens = model.split(".")
    vname = dvar[var]
    DTime = datetime(Year,Mon,1)
    # Search file
    srcdir  = dbasedir[var] + "/%s"%(model)
    ssearch = srcdir + "/%s_Amon*.nc"%(vname)
    lsrcpath = np.sort(glob(ssearch))
    for srcpath in lsrcpath:
        stime = os.path.basename(srcpath).split("_")[6].split(".")[0]
        stime0, stime1 = stime.split("-")

        dtime0 = datetime.strptime(stime0, "%Y%m")
        dtime1 = datetime.strptime(stime1, "%Y%m")
        if (dtime0<=DTime)&(DTime<=dtime1):
            break
    nc = netCDF4.Dataset(srcpath)
    #print(nc.variables)
    #print(srcpath)

    #print(nc.variables["time"][:])
    #print(len(nc.variables["time"][:]))

    # Find time index
    Year0,Mon0 = dtime0.timetuple()[:2]
    Year1,Mon1 = dtime1.timetuple()[:2]
    idxtime = int(Year-Year0)*12 -Mon0 + Mon
    #print(idxtime)
    return nc.variables[vname][idxtime]



def Load_const(model, var):
    vname = dvar[var]
    srcdir  = dbasedir[var] + "/%s"%(model)
    ssearch = srcdir + "/%s_*.nc"%(vname)
    lsrcpath= glob(ssearch)
    srcpath = lsrcpath[0]
    nc = netCDF4.Dataset(srcpath)
    #print(nc.variables)
    return nc.variables[vname][:]



