from numpy import *
from datetime import datetime, timedelta
import netCDF4
import dataloader_CMIP6 as dl
#model = "MIROC6.piControl.r1i1p1f1"
#model = "MRI-ESM2-0.piControl.r1i1p1f1"
#model = "MRI-ESM2-0.piControl.r1i1p1f1"
model = "MRI-ESM2-0.historical.r1i1p1f1"
#model = "MPI-ESM1-2-HR.piControl.r1i1p1f1"

#var   = "va"
#plev  = 850

#dtime = datetime(3376, 2, 1, 0)   #359701010600
dtime = datetime(1980, 2, 1, 0)   #359701010600


#a=dl.Load_6hrPlev(model, var, dtime, plev)
#print(a.shape)
#
var   = "slp"
plev  = 850
b=dl.Load_6hrSfc(model, var, dtime)
print(b)
#
#var   = "topo"
#c=dl.Load_const(model, var)
#print(c)

#var = "sst"
##Year, Mon = 2550,12
#Year, Mon = 1950,1
#a=dl.Load_monSfc(model, var, Year,Mon)
#print(a)
#
#
#print(dl.ret_lats(model))
#print(dl.ret_lons(model))
