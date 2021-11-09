# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import sys, os
from   numpy import *
import numpy as np
#from   detect_fsub import *
from   datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#import detect_func
import util
import Cyclone
import ConstCyclone
import dataloader as dl
#--------------------------------------
prj     = "GANAL_JRA55" # project name. Used to read Const files
model   = "__"          # model name.   Used to read Const files
run     = 'test'        # run name.     Used to name output directories
wsbaseDir = '/home/utsumi/temp/ws'   # Output weather system directory
wsDir     = wsbaseDir + '/%s'%(run)

#prj     = "CMIP6"
#model   = "MRI-ESM2-0.historical.r1i1p1f1"
#run     = "test"   # {expr}-{scen}-{ens}
#wsbaseDir = '/tank/out/ws/CMIP6/%s'%(model)
#wsDir     = wsbaseDir + '/%s'%(run)


iDTime = datetime(2014,7,1,0)
eDTime = datetime(2014,7,31,18)

lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(hours=6))

a1lat   = dl.ret_lats(model)
a1lon   = dl.ret_lons(model)
ny      = len(a1lat)
nx      = len(a1lon)


const  = ConstCyclone.Const(prj=prj, model=model)
#thsst = 27
thsst = 27
exrvort= 3.0*1e-5
tcrvort= 3.0*1e-5
thwcore= 0.0
thdura = 36
thwind = 14.
thwdif = -9999.

const['Lat']     = a1lat
const['Lon']     = a1lon
const['thsst']   = thsst + 273.15   # K
const['exrvort'] = exrvort
const['tcrvort'] = tcrvort
const['thwcore'] = thwcore
const['thdura']  = thdura
const['thwind']  = thwind
const['thwdif']  = thwdif

exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5
slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

cy     = Cyclone.Cyclone(baseDir=wsDir, const=const, tc=True)

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
[[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[47,150]]

#----------------------------------

lonlatfontsize = 10.0
#lonrotation    = 90
lonrotation    = 0
miss_int= -9999

#------------------------
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
_, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
_, dtcppos= cy.mkInstDictC_objTC(iYM,eYM,varname='prepos')

print(_)
print(dtcppos)
sys.exit()
#------------------------
figmap   = plt.figure(figsize=(6,4))
axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())

gl        = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,901, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)

axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.coastlines(color="k")

for DTime in lDTime:
  lxyz = dtcxy[DTime]
  lprepos = dtcppos[DTime]

  if len(lxyz)==0: continue

  for i,(x,y,z) in enumerate(lxyz):
    lon = a1lon[x]
    lat = a1lat[y]

    xpre,ypre = Cyclone.fortpos2pyxy(lprepos[i][2], nx, miss_int=-9999)
    if (xpre>0):
      lonpre,latpre = a1lon[xpre], a1lat[ypre]
    else:
      lonpre,latpre = lon, lat

    #if ((lllon<lon)&(lon<urlon)&(lllat<lat)&(lat<urlat)):
    #  print x,y, "***",lon, lat
    #  M.plot( lon, lat, "o")

    #------------------------------------
    lon1, lat1 = lon, lat
    lon2, lat2 = lonpre, latpre

    scol = 'r'
    if abs(lon1 - lon2) >= 180.0:
      #--------------
      if (lon1 > lon2):
        lon05_1  = 360.0
        lon05_2  = 0.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      elif (lon1 < lon2):
        lon05_1  = 0.0
        lon05_2  = 360.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      #--------------
      axmap.plot( (lon1, lon05_1), (lat1, lat05), linewidth=1, color=scol)
      axmap.plot( (lon05_2, lon2), (lat05, lat2), linewidth=1, color=scol)

      #--------------
    else:
      axmap.plot( (lon1, lon2), (lat1, lat2), linewidth=1, color=scol)


#-- coastline ---------------
print("coastlines")
axmap.coastlines(color="k")

#stitle  = '%s %03d %s-%s'%(scen, ens, iDTime,eDTime) + '\n' + slabel
stitle  = '%s %s-%s'%(run, iDTime,eDTime) + '\n' + slabel
plt.title(stitle)
figdir  = '/home/utsumi/temp/ws/fig'
util.mk_dir(figdir)
figpath = figdir + '/tc-track.%s.png'%(run)
#figpath = figdir + '/temp-hk.png'
plt.savefig(figpath)
print(figpath)
#plt.show()
sys.exit()



# %%

# %%
