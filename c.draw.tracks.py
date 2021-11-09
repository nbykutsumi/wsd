import matplotlib
matplotlib.use('Agg')
import sys, os
import numpy as np
from   numpy import *
import dataloader as dl
from   datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import detect_func
import util
import Cyclone
import ConstCyclone
#--------------------------------------

prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
res     = "145x288"
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'
figdir   = '/home/utsumi/temp/ws/fig'

iDTime = datetime(2014,7,1,0)
eDTime = datetime(2014,7,31,18)
lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(hours=6))

wsDir  = wsbaseDir + '/%s'%(run)
const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = dl.ret_lats()
const['Lon'] = dl.ret_lons()
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,360]]
#----------------------------------
a1lat  = dl.ret_lats()
a1lon  = dl.ret_lons()
ny     = dl.ret_ny()
nx     = dl.ret_nx()


lonlatfontsize = 10.0
#lonrotation    = 90
lonrotation    = 0
miss_int= -9999

#------------------------
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
#_, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
#_, dtcppos= cy.mkInstDictC_objTC(iYM,eYM,varname='prepos')

dcxy ,_= cy.mkInstDictC_noTC(iYM,eYM,varname='vortlw')
dcppos,_= cy.mkInstDictC_noTC(iYM,eYM,varname='prepos')

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
  lxyz = dcxy[DTime]
  lprepos = dcppos[DTime]


  if len(lxyz)==0: continue


  for i,(x,y,z) in enumerate(lxyz):
    lon = a1lon[x]
    lat = a1lat[y]

    xpre,ypre = Cyclone.fortpos2pyxy(lprepos[i][2], nx, miss_int=-9999)
    if (xpre>0):
      lonpre,latpre = a1lon[xpre], a1lat[ypre]
    else:
      lonpre,latpre = lon, lat

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

util.mk_dir(figdir)
figpath = figdir + '/exc.track.png'
plt.savefig(figpath)
print(figpath)

sys.exit()


