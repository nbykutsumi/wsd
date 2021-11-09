import sys, os
import numpy as np
import dataloader as dl
from   datetime import datetime, timedelta
import detect_func
import util
import ConstFront
import Front
import calendar
#--------------------------------------

prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
res     = "145x288"
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'

iDTime = datetime(2014,7,1,6)
eDTime = datetime(2014,7,30,0)

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, noleap, wsbaseDir = largv[1:1+6]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()

  iYear,iMon, eYear, eMon = list(map(int,largv[1+6:]))
  eDay   = calendar.monthrange(eYear,eMon)[1]
  iDTime = datetime(iYear,iMon,1,6)
  eDTime = datetime(eYear,eMon,eDay,18)
#------------------------
wsDir  = wsbaseDir + '/%s'%(run)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime = ret_lDTime(iDTime,eDTime, timedelta(hours=6))

const  = ConstFront.Const(prj=prj, model=model)
const['Lat'] = dl.ret_lats()
const['Lon'] = dl.ret_lons()
f      = Front.Front(baseDir=wsDir, const=const)

outbaseDir = wsbaseDir + "/%s/front.out"%(run)
#----------------------------------
a1lat  = dl.ret_lats()
a1lon  = dl.ret_lons()
ny     = dl.ret_ny()
nx     = dl.ret_nx()

#------------------------
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
#------------------------
for DTime in lDTime:
    year,mon,day,hour = DTime.timetuple()[:4]
    M1 = const["thMt1"]
    M2 = const["thMt2"]
    a2front = f.mk_tfront(DTime, M1=M1, M2=M2)

    outDir  = wsbaseDir + "/%s/6hr/front.out/%04d/%02d"%(run,year,mon)
    outPath = outDir + "/front.%04d.%02d.%02d.%02d.%dx%d"%(year,mon,day,hour,ny,nx)
    util.mk_dir(outDir)
    a2front.astype("float32").tofile(outPath)
    print(outPath)
    ##-- figure ---
    #a2fig = a2front
    #a2fig = np.roll(a2fig, int(nx/2), axis=1)

    #projection=ccrs.PlateCarree()
    #figmap   = plt.figure(figsize=(6,4))
    #axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=projection)
    #
    #dlat = 1.25
    #dlon = 1.25
    #a1latBnd = [-90] + np.linspace(-90+dlat*0.5, 90-dlat*0.5, ny-1).tolist() + [90]
    ##a1lonBnd = [0] + np.linspace(0+dlon*0.5, 360-dlon*0.5, nx-1).tolist() + [360]
    #a1lonBnd = np.linspace(-180-dlon*0.5, 180+dlon*0.5, nx+1).tolist()
    #X,Y = np.meshgrid(a1lonBnd, a1latBnd)

    #im = axmap.pcolormesh(X, Y, a2fig)
    #
    #axmap.coastlines(color="k")
    #plt.colorbar(im)

    #figdir  = '/home/utsumi/temp/ws/fig'
    #util.mk_dir(figdir)
    #figpath = figdir + '/test.front.%04d%02d%02d%02d.png'%(year,mon,day,hour)
    #plt.savefig(figpath)
    #print(figpath)

sys.exit()


