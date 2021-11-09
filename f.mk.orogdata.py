#from numpy import *
from detect_fsub import *
import numpy as np
import util, sys
import dataloader as dl

calcflag  = True
prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
res     = "145x288"
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, wsbaseDir = largv[1:1+5]
#-------------------------
radkm = 300.  # (km)

a1lat  = dl.ret_lats()
a1lon  = dl.ret_lons()
ny     = dl.ret_ny()
nx     = dl.ret_nx()
miss  = -9999.
#----------------

oDir       = wsbaseDir + "/%s/const"%(run)
util.mk_dir(oDir)
maxorogname= oDir + "/maxtopo.%04dkm.%s"%(radkm, res)

if calcflag==True:
  a2orog     = dl.Load_const(model=model, var="topo")
  a2maxorog  = detect_fsub.mk_a2max_rad(a2orog.T, a1lon, a1lat, radkm, miss).T
  #--- write to file -------
  a2maxorog.tofile(maxorogname)
  print(maxorogname)
##--- figure: max orog ----
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
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
#
#[[lllat,lllon],[urlat,urlon]] = [[-90,-180],[90,180]]
#a2fig = np.roll(a2maxorog, int(nx/2), axis=1)
#axmap.coastlines(color="k")
#im = axmap.pcolormesh(X, Y, a2fig)
#plt.colorbar(im)
#figname = maxorogname + ".png"
#plt.savefig(figname)
#print(figname)
