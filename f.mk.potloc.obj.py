#from numpy    import *
from datetime import datetime, timedelta
import numpy.ma as ma
import numpy as np
import util
import Front
import ConstFront
import dataloader as dl
import calendar
import detect_func
import sys, os
#from dtanl_fsub import *
from front_fsub import *
#-----------------------
prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
res     = "145x288"
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'

iDTime = datetime(2014,7,1,0)
eDTime = datetime(2014,7,31,18)

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
#-------------------------
#ltq    = ["t","q"]
ltq    = ["t"]
#miss  = -9999.0
miss_in= -9999.
miss_out=-9999.
dvar  = {"t":"ta", "q":"q"}

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, timedelta(hours=6))

wsDir  = wsbaseDir + '/%s'%(run)
const  = ConstFront.Const(prj=prj, model=model)
const['Lat'] = dl.ret_lats()
const['Lon'] = dl.ret_lons()

a1lat  = dl.ret_lats()
a1lon  = dl.ret_lons()
ny     = dl.ret_ny()
nx     = dl.ret_nx()
#----------------
front  = Front.Front(baseDir=wsDir, const=const)
#------------------------
plev     = 850   #(hPa)
#************************
# FUNCTIONS
#************************
#*************************
# front locator :contour
#---------------
def mk_front_loc_contour(a2thermo, a1lon, a1lat, miss_in, miss_out):
  a2fmask1 = front_fsub.mk_a2frontmask1(a2thermo.T, a1lon, a1lat, miss_in).T
  a2fmask2 = front_fsub.mk_a2frontmask2(a2thermo.T, a1lon, a1lat, miss_in).T
  a2fmask1 = a2fmask1 * (1000.0*100.0)**2.0  #[(100km)-2]
  a2fmask2 = a2fmask2 * (1000.0*100.0)       #[(100km)-1]

  a2loc    = front_fsub.mk_a2meanaxisgrad3_h98_eq6(a2thermo.T, a1lon, a1lat, miss_in).T

  a2loc    = front_fsub.mk_a2contour(a2loc.T, 0.0, 0.0, miss_in).T
  a2loc    = ma.masked_equal(a2loc, miss_in)  

  a2loc    = ma.masked_where(a2fmask1 < 0.0, a2loc)
  a2loc    = ma.masked_where(a2fmask2 < 0.0, a2loc)
  a2loc1   = ma.masked_where(a2loc.mask, a2fmask1).filled(miss_in)
  a2loc2   = ma.masked_where(a2loc.mask, a2fmask2).filled(miss_in)

  a2loc1   = ma.masked_equal(a2loc1, miss_in).filled(miss_out)
  a2loc2   = ma.masked_equal(a2loc2, miss_in).filled(miss_out)
  return a2loc1, a2loc2

#******************************************************
##-- orog & grad orog ----

#a2orog  = dl.Load_const(var="topo")

#******************************************************
for tq in ltq:
  var = dvar[tq]
  #-----------
  for DTime in lDTime:
    a2thermo  = dl.Load_6hrPlev(model, var, DTime, plev)
    a2loc1,a2loc2  = mk_front_loc_contour(a2thermo, a1lon, a1lat, miss_in, miss_out)
    sodir, soname1, soname2   = front.path_potloc(DTime, tq)
    detect_func.mk_dir(sodir)
    #------
    a2loc1.tofile(soname1)
    a2loc2.tofile(soname2)
    print(soname1)
    print(soname2)
  
 
