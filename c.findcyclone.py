from numpy import *
from detect_fsub import *
from datetime import datetime, timedelta
from detect_func import box_filtering
import util
import dataloader as dl
import calendar
import os, sys, shutil
import numpy as np
import Cyclone
import ConstCyclone
#------------------------------------
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
  #iDTime = datetime(iYear,iMon,1,6)
  iDTime = datetime(iYear,iMon,1,0)
  eDTime = datetime(eYear,eMon,eDay,18)
#-------------------------

dDTime = timedelta(hours=6)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)

tstp        = "6hr"
wsDir = wsbaseDir + '/%s'%(run)

const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = dl.ret_lats(model)
const['Lon'] = dl.ret_lons(model)
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)


a1lat  = dl.ret_lats(model)
a1lon  = dl.ret_lons(model)
ny     = dl.ret_ny(model)
nx     = dl.ret_nx(model)
miss_in= dl.ret_miss(model)
miss_out= -9999.

dlon = a1lon[1] - a1lon[0]
dlat = (a1lat[1:] - a1lat[:-1]).mean()
dx_4deg = int(4.0/dlon)
dy_4deg = int(4.0/dlat)
#####################################################
def check_file(sname):
  if not os.access(sname, os.F_OK):
    print("no file:",sname)
    sys.exit()
#####################################################
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#################################################
def mk_dir_tail(var, tstp, model, expr, ens):
  odir_tail = var + "/" + tstp + "/" +model + "/" + expr +"/"\
       +ens
  return odir_tail
#####################################################
def mk_namehead(var, tstp, model, expr, ens):
  namehead = var + "_" + tstp + "_" +model + "_" + expr +"_"\
       +ens
  return namehead
#****************************************************
def read_txtlist(iname):
  f = open(iname, "r")
  lines = f.readlines()
  f.close()
  lines = list(map(float, lines))
  aout  = array(lines, float32)
  return aout
##**************************************************
## Mean Sea Level Pressure
##------------------------
for DTime in lDTime:
  print(DTime)
  year = DTime.year
  mon  = DTime.month
  day  = DTime.day
  hour = DTime.hour
  #***************************************
  # pgrad
  #---------------------------------------
  a2slp   = dl.Load_6hrSfc(model, "slp", DTime)
  try:
    a2slp = a2slp.filled(miss_in)
  except:
    pass
 

  a2center = -detect_fsub.find_localmax(-a2slp.T, a1lat, a1lon, miss_in,  0).T

  a1y,a1x = np.where( ma.masked_not_equal(a2center,0))

  a1slp = a2center[a1y,a1x]

  xdir, xname = cy.path_clist_org("x",DTime)
  ydir, yname = cy.path_clist_org("y",DTime)

  slpdir, slpname = cy.path_clist_org("slp",DTime)

  mk_dir(xdir)
  mk_dir(ydir)
  mk_dir(slpdir)
  np.save(xname, a1x)
  np.save(yname, a1y)
  np.save(slpname, a1slp)
  #***************************************
  # Mean slp
  #---------------------------------------
  a1slp_mean_adj = box_filtering(a2slp, a1y, a1x, 1, 1, func='mean', miss_in=miss_in, miss_out=miss_out)
  slpmean_adj_dir, slpmean_adj_name = cy.path_clist_org("slp_mean_adj",DTime)
  mk_dir(slpmean_adj_dir)
  np.save(slpmean_adj_name, a1slp_mean_adj)


  a1slp_mean_box = box_filtering(a2slp, a1y, a1x, dy_4deg, dx_4deg, func='mean', miss_in=miss_in, miss_out=miss_out)
  slpmean_box_dir, slpmean_box_name = cy.path_clist_org("slp_mean_box",DTime)
  mk_dir(slpmean_box_dir)
  np.save(slpmean_box_name, a1slp_mean_box)

  ##***************************************
  ## rvort @ 850
  ##---------------------------------------
  rvortdir, rvortname = cy.path_clist_org("vortlw",DTime)
  mk_dir(rvortdir)

  a2u       = dl.Load_6hrPlev(model, "ua", DTime, 850)
  a2v       = dl.Load_6hrPlev(model, "va", DTime, 850)

  try:
    a2u = a2u.filled(miss_in)
    a2v = a2v.filled(miss_in)
  except:
    pass


  a2rvort   = detect_fsub.mk_a2rvort(a2u.T, a2v.T, a1lon, a1lat, miss_in, miss_out).T

  a2mask    = ma.masked_equal(a2rvort, miss_out).mask
  a2rvort[:int(ny/2)] = -a2rvort[:int(ny/2)]  # The signs of the missing values in the south hemisphere are also fliped
  a2rvort   = ma.masked_where(a2mask, a2rvort).filled(miss_out)

  #- find maximum vorticity in 9x9 box --
  a1maxvort_adj = box_filtering(a2rvort, a1y, a1x, 1, 1, func='max', miss_in=miss_out, miss_out=miss_out)

  ##- find maximum vorticity in a (2*dy+1) x (2*dx+1) box --
  a1maxvort_box = box_filtering(a2rvort, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss_out, miss_out=miss_out)


  #-- Save ----
  vort_max_box_dir, vort_max_box_name = cy.path_clist_org("vortlw_max_box",DTime)
  vort_max_adj_dir, vort_max_adj_name = cy.path_clist_org("vortlw",DTime)


  mk_dir(vort_max_box_dir)
  mk_dir(vort_max_adj_dir)
  np.save(vort_max_box_name, a1maxvort_box)
  np.save(vort_max_adj_name, a1maxvort_adj)

