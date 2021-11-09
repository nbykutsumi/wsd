# %%
from numpy import *
from   datetime import datetime, timedelta
import os, sys
import calendar
import util
import dataloader as dl
import numpy as np

#******************************************************
model   = "_"
run     = "test"   # {expr}-{scen}-{ens}
res     = "145x288"
tstp    = '6hr'
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'

iDTime = datetime(2014,7,1)
eDTime = datetime(2014,7,31)


#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  model, run, res, tstp, noleap, wsbaseDir = largv[1:1+6]

  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()

  iYear,iMon, eYear, eMon = list(map(int,largv[1+6:]))

  eDay1 = calendar.monthrange(eYear,eMon)[1]
  iDTime = datetime(iYear,iMon,1,0)
  eDTime = datetime(eYear,eMon,eDay1,18)

#-------------------------

print("*"*50)
#print(__file__[0])

lvar   = ["ua","va"]
plev   = 500

lhour  = {"6hr": [0,6,12,18]
         ,"day": [0]
         }[tstp]

dDTime = timedelta(days=1)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)

miss   = -9999.

#******************************************************
# set dlyrange
#******************************************************
dnx    = {}
dny    = {}
#****************************************************
wsDir = wsbaseDir + '/%s'%(run)
ny = dl.ret_ny(model)
nx = dl.ret_nx(model)

dw         = 3
ldaydelta  = list(range(-dw, dw+1))

if   tstp=="6hr": Load_Var = dl.Load_6hrPlev
elif tstp=="day": Load_Var = dl.Load_dayPlev
#####################################################
# Function
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
#******************************************************
def date_slide(year,mon,day, daydelta, noleap):
  today       = datetime(year, mon, day)
  target      = today + timedelta(daydelta)
  targetyear  = target.year
  #***********
  if noleap == True:
    if ( calendar.isleap(targetyear) ):
      leapdate   = datetime(targetyear, 2, 29)
      #---------
      if (target <= leapdate) & (leapdate < today):
        target = target + timedelta(days=-1)
      elif (target >= leapdate ) & (leapdate > today):
        target = target + timedelta(days=1)
  #-----------
  return target

#******************************************************
for var in lvar:
  #------
  odir_root = os.path.join(wsDir,"run.mean",var)
  #------------------------------
  # make heads and tails
  #------------------------------
  for DTime in lDTime:
    #print(DTime)
    #*************
    year   = DTime.year
    mon    = DTime.month
    day    = DTime.day
    odir   = odir_root + "/%04d/%02d"%(year, mon)
    mk_dir(odir)
    #*************
    stime  = "%04d%02d%02d%02d"%(year,mon,day, 0)
    #***********
    oname  = odir + "/run.mean.%s.%04dhPa.%s.%s"%(var, plev, stime, res)
    #*********************
    # start running mean
    #*********************
    # container
    #********
    a3out  = zeros([len(ldaydelta)*len(lhour),ny,nx], float32)
    #********
    i=-1
    for daydelta in ldaydelta:
      target     = date_slide( year, mon, day, daydelta, noleap)
      targetyear = target.year
      targetmon  = target.month
      targetday  = target.day
      #-------------------
      for targethour in lhour:
        i = i+1
        tDTime = datetime(targetyear, targetmon, targetday, targethour)
        try:
          ain  = Load_Var(model, var, tDTime, plev)
        except:
          ain  = np.ones([ny,nx],'float32')*miss



        try:
          ain = ain.filled(miss)
        except:
          pass
        #--------------------
        # add 
        #--------------------
        a3out[i] = ain
    #*****************
    aout = ma.masked_equal(a3out, miss).mean(axis=0).astype('float32')
    if ma.isMA(aout):
      aout = aout.filled(miss)
    #*****************
    aout.tofile(oname)
    print(oname)


# %%
