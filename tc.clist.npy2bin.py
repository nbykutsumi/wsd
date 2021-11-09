import numpy as np
import Cyclone
import ConstCyclone
import dataloader as dl
import util
import os, sys

prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
res     = "145x288"
noleap  = False
wsbaseDir= '/home/utsumi/temp/ws'
wsDir    = wsbaseDir + "/%s"%(run)

iYear, iMon = [2014,7] # Start of the detection period
eYear, eMon = [2014,7] # End of the detection period

#--- argv -----------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, wsbaseDir = largv[1:1+5]
  iYear,iMon, eYear, eMon = list(map(int,largv[1+5:]))
  

lYM = util.ret_lYM([iYear,iMon], [eYear,eMon])

const  = ConstCyclone.Const(prj=prj, model=model)
const["Lat"] = dl.ret_lats()
const["Lon"] = dl.ret_lons()
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

#lvar   = ["dura","nowpos","time","vortlw"]
lvar  = ["dtlw", "dtmd", "dtup", "wmaxlw", "wmaxup", "sst", "land", "initsst", "initland"]

for var in lvar:
    dtype = cy.dNumType[var]
    for Year,Mon in lYM:
        ipath = cy.path_clist(var, Year, Mon)[1]
        odir  = wsDir + "/6hr/clist.bin/%04d/%02d"%(Year,Mon)
        opath = odir + "/" + os.path.basename(ipath)[:-4] + ".bin"
        util.mk_dir(odir)

        adat = np.load(ipath)
        adat.astype(dtype).tofile(opath) 
        print(opath)

