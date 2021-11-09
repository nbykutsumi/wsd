import Cyclone
import ConstCyclone
import dataloader as dl

prj     = "GANAL_JRA55"
model   = "__"
run     = "test"
wsbaseDir= '/home/utsumi/temp/ws'

wsDir  = wsbaseDir + '/%s'%(run)
const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = dl.ret_lats()
const['Lon'] = dl.ret_lons()
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

iYM = [2014,7]
eYM = [2014,7]
dcxy ,_= cy.mkInstDictC_noTC(iYM,eYM,varname='vortlw')

ldtime = dcxy.keys()
for dtime in ldtime:
    print(dtime)
    print(dcxy[dtime])
    print()
