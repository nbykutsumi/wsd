import Cyclone
import ConstCyclone

prj     = "GANAL_JRA55" # project name. Used to read Const files
model   = "__"          # model name.   Used to read Const files
run     = 'test'        # run name.     Used to name output directories
wsbaseDir = '/home/utsumi/temp/ws'   # Output weather system directory
wsDir     = wsbaseDir + '/%s'%(run)

Year = 2014
Mon  = 7
iYM = eYM = [Year,Mon]
#-- Set TC threshold parameters
const  = ConstCyclone.Const(prj=prj, model=model)
thsst = 27
exrvort= 3.0*1e-5
tcrvort= 3.0*1e-5
thwcore= 0.0
thdura = 36
thwind = 14.
thwdif = -9999.

#const['Lat']     = a1lat
#const['Lon']     = a1lon

const['thsst']   = thsst + 273.15   # K
const['exrvort'] = exrvort
const['tcrvort'] = tcrvort
const['thwcore'] = thwcore
const['thdura']  = thdura
const['thwind']  = thwind
const['thwdif']  = thwdif
#--

cy     = Cyclone.Cyclone(baseDir=wsDir, const=const, tc=True)
dex, dtc= cy.mkInstDictC_objTC(iYM,eYM,varname='prepos')

print(dex)
print("*********************************")
print(dtc)
