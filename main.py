import subprocess
import os, sys
import util
import socket
import dataloader as dl
#*********************************
# Edit here
#*********************************
prj     = "GANAL_JRA55" # project name. Used to read Const files
model   = "__"          # model name.   Used to read Const files
run     = 'test'        # run name.     Used to name output directories
noleap  = False         # True: with leap years,    False: without leap years
tstp_runmean = "6hr"    # Do not change for default setting
wsbaseDir = '/home/utsumi/temp/ws'   # Output weather system directory

iYear, iMon = [2014,7] # Start of the detection period
eYear, eMon = [2014,7] # End of the detection period

iYear_data  = iYear # First year of the available data
eYear_data  = eYear # Last year of the available data
iMon_data   = iMon # First month of the available data

ny          = dl.ret_ny(model)
nx          = dl.ret_nx(model)
res         = "%dx%d"%(ny,nx)     # "NXxNY"       Used to name output files


#
#
#prj     = "CMIP6" # project name. Used to read Const files
#model   = "MRI-ESM2-0.historical.r1i1p1f1"          # model name.   Used to read Const files
#run     = 'test'        # run name.     Used to name output directories
#noleap  = False         # True: with leap years,    False: without leap years
#tstp_runmean = "6hr"    # Do not change for default setting
#wsbaseDir = '/tank/out/ws/CMIP6/%s'%(model)   # Output weather system directory
#
##iYear, iMon = [1850,1] # Start of the detection period
##eYear, eMon = [2049,12] # End of the detection period 
##iYear, iMon = [1980,1] # Start of the detection period
#iYear, iMon = [1980,1] # Start of the detection period
#eYear, eMon = [2014,12] # End of the detection period 
#
#*********************************
def exec_func(lcmd):
    scmd = " ".join(map(str, lcmd))
    print(scmd)
    lcmd = list(map(str, lcmd))
    #p = subprocess.Popen(lcmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.Popen(lcmd, shell=False)
    p.wait()

#*********************************
# Cyclone (ExC and TC)
#---------------------------------
cmd = ["python","c.runmean.wind.py"
        , model, run, res, tstp_runmean, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","c.findcyclone.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","c.connectc.fwd.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)


flagresume = False  # True or False
cmd = ["python","c.connectc.bwd.py"
        , prj, model, run, res, noleap, flagresume
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","c.mk.clist.obj.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","c.clist.npy2bin.py"
        , prj, model, run, res
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)


#*********************************
# Cyclone (TC)
#---------------------------------
cmd = ["python","tc.mk.clist.obj.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","tc.mk.clist.obj.initState.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
        , iYear_data, iMon_data
    ]
exec_func(cmd)

cmd = ["python","tc.clist.npy2bin.py"
        , prj, model, run, res
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

#
#*********************************
# Front
#---------------------------------
cmd = ["python","f.mk.orogdata.py"
        , prj, model, run, res
        , wsbaseDir
    ]
exec_func(cmd)

cmd = ["python","f.mk.potloc.obj.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)

cmd = ["python","f.mk.front.bin.py"
        , prj, model, run, res, noleap
        , wsbaseDir
        , iYear, iMon, eYear, eMon
    ]
exec_func(cmd)



