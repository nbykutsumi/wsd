import sys
import os

def Const(prj=None,model=None):

    cst = {}  # constants
    cst['thtopo']    = 1500 # m
    cst['thdura']    = 36   # hours

    if prj=="GANAL_JRA55":
        
        cst['thpdif_min'] = 50  # Pa Mean(8deg x 8deg box) - Center for connect.fwd
        cst['rvort_min'] = 3.0*1.0e-5  # s-1, for connect.fwd

    elif prj=="CMIP6":
        cst['thpdif_min'] = 30  # Pa Mean(8deg x 8deg box) - Center for connect.fwd
        cst['rvort_min'] = 2.0*1.0e-5  # s-1, for connect.fwd

    else:
        print(__file__)
        print('check prj',prj)
        sys.exit()

    return cst
