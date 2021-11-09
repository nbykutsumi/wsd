import sys

myname = "ConstFront.py"

def Const(prj=None, model=None):
    cst = {}
    cst["thorog"]     = 1500  # m
    cst["thgradorog"] = 1.0e+8 # m/m
    cst["trace_coef"] = 0.8

    if (prj,model) == ("GANAL_JRA55","__"):
        cst["thMt1"] = 0.30  # K/100km/100km
        cst["thMt2"] = 0.6   # K/100km
        cst["thMq1"] = 2.3*1.0e-4   # test
        cst["thMq2"] = 0.9*1.0e-3   # test
        cst["thgrids"]= 5/1.25
        cst["miss_in"]= -9999.
        cst["miss_out"]= -9999.

    else:
        print(myname,":check prj,model",prj,model)
        sys.exit()

    return cst


