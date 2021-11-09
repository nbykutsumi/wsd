from numpy import *
from datetime import datetime, timedelta
from front_fsub import *
from detect_fsub import *
import os
import socket
#****************************************************
def read_txtlist(iname):
  f = open(iname, "r")
  lines = f.readlines()
  f.close()
  lines = map(float, lines)
  aout  = array(lines, float32)
  return aout


class Front(object):
  #def __init__(self, cfg, miss=-9999.):
  def __init__(self, baseDir, const):
    #----------------
    self.baseDir = baseDir

    #------------
    self.Lat     = const["Lat"]
    self.Lon     = const["Lon"]
    self.ny      = len(self.Lat)
    self.nx      = len(self.Lon)
    self.res     = "%dx%d"%(self.ny, self.nx)
    #------------
    self.miss_in = const["miss_in"]
    self.miss_out= const["miss_out"]
    self.thgrids = const["thgrids"]
    self.Mt1     = const["thMt1"]
    self.Mt2     = const["thMt2"]
    self.Mq1     = const["thMq1"]
    self.Mq2     = const["thMq2"]
    self.thorog  = const["thorog"]
    self.thgradorog = const["thgradorog"]
    self.trace_coef = const["trace_coef"]
    #-- orog ------
    ny,nx      = self.ny, self.nx
    res        = self.res
    thorog     = self.thorog
    thgradorog = self.thgradorog

    self.maxorogname= os.path.join(self.baseDir, "const", "maxtopo.0300km.%s"%(res))
    self.a2maxorog  = fromfile(self.maxorogname, float32).reshape(ny,nx)

    a2orogmask      = zeros([self.ny, self.nx], float32)*self.miss_out
    a2orogmask      = ma.masked_where(self.a2maxorog > thorog, a2orogmask)
    self.a2orogmask = a2orogmask
    #--------------


  def path_potloc(self, DTime, tq="t"):
    """
    returns: srcDir, srcPath1(M1), srcPath2(M2)
    """
    Year,Mon,Day,Hour = DTime.year, DTime.month, DTime.day, DTime.hour
    srcDir    = os.path.join( self.baseDir,"6hr", "front.%s"%(tq),"%04d"%Year,"%02d"%Mon)
    srcPath1  = os.path.join( srcDir, "front.%s.M1.%04d.%02d.%02d.%02d.%s"%(tq,Year,Mon,Day,Hour,self.res))
    srcPath2  = os.path.join( srcDir, "front.%s.M2.%04d.%02d.%02d.%02d.%s"%(tq,Year,Mon,Day,Hour,self.res))

    return srcDir, srcPath1, srcPath2

  def path_finloc(self, DTime, tq="t"):
    Year,Mon,Day,Hour = DTime.year, DTime.month, DTime.day, DTime.hour
    srcDir    = os.path.join( self.baseDir,"6hr", "front.%s.fin"%(tq),"%04d"%Year,"%02d"%Mon)
    srcPath   = os.path.join( srcDir, "front.%s.%04d.%02d.%02d.%02d.%s"%(tq,Year,Mon,Day,Hour,self.res))

    return srcDir, srcPath

  def path_mask(self, DTime, tq="t", radkm=1000.):
    Year,Mon,Day,Hour = DTime.year, DTime.month, DTime.day, DTime.hour
    srcDir    = os.path.join( self.baseDir,"6hr", "mask.front.%s"%(tq),"%04d"%Year,"%02d"%Mon)
    srcPath   = os.path.join( srcDir, "%s.%04dkm.%04d.%02d.%02d.%02d.%s"%(tq,radkm,Year,Mon,Day,Hour,self.res))

    return srcDir, srcPath


  def mk_tfront(self, DTime, M1=False, M2=False):
    if type(M1)==bool:
      M1,M2 = self.Mt1, self.Mt2

    ny,nx     = self.ny, self.nx
    miss_out  = self.miss_out
    trace_coef= self.trace_coef

    srcDir, srcPath1, srcPath2  = self.path_potloc(DTime, "t")

    a2potloc1 = fromfile(srcPath1, float32).reshape(ny,nx)
    a2potloc2 = fromfile(srcPath2, float32).reshape(ny,nx)
    a2loc     = ma.masked_less(a2potloc1, M1)
    a2loc     = ma.masked_where(a2potloc2 < M2, a2loc).filled(miss_out)

    a2loc     = ma.masked_where(self.a2orogmask.mask, a2loc).filled(miss_out)
    #- fill --
    a2trace   = ma.masked_less(a2potloc1, M1* trace_coef)
    a2trace   = ma.masked_where(a2potloc2 < M2* trace_coef, a2trace).filled(miss_out)
    a2loc     = front_fsub.fill_front_gap_trace(a2loc.T, a2trace.T, miss_out).T
    #---------
    a2loc     = front_fsub.del_front_lesseq_ngrids_wgt(a2loc.T, self.Lat,  miss_out, self.thgrids).T
    return a2loc

  def mk_qfront(self, DTime, M1=False, M2=False, Mt1=False, Mt2=False):
    if type(M1)==bool:
      M1,M2 = self.Mq1, self.Mq2

    ny,nx     = self.ny, self.nx
    miss      = self.miss
    trace_coef= self.trace_coef

    srcDir, srcPath1, srcPath2  = self.path_potloc(DTime, "q")

    a2potloc1 = fromfile(srcPath1, float32).reshape(ny,nx)
    a2potloc2 = fromfile(srcPath2, float32).reshape(ny,nx)

    a2loc     = ma.masked_less(a2potloc1, M1)
    a2loc     = ma.masked_where(a2potloc2 < M2, a2loc).filled(miss_out)

    a2loc_t   = self.mk_tfront(DTime, Mt1, Mt2)
    a2loc_t   = detect_fsub.mk_territory_ngrids( a2loc_t.T, 2, miss_out).T
    a2loc     = ma.masked_where( a2loc_t !=miss, a2loc).filled(miss_out)

    a2loc     = ma.masked_where(self.a2orogmask.mask, a2loc).filled(miss_out)
    #- fill --
    a2trace   = ma.masked_less(a2potloc1, M1* trace_coef)
    a2trace   = ma.masked_where(a2potloc2 < M2* trace_coef, a2trace).filled(miss_out)
    a2loc     = front_fsub.fill_front_gap_trace(a2loc.T, a2trace.T, miss_out).T
    #---------
    a2loc     = front_fsub.del_front_lesseq_ngrids_wgt(a2loc.T, self.Lat, miss_out, self.thgrids).T

    return a2loc

  def mkMask_tfront(self, DTime, radkm=500, M1=False, M2=False, miss=False):

    if type(miss) == bool: miss = self.miss

    return detect_fsub.mk_territory(\
              self.mk_tfront(DTime, M1=M1, M2=M2).T, self.Lon, self.Lat, radkm*1000., imiss=self.miss, omiss=miss\
                                   ).T
 
  def mkMask_qfront(self, DTime, radkm=500, M1=False, M2=False, Mt1=False, Mt2=False, miss=False):

    if type(miss) == bool: miss = self.miss

    return detect_fsub.mk_territory(\
              self.mk_qfront(DTime, M1=M1, M2=M2, Mt1=Mt1, Mt2=Mt2).T, self.Lon, self.Lat, radkm*1000., imiss=self.miss, omiss=miss\
                                   ).T
 
    

