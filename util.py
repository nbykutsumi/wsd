from datetime import datetime, timedelta
import os, sys

def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

def ret_lDTime_noleap(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )

  TstpLeap = []
  for Year in range(iDTime.year, eDTime.year+1):
      if calendar.isleap(Year):
          itstp = ((datetime(Year,2,29,0)-iDTime).total_seconds()-1)\
                 /dDTime.total_seconds() + 1
          etstp = ((datetime(Year,3,1,0)-iDTime).total_seconds()-1)\
                 /dDTime.total_seconds()

          if itstp < 0: itstp=0
          if etstp > total_steps-1: etstp = total_steps-1

          TstpLeap = TstpLeap + list(range(int(itstp), int(etstp)+1))

  ltstp = [i for i in range(total_steps) if i not in TstpLeap]

  return [iDTime + dDTime*i for i in ltstp]

def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass


def ret_lYM(iYM, eYM):
  """
  iYM = [iYear, iMon], eYM = [eYear, eMon]
  """
  iYear, iMon = iYM
  eYear, eMon = eYM
  lYM = []
  for Year in range(iYear, eYear+1):
    if iYear == eYear:
      lMon = list(range(iMon,eMon+1))
    elif Year == iYear:
      lMon = list(range(iMon,12+1))
    elif Year == eYear:
      lMon = list(range(1,eMon+1))
    else:
      lMon = list(range(1,12+1))

    for Mon in lMon:
      lYM.append([Year,Mon])
  return lYM

