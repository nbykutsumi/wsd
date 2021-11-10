# WSD: Weather System Detector
"WSD" is an objective detection tool for weather systems
This tool detects weather systems (currently tropical cyclone, extratropical cyclones, and fronts) from atmospheric fields.

# Requirment (python libraries)
  * numpy       # Required
  * f2py        # This library may come with numpy
  * matplotlib  # [Optional] for drawing programs 
  * cartopy     # [Optional] for drawing programs 

# How to use
The detection programs are written in Python3.X (tested in Python3.8) and fortran.


* Edit f2py.make.py  
    Set fortran compiler path

* Compile fortran subprogram to be called from python programs  
    (f2py is used. Install it if it is not in your environment)

```bash
    python f2py.make.py detect_fsub.f90  
    python f2py.make.py front_fsub.f90  
```

   You can see "Removing build directory ..." at the end of the  
   output message if the compile is successful.  

* Edit dataloader.py  
    Set inpput file directories (explained later)

* Edit ConstCyclone.py and ConstFronts.py  
    Detectoion parameters are set in these files (explained later).


* Edit main.py  
    Set output directory (wsbaseDir).
    You can also change experiment names (explained later).

* Run main.py  
```bash
    python main.py
```

# Input
    slp     :Sea level pressure [Pa]
    ta      :Air temperature [K] at low, middle, and upper levels (normally at 850, 500, 250hPa).
    ua, va  :Meridional and zonal wind speed [m/s] at low and middle levels (normally at 850 and 500hPa).
    sst     :Monthly sea surface temperature [K] (for TC).
    topo    :Surface heigt [m]
    land    :Land/Sea mask: (For TC) Sea is 0. Land >0. 

##  Input file format:
    Big endian plain binary (no header), 32bit floating point number.
    Record grid location changes in the x (=longitudinal, 0~360, from west to east) direction faster
    and in the y (=latitudinal, -90~90, from south to north) direction next.
    e.g.,
    record#  1      2       ... nx        nx+1      nx+2    ... 2nx     2nx+1   ...
    (x,y)   (1,1)   (2,1)   ... (nx,1)    (1,2)     (2,2)   ... (nx,2)  (1,3)   ...


dataloader.py
    This program is used to read input data.
    Edit this program to set input file directories.

##   Naming convention for input binary data files 
###  6-hourly variables at pressure levels (6-hourly)
    {input-dir}/yyyymm/{var}_{plev}.yyyymmddhh.bin     # var: slp, ta, ua, va, topo
                                                       # plev: pressure level (hPa) in 4 digits (e.g., 0850, 0500, 0250)

###  6-hourly surface variables (slp)
    {input-dir}/{var}.yyyymmddhh.bin

###  Monthly surface variables (SST)
    {input-dir}/{var}.yyyymm.bin

###  Surface height
    {input-dir}/topo.bin

###  Land-sea mask
    {input-dir}/land.bin



# Main program
### main.py  
    Thip program calls subprograms.  
    You can comment out subprograms if you wish to turn them off.  

 
# Cyclone programs
## Commen programs for extratropica cyclones and tropical cyclones

### ConstCyclone.py  
    Edit this program to change parameters.  

    #Parameters
    thtopo      : Grid points lower than this height [m] are not used for cyclone detection.  
    thdura      : Minimum lifetime duration [hours] required for cyclones.  
    thpdif_min  : Minimum central pressure difference [Pa] # Surrounding average - center  
    rvort_min   : Minimum relative vorticy required for cyclones to be tracked [s-1].  

### Cyclone.py  
    This program provides functions and python class handling cyclones.  
    e.g., You can obtain a dictionary of cyclone centers (key=datetime) using mkInstDictC_noTC function.   
    See c.sample.py  
     

### c.runmean.py  
    This program makes running mean stearing wind data.  
    No parameters from ConstCyclone.py are used.  

### c.find.cyclone.py  
    This program finds candidates of cyclone centers from sea level pressure (SLP) field.  
    Average SLP around cyclone center (to calculate central pressure difference) and maximum vorticity around cyclone centers are also calculated.  
    Lats, lons from ConstCyclone.py are used.  

### c.connect.fwd.py  
    This program tracks the movement of cyclone centers.  
    Genesis position, genesis time, and previous position are identified.  
    Lats, lons, thpdif_min, rvort_min, thdist_search, thtopo from ConstCyclone.py are used.  

### c.connect.bwd.py  
    This program tracks the movent of cyclone centers in a backward direction.  
    Next position, lifetime duration, and final location are identified.  
    Lats, lons, and thtopo from ConstCyclone.py are used.  

### c.mk.clist.obj.py  
    This program joins 6-hourly data file (.npy) to make monthly file (.npy).  
    Lats, Lons from ConstCyclone.py are used (only to load Cyclone.py and to use load_clist_org)  

### c.clist.npy2bin.py  
    This program converts numpy format data (.npy) to plain binary data (.bin).  

## Programs for tropical cyclones  
### tc.mk.clist.obj.py  
    This program joins TC-related 6-hourly data file (.npy) to make monthly file (.npy).  
    Lats, Lons from ConstCyclone.py are used  


### tc.mk.clist.obj.initState.py  
    This program checks SST and land/sea mask at the genesis location/time of TCs  
    Lats, Lons from ConstCyclone.py are used (only to load Cyclone.py and to use load_clist_org)  

### tc.clist.npy2bin.py  
    This program converts TC-related numpy format data (.npy) to plain binary data (.bin).  

#

## Output variables  
lat             : Latitude  
lon             : Longitude  
nowpos          : Current position  
age             : Hours from genesis time  
dura            : Lifetime hours from genesis to the end  
epos            : Final (end) position  
idate           : Genesis (initial) time  
ipos            : Genesis (initial) position  
nextpos         : Next timestep position  
prepos          : Previous timestep position  
slp_mean_adj    : SLP averaged over adjacent 8 gird boxes [Pa]  
slp_mean_box    : SLP averaged over 8deg x 8deg box centered on the cyclone [Pa]  
vortlw          : Low level (default:850hPa) vorticity [s-1]  
vortlw_max_box  : Maximum low level (default:850hPa) vorticity [s-1] in the 8deg x 8deg box centered on the cyclone  
x               : x (starts from 0) position of the cyclone center  
y               : y (starts from 0) position of the cyclone center  


## Output data format  
### ".py"  
    numpy data format  
    This data can be read using Python and numpy  
    With python,  
```
    import numpy  
    dat = np.load(filename)  
```

### ".bin"  
    Plain binary (no header) file. Big endian.  
    See Cyclone.py for the data type for each variable (copied below).  
```
    self.dNumType= {  
                    "x"       :int32,  
                    "y"       :int32,  
                    "life"    :int32,  
                    "dura"    :int32,  
                    "ipos"    :int32,  
                    "epos"    :int32,  
                    "idate"   :int32,  
                    "time"    :int32,  
                    "age"     :int32,  
                    "nowpos"  :int32,  
                    "lastpos" :int32,  
                    "prepos"  :int32,  
                    "nextpos" :int32,  
                    "pgrad"   :float32,  
                    "vortlw"  :float32,  
                    "vortlw_max_box":float32,  
                    "pgmax"   :float32,   
                    "pmean"   :float32,   
                    "iedist"  :float32,   
                    "dtlw"    :float32,   
                    "dtmd"    :float32,   
                    "dtup"    :float32,   
                    "wmeanlow":float32,   
                    "wmeanup" :float32,   
                    "wmaxlw"  :float32,   
                    "wmaxup"  :float32,   
                    "sst"     :float32,   
                    "land"    :float32,   
                    "initsst" :float32,   
                    "initland":float32,   
                    "iedist"  :float32,   
                    "lat"     :float32,   
                    "lon"     :float32,   
                    "slp"     :float32,   
                    "slp_mean_adj":float32,  
                    "slp_mean_box":float32,  
                   }  
```


    The record length of the output variable files is the same for month. Records are saved in the same order.  
    e.g., N-th record of "dura" is the life time duration of the cyclone at N-th record of "nowpos" location at the time of N-th record of "time".  

## Decording "time"  
```
    def solve_time(stime):  
      year = int( stime/10**6 )  
      mon  = int( (stime - year*10**6)/10**4 )  
      day  = int( (stime - year*10**6 - mon*10**4)/10**2)  
      hour = int( (stime - year*10**6 - mon*10**4 - day*10**2) )  
      return year, mon, day, hour  
```

## Decording position (ipos, epos, prepos, nowpos, nextpos)
```
    def fortpos2pyxy(number, nx, miss_int):
      if (number == miss_int):
        iy_py = miss_int
        ix_py = miss_int
      else:
        iy_py = int((number-1.0)/nx)      # iy_py = 0,1,2,..
        ix_py = int(number - nx*iy_py-1)   # ix_py = 0,1,2,..
      #----
      return ix_py, iy_py
```


### c.draw.tracks.py  
    Sample python program for drawing all cyclone tracks (numpy, matplotlib, cartopy are required)  

### tc.draw.tracks.py  
    Sample python program for drawing tropical cyclone tracks (numpy, matplotlib, cartopy are required)  


### Finding tropical cyclones  
    Tropical cyclones are detected if a cyclone center satisfies the following conditions
```
    "dura"                  >= thdura  
    "vortlw"                >= thrvortd  
    "dtlw"+"dtmd"+"dtup"    >= thwcore  
    "wmaxlw"                >= thwind  
    "wmaxlw" - "wmaxup"     >= thwdif
    "initsst"               >= thinitsst
    "initland"              == 0
```

    You can also use "mkInstDictC_objTC" method in Cyclone.py  
    This method generates python dictionary for extratropical cyclones and tropical cyclones.  
    Dictionary key: datetime  
    Dictionary value: specified variable (e.g., "nextpos","prepos", "vortlw", etc.)  
    See **tc.test.py** for detail.  


# Front programs  
### ConstFront.py  
    Edit this program to change parameters.  

    #Parameters
    thorog      : Grid points lower than this height [m] close to such grid points are not used for front detection.  
    thgradorog  : Threshold of average slope between adjacent grid points. Grid points with large slopes are excluded from front detection.  
    trace_coef  : Parameter to be used for filling gaps of fronts. Nomally, this parameter is not tuned.  
    thMt1       : Threshold for M1 for front detection   # K/100km/100km  
    thMt2       : Threshold for M2 for front detection   # K/100km  
    thMt1       : [Option] Parameter 1 for front detection. Not used for the default setting.   # K/100km/100km  
    thMt2       : [Option] Parameter 2 for front detection. Not used for the default setting.   # K/100km  
    thgrids     : This parameter is used to delete too short fronts. Fronts less than this number of grids are removed.  
    miss_in     : Missing value for input data  
    miss_out    : Missing value for output data  


### Front.py  
    This program provides functions and python class handling fronts.  

### f.mk.orogdata.py   
    This program finds maximum surface height within a given radius at each grid point.  

### f.mk.potloc.obj.py  
    This program calculates parameter M1 and M2.  

### f.mk.front.bin.py  
    This program detects fronts from parameter M1 and M2.  
    Output: Plain binary file.  

### f.draw.front.py   
    Sample python program for drawing fronts (numpy, matplotlib, cartopy are required)  



# Author
Nobuyuki UTSUMI @KUAS, Japan, Sep 13, 2021.





