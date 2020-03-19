#!/usr/bin/env python
#
"""
This script filters a 2D field of data with Lanczos windowing method and writes the filtered signal in a new file that has the same structure than the input file.
It is possible to filter one or multiple files (same variable) and maps of total and filtered signal can also be produced.
External module needed : 
  - WavenumberSpectrum : https://github.com/lesommer/codes/blob/master/WavenumberSpectrum.py
  - oocgcm filtering module : https://github.com/lesommer/oocgcm
"""

## path for mdules

import sys
sys.path.insert(0,"/home/albert7a/lib/python")

## imports

import numpy as np
import xarray as xr
import GriddedData
sys.path.insert(0,"/home/albert7a/git/xscale")
import xscale

#- Other modules
import numpy.ma as ma
from netCDF4 import Dataset

### local/specific imports
import oocgcm
import oocgcm.filtering
import oocgcm.filtering.linearfilters as tf

from datetime import date
import time

## read the data

def read(filename,varname):
   """Return navlon,navlat,data.
   """
   navlon = xr.open_dataset(filename,chunks={'x':500,'y':500})['nav_lon']
   navlat = xr.open_dataset(filename,chunks={'x':500,'y':500})['nav_lat']
   data = xr.open_dataset(filename,chunks={'x':500,'y':500})[varname]
   return navlon,navlat,data


## filter the data

def filt(data,nwin,fcut):
    """ Filter the data with Lanczos window of size nwin and fcut cut-off frequency
        Return signal_LS[Large scale] and signal_SS[Small scale]
    """
    win_box = data.window
    win_box.set(n=nwin, dim=['time_counter'], cutoff=fcut)
    signal_LS = win_box.convolve()
    signal_SS=data-signal_LS
    return signal_LS,signal_SS

## write output file
def write(filein,signal_SS,signal_LS,nwin,fcut,varname):
    """Write the output file with the same structure than the input file
       In the variable varname_filt, the fine scale signal is written
    """
    
    outname=filein[0:len(filein)-3]+'_filt-n'+str(nwin)+'-f'+str(fcut)+'.nc'
    print('output file is '+outname)
    dstin=Dataset(filein,'r')
    dstout=Dataset(outname,'w')

    today=date.today()
    dstout.description = "Data time-filtered with window size of "+str(nwin)+" and cut-off frequency of "+str(fcut)+" obtained with TimeFilter.py script "+str(today.day)+"/"+str(today.month)+"/"+str(today.year)

    #Copy the structure of the input file
    for dname, the_dim in dstin.dimensions.iteritems():
      dstout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    for v_name, varin in dstin.variables.iteritems(): 
      if v_name == varname: 
        continue
      outVar = dstout.createVariable(v_name, varin.datatype, varin.dimensions)
      outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
      outVar[:] = varin[:]

    #Create output variable
    datain=dstin[varname]
    dataout1=dstout.createVariable(varname+'_wave',datain.datatype,datain.dimensions)
    dataout2=dstout.createVariable(varname+'_nowave',datain.datatype,datain.dimensions)
    dataout1.setncatts({k: datain.getncattr(k) for k in datain.ncattrs()})
    dataout2.setncatts({k: datain.getncattr(k) for k in datain.ncattrs()})

    dataout1[:] = signal_SS[:]
    dataout2[:] = signal_LS[:]
    dstout.close()    


## parser and main
def script_parser():
    """Customized parser.
    """
    from optparse import OptionParser
    usage = "usage: %prog [options] file.nc varname n[size of Lanczos window] f[cut-off frequency]"
    parser = OptionParser(usage=usage)
    return parser


def main():
    parser = script_parser()
    (options, args) = parser.parse_args()
    if len(args) < 4: # print the help message if number of args is not 3.
        parser.print_help()
        sys.exit()
    optdic = vars(options)
    ## One file in input
    print time.strftime('%d/%m/%y %H:%M',time.localtime())
    if len(args) == 4:
      filein = args[0]
      varname = args[1]
      nwin=int(args[2])
      fcut=float(args[3])
      navlon,navlat,data = read(filein,varname)
      signal_LS,signal_SS = filt(data,nwin,fcut)
      write(filein,signal_SS,signal_LS,nwin,fcut,varname)
    print time.strftime('%d/%m/%y %H:%M',time.localtime())

if __name__ == '__main__':
    sys.exit(main() or 0)
