#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import time
import json
import pathlib
from netCDF4 import Dataset
from datetime import datetime, timedelta
#from netcdftime import datetime as datetimex, date2num
import glob
from socket import gethostname


"""
aggregate files
"""
__svnHeader__ = "$Header: https://svn-kerberos.ecdf.ed.ac.uk/repo/geos/CESD/UM_utilities/aggregate/splitAgg.py 736 2019-09-05 16:11:57Z mjm@EASE.ED.AC.UK $"
__svnRevision__ = "$Rev: 736 $"
__svnId = "$Id: splitAgg.py 736 2019-09-05 16:11:57Z mjm@EASE.ED.AC.UK $"
__version__="0.1"

TIME_NAME="time"
TBOUNDS_NAME="time_bounds"

class MyNCvar:

    def __init__(self,vname, indata,outfmt,outdir=None,prefix_name=None):
        global TIME_NAME
        global TBOUNDS_NAME
        self.tname=TIME_NAME
        self.tb_name=TBOUNDS_NAME
        self.firstFile=True
        self.gotDimensions={}
        self.coords={}
        outFile="%s.nc"%vname
        self.outvar={}
        if prefix_name:
            outFile="%s%s"%(prefix_name, outFile)
        if outdir:
           outFile = outdir/outFile
        self.dout=Dataset(outFile,"w",format=outfmt)
#        copy_dimensions(indata,vname,self.dout, self.gotDimensions,self.coords)
        self.prepare_nc(indata, vname)

    def add_data(self,outname,varin,this_time,this_tbounds, tindex):
        if outname  in self.outvar and self.dout[outname].dimensions[0]==self.tname:
            self.outvar[outname][tindex,:]=varin[0,:]
            self.tout[tindex]=this_time
            self.outvar[self.tb_name][tindex,:]=this_tbounds[:]
        else:
            print( "failed on ",outname)

    def write_file(self,vv):
        history="%s on %s in %s using %s\n%s"%(time.asctime( time.localtime(time.time())),gethostname(),
                      os.getcwd(), __svnHeader__,my_args())
        self.dout.history=history
        self.dout.source="Univ. of Edinburgh, GeoSciences "
        self.dout.use="aggregation of climate data, reordering size1 dimensions where time is not first"
        self.dout.contact=os.environ['USER']

                # close the dataset
        # print "closing ",vv
        self.dout.close()

    def copy_dimensions(self,din, vname):
       # copy dimensions for a variable called vname
    #print "copy dims for var ",vname
        reordered=None
        varin=din.variables[vname]
        outname=vname
                          # check desired order
        mydims=[None]*len(varin.dimensions)
        if len(varin.dimensions) > 2 and \
                    varin.dimensions[0] != self.tname and  \
                    din.dimensions[varin.dimensions[0]].size ==1 and \
                    varin.dimensions[1]==self.tname:
            mydims[0]=varin.dimensions[1]
            mydims[1]=varin.dimensions[0]
            mydims[2:]=varin.dimensions[2:]
            reordered=varin.dimensions[0]
        elif varin.dimensions[0] == self.tname or len(varin.dimensions) <= 2:
                    mydims=varin.dimensions[:]
        else:
                    print("cant do %s "%vname)
                    print(varin.dimensions)
                    print("cannot manage varaibles with more than size 1 in time")
                    print("cannot manage a file with 2 time dimensions")
                    exit(1)
        for dname in mydims: #varin.dimensions[:]:
           if not dname in self.gotDimensions:
               self.gotDimensions[dname]=True
#               if dname == "longitude":
                        #print "line for debug break"
               if dname in din.variables.keys():
                   the_dim=din.variables[dname]

                   # print "copy dimension variable",dname
                   self.dout.createDimension(dname, len(the_dim))
        
                   self.outvar[dname] = self.dout.createVariable(dname, the_dim.datatype, the_dim.dimensions)
                                    # Copy attributes
                   self.outvar[dname].setncatts({k: the_dim.getncattr(k) for k in the_dim.ncattrs()})
                                     # copy array for dimension
                   self.outvar[dname][:]=din.variables[dname][:]
                          # note a diension cannot have time dimensions, so just copy it
                   if "bounds" in the_dim.ncattrs():
                       bounder=the_dim.getncattr("bounds")
                       #print " bounds found",bounder
                                    # call is recursive...
                       (breordered, boutname, bmydims)=self.copy_dimensions(din,bounder)
                       self.create_agg_variable(din, bounder,bounder, bmydims,breordered)
               else:
                                     # scalar dimension, no array to copy
                   #print "scalardimension?",dname
                   the_dim=din.dimensions[dname]
                   self.dout.createDimension(dname, the_dim.size)
        if "bounds" in varin.ncattrs():
             bounder=varin.getncattr("bounds")
             #print " bounds in var found",bounder
                                    # call is recursive...
             (breordered, boutname, bmydims)=self.copy_dimensions(din, bounder)
             self.create_agg_variable(din, bounder,bounder, bmydims,breordered)

        return reordered, outname, mydims

    def create_agg_variable(self,din,  vname, outname, mydims, reordered):
    # create variable, and copy if its not time dependent.
                      
        varin=din.variables[vname]
      #  print "creating var",outname, "with ",mydims
        if "coordinates" in varin.ncattrs():
           coList=varin.getncattr("coordinates")
           for cc in coList.split():
               self.coords[cc]=True
        if '_FillValue' in varin.ncattrs():
             fv=varin.getncattr('_FillValue')
             self.outvar[outname] = self.dout.createVariable(outname, varin.datatype, mydims, fill_value=fv)
        else:
             self.outvar[outname] = self.dout.createVariable(outname, varin.datatype, mydims)
       
        for k in varin.ncattrs():
            if k != "_FillValue":
             self.outvar[outname].setncatts({k: varin.getncattr(k)})
                #    outvar[outname].setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
       
        if reordered:            # add another attribute to flag this.... for user's testing only
                     self.outvar[outname].setncatts({'reordered':reordered})


        if mydims[0] != self.tname:
                        # eg. air_pressure
                      # print "t invariant",outname
                      self.outvar[outname][:]=varin[:]
        return 

    def prepare_nc (self,din, vname):
        self.dout.createDimension(self.tname, None)
        self.tout= self.dout.createVariable(self.tname, 'f', [self.tname])
        the_dim=din[self.tname]
        self.tout.setncatts({k: the_dim.getncattr(k) for k in the_dim.ncattrs()})
        self.gotDimensions[self.tname]=True
        self.coords={} # to know what is asked for via coordinates attribute.

        if self.tb_name in din.variables:
           (reordered, self.tb_name, mydims) =self.copy_dimensions(din, self.tb_name)
           self.create_agg_variable(din,  self.tb_name, self.tb_name, mydims, reordered)
           self.gotDimensions[self.tb_name]=True

               # copy each required variable in turn
               # copy other dimension definitions as meet them
        vlist=[vname]  # quick test -t idy up if ok!
        for vname in vlist:
           varin=din.variables[vname]

           (reordered, outname, mydims) =self.copy_dimensions(din, vname)
           self.create_agg_variable(din,  vname, outname, mydims, reordered)
       
        for vname in self.coords:
           varin=din.variables[vname]
           outname=vname
           if not vname in self.gotDimensions:
               if vname in din.variables.keys():

                   (reordered, outname, mydims) =self.copy_dimensions(din,vname)
                   self.outvar[vname]=self.create_agg_variable(din,  vname, outname, mydims, reordered)
               else:
                   print("puzzle over vname",vname)
        
        return


def my_args():
    '''
    Build a string without a few characters: '[],
    For inclusion in nc file history
    '''
    sarg=str(sys.argv[1:]).replace("'"," ").replace(' , ','').replace('[','').replace(']','')
    return sarg
if __name__ == '__main__':

    """ Read in files and aggregate each variable into a separate file
    :param arg1 input NetCDF file(s)
    :param arg2 option choice of NetCDF3_CLASSIC or NetCD4_CLASSIC
    :param arg4 option choice of variable
    :param arg5 option: filename with list of required variables
    """
    parser = argparse.ArgumentParser(description=
             "Aggregating netcdf files, so each output file holds one variable at all modelled times")
    parser.add_argument("--input","-i",help="input NC file(s)",nargs='+', required=True)
    parser.add_argument("--listVars","-l", help="file with list of nc variable names to keep",
                        required=False)
    parser.add_argument("--var","-v", help="nc variable name to select",required=False)
    parser.add_argument("--prefix","-p", help="prefix for each file, e.g. runid-apm-...",required=False)
    parser.add_argument("--format","-f",help="NetCDF format, NETCDF4_CLASSIC is default",
                        choices=['3','4'], default="4")

    parser.add_argument("--outdir","-o",help="Output directory. Default is current working dir",default='.')

    args = parser.parse_args()
    inFiles=sorted(args.input) # sort files in lexographic order. *Likely* in time order for FAMOUS sims

    if args.listVars:
          vlist=[line.strip() for line in open(args.listVars)]
          allVars=False
    elif args.var:
          vlist=[args.var]
          allVars=False
    else:
          vlist=[]
          allVars=True

    if args.format =='3':
          outfmt="NETCDF3_CLASSIC"
    else:
         outfmt="NETCDF4_CLASSIC"
#    print "Writing files using ",outfmt

    outdir= pathlib.Path(args.outdir)
    outdir.mkdir(parents=True,exist_ok=True) # create directory (and path to it) if necessary
    
    doutAll={}
    tindex=0

    file_prefix=None
    if  args.prefix:
       file_prefix=args.prefix 
    din=Dataset(inFiles[0],"r")
    if allVars:     # get all variables with STASH codes
        vvlist=din.get_variables_by_attributes(um_stash_source= lambda v:v is not None)
        for v in vvlist:
            vlist.append(v.name)

         # set up dimensions and output data structures to receive copies
    for vv in vlist:
        doutAll[vv]=MyNCvar(vv, din, outfmt,outdir=outdir,prefix_name=file_prefix)
    din.close()   

         # read each file in turn - assuming time ordered!


    for infile in inFiles:
        # not showing all....ed withut LF. print '\r{0}'.format( infile),
        din=Dataset(infile,"r")

        my_time=din.variables[TIME_NAME][0]
        my_tbounds=din.variables[TBOUNDS_NAME][0,:]
        if din[TIME_NAME].size!=1:
            print("puzzled with",infile, "time size is",din[TIME_NAME].size)
                                          #if ever generalised to multiple times in input,
                                          # need to modify code
            sys.exit(1)
                        
        for vname in vlist:
#            print("        variable ",vname)
            varin=din.variables[vname]
            outname=vname
            doutAll[vname].add_data(outname,varin,my_time, my_tbounds,tindex)

        tindex+=1
        din.close()

    for vv in vlist:
        doutAll[vv].write_file(vv)

    print("Done!")




