#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
import numpy as np
import time
import json
from netCDF4 import Dataset
from datetime import datetime, timedelta
#from netcdftime import datetime as datetimex, date2num
from netCDF4  import  date2num, num2date



"""
utility building netcdf for an ancil file
"""

MASK_FILE="/exports/csce/eddie/geos/groups/cesd/UM/FAMOUS/ancil/qrparm.mask.nc"

LAND_AREA=1.475354e+14
GTTOKG=1.e12
CTOCO2=44./12.
SCALING=CTOCO2*GTTOKG/(360*24*60*60*LAND_AREA)
# desired emission=SCALING*input emission where input is GT of C per yr

USECOLS=['v YEARS/GAS >', 'FossilCO2', 'OtherCO2']

             # time: is hardcoded to 360 day year, as 't' 
            # define reference date  for use in t (time) attributes
            # and make time for the ancil data be 1 July in the input year.
YR1=1850 
ANCIL_DAY_IN_YEAR=180 # data are to be timed as for mid-year, 1 July.



__svnHeader__ = "$Header: https://svn-kerberos.ecdf.ed.ac.uk/repo/geos/CESD/UM_utilities/csvAncil/ncAncil.py 689 2017-11-01 13:02:49Z mjm@EASE.ED.AC.UK $"
__svnRevision__ = "$Rev: 689 $"
__svnId = "$Id: ncAncil.py 689 2017-11-01 13:02:49Z mjm@EASE.ED.AC.UK $"
__version__="0.1"


def myArgs():
    '''
    Build a string without a few characters: '[],
    For inclusion in nc file history
    '''
    sarg=str(sys.argv[1:]).replace("'"," ").replace(' , ','').replace('[','').replace(']','')
    return sarg

def findHdr(filename):
#    lookup="UNITS,GtC,GtC"
    lookup=USECOLS[0]
    with open(filename) as myFile:
      for num, line in enumerate(myFile, 0):
        if lookup in line:
#            print('found at line:', num)
            myFile.close()
            return num
    myFile.close()
    print("not found header in csv file, seeking $s\n"%lookup)
    exit()

def checkDF(df):
    if df.keys()[1]==USECOLS[1] and \
        df.keys()[2]==USECOLS[2]:
         return
    else:
       print("invalid csv? got wrong keys of")
       df.keys()
       exit()
     
    


#def extendHistory(newDS, 
if __name__ == '__main__':

    """ Read in template nc file and csv file with global emissions
        Derive global grid and write netcdf
        Note: not general purpose, but we are creating fields with
        coords={'t':timList,'depth':template.coords['depth'], \
               'latitude':template.coords['latitude'], \
               'longitude':template.coords['longitude']}
    nlat=template.dims['latitude']
    :param arg1 template NetCDF file
    :param arg2 csv file: of dates and global emissions (2 columns with GT of C per year) 
    :param arg3 json file with attributes for ancil field
    :param arg4 input field name (field1561)
    :param arg5 output field name
    :param arg6 output NetCDF file
    """
    parser = argparse.ArgumentParser(description="building netcdf for an ancil file")
    parser.add_argument("templateNCFile",help="templateNC filename")
    parser.add_argument("csvEmissions",help="csv export of the RCPemissions")
    parser.add_argument("jsonAtts",help="json file with attributes for new nc file")
    parser.add_argument("fieldRead",help="field name for template input field - sets dimension names")
    parser.add_argument("fieldWrite",help="field name for output data")
    parser.add_argument("outNCFile",help="output NC filename")

    args = parser.parse_args()
    jsonFile=args.jsonAtts
    inFieldName=args.fieldRead
    outFieldName=args.fieldWrite
    outFile=args.outNCFile
    emFile=args.csvEmissions

          # open template nc file of ancils 
    din=Dataset(args.templateNCFile,"r")
    maskData=Dataset(MASK_FILE,"r")
    mask=maskData['temp']
    dout=Dataset(outFile,"w",format='NETCDF3_CLASSIC')

          # read csv file of emissions: find where the header line is and use it:
    hdrLine=findHdr(emFile)
    df=pd.read_csv(emFile,header=hdrLine,dtype={USECOLS[0]:np.int64}, \
                   usecols=USECOLS)

         # check the header seems ok - hardcoded strings as this is not a general utility.
    checkDF(df)

    yrLabel=df.keys()[0]
    fossilLabel=df.keys()[1]
    otherCLabel=df.keys()[2]

    nAncilFields=len(df[yrLabel])
    #print 'nAncilFields', nAncilFields

                     # build list of time and emissions in format needed
          # sum and scale emissions
    df['scaled']=SCALING*(df[fossilLabel]+df[otherCLabel])
         # to help debug, write out the data in the pandas  dataframe
    df.to_csv("dump_emissions_%s"%emFile)
    tList=[]
    emissList=[]

    emissMax=max(df['scaled'])
    emissMin=min(df['scaled'])
    emissMin=min(emissMin,0.0) # masking later on so some cells are 0.

        # build list of times in terms of days since ref time defined above 

    for yr,emiss in zip(df[yrLabel], df['scaled']):

        emissList.append(emiss)

               # use 1 July as the emissions are mid-year
        #tList.append((yr-YR1)*360+ANCIL_DAY_IN_YEAR)
        tList.append(datetime(yr,7,1)) #ANCIL_DAY_IN_YEAR)

    with open(jsonFile) as jdata:
         jsonattr=json.load(jdata)
         jdata.close()

    if outFieldName in jsonattr:
         outattr= jsonattr[outFieldName]
    else:
         outattr['name']=outFieldName

              # define time as 't', dimension and then variable
    dout.createDimension('t', len(tList))
    tout= dout.createVariable('t', 'f', ['t'])
    tout.setncatts({"units":"days since %4d-01-01 00:00:00"%YR1,
                           "long_name":"t", "calendar":"360_day",
                           "standard_name":"time",
                           "time_origin":"01-JAN-%4d:00:00:00"%YR1})
    print (".... created time")

    tout[:]=date2num(tList, units=tout.units, calendar=tout.calendar)

               # copy other dimension definitions from template

    for dname in din.variables[inFieldName].dimensions[1:]:
               # dname is dimension's name - time was assumed to be 1st, and skipped.         
       the_dim=din.variables[dname]
       dout.createDimension(dname, len(the_dim))
        
       outvar = dout.createVariable(dname, the_dim.datatype, the_dim.dimensions)

                         # copy data 
       outvar[:] = the_dim[:] 
       print(".... copied %s"%dname)

                # Copy attributes
       outvar.setncatts({k: the_dim.getncattr(k) for k in the_dim.ncattrs()})

                  # create variable for new ancil

    varin=din.variables[inFieldName]
    outdata = dout.createVariable(outFieldName, varin.datatype, ('t',)+varin.dimensions[1:])
                                # bit pedantic, but allows first dimension to be called
                                # other than 't' in the input template.
    print(".... created %s"%outFieldName )
                # get sizes of dimensions

    nlat=dout.dimensions['latitude'].size
    nlon=dout.dimensions['longitude'].size
    ndepth=dout.dimensions['depth'].size

    ninfield=nlat*nlon*ndepth

                  # generate a list of fields from the emission values

    for i,(tim,em) in  enumerate(zip(tList,emissList)):
        outdata[i]=(np.array([em]*ninfield).reshape((ndepth,nlat,nlon))*mask[:])
        
        #tout[i]=tim

                   # set extra attributes: note if we hada  sparse matrix we would have
                   # filled the array with the _FillValue then added the  data
    outdata.setncatts({'valid_max':emissMax,'valid_min':emissMin,'_FillValue':2.e+20, 
                       'missing_value':2.e+20})
    outdata.setncatts({k: outattr[k] for k in outattr})

                # define global attributes

    history="%s using %s\n%s"%(time.asctime( time.localtime(time.time())),__svnHeader__,myArgs())
    dout.history=history
    dout.source="Univ. of Edinburgh, GeoSciences  MERLiN project"
    dout.use="in FAMOUS with mods to permit emissions scenarios"
    dout.contact="Viv Scott"

                # close the dataset

    dout.close()

    print("Done!")




