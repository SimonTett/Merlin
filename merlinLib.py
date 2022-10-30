#
"""
useful functions and variables for MERLIN analysis


See options below -- but in particular _merlin_root_dir will need to be setup.
see file_path (below) for how this is done and read_data/aggload which does the data processing/caching.
They will likely need modification for your setup.
In using this it is also worth modifying lookup.csv & references (variable defined in here).
For some cases you will also need a land/sea mask -- see land_frac_path (below) for where this is.

"""
import functools
import pathlib
import copy
import cf_units
import iris
import iris.analysis
import iris.coord_categorisation
import iris.exceptions
import iris.coords
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import xarray # to support conversion from iris cubes to xarray DataArray's
# and turn off all the annoying warning messages from iris.
import warnings
warnings.filterwarnings("ignore", message="Collapsing a non-contiguous coordinate.")
warnings.filterwarnings("ignore",message = "Collapsing spatial coordinate")
warnings.filterwarnings("ignore","Using DEFAULT_SPHERICAL_EARTH_RADIUS.")

clean_cache = False  # set True should you want to clean the cache out -- data will be overwritten. Seems to cause problems...
debug = True  # set True for messages
return_xarray = False # if True return xarray dataArrays using DataAray.from_iris()


# Values below are basic configuration variables. Set values to what needed for your data.
# these might be better in a separate file.. though then would need to modify their names!

_merlin_root_dir = pathlib.Path(r"\\csce.datastore.ed.ac.uk\csce\geos\groups\merlin\top\stett2\AggData") # agg data on datastore. Note windows path.
# cache path root
_merlin_cache_dir = pathlib.Path(r"C:\Users\stett2\OneDrive - University of Edinburgh\data\MERLIN_data") # local windows path (data can be copied here from linux)
land_frac_path = _merlin_root_dir / 'ancil' / 'qrparm.mask_frac.nc'
try:
    land_frac = iris.load_cube(str(land_frac_path)) #TODO addd local cache for when on the train...
except OSError:
    print("File not found. Using local version")
    land_frac_path = _merlin_cache_dir / 'ancil' / 'qrparm.mask_frac.nc'

# variables for human editing..
info = pd.read_csv('lookup.csv', header=0, index_col=0)
lookup = info.query("FixCO2==False") # for compatibility with old code
lookup_fix = info.query("FixCO2==True")
#references = {'Historical': 'xnyvh', 'Control': 'xoauj', 'Spinup': 'xnnqm'}  # reference experiments
references = {'Historical':'xovdz','Historical_bad': 'xnyvh','Historical_fixC':'xovdx','Control':'xovfj',
              'Control_bad': 'xoauj', 'Control_fixC':'xovfr','Spinup': 'xnnqm','Historical_bad_1750':'xnphd','Spinup_old':'xnphc'}  # reference experiments
# NOTE that Historical_bad_1750 & Spinup_old only have decadal mean data
# _fixC cases are with fixed CO2 (4.4e-4 mmr)  in the radiation scheme

start_emission = 2020 # when emissions start

fsize = [11.7, 8.3]  # A4 in inches -- curses to matplotlib for using inches..
def file_path(experiment, mean_dir, file):
    """
    Work out path to variable given experiment, variable and mean_dir
    """

    path = _merlin_root_dir / experiment / 'agg' / mean_dir / (experiment + '-' + mean_dir + '-' + file)

    if not path.exists():
        print("Path ", path, " Does not exist")
    return path


def ocn_vol(file):
    """
    Compute the volume of the ocean from first time in a file
    """
    cube = iris.load_cube(str(file))[0, ...]  # get in file and extract first time.
    cube.data = cube.data / cube.data  # make it all 1 (where not masked)
    layer_thickness = cube.coord('depth').bounds[:, 1] - cube.coord('depth').bounds[:, 0]
    grid_areas = iris.analysis.cartography.area_weights(cube)
    volume_weights = layer_thickness[:, np.newaxis, np.newaxis] * grid_areas
    vol = cube.collapsed(['depth', 'latitude', 'longitude'], iris.analysis.SUM, weights=volume_weights)
    return vol.data


def netflux(var,exper,**kwargs):
    """
    Compute the netflux
    :param var: Should be related to netflux!
    :param exper: name of experiment to use
    :param kwargs: All passed to read_data
    :return: netflux timeseries
    """
    vars = ['OLR','InSW','OutSW']
    if var == 'NetFlux':
        pass # vars are fine
    elif var == 'NetFlux_ZM':
        vars = [v+"_ZM" for v in vars]
    else:
        raise NotImplementedError(f"Not implemented {var}")

    netFlux=0
    if kwargs.get('readCube',False):
        read_fn = read_file
    else:
        read_fn = read_data
    for v in vars:
        ts=read_fn(v,exper,**kwargs)
        if v == 'InSW':
            ts = -ts

        netFlux -= ts

    netFlux.rename(var + '-' + exper)  # not a good idea to do this when caching. Not sure why
    return netFlux


def read_land(var, exper, use_cache=True,readCube=False,decadal=False):
    if var == 'LandC':
        vars = ['SoilC', 'VegC']
    elif var == 'LandC_ZM':
        vars = ['SoilC_ZM', 'VegC_ZM']
    else:
        raise NotImplementedError(f"Not implemented {var}")

    if readCube:
        read_fn = read_file
    else:
        read_fn = read_data

    cube = 0
    for v in vars:
        cube += read_fn(v, exper, use_cache=use_cache,decadal=decadal)

    cube.rename(var + '-' + exper)  # not a good idea to do this when caching. Not sure why
    return cube

def read_totalC(var, exper, use_cache=True,readCube=False,decadal=False):
    vars = ['AtmC','OcC','SoilC', 'VegC'] # vars wanted
    if var == 'TotalC':
        pass # nothing to do!
    elif var == 'TotalC_ZM':
        vars = [v+"_ZM" for v in vars] # zonal means
    else:
        raise NotImplementedError(f"Not implemented {var}")

    if readCube:
        read_fn = read_file
    else:
        read_fn = read_data

    cube = 0
    for v in vars:
        cube += read_fn(v, exper, use_cache=use_cache,decadal=decadal)

    cube.rename(var + '-' + exper)  # not a good idea to do this when caching. Not sure why
    return cube


@functools.lru_cache(maxsize=6000)
def read_forcing(var, exper, use_cache=True,readCube=False,decadal=False):
    """
    Compute forcing from CO2 using Myhre formula.
    :param var: must be forcing -- done this way to have signature that matches other reads
    :param exper: experiment name
    :param use_cache: not used if  True use the cached values
    :param readCube: If True read a cube using read_file rather than uisng read_data. later does processing.
    :return:
    """
    if var != 'Forcing':
        raise NotImplementedError(f"Not implemented {var}")

    if readCube:
        read_fn = read_file
    else:
        read_fn = read_data


    cube = read_fn('AtmCO2', exper,use_cache=use_cache,decadal=decadal)
    # get in the control values (which we are using as a reference)
    ref_atm_co2  = 4.4e-4 # PI value
    f = cube/ref_atm_co2
    f = 5.35 * iris.analysis.maths.log(f)  # Myhre formula.
    f.units=cf_units.Unit('W m-2') # forcing is W/m^2
    f.rename(var + '-' + exper)
    return f

_ocn_file = _merlin_cache_dir/'ancil'/ 'sea_water_potential_temperature.nc'
if not _ocn_file.exists(): # no local copy so copy it!
    import shutil

    shutil.copy(file_path(references['Historical'], 'opx', 'sea_water_potential_temperature.nc'),_ocn_file)

ocnVol = ocn_vol(str(_ocn_file))  # volume of ocean

mAtmCol = 98630.1 / 9.8  # mass of column (Kg) of atmosphere/m^2 -- int of pstar for one exper/time -- but model should conserve
areaEarth = 4 * np.pi * 6371229.0 ** 2  # area of earth in m^2
MoAtm = areaEarth * mAtmCol * 1e3
# Mass of atmosphere in grams - 6371229.0 = radius earth in m, 98630.1 = surface pressure in hPa
scale_Pg_ppm = (44./12. * 1e15/MoAtm)* (29.0/44)*1e6 # conversion factor for atmospheric Pg C to ppm CO2.
# first convert from Pg C to mass mixing ratio of CO2 (kg CO2/kg Air) then convert to part ratio and then to ppm.
# if var_lookup is changed then need to remove the cached data in the data directory (or at least the ones that were affected by the change
PgC = cf_units.Unit('Pg C')  # units for Pg Carbon
PgC_per_yr = cf_units.Unit('Pg/year C')
## THINK ocean carbon is mol/liter of water = 12.0107 g/liter = 12e-3 g/m^3 = 12e-18 Pg/m^2
ocn_C_scale = 12.0107e-18  # mass of one mole of C in Pg...
CO2_to_C = 12.0107 / 44.0095  # convert mass CO2 to mass C.
secs_year = 360. * 24. * 60. * 60.  # number of seconds in UM year.
ocn_dims = ['depth','longitude','latitude']

var_lookup = {
    'VAOT': dict(file='sea_water_potential_temperature.nc', collapse_dims=ocn_dims, Ocean=True),
    'OcT': dict(file='sea_water_potential_temperature.nc', Ocean=True,  method=iris.analysis.MEAN,
                collapse_dims=ocn_dims),
    'OcT_CS': dict(file='sea_water_potential_temperature.nc', Ocean=True,  method=iris.analysis.MEAN,
                  collapse_dims=['longitude']),
    'OcT_Profile': dict(file='sea_water_potential_temperature.nc', Ocean=True, method=iris.analysis.MEAN,
                    collapse_dims=['longitude','latitude']),
    'SST': dict(file='sea_water_potential_temperature.nc', Ocean=True,
                Constraint=iris.Constraint(model_level_number=1)),
    # CO2 emissions -- Kg CO2/m^2/second.
    'CO2emis': dict(file='UM_m01s00i251_vn405.0.nc', Land=True,
                    method=iris.analysis.SUM, scale=CO2_to_C * secs_year * 1e-12, units=PgC_per_yr),
    'CO2fluxAtm': dict(file='UM_m01s00i250_vn405.0.nc',
                    method=iris.analysis.SUM, scale=CO2_to_C * secs_year * 1e-12, units=PgC_per_yr),
    'AtmC': dict(file='mass_fraction_of_carbon_dioxide_in_air.nc',
                   Constraint=iris.Constraint(model_level_number=11), method=iris.analysis.SUM,
                   scale=CO2_to_C * mAtmCol * 1e-12, units=PgC),
    'AtmCO2': dict(file='mass_fraction_of_carbon_dioxide_in_air.nc',
                 Constraint=iris.Constraint(model_level_number=11)),
    'AtmCO2_Profile': dict(file='mass_fraction_of_carbon_dioxide_in_air.nc',
                  method=iris.analysis.MEAN),
    'AtmT_Profile': dict(file='air_temperature_3.nc', method=iris.analysis.MEAN,
                    collapse_dims=['longitude','latitude']),
    'Forcing':dict(func=read_forcing),
    'MLD': {'file': 'UM_m02s00i137_vn405.0.nc', 'Ocean': True},
    'SAT': dict(file='air_temperature.nc',squeeze=True),
    'Land_SAT': dict(file='air_temperature.nc',  Land=True),
    'precip': dict(file='precipitation_flux.nc'),
    'Land_precip': dict(file='precipitation_flux.nc',  Land=True),
    'Snow':dict(file='snowfall_amount.nc',land=True), # snow amount (km/m^2
    'OcC': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True, scale=ocn_C_scale,
                     method=iris.analysis.SUM, units=PgC, collapse_dims=ocn_dims),
    # data is CO2 in moles per liter??
    'OcC_CS': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True,  scale=ocn_C_scale,
                method=iris.analysis.SUM, units=PgC, collapse_dims=['depth','longitude']),
    'OcC_Profile': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True, scale=ocn_C_scale,
                   method=iris.analysis.SUM, units=PgC, collapse_dims=['latitude', 'longitude'],scale_dims=['depth']),
# for PFT how to deal with area of box that is not veg... i.e. PFTs 6-9.
# Think need to upscale by the plants in someway. (All fluxes are per m -- we want totals.)
    'VegC': dict(file='vegetation_carbon_content.nc', Land=True,
                      method=iris.analysis.SUM, scale=1e-12, units=PgC),
    'VegC_PFT': dict(file='UM_m01s19i001_vn405.0.nc',  fileScale='UM_m01s03i317_vn405.0.nc',
                     method=iris.analysis.SUM, scale=1e-12, units=PgC),
    # total veg carbon as Pg. Land Sfc values are kg/m^2
    'SoilC': dict(file='soil_carbon_content.nc', Land=True,
                       method=iris.analysis.SUM, scale=1e-12, units=PgC),

    'LandC': dict(func=read_land),
    # total veg carbon as Pg
    'TotalC':dict(func=read_totalC), # total carbon (atmosphere + ocean + land)
    'Litter': dict(file='UM_m01s19i005_vn405.0.nc', Land=True, # units fro diagnostic are Kg C/m2/year (year assumes 360 day year)
                   method=iris.analysis.SUM, scale=1e-12, units=PgC_per_yr),
    'NPP': dict(file='net_primary_productivity_of_carbon.nc', Land=True,
                scale=1e-12 * secs_year,
                method=iris.analysis.SUM, units=PgC_per_yr),  # convert to PgC/year

    'NPP_PFT': dict(file='UM_m01s19i009_vn405.0.nc', Land=True, # raw units are Pg C/year
                    scale=1e-12 ,fileScale='UM_m01s03i317_vn405.0.nc',
                    method=iris.analysis.SUM, units=PgC_per_yr),  # convert to PgC/year
    'PFT': dict(file='UM_m01s03i317_vn405.0.nc', Land=True, method=iris.analysis.MEAN,units='m^2'), # area (of land box) for each PFT

    'Litter_PFT': dict(file='UM_m01s19i004_vn405.0.nc', Land=True,
                   method=iris.analysis.SUM, scale=1e-12, units=PgC_per_yr, fileScale='UM_m01s03i317_vn405.0.nc'),


    'PlantResp': dict(file='plant_respiration_carbon_flux.nc', Land=True,
                      method=iris.analysis.SUM, scale=1e-12 * secs_year, units=PgC_per_yr),
    'SoilResp': dict(file='soil_respiration_carbon_flux.nc', Land=True,
                     method=iris.analysis.SUM, scale=1e-12 * secs_year, units=PgC_per_yr),

    'SoilResp_PFT': dict(file='soil_respiration_carbon_flux.nc', Land=True,fileScale='UM_m01s03i317_vn405.0.nc',
                     method=iris.analysis.SUM, scale=1e-12 * secs_year, units=PgC_per_yr),
    # Resp partitioned over PFTs

    'icef': dict(file='sea_ice_area_fraction.nc', Ocean=True, method=iris.analysis.SUM, scale=1e-12,units=cf_units.Unit('Gm^2')),  # as million km^2
    # ice area in milion km^2
    # ocean carbon single level diags - think they are in K mole/m2/year..
    'OcSfcFluxC': dict(file='UM_m02s30i249_vn405.0.nc',Ocean=True,squeeze=True,
                   method=iris.analysis.SUM, scale = ocn_C_scale*1000, units=PgC_per_yr), # units for fluxes appear to be Kmoles/year...
    'OcPrimFluxC': dict(file='UM_m02s30i252_vn405.0.nc', Ocean = True, squeeze=True,
                    method=iris.analysis.SUM,scale=1e-15*360,units=PgC_per_yr), # diag claims to be  g C m^{-2} day^-1

    # diag claims to be  g C m^{-2} day^-12
    'OcZooFluxC': dict(file='UM_m02s30i253_vn405.0.nc', Ocean=True, squeeze=True, # not convinced scaling is correct..
                    method=iris.analysis.SUM, scale=ocn_C_scale*1000, units=PgC_per_yr),  # diag is g C m^{-2} day^-1
    # Radiation fluxes
    'OLR': dict(file='toa_outgoing_longwave_flux.nc'),
    'InSW':dict(file='toa_incoming_shortwave_flux.nc'),
    'OutSW':dict(file='toa_outgoing_shortwave_flux.nc'),
    'NetFlux':dict(func=netflux),

    ## non conservation term.
    # Unit is mili-moles per liter per second (or 10 moles over the whole of the 10m top layer)
    'nonCons': dict(file='UM_m02s30i292_vn405.0.nc', Ocean=True, method=iris.analysis.SUM,
                    scale=ocn_C_scale * secs_year * 10, units=PgC_per_yr) # diagnostics analysis suggests that this should be scaled by around 2.5-2.7 to close the budget...
    # can not see any obvious explanation for that...

}

## add ZM to data -- this involves ading longitude to the collapse dimensions.
new_vars  ={}
# compute dlat for land and ocn. So we can convert zontal totlas to amount/degree.
_lat_land= land_frac.coord('latitude').points
dlat_land = _lat_land[0]-_lat_land[1]
_lat_ocn =iris.load_cube(str(_ocn_file)).coord('latitude').points
dlat_ocn = _lat_ocn[1]-_lat_ocn[0]

for k,v in var_lookup.items():
    if k.endswith('ZM'): # skip zonal mean keys.
        continue

    collapse = copy.copy(v.get('collapse_dims',['longitude']))
    if 'longitude' not in collapse:
        collapse.append('longitude')
    try: # remove latitude
        collapse.remove('latitude')
    except ValueError:
        pass
    vv = copy.deepcopy(v) # copy the information.
    vv['collapse_dims'] = collapse # add modify the collapse info
    # and if method is sum convert to X/degree.
    if v.get('method') is iris.analysis.SUM:
        if v.get('Ocean',False): # land
            dy=dlat_land
        else:
            dy=dlat_ocn
        scale = v.get('scale',1.0)/dy
        vv['scale']=scale
        vv['units'] = v.get('unit',cf_units.Unit(''))/'degree'

    new_vars[k+'_ZM'] = vv


var_lookup.update(new_vars)

def gen_name(exper,series=None):
    """
    Generate standard name from experiment name using lookup
    :param exper: name of experiment
    :return: standard name
    """
    if series is None:
        series = lookup.loc[exper]

    if series.Reference == 'Control' or series.Reference == 'Control_fixC':
        name = 'C'
    else:
        name = 'H'
    if series.FixCO2:
        end_name='F'
    else:
        end_name=''
    name += f"_{series.Carbon:d}C_{series.Time:d}y{end_name}"

    return name

def named_periods(series, pnames=None, len=30, emission_stop=None,constraint=False,additional=None):
    """
    Return dict of named periods
    :param series information on experiment. Only need time
    :param len Length of averaging period.
    :param emission_stop year at  which emissions stop. If None default value of start_emission+10 is used
    :param pnames List of names to be returned.
    :param constraint -- if True return iris.constraints rather than slices.
    :return: dict of names and years as slices (or as iris.constraints)
       +E len years after emissions stop
       -D len years before drawdown
       +D len years after drawdown finishes
       R: Recovery 2450 to 2450+len year.
    """
    time=series.Time
    if emission_stop is None:
        emission_stop = start_emission + 10

    periods = {'+E':slice(emission_stop, emission_stop + len), '-D':slice(emission_stop + time - (len+10), emission_stop -10 + time),
              '+D':slice(emission_stop + time , emission_stop + time  + len),
              "R":slice(2450,2450+len)}
    if additional is not None:
        periods.update(additional) # add in additional names/slices.
    if pnames is not None:
        periods = {key:value for (key,value) in periods.items() if key in pnames} # only keep pnames
    else:
        pnames = list(periods.keys())

    if constraint:
        for name in pnames:
            def gen_constraint_fn(sl): # making python closure work for us. "Obvious" way means constrain has same slice..
                fn = lambda cell: sl.start <= cell < sl.stop
                return fn
            periods[name]=iris.Constraint(year=gen_constraint_fn(periods[name]))
    return periods
def plot_emission(times=None, ax=None,overshoot=None):
    """
    plot grey shading region where emissions are happening
    :param time -- time (in years)

    """
    if ax is None:
        ax = plt.gca()

    if overshoot is None:
        overshoot = lookup.Time.unique()

    if times is None:
        times = [start_emission]
        times.extend([t + start_emission for t in overshoot])

    y = ax.get_ylim()

    if isinstance(times,list):
        for time in times:
            ax.fill_betweenx(y, time, time + 10, color='slategrey',zorder=0) # make sure it is at the back
    else:
        ax.fill_betweenx(y, times, times + 10, color='slategrey',zorder=0) # make sure it is at the back
    return

# standard linestyle, colors and labels for PFT
pft_linestyles = dict(BLT='solid', NLT='solid', C3='dashed', C4='dashed', Shrub='dotted', Urb='dotted',
                         Wat='dotted', Soil='solid', ice='dotted',Total='solid')
pft_colors = dict(BLT='lightgreen', NLT='darkgreen', C3='goldenrod', C4='orange', Shrub='limegreen', Urb='gray',
                      Wat='blue', Soil='brown', ice='skyblue',Total='black')
pft_labels = ['BLT', 'NLT', 'C3', 'C4', 'Shrub', 'Urb', 'Wat', 'Soil', 'Ice','Total']
def plot_pft(pft, ax=None,figLabel=None, labels=None,colors=None,linestyles=None,xcoord=None,pftindices=None,**kwargs):
    """
    PLot timeseries of PFT on axis (or current axis)
    pft -- timeseries of pfts
    ax (default None) axis on which to plot

    """
    # set up defaults
    if ax is None:
        ax = plt.gca()
    if labels is None:
        labels = pft_labels
    if colors is None:
        colors = pft_colors
    if linestyles is None:
        linestyles = pft_linestyles
    if xcoord is None:
        xcoord='year'
    if pftindices is None:
        pftindices = pft.coord('pseudolevel').points

    for indx in pftindices: # loop over co-ords.
        label =labels[indx-1]
        style=dict(marker='o', color=colors.get(label), label=label,linestyle=linestyles.get(label, 'solid'))
        style.update(kwargs)
        iris.plot.plot(pft.coord(xcoord), pft.extract(iris.Constraint(pseudolevel=indx)), axes=ax,**style)

    # plot the total
    """    total = pft.collapsed('pseudolevel',iris.analysis.SUM)
    label='Total'
    iris.plot.plot(pft.coord(xcoord), total, label=label,
                   color=colors.get(label),
                   linestyle=linestyle.get(label), axes=ax, **kwargs)"""

    return






def var_properties(var,exper=None,decadal=False):
    """
    Return properties for specified climate variable. Uses var_lookup
    """

    if var not in var_lookup:
        raise Exception(f"Do not know about {var}")

    lookup= var_lookup[var]
    if lookup.get('Ocean', False):
        direct = 'opy'
    else:
        direct = 'apy'

    if decadal:
        direct = direct[0:2]+'x'

    scale = lookup.get('scale', 1.0)
    volAvg = lookup.get('volAvg', False)
    file = lookup.get('file')
    constraint = lookup.get('Constraint')
    method = lookup.get('method', iris.analysis.MEAN)
    thick_scale = lookup.get('thick_scale', False)
    collapse_dims = lookup.get('collapse_dims', ['longitude', 'latitude'])  # TODO merge volAvg into this logic
    units = lookup.get('units')
    if lookup.get('Land', False):
        lf = np.squeeze(land_frac.data)
    else:
        lf = 1.0

    if lookup.get('fileScale',False):
        # decode path
        path = file_path(exper,direct,lookup['fileScale'])

        fileScale = iris.load_cube(str(path))
        addAuxCoords(fileScale) # fix co-ords

    else:
        fileScale=None

    func = var_lookup[var].get('func', None)
    squeeze = var_lookup[var].get('squeeze',False)
    scale_dims = var_lookup[var].get('scale_dims',[])
    return dict(dir=direct, volAvg=volAvg, file=file, constraint=constraint, scale=scale,
                method=method, thick_scale=thick_scale, lf=lf, func=func,
                collapse_dims=collapse_dims, units=units,scale_dims=scale_dims,fileScale=fileScale,squeeze=squeeze)


def addAuxCoords(cube):
    """
    Add useful aux coords to a cube and fix various problems
    :param cube: cube to be modified
    :return: nada as cube modified in place
    """
    cube.coord('longitude').circular = True  # make longitude circular
    try:
        iris.coord_categorisation.add_year(cube, 'time')  # add year
        cube.coord('year').var_name='year' # sometimes we have var_name and iris complains if co-ords differ in meta-data...
        # iris.coord_categorisation.add_month(cube, 'time')  # add month
        iris.coord_categorisation.add_month_number(cube, 'time')  # add month
        cube.coord('month_number').var_name='month_number' # sometimes we have var_name and iris complains if co-ords differ in meta-data...
    except (iris.exceptions.CoordinateNotFoundError, ValueError):
        pass
    for bndCoord in ['time', 'longitude', 'latitude']:
        try:
            cube.coord(bndCoord).guess_bounds()
        except (ValueError, iris.exceptions.CoordinateNotFoundError):
            pass

    # fix latitude bounds.
    try:
        cube.coord('latitude').bounds = np.clip(cube.coord('latitude').bounds, -90, 90)

    except (ValueError, iris.exceptions.CoordinateNotFoundError):
        pass
    try:
        cube.coord('pseudolevel').var_name='UM_pseudolevel'
    except (ValueError, iris.exceptions.CoordinateNotFoundError):
        pass


# ------------
def exper_properties(exper,**override):
    """
    :param exper  -- experiment  for which properties are required.
      properties (colour& marker defined from Reference, Carbon and Time
      alpha defned as 0.8 & linewidth a 2.
    kwargs override defined values.

    """

    series = lookup.loc[exper]
    p = properties(series,**override)
    return p
def properties(series, **override):
    """
    :param series -- series for which properties are required.
      properties (colour& marker defined from Reference, Carbon and Time
      alpha defined as 1 & linewidth  2.
    kwargs override defined values.

    """

    hist = {1000: 'firebrick', 500: 'red', 250: 'coral'}
    # add in the carbon values relative to the Control.
    for k,v in hist.copy().items():
        hist[k+552.5]=v
    colour_lookup = dict(Historical=hist,
                         Special={552.5:'magenta'}, # Historical Simulation.
                         Control={1000: 'blue', 500: 'cornflowerblue', 250: 'cyan'})

    marker_lookup = {50: 'd', 100: 'h', 200: '*'}  # indexed by time
    if series.FixCO2:
        linestyle='dashed'
    else:
        linestyle='solid'
    colour = colour_lookup[series.Reference.replace('_fixC','')][series.Carbon]
    # alpha = alpha_lookup[series.Carbon]
    alpha = 1
    marker = marker_lookup[series.Time]
    label = gen_name('',series=series)
    result = dict(color=colour, alpha=alpha, marker=marker, linewidth=2,linestyle=linestyle,label=label)
    result.update(**override)
    return result



def proc_all(var, plot_constraint=None, ratio=False):
    """
    Process all data for a specified variable.

    """

    refs = {}
    diffs = {}
    timeseries = {}
    for exper, series in lookup.iterrows():
        ts = read_data(var, exper)
        ref = read_data(var,references[series.Reference])
        refs[series.Reference] = ref
        t = diff(ts, ref,ratio=ratio)
        diffs[exper] = t.extract(plot_constraint)
        timeseries[exper] = ts

    return (refs, diffs, timeseries)

@functools.lru_cache()
def delta(variable,experiment,refName=None,ratio=None,meanValue=None):
    """
    Returns processed variable for specified experiment and then differenced against reference simulation.
    Caches difference.
    :arg variable -- variable you want.
    :arg experiment -- experiment wanted. If list like then will interate over
    :arg refName -- name of Reference (Control or Historical normally). If not specified worked out from lookup table
    :arg ratio -- If true compute ratio between experiment and reference.

    """

    if isinstance(experiment,str):
        if refName is None:
            refName = lookup.loc[experiment].Reference
        ref_exper = references[refName]
        ts_sim = read_data(variable,experiment)
        ts_ref = read_data(variable,ref_exper)
        if meanValue: # use mean from reference value
            ts_ref = ts_ref.collapsed('time',iris.analysis.MEAN)
            if ratio:
                delta_ts = ts_sim/ts_ref
            else:
                delta_ts = ts_sim - ts_ref
        else:
            delta_ts = diff(ts_sim,ts_ref,ratio=ratio)
        delta_ts.rename(f"{ts_sim.name()} minus {ts_ref.name().split('-')[-1]}")
    else:
        delta_ts = dict()
        for exper in experiment: # a list so iterate over it.
            delta_ts[exper]= delta(variable,exper,refName=refName)

    return delta_ts
@functools.lru_cache(maxsize=6000)
def read_file(var, exper,use_cache=False,decadal=False):
    """
    Read file for given variable and experiment. Does no processing except adding aux coords
    Takes in use_cache but ignores it.
    """

    varprops = var_properties(var,exper=exper,decadal=decadal)  # properties for variable
    func = varprops.pop('func')
    if func is not None:
        cube = func(var, exper, readCube=True,decadal=decadal)
    else:
        dir = varprops.pop('dir')
        file = varprops.pop('file')
        path = file_path(exper, dir, file)  # path
        cube = iris.load_cube(str(path))  # read the  data

    addAuxCoords(cube)
    if return_xarray:
        cube = xarray.DataArray.from_iris(cube)
    return cube

topics=lambda  cell: -30 <= cell < 30.
named_regions=dict(Amazonia=iris.Constraint(longitude=lambda  cell: 210. <= cell < 330. , latitude=topics),
                   Africa=iris.Constraint(longitude=lambda  cell: (0 <= cell < 90.) | (330 < cell <= 360), latitude=topics),
                   Indonesia=iris.Constraint(longitude=lambda  cell: 90 <= cell < 210, latitude=topics)
                   )
@functools.lru_cache(maxsize=6000)
def read_data(var, exper, use_cache=None,decadal=False,region=None):
    """
    Read and process data for variable and experiment. Note value is cached which means all parameterers must be
       allowed keys for dict.
    :param var -- variable wanted. See var_lookup for what is allowed
    :param exper -- UM experiment name
    :param use_cache (optional -- default is None) use the cache -- though overwritten by clean_cache
       If True -- use the cache
       if False -- do not use the cache
       If None -- set True if cache newer than file; False if not
    :param regions -- list of names regions or single string. regions must exist in merlinLib.named_regions
    example: zm = merlinLib.read_data('SAT_ZM','xovdz')
    """
    # various magic things..


    varprops = var_properties(var,exper=exper,decadal=decadal)  # properties for variable
    func = varprops.pop('func')
    if func is not None:
        cube = func(var, exper, use_cache=use_cache,decadal=decadal)
        if return_xarray and (type(cube) == iris.cube.Cube):
            cube = xarray.DataArray.from_iris(cube)
        return cube

    direct = varprops.pop('dir')
    file = varprops.pop('file')
    path = file_path(exper, direct, file)  # path to raw file.
    proc_dir = _merlin_cache_dir/'data_cache' # where cached data will live
    if decadal:
        proc_dir = proc_dir /'decadal'

    proc_dir.mkdir(parents=True,exist_ok=True)
    if region is not None:
        region_str = region+"_"
        varprops['region']=region
    else:
        region_str = ''

    proc_path = proc_dir / (region_str+"_".join((var, exper)) + ".nc")
    if use_cache is None:
        use_cache = proc_path.exists() and (proc_path.stat().st_mtime > path.stat().st_mtime)

    if (not clean_cache) and use_cache and proc_path.exists():
        if debug:  print("Using cache at ", proc_path)
        cube = iris.load_cube(str(proc_path))  # read the processed data

    else:
        if debug:  print("Computing using ", path, varprops)
        if not path.exists():
            raise OSError(f"{path} does not exist for {var} in {exper}")
        cube = aggLoad(str(path), **varprops)

        iris.save(cube, str(proc_path))  # save the data

    cube.rename(region_str+var + '-' + exper) # not a good idea to do this when caching. Not sure why
    if return_xarray:
        cube = xarray.DataArray.from_iris(cube)
    return cube


def aggLoad(filepath, volAvg=False, method=iris.analysis.MEAN, constraint=None, lf=1,
            scale=None, thick_scale=False, collapse_dims=['longitude', 'latitude'], units=None,
            collapse=True,fileScale=None, squeeze=False,scale_dims=[],weighting=True,return_weights=False,
            region=None):  # default method is global mean

    """
    load and do some processing on data
    """
    import iris
    if thick_scale:
        raise NotImplementedError("thick_Scale is no longer implemented -- use scale_dims")

    cube = iris.load_cube(filepath)
    addAuxCoords(cube)
    if fileScale:

        ## evil hack as PFT areas all types but other variables not nesc.. so work out intersection on coord 1..
        cnames = [c.name() for c in cube.coords()]
        fnames = [c.name() for c in fileScale.coords()]
        if 'pseudolevel' in cnames:
            print("Doing evil hack for pseudonames...")
            cube = cube * fileScale[:,0:5,:,:] # should fail if co-ords do not match.

        elif 'pseudolevel' in fnames:
            # iris does not  broadcast while xarray does... So will use that for the hard work!
            import xarray
            cube = (xarray.DataArray.from_iris(cube)*xarray.DataArray.from_iris(fileScale)).to_iris()
            # need to fix things so that order is long/lat/time/ps
            cube.transpose([0,3,1,2])
            addAuxCoords(cube) # put the aux cords back in.,.
            print("used xarray to compute product... ")
        else:
            cube = cube*fileScale

    if region is not None:
        import iris.coords
        try:
            cube=cube.extract(named_regions[region])
        except KeyError:
            print("Allowed regions are ",list(named_regions.keys()))
            raise

        # add singular co-ordinate.
        coord = iris.coords.DimCoord([list(named_regions.keys()).index(region)],standard_name='region')
        cube=iris.util.new_axis(cube)
        cube.add_dim_coord(coord,0)

    cube = cube.extract(constraint) # then apply the constraint for the generic type


    try:
        weights = iris.analysis.cartography.area_weights(cube) * lf.filled(0.0)  # assuming long/lat here. Using lf.data to get rid of missing data
    except AttributeError: # not numpy masked array
        weights = iris.analysis.cartography.area_weights(cube) * lf
    dims = collapse_dims[:]
    if volAvg: # add depth to the dimensions to collapse over.
        dims.append('depth')

    if ('depth' in dims) or ('depth' in scale_dims):
        layer_thickness = cube.coord('depth').bounds[:, 1] - cube.coord('depth').bounds[:, 0]
        weights = layer_thickness[:, np.newaxis, np.newaxis] * weights


    if collapse:
        mod_global = cube.collapsed(dims, method, weights=weights)
    else:
        mod_global = cube.copy()
        if weighting:
            mod_global.data *= weights

    if scale is not None:
        mod_global *= scale

    if units is not None:
        mod_global.units = units

    if squeeze:
        mod_global = iris.util.squeeze(mod_global)

    if return_weights: # be good to conver weights to a cube...
        # convert weights into a cube
        wt = mod_global[0,...]
        wt.units= None
        wt.data = weights[0,...]
        # drop the unneded co-ords
        for c in ['time','year','month_number']:
            try:
                wt.remove_coord(c)
            except iris.exceptions.CoordinateNotFoundError:
                print(f"Failed to remove coord: {c}")
        return mod_global, wt
    else:
        return mod_global

# interpolate values
def interp_by_coord(var,ref,coord='year'):
    """Interpolate var to ref times and add (back in ) bounds as iris interpolate does not do this. Sigh
    :param var -- variable to interpolated
    :param ref -- reference variable that provides co-ord values
    :param coord (optional -- default 'time') coord to interpoalte
    """

    # deal with dataset/iris cube.

    try:
        try:
            result = var.interpolate([(coord, ref.coord(coord).points)], iris.analysis.Linear())

        except ValueError:
            return var # no appropriate dim so just return it.
        result.coord(coord).guess_bounds()
    except AttributeError: # hope this is an xarray

        result = var.interp({coord:ref[coord]})

    return result

def common_period(cube_list,coord='year'):
    """
    Constrain everything to the same common period
    :param ts_list: list of timeseries
    :return: list of timeseries constrained to max of min and min of max
    """
    minC=[]
    maxC=[]
    for cube in cube_list:
        coordV = cube.coord(coord).points
        minC.append(coordV.min())
        maxC.append(coordV.max())

    minC = np.max(minC)
    maxC = np.min(maxC)
    cfunc = lambda pt: minC <= pt <= maxC
    constraint = iris.Constraint(coord_values={coord:cfunc})
    cons_cubes = [cube.extract(constraint) for cube in cube_list]
    return cons_cubes


def diff(ts, ref, ratio=False):
    """
    Compute difference between ts and reference field.
    Both fields are extracted to common minimum period
    Reference field is interpolated to time coords of ts field
    :param ts -- field 1
    :param ref -- field 2

    """

    ts_cons,ref_cons = common_period([ts,ref])
    interp = interp_by_coord(ref_cons, ts_cons)
    try:
        if ratio:
            diff = ts_cons/interp
        else:
            diff = ts_cons-interp
    except ValueError:  # iris sucks and throws error for this with no useful explanation...
        diff = ts_cons.copy()
        if ratio:
            diff.data = ts_cons.data / interp.data
        else:
            diff.data = ts_cons.data - interp.data

    return diff


def std_error(cube, window=None):
    """
    Compute SD for time-filtered data
    :param cube -- cube (of reference data) to compute sd
    :param window -- window length -- if specified
    """

    c = cube
    if window is not None:
        c = c.rolling_window('time', iris.analysis.MEAN, window)
    # fit 2nd order polynomial.
    ZeroRef_poly, resid_var, rank, singular_values, rcond = np.polyfit(c.coord('time').points, c.data, 2, full=True)
    # resid_var is sum of squares. Need it as an average...
    var = resid_var / c.data.shape[0]  # assuming time is 1st dim
    sd = np.sqrt(var)
    return sd


@functools.lru_cache(maxsize=6000)
def read_cube(var, exper,  weight=True,raw=False,return_weights=False,**kwargs):
    """
    Read in specified cube and apply standard things to it but no collapse
    :param var -- variable to read in,
    :param exper -- name of experiment to read in
    :param weight (default True) -- I don't think does anything... so weighting applied...
    :param raw (default False) -- if True then raw data is returned,
    any other kwargs are ignored

    """
    varprops = var_properties(var,exper=exper)  # properties for variable
    func = varprops.pop('func')
    if func is not None:
        cube = func(var,exper,readCube=True)
    else:
        dir = varprops.pop('dir')
        file = varprops.pop('file')
        path = file_path(exper, dir, file)  # path
        if raw:
            if return_weights:
                raise ValueError("Can not set return_weights with raw")
            cube = iris.load_cube(str(path)) # just read the cube.

        else:
            cube = aggLoad(str(path),collapse=False,weighting=weight,return_weights=return_weights,**varprops)
            if return_weights:
                (cube, weights) = cube
    if return_xarray:
        cube = xarray.DataArray.from_iris(cube)
    if return_weights:
        return cube,weights
    else:
        return cube


## make a label object
class plotLabel:
    """
    Class for plotting labels on sub-plots
    """

    def __init__(self, upper=False, roman=False):
        """
        Make instance of plotLabel class
        parameters:
        :param upper -- labels in upper case if True
        :param roman -- labels use roman numbers if True
        """

        import string
        if roman:  # roman numerals
            strings = ['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii', 'ix', 'x', 'xi', 'xii']
        else:
            strings = [x for x in string.ascii_lowercase]

        if upper:  # upper case if requested
            strings = [x.upper() for x in strings]

        self.strings = strings
        self.num = 0

    def label_str(self):
        """
        Return the next label
        """
        string = self.strings[self.num] + " )"
        self.num += 1
        self.num = self.num % len(self.strings)
        return string

    def plot(self, ax=None, where=None):
        """
        Plot the label on the current axis
        """

        if ax is None:
            plt_axis = plt.gca()
        else:
            plt_axis = ax


        if where is None:
            x = -0.03
            y = 1.03
        else:
            (x, y) = where

        # try and iterate over axis
        try:
            axis_to_use = plt_axis.flatten()
        except  AttributeError: # no flatten method so assume single axis
            axis_to_use = [plt_axis]
        for ax in axis_to_use:
            text = self.label_str()
            ax.text(x, y, text, transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='bottom')


def saveFig(fig, name=None, savedir=pathlib.Path("figures"), figtype=None, dpi=None):
    """
    :param fig -- figure to save
    :param name (optional) set to None if undefined
    :param savedir (optional) directory as a pathlib.Path to save figure to . Default is figures path
    :param figtype (optional) type of figure. (If nto specified then png will
    """

    defFigType = '.png'
    if dpi is None:
        dpi = 300
    # set up defaults
    if figtype is None:
        figtype = defFigType
    # work out sub_plot_name.
    if name is None:
        fig_name = fig.get_label()
    else:
        fig_name = name

    outFileName = savedir / (fig_name + figtype)
    fig.savefig(outFileName, dpi=dpi)


def latitudeLabel(value, pos):
    """
    :param values -- value of label
    """

    deg = r'$^\circ$'  # what we need for a degree symbol
    if value < 0:
        end = deg + 'S'
    elif value > 0:
        end = deg + 'N'
    else:
        end = ''

    c = mpl.ticker.ScalarFormatter()
    c.set_scientific(False)
    str = c.format_data_short(abs(value)).strip()  # trailing strip removes whitespace
    str += end
    if abs(value) < 1e-6:
        str = 'Eq'

    return str


lat_format = FuncFormatter(latitudeLabel)
