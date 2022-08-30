"""
Check C budgets
"""
import iris
import iris.plot
import numpy as np
import matplotlib.pyplot as plt
import merlinLib
inv_references = {v: k for k, v in merlinLib.references.items()}
total_carbon = dict()
c_emit = dict()
for exper in ['xovfd','xoaue','xoauh','xoauj','xovdc','xocyb','xocyf','xovdz','xnnqm']:
    total_c =0.0
    for var in ['LandC','OcC', 'CO2atm']: # add up the stocks
        total_c += merlinLib.read_data(var,exper)

    # add in the non-cons term
    ts= merlinLib.read_data('nonCons',exper)
    #total_c.data += np.cumsum(ts.data)

    total_carbon[exper] = total_c.copy()

    ts =(merlinLib.read_data('CO2emis',exper)).copy()
    ts.data = np.cumsum(ts.data)[:]
    c_emit[exper] = ts

## plot the data
fig,(axT,axE) =plt.subplots(nrows=2,ncols=1,sharex=True,num='total_c',clear=True,figsize=merlinLib.fsize)
for exper,total_c in total_carbon.items():
    try:
        series = merlinLib.lookup.loc[exper]
        label = f'{series.Reference} C:{series.Carbon} Sep:{series.Time} FC {total_c.data[-1]:5.1f}'
    except KeyError:
        label = f'{inv_references[exper]} SC {total_c.data[0]:5.1f} FC {total_c.data[-1]:5.1f}'

    iris.plot.plot(total_c,label=label,axes=axT)
    iris.plot.plot(c_emit[exper], label=label, axes=axE)

for ax in [axT,axE]:
    ax.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('PgC')
axT.set_title('Total C')
axE.set_title('Cum C emissions')
fig.show()
merlinLib.saveFig(fig)
## some handy fns
import cf_units

def cube_cumsum(cube):
    """
    compute cumsum for cube data
    :param cube -- cube to have cumsum applied to
    """
    result = cube.copy()
    result.data = np.cumsum(cube.data)
    # coudl change the units...
    result.units='1'
    return result

def cube_delta(cube):

    """
    Compute the difference for a cube from the start
      :param cube
    """
    delta = cube.copy()
    delta.data = cube.data-cube.data[0]
    delta.units='1'
    return delta

## check each cpt seperately.
import scipy.stats
#name='Spinup'
#name='Historical'
#name='Control'
name='new1000G'
#exper = merlinLib.references[name]
exper='xovdc'

fig_num='Check'+name
# 1) The ocean first
OcC = merlinLib.read_data('OcC', exper)
AtmOcCflux = merlinLib.read_data('CO2fluxAtm',exper)
#OcCflux = merlinLib.read_data('OcSfcFluxC',exper)
#nonCons = merlinLib.read_data('nonCons',exper)
IntAtm = cube_cumsum(AtmOcCflux) # Int flux
#IntOcn = cube_cumsum(OcCflux) # Int flux
#IntNonCons = cube_cumsum(nonCons)
DeltaOc = cube_delta(OcC) # change in C
delta_Oc = DeltaOc-IntAtm
resid = DeltaOc-IntAtm

# Look like we need to scale nonCons by 2.5ish  But why 2.5ish?
fig,(ax_ocn,ax_land,ax_soil)=plt.subplots(ncols=3,nrows=1,num=fig_num,clear=True,figsize=merlinLib.fsize)
#iris.plot.scatter(IntNonCons,delta_Oc) # this is not constant when CO2 emitted..  Suggests some small error in my Ocean C calculation.
#breakpoint()


iris.plot.plot(cube_delta(OcC), label='Delta OcC',axes=ax_ocn)
iris.plot.plot(cube_cumsum(AtmOcCflux),label='Int Atmos C flux',axes=ax_ocn)
#iris.plot.plot(cube_cumsum(OcCflux),linestyle='dashed',label='Int OcC flux',axes=ax_ocn)
#iris.plot.plot(cube_cumsum(-1*nonCons),label='Int Non Cons flux',axes=ax_ocn)
#iris.plot.plot(cube_cumsum(AtmOcCflux-nonCons),label='Int net flux',axes=ax_ocn)
iris.plot.plot(resid,label='Diag Int Non Cons',axes=ax_ocn)



# next plot net flux into land.


LandC = merlinLib.read_data('LandC', exper)
NPP = merlinLib.read_data('NPP',exper)
SoilResp = merlinLib.read_data('SoilResp',exper)
netLand = NPP - SoilResp
IntNet = cube_cumsum(netLand) # Int flux

DeltaLandC = cube_delta(LandC) # change in C

resid = DeltaLandC-IntNet

iris.plot.plot(DeltaLandC, label='Delta LandC',axes=ax_land)
iris.plot.plot(IntNet,label='Int Net Land C flux',axes=ax_land)
iris.plot.plot(resid,label='Diag Int Non Cons',axes=ax_land)
# and now the soil

SoilC=merlinLib.read_data('SoilC', exper)
Litter = merlinLib.read_data('Litter',exper)
NetSoil = Litter-SoilResp
IntNetSoil = cube_cumsum(NetSoil) # Int flux

DeltaSoilC = cube_delta(SoilC) # change in C

residSoil = DeltaSoilC-IntNetSoil
iris.plot.plot(DeltaSoilC, label='Delta SoilC',axes=ax_soil)
iris.plot.plot(IntNetSoil,label='Int Net Land C flux',axes=ax_soil)
iris.plot.plot(residSoil,label='Diag Int Non Cons',axes=ax_soil)
for ax,title in zip([ax_ocn,ax_land,ax_soil],['Ocean','Land','Soil']):
    ax.axhline()
    ax.legend()
    ax.set_title(title)

fig.show()

# having problem with Litter looking to be to about 3% small compared to NPP & SOil Resp in spinup
# and other simulations.
# Litter (on average) is too small but there is a lot of geographical varn

import cartopy.crs as ccrs
fig,axis = plt.subplots(nrows=2,ncols=3,num=f'Check_Soil_{name}',figsize=merlinLib.fsize,clear=True,subplot_kw=dict(projection=ccrs.PlateCarree()))
fig2,axis2 = plt.subplots(nrows=1,ncols=1,num=f'Check_Soil_ZM_{name}',figsize=merlinLib.fsize,clear=True)
SoilC=merlinLib.read_cube('SoilC', exper,weight=True)
deltaSoilC = (SoilC[-1,:,:]-SoilC[0,:,:])/(SoilC.shape[0]-1)
deltaSoilC.units=''
Litter = merlinLib.read_cube('Litter',exper,weight=True).collapsed(['time'],iris.analysis.MEAN)
#Litter = merlinLib.read_cube('Litter_PFT',exper).collapsed('pseudolevel',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
SoilResp = merlinLib.read_cube('SoilResp',exper,weight=True).collapsed('time',iris.analysis.MEAN)
deltaFlux = Litter - SoilResp
deltaFlux.units=''
deltaS = Litter-SoilResp
SoilResp.units=''
diagLitter = deltaSoilC+SoilResp
for ax,var,title in zip(axis.flatten(),[deltaFlux,deltaSoilC,(deltaSoilC-deltaFlux),Litter,diagLitter,SoilResp],
                        ['Int Flux','Delta Soil C','Delta','Litter','Diag Litter','SoilResp']):
    cm=iris.plot.pcolormesh(var,axes=ax)
    fig.colorbar(cm,ax=ax,orientation='horizontal',fraction=0.1)
    ax.set_title(title)
    # plot ZM
    if title in ['Litter','Diag Litter','SoilResp','Int Flux','Delta Soil C']:
        zm = var.collapsed('longitude',iris.analysis.SUM)
        lines='solid'
        if title is 'Diag Litter':
            lines='dashed'
        iris.plot.plot(zm,label=title,axes=axis2,linestyle=lines)
        print(f"{title}:  {zm.data.sum()} Pg")

axis2.legend()
fig.show()
fig2.show()
## produce scatter plot
import iris.analysis.calculus
Litter = merlinLib.read_cube('Litter',exper)
#Litter = merlinLib.read_cube('Litter_PFT',exper).collapsed('pseudolevel',iris.analysis.SUM)
SoilResp = merlinLib.read_cube('SoilResp',exper)
SoilC = merlinLib.read_cube('SoilC',exper)
SoilCD = iris.analysis.calculus.cube_delta(SoilC,'time')
fig3,(ax,ax_ts)=plt.subplots(1,2,num=f'Scatter_{name}',figsize=merlinLib.fsize,clear=True)

ax.scatter(SoilResp.data.flatten(), Litter.data.flatten(),marker='x',s=10)
ax.set_xlabel('Soil Resp (Pg/year)')
ax.set_ylabel('Litter (Pg/Year)')
ax.axhline(color='black',linestyle='dashed')
ax.axvline(color='black',linestyle='dashed')
minv = Litter.data.min()
pos = np.where(Litter.data==minv)
for var,title in zip([Litter,SoilResp,SoilCD],['litter','Soil Resp','Soil CD']):
    ax_ts.plot(var.data[:,pos[1],pos[2]],label=title)

ax_ts.axhline(color='black',linestyle='dashed')
ax_ts.legend()

fig3.show()


## thoughts

"""
Thoughts:
1) Triffid uses 365 day => might explain the global avg differences/. Nope code suggests it is assumign a 360 day.
2) Missing distrubance -- but that would increase litter - does not explain geographic variation.
3) Have -ve values of Soil Respiration & Litter. 
4) Wonder if TRIFFID has gone unstable and that is what we are seeing. If so are the results reliable?
5) Wonder if because sampling TRIFFID every 480 timesteps (starting at 480) which is I think once every 20 days... 
6) But TRIFFID runs every 10 days so not sampling all the time. In particular missing first chunk. 
TDUMPMN diagnostic -- samples every dump period and means every timestep with no offset
TTRIFID diagnostic -- samples every 480 timesteps starting at timestep 480.
Trial -- rerun one case with diagnostics every 240 timesteps and check for conservation... 
"""





