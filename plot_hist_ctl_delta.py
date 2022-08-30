"""
Plot relevant differences between hist and ctl
"""
import merlinLib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import iris.quickplot as qplt
import iris
import numpy as np
import iris.plot # gives plotting and brewer colour tables


def delta_h_c(var,exper1,exper2,cumsum=False,smooth=None):
    cubes = []

    for exper in [exper1,exper2]:

        series = merlinLib.lookup.loc[exper]
        ref=merlinLib.references[series.Reference]

        ref_cube = merlinLib.read_data(var,ref)
        cube = merlinLib.read_data(var,exper)

        diff = merlinLib.diff(cube, ref_cube)
        cubes.append(diff)
    delta = merlinLib.diff(cubes[0],cubes[1])
    if smooth is not None:
        delta = delta.rolling_window('time',iris.analysis.MEAN,smooth)
    if cumsum:
        d2 = delta.copy()
        d2.data = np.cumsum(d2.data,axis=0)
        delta = d2
    delta.long_name='H-R diff '+var
    print("delta_h_c ",var)
    return delta

## read in the data
levels=np.array([-10,-5,-2,-1,1,2,5,10])/2.
levels=[-60,-40,-20,-10,-5,5,10,20,40,60]
levels=np.linspace(-6,6,13)
levels = np.array([-6,-4,-2,-1,-0.5,-0.25,0,0.25, 0.5,1,2,4,6])
smooth=None
cumsum=False

vars = [ 'VegC_ZM','SoilC_ZM','LandC_ZM','OcC_ZM','OcC_CS','OcnT_CS','Land_SAT_ZM']

titles=['Vegetation Carbon (PgC/deg)','Soil Carbon (PgC/deg)','Ocean Carbon (Pg/C)']

vars_fluxes = [ 'NPP_ZM','Litter_ZM','SoilResp_ZM','OcSfcFluxC_ZM','OcPrimFluxC_ZM']
vars_PFT = ['PFT_ZM','NPP_PFT_ZM','Litter_PFT_ZM','SoilResp_PFT_ZM','VegC_PFT_ZM'] #
c = ['soil','veg','land']
carbon_colors=dict(soil='brown',veg='green',land='red',ocean='blue') # color lookups
#titles=['NPP (PgC/yr deg)','Litter (PgC/yr deg)','SoilResp(Pg/ yrC)']
hist=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Historical"').index[0]
ctl=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Control"').index[0]
delta=dict()

def add_soil(cube):
    """
    add soil to a PFT cube as pseudolevel=6 with value = 0.0.

    :param cube: Cube for which info to be added
    :return: modified cube
    """
    pl = cube.coord('pseudolevel')
    if 6 not in pl.points:  # need to add in soil and above... (6, 7,8, 9)
        extra_pft = cube.extract(iris.Constraint(pseudolevel= lambda l: l in [1,2,3,4]))  # grab four values
        extra_pft.data = extra_pft.data * 0.0  # set all values to zero but preserve mask
        extra_pft.coord('pseudolevel').points = np.arange(6,10)

        result = iris.cube.CubeList([cube, extra_pft]).concatenate_cube()
        return  result

    else:
        return cube # no change
for var in vars + vars_fluxes+vars_PFT:
    c=delta_h_c(var,hist,ctl,smooth=smooth,cumsum=cumsum)
    if 'PFT' in var: # A PFT variable then add in values for non plant sfc types if not present
        try:
           c = add_soil(c)
        except iris.exceptions.CoordinateNotFoundError:
            pass
    delta[var] =c
    if var not in ['SAT_ZM','PFT_ZM','PFT']:
        y = c.coord('latitude').points
        dy = (y.max() - y.min()) / (len(y) - 1)
        delta[var].data = delta[var].data/dy # Pg/degree.

    # now add the t/s to the data.
    if var.endswith('_ZM'):
        var2 = '_'.join(var.split('_')[0:-1])
        c2 =  delta_h_c(var2,hist,ctl,smooth=smooth,cumsum=cumsum)
        if 'PFT' in var:  # A PFT variable then add in values for bare soil if not present
            try:
                c2 = add_soil(c2)
            except iris.exceptions.CoordinateNotFoundError:
                pass
        delta[var2]=c2



varsRead = ['VegC','SoilC','OcC', 'LandC','CO2atm','nonCons','CO2emis']
ts_diff=dict()
for var in varsRead:

    ts_diff[var]=delta_h_c(var,hist,ctl)
# budget delta
vars_sum = ['VegC','SoilC','OcC', 'CO2atm']
total_c=ts_diff[vars_sum[0]].copy()
for v in vars_sum:
    total_c.data += ts_diff[v].data
total_c.data += np.cumsum(ts_diff['nonCons'].data)
## plot everything
cmap = plt.get_cmap('brewer_Spectral_11').reversed()# want reversed version 'brewer_Spectral_11_r'

levels=levels[levels !=0]
nvars = len(titles)
fig, ax_time_flux=plt.subplots(nvars + 1, 1, num='Change', figsize=merlinLib.fsize, clear=True)

#fig,ax_time_flux = plt.subplots(nrows=2,ncols=1,num='Cchange',figsize=[8.3,11.7],clear=True)

c=['year','latitude']

for var in varsRead:
    t=ts_diff[var]
    iris.plot.plot(t.coord('year'), t, axes=ax_time_flux[0], linewidth=2, label=var)
ax_time_flux[0].legend(loc='upper right', ncol=2)
ax_time_flux[0].axhline(0.0, color='black', linestyle='dashed')
for ax,var,title in zip(ax_time_flux[1:], vars, titles):
    d=delta[var]
    cf=iris.plot.contourf(d,coords=c,axes=ax,levels=levels,cmap=cmap)
    cs=iris.plot.contour(d,coords=c,axes=ax,levels=levels,colors='black',linewidths=1)
    ax.clabel(cs,colors='black',fmt='%3.0f')
    #ax.coastlines()
    ax.set_title(title)
    ax.set_xlabel('Year')
    ax.set_ylabel('Latitude')
    ax.set_yticks([-60,-30,0,30,60])

    ax.yaxis.set_major_formatter(merlinLib.lat_format)
# put color bar on
label = merlinLib.plotLabel()
for ax in ax_time_flux:
    ax.get_xaxis().set_visible(False)
    label.plot(ax)
    #ax.label_outer()
    for yr in [2020,2030,2220,2230]:
        ax.axvline(yr,linestyle='dashed',color='black')

cbar = fig.colorbar(cf, ax=ax_time_flux, ticks=levels, boundaries=levels, orientation='horizontal',
                    extend='both', extendfrac='auto', drawedges=True, fraction=0.1)
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)
##






# Compute fluxes for land, soil & veg
def resid(PFT,box):
    total = PFT.collapsed('pseudolevel',iris.analysis.SUM)
    resid = total - box
    return resid



for end in ['','_PFT','_ZM','_PFT_ZM']:
    delta['veg'+end]= delta['NPP'+end]-delta['Litter'+end]
    delta['soil'+end]= delta['Litter'+end]-delta['SoilResp'+end]
    delta['land'+end] = delta['NPP'+end]-delta['SoilResp'+end]

min_yr = 2030
max_yr = 2080
cons = iris.Constraint(year=lambda yr: min_yr < yr < max_yr)

def comp_derived(var,delta,cons):
    """
    Compute derived quantities for PFT variables.
    :param var -- variable name (i.e. Litter)
    :param delta -- dict of deltas
    :param cons -- time to extract to
    """

    def comp_resid(PFT, box):
        total = PFT.collapsed('pseudolevel', iris.analysis.SUM)
        resid = total - box
        return resid

    gm = delta[var]
    zm = delta[var+'_ZM']
    gm_PFT = delta[var+'_PFT']
    zm_PFT = delta[var+'_PFT_ZM']
    resid = comp_resid(gm_PFT,gm)
    resid_zm = comp_resid(zm_PFT,zm)
    values = gm, zm,gm_PFT,zm_PFT,resid,resid_zm
    result=[]
    for v in values:
        result.append(v.extract(cons).collapsed('time',iris.analysis.MEAN))
    return result

PFT = delta['PFT'].extract(cons).collapsed('time',iris.analysis.MEAN)*100
PFT_ZM = delta['PFT_ZM'].extract(cons).collapsed('time',iris.analysis.MEAN)

NPP, NPP_ZM, NPP_PFT, NPP_PFT_ZM, NPP_resid, NPP_resid_ZM=comp_derived('NPP',delta,cons)
Litter, Litter_ZM, Litter_PFT, Litter_PFT_ZM, Litter_resid, Litter_resid_ZM=comp_derived('Litter',delta,cons)
veg, veg_ZM, veg_PFT, veg_PFT_ZM, veg_resid, veg_resid_ZM=comp_derived('veg',delta,cons)
soil, soil_ZM, soil_PFT, soil_PFT_ZM, soil_resid, soil_resid_ZM=comp_derived('soil',delta,cons)


# now to plot.
## plot zonal means from 2020-2100
colors = dict(BLT='green', NLT='darkgreen', C3='goldenrod', C4='orange', Shrub='palegreen', Urb='gray',
              Wat='blue', Soil='brown', ice='skyblue',Resid='blue',Total='grey')
linestyles = dict(BLT='solid', NLT='dashed', C3='solid', C4='dashed', Shrub='solid', Urb='dotted',
                 Wat='dotted', Soil='solid', ice='dotted')
pft_indx = {k:ind for ind,k in enumerate(colors.keys())}

figZM, axes = plt.subplots(3,1, num='ZM_delta_sig', squeeze=False, figsize=merlinLib.fsize[::-1],clear=True,sharex=True,constrained_layout=True)
(ax_land, ax_veg_pft, ax_soil_pft)=axes.flatten()
# show
# 1) Land C fluxes -- NPP, Litter & Soil resp
# 2) Veg flux (NPP-Litter) per PFT
# 3 Soil Flux (Litter - SoilResp) per PFT

# plot Oceanm, Resp, Litter and NPP zonal means.
for var in ['OcSfcFluxC_ZM','NPP_ZM','Litter_ZM','SoilResp_ZM']:
    zm = delta[var].extract(cons).collapsed('time', iris.analysis.MEAN)
    var2 = '_'.join(var.split('_')[0:-1])
    t = delta[var2].extract(cons).collapsed('time', iris.analysis.MEAN)
    label = f"{var.split('_')[0]} {t.data:3.1f} Pg/yr"

    iris.quickplot.plot(zm, axes=ax_land, linewidth=2, marker='o', label=label)

for name in ['BLT','NLT','C3','C4','Shrub','Soil']:
    ind = pft_indx[name]
    style=dict(linewidth=2, marker='o', color=colors[name],linestyle=linestyles.get(name,'solid'))

    label = f'{veg_PFT[ind].data:3.1f} Pg/yr'
    label_litter = f'{soil_PFT[ind].data:3.1f} Pg/yr'
    iris.quickplot.plot(veg_PFT_ZM[ind, :], label=label, axes=ax_veg_pft, **style)
    iris.quickplot.plot(soil_PFT_ZM[ind, :], label=label_litter, axes=ax_soil_pft, **style)

# plot resids
style['color'] = colors['Resid']
label = f'{veg_resid.data:3.1f} Pg/yr'
label2 = f'{soil_resid.data:3.1f} Pg/yr'
iris.quickplot.plot(veg_resid_ZM, axes=ax_veg_pft, label=label, **style)
iris.quickplot.plot(soil_resid_ZM, axes=ax_soil_pft, label=label2, **style)

# plot total
style['color'] = colors['Total']
label = f' {veg.data:3.1f} Pg/yr'
label2 = f' {soil.data:3.1f} Pg/yr'
iris.quickplot.plot(veg_ZM, axes=ax_veg_pft, label=label, **style,zorder=-10)
iris.quickplot.plot(soil_ZM, axes=ax_soil_pft, label=label2, **style,zorder=-10)

plt_label = merlinLib.plotLabel()
for ax,title in zip((ax_land, ax_veg_pft, ax_soil_pft), ['C fluxes', 'Veg C Flux by PFT', 'Soil C Flux by PFT']):
    plt_label.plot(ax)
    ax.legend(ncol=4,fontsize='small',handlelength=3,loc='lower left')
    ax.axhline(color='black', linestyle='dashed')
    ax.set_title(title)
    ax.set_ylabel('Pg C/(yr deg)')
    ax.set_xlabel('Latitude')
    ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.xaxis.set_major_formatter(merlinLib.lat_format)

ax_veg_pft.set_ylim(ax_soil_pft.get_ylim())


# add legends to the whole figure.
import matplotlib.lines as mlines
lines=[]
for name in ['BLT','NLT','C3','C4','Shrub','Soil','Resid','Total']:
    style = dict(linewidth=2, marker='o', color=colors[name], label=name,linestyle=linestyles.get(name, 'solid'))
    lines.append(mlines.Line2D([],[],**style))
legend = figZM.legend(ncol=5,handles=lines,handlelength=3,bbox_to_anchor=[0.2,0.6],loc='lower left') # hand placed -- probably will need adjustment

figZM.tight_layout()
figZM.show()

## some extra plots.
figExtra, axes = plt.subplots(2, 3, num='extraZM', squeeze=False, figsize=merlinLib.fsize,clear=True)
(ax_time_flux, ax_fluxes,ax_pft,ax_temp,ax_landC,ax_spare)=axes.flatten()
# plot timeseries of NPP, Litter & Soil Resp for two pulse cases
plot_constraint = iris.Constraint(year=lambda cell: cell >= 2010)
subSet = merlinLib.lookup.query('Carbon==1000 & Time==200')
vars = ['OcSfcFluxC','Litter','NPP','SoilResp']
flux_linestyles = ['dotted','solid','dashed','-.','--']
for var,linestyle in zip(vars, flux_linestyles):
    for exper,series in subSet.iterrows():
        ts=merlinLib.read_data(var,exper).extract(plot_constraint)
        ref = merlinLib.read_data(var,merlinLib.references[series.Reference])
        diff = merlinLib.diff(ts,ref)
        #diff.data = np.cumsum(diff.data)
        diff = diff.rolling_window('year',iris.analysis.MEAN,21)
        style = merlinLib.properties(series,linestyle=linestyle,marker='None')
        if series.Reference=='Historical':
            label = var
        else:
            label = None
        iris.plot.plot(diff.coord('year'),diff,axes=ax_time_flux, label=label,markevery=20,**style)
# add on the min and max yrs
for y in [min_yr,max_yr]:
    ax_time_flux.axvline(y,linestyle='dashed',linewidth=2,color='black')
# plot start and end of emissions
start_emission = 2020
y=ax_time_flux.get_ylim()
ax_time_flux.fill_betweenx(y, start_emission, start_emission + 10, alpha=0.6, color='grey')
ax_time_flux.fill_betweenx(y, 200 + start_emission, 200 + start_emission + 10, alpha=0.4, color='grey')
ax_time_flux.set_title('GM C Fluxes')
ax_time_flux.set_xlabel('Time')
ax_time_flux.set_ylabel('Pg C/Year')
ax_time_flux.legend(ncol=1)
ax_time_flux.axhline(0,linestyle='dashed',color='black',linewidth=2)

colours = [carbon_colors[name] for name in ['ocean','soil','veg']]
for var,colour in zip(['OcSfcFluxC_ZM','soil_ZM','veg_ZM'],colours):
    if var is 'land_ZM':
        linewidth=3
    else:
        linewidth =1.5
    zm=delta[var].extract(cons).collapsed('time',iris.analysis.MEAN)
    var2 = '_'.join(var.split('_')[0:-1])
    t = delta[var2].extract(cons).collapsed('time',iris.analysis.MEAN)
    label = f"{var.split('_')[0]} {t.data:3.1f} Pg/yr"

    iris.quickplot.plot(zm,axes=ax_fluxes,linewidth=linewidth,marker='o',label=label)

for name in ['BLT','NLT','C3','C4','Shrub','Soil']:
    ind = pft_indx[name]
    style=dict(linewidth=2, marker='o', color=colors[name],linestyle=linestyles.get(name,'solid'))
    label = f'{PFT[ind].data:3.1f} '
    iris.quickplot.plot(PFT_ZM[ind,:],label=label,axes=ax_pft,**style)

# plot the land temp diff
land_SAT_ZM= delta['Land_SAT_ZM'].extract(cons).collapsed('time',iris.analysis.MEAN)
iris.quickplot.plot(land_SAT_ZM,axes=ax_temp)
# plot the land C for each case
labels = ['Soil','Veg','Land']
for var,colour,label_var in zip(['SoilC_ZM','VegC_ZM','LandC_ZM'],[carbon_colors[name] for name in ['soil','veg','land']],labels):
    for exper,label,linestyle in zip([hist,ctl],['H_'+label_var,'C_'+label_var],['solid','dashed']):
        zm = merlinLib.read_data(var,exper).extract(cons).collapsed('time',iris.analysis.MEAN)/5.
        ts_var = var.split('_')[0]
        ts = merlinLib.read_data(ts_var,exper).extract(cons).collapsed('time',iris.analysis.MEAN)
        label+= f" {ts.data:4.0f}"
        iris.quickplot.plot(zm,axes=ax_landC,label=label,linestyle=linestyle,linewidth=linewidth,color=colour)
ax_landC.legend()
ax_landC.set_title('Land C')
ax_landC.set_ylabel('Pg C/degree')
ax_landC.set_xlabel('Latitude')
ax_landC.set_xticks([-90, -60, -30, 0, 30, 60, 90])
ax_landC.xaxis.set_major_formatter(merlinLib.lat_format)

for ax,title in zip((ax_time_flux, ax_pft),['C Fluxes','PFT']):
    ax.legend(ncol=1,fontsize='small',handlelength=3,loc='lower left')
    ax.axhline(color='black', linestyle='dashed')
    ax.set_title(title)
    if ax is not ax_pft:
        ax.set_ylabel('Pg C/(yr deg)')
    else:
        ax.set_ylabel('PFT (%)')






figExtra.tight_layout()
figExtra.show()