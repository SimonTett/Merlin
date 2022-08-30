"""
Plot timeseries of PFT for H_1000C_200y
"""
import iris.analysis
import numpy as np
import iris.plot
import merlinLib
import matplotlib.pyplot as plt
import matplotlib.colors
import functools
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

@functools.lru_cache()
def delta_var(var, exper, smooth=None,region=None):
    """
    Compute difference from reference for simulation
    :param var -- name of variable wanted
    :param exper:
    :return: difference timeseries.
    """

    series = merlinLib.lookup.loc[exper,:]
    ref= merlinLib.references[series.Reference]
    ts= merlinLib.read_data(var,exper,region=region)
    ref_ts = merlinLib.read_data(var,ref,region=region)
    delta = merlinLib.diff(ts,ref_ts)

    if smooth is not None:
        delta = delta.rolling_window('time',iris.analysis.MEAN,smooth)
    return  delta

hist = merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Historical"').iloc[0]
ctl = merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Control"').iloc[0]



#figure = plt.figure(num='pft_plot',clear=True,figsize=[7.5,7])
#subfigs = figure.subfigures(ncols=1,nrows=2)
#axes = subfigs[0].subplots(nrows=1,ncols=2,sharex=True,sharey=True)
figure,axes = plt.subplots(nrows=2,ncols=2,clear=True,figsize=[7.5,4],sharey='row',sharex='col',num='pft_plot')

label=merlinLib.plotLabel()
colors=dict(Amazonia='blue',Africa='purple',Indonesia='magenta')
for series,ax,ax2 in zip([ctl,hist],axes[0],axes[1]):
    pftC_delta = delta_var('VegC_PFT', series.name)  # *100
    rgn_pft_BLT=dict()
    sum=0
    p = merlinLib.named_periods(series, constraint=True, pnames=['R'])['R']
    for reg in merlinLib.named_regions.keys():
        ts = iris.util.squeeze(delta_var('VegC_PFT', series.name, region=reg).extract(iris.Constraint(pseudolevel=1)))
        mn = ts.extract(p).collapsed('time', iris.analysis.MEAN)
        sum += ts
        iris.plot.plot(ts.coord('year'),ts,axes=ax2,color=colors[reg],label=reg,linewidth=1.5)
        #ax2.plot([2250,2500],[float(mn.data)]*2,linestyle='dotted',color=colors[reg])
    soilC_delta = delta_var('SoilC', series.name)
    LandC_delta = delta_var('LandC',series.name)
    merlinLib.plot_pft(pftC_delta,ax=ax,pftindices=[1,2,3,4,5],marker='')
    # decomposition of BLT carbon
    merlinLib.plot_pft(pftC_delta, ax=ax2, pftindices=[1], marker='') # just plot the BLT on 2nd row
    resid =pftC_delta.extract(iris.Constraint(pseudolevel=1)) - sum
    iris.plot.plot(resid.coord('year'),resid,axes=ax2,color='darkgreen',label='Ex-T')
    iris.plot.plot(soilC_delta.coord('year'),soilC_delta,axes=ax,color='brown',label='Soil')
    iris.plot.plot(LandC_delta.coord('year'), LandC_delta, axes=ax, color='black', label='Land',linewidth=2)


    title=merlinLib.gen_name(series.name)
    for a in [ax,ax2]:
        a.set_title(title,fontstyle='italic')


for ax in axes.T[1][0:2]:
    ax.legend(ncol=3,fontsize='small',handletextpad=0.1,columnspacing=0.5,labelspacing=0.25,handlelength=1.2,borderaxespad=0)
for ax in axes.T[0][0:2]:
    ax.set_ylabel("$\Delta$ Reference (Pg C)")
for ax in axes[:][0:2].flatten():
    label.plot(ax)
    ax.axhline(linestyle='dashed', color='black')
    ax.locator_params(axis='y', nbins=11)
    merlinLib.plot_emission(ax=ax,overshoot=[200])

for ax in axes[0:2][1]:
    ax.set_xlabel("Year")
    ax.set_ylim(-125,175)
figure.tight_layout()
figure.show()
merlinLib.saveFig(figure)
## now plot Temp & land precip from C_1000C_200y relative to mean from start of +D.
label = merlinLib.plotLabel()
figure_zm,axes = plt.subplots(nrows=1,ncols=3,clear=True,figsize=[7.5,3],sharey='row',sharex='col',num='zm_C_1000C_200y')
colorbar_kws = dict(orientation='horizontal',fraction=0.15,pad=0.1,format='%.1f')
start_yr=merlinLib.named_periods(ctl,pnames='+D')['+D'].start+15
land_precip = delta_var('Land_precip_ZM',ctl.name,smooth=30).extract(iris.Constraint(year=lambda cell: cell >= start_yr))
land_precip *= 24*60*60 # mm/day
precip = delta_var('precip_ZM',ctl.name,smooth=30).extract(iris.Constraint(year=lambda cell: cell >= start_yr))
precip *= 24*60*60 # mm/day
sat = delta_var('SAT_ZM',ctl.name,smooth=30).extract(iris.Constraint(year=lambda cell: cell >= start_yr))

sat_lev = np.linspace(-1,1,5)
sat_lev = [-1.5,-1,-0.6,-0.2,0.2,0.6,1,1.5]
precip_lev = [-0.5,-0.3,-0.1,0.1,0.3,0.5]
for var,ax,cmap,levels in zip([sat,precip,land_precip],axes.flatten(),['RdYlBu_r','RdYlBu','RdYlBu'],
                              [sat_lev,precip_lev,precip_lev]):
    cm=iris.plot.contourf(var,coords=['year','latitude'],cmap=cmap,axes=ax,levels=levels,extend='both')
    figure_zm.colorbar(cm,ax=ax,**colorbar_kws)


# add lat stuff
from cartopy.mpl.ticker import LatitudeFormatter
for ax,title in zip(axes,['Temperature (K)','Precip (mm/day)','Land P. (mm/day)']):
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    label.plot(ax)
    ax.set_title(title)
    ax.locator_params(axis='x', nbins=7)



#
figure_zm.tight_layout()
figure_zm.show()
merlinLib.saveFig(figure_zm)



