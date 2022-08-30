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



figure = plt.figure(num='pft_plot',clear=True,figsize=[7.5,7])
subfigs = figure.subfigures(ncols=1,nrows=2)
axes = subfigs[0].subplots(nrows=1,ncols=2,sharex=True,sharey=True)

label=merlinLib.plotLabel()
for series,ax in zip([ctl,hist],axes.flatten()):
    pftC_delta = delta_var('VegC_PFT', series.name)  # *100
    rgn_pft_BLT=dict()
    for reg in merlinLib.named_regions.keys():
        rgn_pft_BLT[reg]=iris.util.squeeze(delta_var('VegC_PFT',series.name,region=reg))
    soilC_delta = delta_var('SoilC', series.name)
    LandC_delta = delta_var('LandC',series.name)
    merlinLib.plot_pft(pftC_delta,ax=ax,pftindices=[1,2,3,4,5],marker='')
    for (key,ts),col in zip(rgn_pft_BLT.items(),['red','orange','blue']): ## TODO set up named colours and hten plot
        merlinLib.plot_pft(ts, ax=ax, pftindices=[1], colors=dict(BLT=col), marker='')
    iris.plot.plot(soilC_delta.coord('year'),soilC_delta,axes=ax,color='brown',label='Soil')
    iris.plot.plot(LandC_delta.coord('year'), LandC_delta, axes=ax, color='black', label='Land',linewidth=2)
    merlinLib.plot_emission(ax=ax,overshoot=[200])
    title=merlinLib.gen_name(series.name)
    ax.set_title(title)
    ax.axhline(linestyle='dashed', color='black')
    label.plot(ax)
    ax.set_xlabel("Year")
    ax.set_ylabel("$\Delta$ Reference (Pg C)")
    ax.locator_params(axis='y', nbins=11)

axes[0].legend(ncol=4,fontsize='small',handletextpad=0.2,columnspacing=1)

# plot maps of blt C
blt_ctl_pulse = merlinLib.read_cube('VegC_PFT',ctl.name).\
                    extract(iris.Constraint(year=lambda cell: 2300 <= cell < 2400)).\
                    collapsed('time',iris.analysis.MEAN)[0,:,:]
blt_ref = merlinLib.read_cube('VegC_PFT',merlinLib.references[ctl.Reference]).\
                    collapsed('time',iris.analysis.MEAN)[0,:,:]
axes_ll = subfigs[1].subplots(nrows=1,ncols=2,sharex=True,sharey=True,subplot_kw=dict(projection=ccrs.PlateCarree()))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
c_level = [1,2,3,4,5,6,8]
for blt,title,ax in zip([blt_ref,blt_ctl_pulse],['Control','C_1000C_200y'],axes_ll.flatten()):
    norm = matplotlib.colors.BoundaryNorm(c_level, 255)
    cf = iris.plot.pcolormesh(blt,cmap='YlGn',norm=norm,axes=ax)
    ax.set_xticks([-180, -90, 0,90, 180], crs=ccrs.PlateCarree())  # from cartopy example
    # but note the co-ords need to be given in the right order...
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    ax.coastlines()
    unit = 'kg C m$^{-2}$'
    ax.set_title(rf'{title} ({unit})')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

cbar = subfigs[1].colorbar(cf, ax=axes_ll, ticks=c_level, orientation='horizontal', drawedges=True,fraction=0.1,pad=0.2,extend='both')

#figure.tight_layout()
figure.show()
merlinLib.saveFig(figure)



