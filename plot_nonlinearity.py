"""
Plot to show non-linearity (or lack off)
Things to show:
1) SST delta at +30-40 years after emissions start.
2) SST delta at  10-0 years before drawdown.
3) Zonal-mean SAT & NPP at 500 years or so but only for the two +1000 Gtonne cases
"""
import cartopy.crs
import iris
import iris.plot
import matplotlib.colors
import matplotlib.pyplot as plt
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import merlinLib
import numpy as np
import functools
import cartopy.crs as ccrs
import matplotlib.colors

plot_constraint = iris.Constraint(year=lambda cell: cell >= 2010)
#start_emission = 2020
units={'SAT':'K', 'Land_precip':'mm/day','precip':'mm/day','Forcing':'Wm$^{-2}$'}
titles = {'Land_precip': 'Land P', 'LandC': 'Land C', 'SoilC': 'Soil C','VegC': 'Veg. C', 'OcC': 'Ocn. C','AtmC':'Atmos. C'}
varsRead = dict(SAT_ZM=30, Land_precip_ZM=30, precip_ZM=30, VegC_ZM=30, SoilC_ZM=30, LandC_ZM=30,SST=30, SAT=30, VegC=30,AtmC=30,
                SoilC=30,Land_precip=30, precip=30,LandC=30,OcC=30,Forcing=30)  # variables we want and averaging period for sd
refs = dict()
diffs_ctl=dict()
diffs = dict()
timeseries = dict()
sd = dict()
for var, sd_time in varsRead.items():
    # ratio = var in ['Land_precip_ZM','NPP_ZM']
    rr, dd, tt = merlinLib.proc_all(var, plot_constraint=plot_constraint)
    refs[var] = rr
    diffs[var] = dd
    timeseries[var] = tt
    delta=dict()

    # for forcing need to recompute rel to control. Done below...

    diffs_ctl[var] = {name: merlinLib.diff(tt_ind, rr['Control']) for name, tt_ind in
                      tt.items()}  # difference between timeseries and control
    sd[var] = merlinLib.std_error(refs[var]['Control'], window=sd_time)


##

def comp_prd_means(variab,series,pnames=None,emission_stop=None):
    periods = merlinLib.named_periods(series,pnames=pnames,emission_stop=emission_stop,constraint=True)
    summary = dict()
    for name,constrain in periods.items():
        v = variab.extract(constrain)  #
        v = float((v.collapsed('time', iris.analysis.MEAN) * 1000 / series.Carbon).data)
        summary[name] = v

    summary = pd.Series(summary).rename(exper).reindex(pnames) # giving null for anything missing.
    return summary


scales = dict(Land_precip=24 * 60 * 60)

# do the data processing.
dataframes=dict()
dataframes_ctl=dict()
hist_ctl_carbon=552.5 # rough amount of Pg emitted for hist_ctl corrected for initial error which reset carbon and removed 38.7 PgC from system.
for var in ['SAT', 'Land_precip','LandC','SoilC', 'VegC','OcC','AtmC','Forcing']:
    print(var,end='.. ')
    df = []
    df_ctl = []
    for exper, series in merlinLib.lookup.iterrows():
        summary = comp_prd_means(diffs[var][exper],series)
        series_ref_ctl = series.copy()
        if series_ref_ctl.Reference == 'Historical':
            series_ref_ctl.Carbon += hist_ctl_carbon
        summary_ctl = comp_prd_means(diffs_ctl[var][exper],series_ref_ctl)
        df.append(summary)
        df_ctl.append(summary_ctl)
    dataframes[var]=merlinLib.lookup.merge(pd.DataFrame(df),left_index=True,right_index=True).sort_values(by=['Carbon','Time','Reference'])
    df_ctl = merlinLib.lookup.merge(pd.DataFrame(df_ctl), left_index=True, right_index=True)
    L=(df_ctl.Reference == 'Historical')
    df_ctl.loc[L,'Carbon'] += hist_ctl_carbon # Carbon emission rel to PI
    df_ctl = df_ctl.sort_values(by=['Carbon','Time','Reference'])
    # add on the Historical Simulation information to the ctl.
    hist_series=[]
    for time in [50,100,200]:
        series=pd.Series(dict(Carbon=hist_ctl_carbon,Time=time,Reference='Special'))
        svalues = comp_prd_means(refs[var]['Historical']-float(refs[var]['Control'].data.mean()),series,emission_stop=2020,pnames=['+E','-D'])
        series = pd.concat([series,svalues]).rename(merlinLib.references['Historical'])
        hist_series.append(series)
    df_ctl = pd.concat([df_ctl,pd.DataFrame(hist_series)])

    dataframes_ctl[var] = df_ctl
print("") # print a newline!
##
def plot_means(df, sd, pnames, plotNames=False,xlabel=None,ylabel=None,ax=None, title=None,errorBar=None,**kwargs):
    """
    Plot changes in 30 year means wrt timescale and Carbon. Timescale grouped.
    :param ax:
    :param df:
    :param sd:
    :param pnames:
    :param errorBar:
    :return:
    """
    if ax is None:
        ax=plt.gca()

    ucarbon = np.sort(df.Carbon.unique())
    carbon_offsets = dict(zip(ucarbon,np.arange(len(ucarbon)))) # offsets for Carbon
    dp=0.075
    utime = np.sort(df.Time.unique())
    time_offsets=dict(zip(utime,np.arange(len(utime))*dp))
    nnames=len(pnames)
    dname = 1.0/nnames
    offsets=np.arange(nnames)*dname-0.45 # plot values for four periods.

    ms = kwargs.pop('ms',6) # extract the marker size as used in the error plot.
    linestyle=kwargs.pop('linestyle','none')
    for exper, series in df.iterrows():  # iterate over data
        prop = merlinLib.properties(series, linestyle=linestyle, ms=ms,**kwargs)
        x = carbon_offsets[series.Carbon] + offsets + time_offsets[series.Time]
        y = series.loc[pnames]
        if series.Reference == errorBar:  # plot error bars
            sd_var = 2 * sd  # double difference so SD doubled. (x-ref)-(y-ref2)
            err = float(2 * sd_var * 1000 / series.Carbon)
            ax.errorbar(x, y, yerr=err, elinewidth=1, capsize=ms / 2, capthick=1, **prop)
        else:
            ax.plot(x, y, **prop)
    # Put dashed lines to separate periods for all C values
    for k, x in carbon_offsets.items():  # iterate over carbon
        for offset in offsets:  # and over periods.
            ax.axvline(offset + x, linestyle='dotted', color='grey')

    # decorate  first put grey error ranges on
    sd_var = np.sqrt(2.) * sd# increase SD as difference between model and reference.
    for carbon in ucarbon:
        err = 2 * sd_var * 1000 / carbon
        ax.fill_between([carbon_offsets[carbon] - 0.55, carbon_offsets[carbon] + 0.45], err, -err, color='grey', alpha=0.4)
        ax.axvline(carbon_offsets[carbon] - 0.55, color='black', linestyle='dashed')

    tv = list(carbon_offsets.values())
    tl = list(carbon_offsets.keys())
    ax.set_xticks(tv)
    ax.set_xticklabels(tl)
    ax.axhline(0.0, color='black', linestyle='dashed')
    ax.set_xlim(np.min(tv) - 0.5, np.max(tv) + 0.5)
    if plotNames:
        text_kw=dict(ha='left',fontsize='small',fontweight='bold')
        ylim = ax.get_ylim()[0]
        x = list(carbon_offsets.values())[0]
        for offset,txt in zip(offsets,pnames):
            ax.text(offset + x, ylim, txt, **text_kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)




## plot +E and -D for C difference to ctl (rather than historical/Control)
def plot_df(df,sd,ax=None,offset=0,alpha=None):
    """
    Plot a dataframe
    :param ax:
    :param df:
    :param offset:
    :return:
    """

    if ax is None:
        ax=plt.gca()
    if alpha is None:
        alpha=1
    df100 = df.query('Time==100')
    ax.scatter(df100.Carbon+offset, y=df100.loc[:,'+E'], c=df100.loc[:,'color'], marker='x',  s=20,alpha=alpha,zorder=20)  # marker_sizes[0])
    err = 2 * sd * 1000 / df100.Carbon  # 2 sigma.
    ax.errorbar(df100.Carbon+offset, df100.loc[:, '+E'], yerr=err, elinewidth=2, capsize=7, capthick=2, color='lime',
                linestyle='none',alpha=alpha,zorder=2)

    for t in df.Time.unique():
        dd=df.query(f"Time=={t}")
        ax.scatter(dd.Carbon+offset, y=dd.loc[:,'-D'], c=dd.color, marker=dd.marker.iloc[0], s=20,alpha=alpha,zorder=20)
        if t == 200: # plot errorbar
            err = 2 * sd * 1000 / dd.Carbon  # 2 sigma.
            ax.errorbar(dd.Carbon+offset, dd.loc[:, '-D'], yerr=err, elinewidth=2, capsize=4, capthick=2, color='black',
                        linestyle='none',alpha=alpha,zorder=2)


ref_lookup = {'Historical': {'colour': 'red', 'linestyle': 'solid'},
              'Control': {'colour': 'blue', 'linestyle': 'solid'},
              }
fig_ref_ctl, axes  =plt.subplots(num='nonlinearity', clear=True, ncols=2, nrows=6,sharex='col')
fig_ref_ctl.set_size_inches(7.5,8)
marker_sizes={0:40,50:20,100:8,200:12}
label=merlinLib.plotLabel()



pnames=['+E','-D','+D','R']
for ax,ax_post,var in zip(axes[:,0],axes[:,1],['SAT', 'Land_precip','Forcing','AtmC','OcC','LandC'] ):
    title = titles.get(var,var)
    scale= scales.get(var, 1.0)
    df_prop=[]
    for name, series in dataframes_ctl[var].iterrows():
        prop = merlinLib.properties(series)#,color=ref_lookup[series.Reference]['colour'])
        df_prop.append(pd.Series(prop).rename(name))
    df_prop = pd.DataFrame(df_prop)
    df=dataframes_ctl[var].merge(df_prop,left_index=True,right_index=True).copy()
    df.loc[:,pnames] *= scale
    df_std = dataframes[var].query("Reference=='Historical'").merge(df_prop,left_index=True,right_index=True).copy()
    df_std.loc[:,pnames]*= scale
    sd_var = np.sqrt(2) * float(sd[var])*scale # difference between two simulations.
    plot_df(df,sd_var,ax=ax)
    plot_df(df_std, sd_var, offset=-50,ax=ax,alpha=0.2)

    ax.axhline(float(df.query(f"Carbon=={df.Carbon.max()} & Time==200").loc[:,'-D']),linestyle='dashed',color='black')
    ax.set_title(title)
    pos = ax.get_position() # position of axis in figure
    xlabel=None
    if pos.ymin < 0.1: # bottom
        ax.set_xlabel('Total Carbon Emitted (Pg C)')
        xlabel = 'Carbon Pulse (Pg C)'
    if pos.xmin <0.15:# left
        ax.set_ylabel(f"{units.get(var,'Pg C')}/Eg C")
    label.plot(ax)
    df = dataframes[var].copy()
    df.loc[:,pnames]*= scale
    plot_means(df, sd_var, ['+D','R'], plotNames=True, xlabel=xlabel, ylabel=None, ax=ax_post, title=title, errorBar=None)
    label.plot(ax_post)

fig_ref_ctl.tight_layout()
fig_ref_ctl.show()
merlinLib.saveFig(fig_ref_ctl)

## work out regression between SAT and forcing for +E & -D periods rel to control
import statsmodels.api as sm
def comp_reg(dataframes):

    forcing = dataframes['Forcing'].loc[:,['+E','-D']].values.flatten()
    sat = dataframes['SAT'].loc[:,['+E','-D']].values.flatten()
    sat = sat[~np.isnan(sat)]

    forcing = forcing[~np.isnan(forcing)]
    #forcing = sm.add_constant(forcing)
    fit = sm.OLS(sat,forcing).fit()
    return fit
print("Forcing Regressions")
print(comp_reg(dataframes_ctl).summary())
print(comp_reg(dataframes).summary())
print("Done Forcing reg")


##plot Zonal means of two *_1000C_200y cases at R period & long lat plots ofor C_1000C_200y delta to Control
print("Plotting Zonal Means")
df=merlinLib.lookup.query('Carbon == 1000 & Time==200')
figZm,axis_zm=plt.subplots(nrows=2,ncols=2, num='nonlinearity_zm', clear=True,sharey=True)
figZm.set_size_inches(7.5,3.5)
#
def plot_zm(ax,var,df,title,linestyle=None,decorate=True):
    """Plot ZM values"""
    for name, series in df.iterrows():
        scale = 1.0
        if var.find('precip') >= 0:
            scale = 24 * 60 * 60.  # conversion to mm/day
        diff = diffs[var][name] * scale
        prop = merlinLib.properties(series, marker=None,linestyle=linestyle,ms=10)
        period = merlinLib.named_periods(series,pnames=['R'],constraint=True)
        c = diff.extract(period['R']).collapsed('time',iris.analysis.MEAN)
        y = c.coord('latitude').points
        ax.plot(c.data, y, markevery=3, **prop)
    if decorate:
        # decorate plot
        ax.fill_betweenx(y, 2 * sdev, -2 * sdev, color='grey', alpha=0.3)
        ax.axvline(0.0, linestyle='dashed', color='black')
        ax.set_xlabel(r"$\Delta$ " + xtitle)
        # axZM.set_ylabel('Latitude')
        ax.set_title(title)
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.set_yticks([-90, 60, -30, 0, 30, 60, 90])

for axZM, var, title, xtitle in zip(axis_zm.flatten() ,
                                    ['SAT_ZM','precip_ZM','SoilC_ZM', 'VegC_ZM'],
                                    ['Temperature', 'Precipitation','Soil C', 'Veg. C'],
                                    ['T (K)','P (mm/day)', 'Soil C (Pg/deg.)','Veg. C (Pg/deg.)']):
    sdev = sd[var] * np.sqrt(2.)
    scale = 1.0

    if var.find('precip') >= 0:
        scale = 24 * 60 * 60.  # conversion to mm/day
    sdev *= scale
    linestyle = 'solid'
    plot_zm(axZM,var,df,title)
    if var == 'precip_ZM':
        plot_zm(axZM, 'Land_precip_ZM', df,title,linestyle='dashed',decorate=False)


# put labels on
label = merlinLib.plotLabel()
for ax in axis_zm.flatten():
    label.plot(ax)
figZm.tight_layout()
figZm.show()
merlinLib.saveFig(figZm)


# plot long/lat plots

@functools.lru_cache(maxsize=6000)
def mean_dif(var, exper, ratio=False, crit=None):
    series = merlinLib.lookup.loc[exper]
    ref_exper = merlinLib.references[series.Reference]
    cubes=[]
    period = merlinLib.named_periods(series)['R']
    for e in [exper,ref_exper]:
        cube = merlinLib.read_file(var, e)
        cube = cube.extract(iris.Constraint(year=lambda cell: period.start <= cell < period.stop))
        cube = cube.collapsed('time', iris.analysis.MEAN)
        cubes.append(cube)
    delta = merlinLib.diff(cubes[0],cubes[1],ratio=ratio)
    if ratio and (crit is not None): # mask below critical values.
        L = np.abs(cubes[1].data) < crit
        delta.data.mask[L] = True


    return delta





diffs_ll = dict()

for var in ['precip', 'LandC','PFT']:
    d = {}
    sub = merlinLib.lookup.query('Time ==200 & Carbon==1000')
    scale = 1.0
    if var.find('precip') >= 0:
        scale = 24 * 60 * 60.  # conversion to mm/day
    for exper, series in sub.iterrows():
        #ratio = (var == 'LandC')
        d[exper] = mean_dif(var, exper) * scale
    diffs_ll[var] = d

## plot now
print("Plotting Long lat maps")
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
proj = ccrs.PlateCarree()
fig, axis_ll = plt.subplots(ncols=2, nrows=1, sharex=True, sharey=True,num='delta_ll', clear=True,
                         squeeze=False,subplot_kw=dict(projection=proj))
fig.set_size_inches(7.5,2.75) # shame that plt.subplots does not just do this!

#label = merlinLib.plotLabel()  -- keep on with the previous label. Then we can just merge the figures in the text!
levels_t = [-3, -2, -1, -0.5,  0.5, 1, 2, 3]
levels_precip = np.array([-2, -1.5, -1, -0.5,  0.5, 1, 1.5, 2])
levels_C = np.arange(-4,5)
for ax, var, levels, title,cmap in zip(axis_ll.flatten(), ['precip', 'LandC'], [levels_precip, levels_C],
                                  [ 'Precip', 'Land C'],
                                       ['RdBu','RdYlGn']):
    diff = diffs_ll[var]['xovfd']
    norm = matplotlib.colors.BoundaryNorm(levels, 255)
    cf = iris.plot.pcolormesh(diff, axes=ax, norm=norm,cmap=cmap)
    ax.set_xticks([-180, -90, 0,90, 180], crs=ccrs.PlateCarree())  # from cartopy example
    # but note the co-ords need to be given in the right order...
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    ax.coastlines()
    unit = units.get(var,'kg C m$^{-2}$')
    ax.set_title(rf'$\Delta$ {title} ({unit})')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    cbar = fig.colorbar(cf, ax=ax, ticks=levels, orientation='horizontal', drawedges=True,fraction=0.1,pad=0.2,extend='both')
    label.plot(ax)
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)
