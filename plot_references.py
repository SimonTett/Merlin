"""
Plot the reference simulations
Will plot Atmos C & SAT as key variables.
"""
import merlinLib
import matplotlib.pyplot as plt
import iris.plot
import  numpy as np
import statsmodels.api as sm
varsPlot = ['AtmC', 'OcC', 'LandC','Forcing','SAT', 'icef']


titles = dict(AtmC='Atmosphere C',Forcing='CO$_2$ Forcing',SAT='Sfc. Air Temp.',OcC='Ocean C',LandC='Land C',icef='Ice Area')
ytitles = dict(AtmC='Pg C', OcC='Pg C',LandC='Pg C',Forcing='Wm$^{-2}$',SAT='K',icef=r'$10^{12}$m$^2$')
refs=['Spinup','Historical_bad_1750','Historical','Control']
labels = dict(Historical_bad_1750='Historical_1750')

cols=dict(Spinup='black',Spinup_old='grey',Historical_bad_1750='brown',Historical='red',Control='blue')
timeseries=dict()

for var in varsPlot:
    data=dict()
    for ref in refs:
        if ref == 'Historical_bad_1750':
            t = merlinLib.read_data(var,merlinLib.references[ref],decadal=True) # only have decadal mean data.
            data[ref]= t.extract(iris.Constraint(year =lambda cell: cell <= 1850))
        else:
            data[ref]=merlinLib.read_data(var,merlinLib.references[ref]) #

    timeseries[var]=data

fig_ref,axis = plt.subplots(nrows=3,ncols=2,clear=True,sharex=True,num='reference',squeeze=False)
fig_ref.set_size_inches(7.5,6)
plabel = merlinLib.plotLabel()
for (var, ts),ax in zip(timeseries.items(),axis.T.flatten()):
    plabel.plot(ax)
    for name,t in ts.items():
        plot_params = dict(axes=ax,color=cols[name],label=labels.get(name,name),linewidth=2)
        iris.plot.plot(t.coord('year'),t,**plot_params)
    ax.set_title(titles[var])
    ax.set_ylabel(ytitles[var])
    # add horizontal line corresponding to *mean* Control values
    mn = float(ts['Control'].collapsed('time',iris.analysis.MEAN).data)
    ax.axhline(mn,linestyle='dashed',color='black')
    ax.axvline(2010,linestyle='solid',color='grey',linewidth=2,zorder=0)
    ax.locator_params(axis='y',nbins=5)
    # add a line connecting the spinup to the startup. (1750 - 2000)
    spinupv = float(ts['Spinup'].extract(iris.Constraint(year=1750)).data) #
    controlv = float(ts['Control'].extract(iris.Constraint(year=2001)).data)
    ax.plot([1750,2001],[spinupv,controlv],color=cols['Control'],linestyle='dotted')
    # add labels for AtmC & ppm
    if var == 'AtmC':
        ax.legend(fontsize='small')
        ax2 = ax.twinx()
        ax2.set_ylim(np.array(ax.get_ylim()) * merlinLib.scale_Pg_ppm)
        ax2.set_ylabel("ppm")#, labelpad=-10, y=50)
        #ax2.tick_params('y', direction="out", pad=-30)
        ax2.locator_params(axis='y',nbins=5)
    #ax.margins()

for ax in axis[-1][:]:
    ax.set_xlabel("Year")
fig_ref.tight_layout()
fig_ref.show()
merlinLib.saveFig(fig_ref)

## work out error from resetting CO2 in Historical case at 1850
totalC=dict()
emis = dict()
for ref in ['Historical','Control']:
    exper = merlinLib.references[ref]
    totalC[ref]=merlinLib.read_data('TotalC',exper)
    emis[ref]=merlinLib.read_data('CO2emis',exper)

error = float((totalC['Historical'][0]-totalC['Control'][0]).data)
print("Carbon change from error in Historical initialisation ",error)
# total carbon emitted and then corrected for error.
totEmit=emis['Historical'].data.sum()
print(f"Hist emis {totEmit:3.1f}  and corrected {totEmit+error:3.1f}")

# CO2 varn for comparision with Law Dome estimate
tsmooth = timeseries['AtmC']['Spinup'].rolling_window('time', iris.analysis.MEAN, 30) # 30 year running mean for quick and dirty comparision with ice core data
Spinup_co2_ppm = tsmooth*merlinLib.scale_Pg_ppm
print(
    f"Spinup CO2 (ppm) sd:{Spinup_co2_ppm.data.std():3.1f} max:{Spinup_co2_ppm.data.max():3.1f} min:{Spinup_co2_ppm.data.min():3.1f}")

# and the linear trend.
def lin_trend(ts):
    x = ts.coord('year').points / 100  # per century
    x -= x[0]  # to start.
    x = sm.add_constant(x)
    fit = sm.OLS(ts.data, x).fit()
    return fit

ts= timeseries['AtmC']['Spinup']*merlinLib.scale_Pg_ppm
print(f"Mean CO2 {ts.data.mean():3.1f} ppm")
fit = lin_trend(ts)
print(fit.summary())
# and the same for Ocean & Land C
for var in ['OcC','LandC','AtmC']:

    vfit = lin_trend(timeseries[var]['Spinup'])
    print("="*10,'  ',(var+" ")*2,"="*10)
    print(vfit.summary())


# and the forcing relationship for Historical.
print("SAT = lambda Forcing")
force = timeseries['Forcing']['Historical']-timeseries['Forcing']['Control'].data.mean()
sat = timeseries['SAT']['Historical']-timeseries['SAT']['Control'].data.mean()
fit = sm.OLS(sat.data,force.data).fit()
print(fit.summary())

##






