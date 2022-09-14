"""
Reimplementation of vivs code to plot sfc temp, CO2 & Ice.
"""
import iris
import iris.plot
import matplotlib.pyplot as plt

import merlinLib
import numpy as np


plot_constraint = iris.Constraint(year=lambda cell: cell >= 2010)
start_emission = 2020
## read data

varsRead = ['SAT', 'icef', 'NPP','AtmC','Forcing', 'OcC', 'LandC', 'VegC','SoilC']
smoothVars = [var for var in varsRead if var[-1] != 'C' ]
refs = dict()
diffs = dict()
timeseries = dict()
sd=dict()
for var in varsRead:

    rr, dd, tt = merlinLib.proc_all(var, plot_constraint=plot_constraint)
    refs[var] = rr
    diffs[var] = dd
    timeseries[var] = tt
    sd[var] = merlinLib.std_error(refs[var]['Control'],window=11)

## plot data

titleDiff='Historical - Control'
varsPlot = ['AtmC', 'OcC', 'LandC','Forcing','SAT', 'icef']
titles = ['Atmosphere C','Ocean C', 'Land C', 'CO$_2$ Forcing','Sfc. Air Temp.', 'Ice Area']
xtitles = ['Pg C', 'Pg C','Pg C','Wm$^{-2}$','K', r'$10^{12}$m$^2$']
fig, axes = plt.subplots(nrows=3, ncols=2, num="sfc_vars", clear=True, figsize=[7.5, 8],sharex=True,gridspec_kw=dict(wspace=0.2))#,top=0.95,bottom=0.05))
hist = merlinLib.lookup.query('Reference=="Historical" & Time==200 & Carbon==1000').index[0]
ctl = merlinLib.lookup.query('Reference=="Control" & Time==200 & Carbon==1000').index[0]

for var, ax, ytitle, title in zip(varsPlot, axes.T.flatten(), xtitles, titles):
    print(var)

    if var == 'Budget': # do the budget
        # plot the difference in the C budget
        for v, t,col, in zip(['AtmC','OcC', 'LandC','VegC','SoilC'],
                             ['Atmos','Ocean','Land','Veg','Soil'],
                             ['red','blue','green','olive','brown']):
            ts = merlinLib.diff(diffs[v][hist],diffs[v][ctl])
            if var in smoothVars:
                ts = ts.rolling_window('time', iris.analysis.MEAN, 11)
            ax.plot(ts.coord('year').points,ts.data,label=t,linewidth=1,color=col)
        ax.legend(loc='lower right')
        y = ax.get_ylim()

    else:
        startIndx = 0
        for k, ts in diffs[var].items():
            series = merlinLib.lookup.loc[k]
            prop = merlinLib.properties(series,linewidth=1)
            t= ts.copy()
            if var in smoothVars:
                t = ts.rolling_window('time', iris.analysis.MEAN, 11)
            # generate name
            name = merlinLib.gen_name(k)
            ax.plot(t.coord('year').points, t.data, label=name, markevery=(startIndx,105), **prop,zorder=100)
            startIndx += 10
        ax.fill_between(ax.get_xlim(), 2 * sd[var], -2 * sd[var], color='grey', alpha=0.4)



    ax.set_title(title)
    ax.set_ylabel(ytitle,labelpad=0)
    ax.axhline(0.0, linestyle='dashed', color='black')
    ax.set_xlim(2000,2500)

    # plot start of emissions
    y = ax.get_ylim()
    ax.fill_betweenx(y, start_emission, start_emission + 10, color='slategrey',zorder=0)
    for t in merlinLib.lookup.Time.unique():# Plot grey bars to show  end of emissions.
        ax.fill_betweenx(y, t + start_emission, t + start_emission+10,  color='slategrey',zorder=0)
    # do some stuff for CO2
    if var == 'AtmC':
        ax2=ax.twinx()
        ax2.set_ylim(np.array(ax.get_ylim())*merlinLib.scale_Pg_ppm)
        ax2.set_ylabel("ppm",labelpad=-10,y=50)
        ax2.tick_params('y',direction="in",pad=-30)
    ax.margins()
    # done do final figure stuff.
# put year on the bottom two axis
for ax in axes[-1,:]:
    ax.set_xlabel('Year')
label = merlinLib.plotLabel()  # new set of labels
for ax in axes.T.flatten():
    label.plot(ax)

#fig.tight_layout(pad=0.,h_pad=0.0,w_pad=0.)
h, l = axes[0][0].get_legend_handles_labels()
fig.legend(h,l,loc='lower left',ncol=6,columnspacing=0.5,fontsize='small')
fig.show()
merlinLib.saveFig(fig)

## fit some data
import scipy.optimize
import pandas as pd
def expfn(time,v0,vinf,timescale):
    return v0+(vinf-v0)*(1-np.exp(-(time-time[0])/timescale))

print(f"{'variable':12s}" ,f"{'name':12s} {'v0':6s} {'vinf':6s} {'tau':6s}"*2)
lst_series=[]
for var,diff in diffs.items():
    print(f"{var:12s}",end=" ")
    for exper,info in merlinLib.lookup.query('Carbon==1000 & Time==200').iterrows():
        name = merlinLib.gen_name(exper)
        ts=diff[exper].extract(iris.Constraint(year  = lambda cell: 2050 < cell < 2220 ))
        optParam = scipy.optimize.curve_fit(expfn,ts.coord('year')._values,ts.data,[ts.data[0],2*ts.data.max(),100])
        print(f"{name:12s} {optParam[0][0]:6.3g} {optParam[0][1]:6.3g} {optParam[0][2]:6.3g}",end=' ')
        series=pd.Series(dict(var=var,name=name,v0=optParam[0][0],vinf=optParam[0][1],tau=optParam[0][2]))
        lst_series.append(series)

    print("")
df = pd.DataFrame(lst_series)
df.to_csv('summary_timescales.csv')


