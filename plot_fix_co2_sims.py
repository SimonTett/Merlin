"""
Plot timeseries for fixed CO2 simulations
"""

import merlinLib
import matplotlib.pyplot as plt
import iris.plot
ppt=False # set False if want figure for paper
varsPlot = ['AtmC', 'OcC', 'LandC','Forcing','SAT', 'icef']
smoothVars = [var for var in varsPlot if var[-1] != 'C' ]
titles = ['Atmosphere C','Ocean C', 'Land C', 'CO$_2$ Forcing','Sfc. Air Temp.', 'Ice Area']
xtitles = ['Pg C', 'Pg C','Pg C','Wm$^{-2}$','K', r'$10^{12}$m$^2$']
if ppt:
    figSize=[11,8]
    nrows=2
    ncols=3
    name='fix_co2_sims_ppt'
else:
    figSize=[7.5,8]
    nrows=3
    ncols=2
    name='fix_co2_sims'

fig,axes = plt.subplots(nrows=nrows,ncols=ncols,clear=True,figsize=figSize,num=name,sharex=True,
                        gridspec_kw=dict(wspace=0.25,top=0.95,bottom=0.1))
for var,title,xtitle,ax in zip(varsPlot,titles,xtitles,axes.T.flatten()):
    for name, experiment in merlinLib.lookup_fix.iterrows():
        ts= merlinLib.delta(var,name,refName=experiment.Reference)
        if var in smoothVars:
            ts = ts.rolling_window('time', iris.analysis.MEAN, 11)
        plotatt=merlinLib.properties(experiment,markevery=50,linestyle='solid')
        iris.plot.plot(ts.coord('year'),ts,axes=ax,**plotatt)
        experiment2=merlinLib.lookup.query(f"Reference==\'{experiment.Reference.replace('_fixC','')}\' and Carbon=={experiment.Carbon} and Time=={experiment.Time}").iloc[0]
        ts2= merlinLib.delta(var,experiment2.name,refName=experiment2.Reference)
        if var in smoothVars:
            ts2 = ts2.rolling_window('time', iris.analysis.MEAN, 11)
        plotatt=merlinLib.properties(experiment2,markevery=50,alpha=0.5)
        iris.plot.plot(ts2.coord('year'),ts2,axes=ax,**plotatt)
    for extra,alpha in zip(['_fixC',''],[1,0.5]): # plot the historical cases
        label = 'Historical'+extra
        exper = merlinLib.references[label]
        ref = 'Control'+extra
        ts_hist = merlinLib.delta(var,exper,refName=ref,meanValue=True)
        if var in smoothVars:
            ts_hist = ts_hist.rolling_window('time', iris.analysis.MEAN, 11)
        iris.plot.plot(ts_hist.coord('year'),ts_hist,axes=ax,color='black',linewidth=2,alpha=alpha,label=label)


    ax.set_title(title)
    ax.set_ylabel(xtitle)
    ax.axhline(color='black',linestyle='dashed')
    merlinLib.plot_emission( ax=ax, overshoot=[200])

# put year on the bottom two axis
for ax in axes[-1,:]:
    ax.set_xlabel('Year')
label = merlinLib.plotLabel()  # new set of labels
for ax in axes.T.flatten():
    label.plot(ax)
h, l = axes[0][0].get_legend_handles_labels()
fig.legend(h,l,loc='lower center',ncol=6,columnspacing=0.5,fontsize='small')
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)



