"""
Plot the land c fluxes
"""
import merlinLib
import matplotlib.pyplot as plt
import iris.plot
import numpy as np

end_time=2500
def delta_ref(var,sim,ref):
    """
    Compute the difference ts for named var and simulation wrt
    appropriate reference simulation.
    """
    ref_name=merlinLib.references[ref]
    ts = merlinLib.read_data(var, sim)
    ref_ts = merlinLib.read_data(var, ref_name)
    delta = merlinLib.diff(ts, ref_ts)#.rolling_window('time', iris.analysis.MEAN, 11)
    if end_time is not None:
        delta=delta.extract(iris.Constraint(year=lambda yr : yr < end_time))
    return delta

experiments = merlinLib.lookup.query('Carbon==1000 and Time==200')
vars_stock=['VegC','SoilC','CO2atm','OcC']
colors=['palegreen','green','red','blue']

fig, axes = plt.subplots(nrows=4,ncols=1,sharex=True,num='land_c_stock_delta',clear=True, squeeze=False,
                      figsize=merlinLib.fsize[::-1])
start_emission=2020
ax=axes[0][0]
ax_ts = axes[1][0]
ax_ref = axes[2][0]
label=merlinLib.plotLabel()
delta = {}
ts = {}
extract =iris.Constraint(year=lambda cell: cell >=2000 )
for var,color in zip(vars_stock,colors):

    for exper,info  in experiments.iterrows():
        delta[(info.Reference,var)] = delta_ref(var, exper,info.Reference).extract(extract)
        ts[(exper,var)] = merlinLib.read_data(var,exper).extract(extract)
        exper_ref = merlinLib.references[info.Reference]
        ts[(info.Reference,var)] = merlinLib.read_data(var,exper_ref).extract(extract)



    delta_stock = merlinLib.diff(delta[('Historical',var)],delta[('Control',var)])
    d_ref = merlinLib.diff(ts[('Historical',var)],ts[('Control',var)])
    iris.plot.plot(delta_stock.coord('year'), delta_stock, axes=ax,
                       markevery=10, color=color,marker='o',
                       label=var,linewidth=2)
    # iris.plot.plot(d_ref.coord('year'), d_ref, axes=ax_ts,
    #                    markevery=10, color=color,marker='o',
    #                    label=var,linewidth=2)

ax.axhline(0.0,color='black',linestyle='dashed')
ax.set_title("Carbon Stocks")
ax.set_ylabel("Pg C")

## plot TotalSoilC
for var,ax in zip(['SoilC','VegC','OcC'],axes.flatten()[1:]):
    for exper,info  in experiments.iterrows():
        props = merlinLib.properties(info)
        times = ts[(exper,var)]
        iris.plot.plot(times.coord('year'),times,axes=ax,
                       label=exper,markevery=20,**props)

        times = ts[(info.Reference,var)]
        iris.plot.plot(times.coord('year'),times,axes=ax,
                       label=info.Reference,linestyle='dotted',markevery=20,**props)
        ax.set_title(var)
for a in axes.flatten():

    a.legend(ncol=2)
    label.plot(a)
    times = experiments.Time.unique()
    y = a.get_ylim()
    a.fill_betweenx(y, start_emission, start_emission+10,
                        alpha=0.4, color='grey')
    for t in times:
        a.fill_betweenx(y, t + start_emission, t + start_emission+10,
                        alpha=0.4, color='grey')




fig.show()
merlinLib.saveFig(fig)