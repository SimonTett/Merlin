"""
Plot the land c fluxes
"""
import merlinLib
import matplotlib.pyplot as plt
import iris.plot
import numpy as np

end_time=2500
def delta_ref(var,info,smooth=11):
    """
    Compute the difference ts for named var and simulation wrt
    appropriate reference simulation.
    """
    sim = info.name
    ref = merlinLib.references[info.Reference]
    ts = merlinLib.read_data(var, sim)
    ref_ts = merlinLib.read_data(var, ref)
    delta = merlinLib.diff(ts, ref_ts)
    if smooth is not None:
        delta = delta.rolling_window('time', iris.analysis.MEAN, smooth)

    return delta

experiments = merlinLib.lookup.query('Carbon==1000 and Time==200')
#cases = ['xocyc','xoaud']
properties = [merlinLib.properties(s) for name,s in experiments.iterrows()]
colours = ['red','blue']
variables = ['NPP','Litter','SoilResp']
#variables = ['Litter','SoilResp','NPP']
vars_stock=['VegC','SoilC','LandC']
markers=['s','*','^']
linestyles=['dotted','dashed','solid']
fig,(ax_stock,ax_flux) = plt.subplots(nrows=2,ncols=1,sharex=True,num='land_c_flux',clear=True,
                      figsize=merlinLib.fsize[::-1])
start_emission=2020
for exper,info  in experiments.iterrows():
    #prop_series=merlinLib.lookup.loc[sim]

    ref=merlinLib.references[info.Reference]

    # plot the stocks

    for var,marker,linestyle in zip(vars_stock,markers,linestyles):
        delta = delta_ref(var, info)
        plt_props = merlinLib.properties(info, alpha=1,
                                         marker=None, linestyle=linestyle)

        # ax_flux.plot(delta.coord('year').points,np.cumsum(delta.data),
        #         label=prop_series.Reference+' '+var,markevery=10,**plt_props)
        iris.plot.plot(delta.coord('year'), delta, axes=ax_stock,
                       markevery=10, **plt_props, ms=4,
                       label=info.Reference + ' ' + var)

    for var,marker,linestyle in zip(variables,markers,linestyles):
        delta = delta_ref(var,info)
        plt_props=merlinLib.properties(info,alpha=1,
                                       marker=None,linestyle=linestyle)

        #ax_flux.plot(delta.coord('year').points,np.cumsum(delta.data),
        #         label=prop_series.Reference+' '+var,markevery=10,**plt_props)
        iris.plot.plot(delta.coord('year'),delta,axes=ax_flux,
                       markevery=10,**plt_props,ms=4,
                       label=info.Reference+' '+var)


label=merlinLib.plotLabel()
for a in [ax_stock,ax_flux]:
    a.axhline(0.0,color='black',linestyle='dashed')
    a.legend(ncol=2)
    label.plot(a)
    t = info.Time
    y = a.get_ylim()
    a.fill_betweenx(y, start_emission, start_emission+10,
                        alpha=0.4, color='grey')
    a.fill_betweenx(y, t + start_emission, t + start_emission+10,
                        alpha=0.4, color='grey')

ax_stock.set_title("Land Carbon Stocks")
ax_stock.set_ylabel("Pg C")

ax_flux.set_title("Land Carbon Fluxes")
ax_flux.set_ylabel("Pg C a$^{-1}$")
fig.show()
merlinLib.saveFig(fig)