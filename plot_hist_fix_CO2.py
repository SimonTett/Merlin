"""
Plot Forcing (from CO2) and SAT from Fixed CO2 cases
"""
import matplotlib.pyplot as plt
import merlinLib
import iris.plot

fig,axes = plt.subplots(nrows=2,ncols=1,clear=True,figsize=[7.5,4],num='hist_fix_CO2',sharex=True)
for  n in ['Historical_fixC','Historical', 'Control_fixC','Control']:
    f = merlinLib.read_data('Forcing',merlinLib.references[n])
    s = merlinLib.read_data('SAT',merlinLib.references[n])
    for ts,ax in zip([f,s],axes.flatten()):
        iris.plot.plot(ts.coord('year'), ts,label=n,axes=ax)

# common decorations
label=merlinLib.plotLabel()
for ax,title in zip(axes.flatten(),['Diag CO2 Forcing','SAT']):
    label.plot(ax)
    ax.set_title(title)
axes[0].legend()
fig.tight_layout()
fig.show()
