"""
Plot Sensitivity cases.
1) Drawdown over 100 years
"""
import merlinLib
import matplotlib.pyplot as plt
import iris.plot
ref_cases = merlinLib.lookup.query("Reference=='Control' & Carbon==1000")

fig,axes = plt.subplots(nrows=1,ncols=1,squeeze=False,clear=True,figsize=[7,3],num='sens')
ax = axes[0][0]
colours={200: 'firebrick', 100: 'red', 50: 'coral'} # use same colours.,
for name, series in ref_cases.iterrows():
    prop=merlinLib.properties(series,linewidth=2,color=colours[series.Time],marker='None')
    reference = merlinLib.references[series.Reference]
    ref_ts = merlinLib.read_data('SAT',reference)
    ts = merlinLib.read_data('SAT',name)
    ts = merlinLib.diff(ts,ref_ts)
    label= merlinLib.gen_name(name)
    iris.plot.plot(ts.coord('year'),ts,axes=ax,label=label,**prop)
# plot the sens case...
sens = merlinLib.read_data('SAT','xovfq')
sens = merlinLib.diff(sens,ref_ts)
label = 'C_2000C_100y_s'
prop['color']='black'
iris.plot.plot(sens.coord('year'),sens,axes=ax,label=label,**prop)
# decorate axis
merlinLib.plot_emission(ax=ax)
ax.axhline(0,color='black',linestyle='dashed')
ax.legend(ncol=2,fontsize='small',handletextpad=0.1,columnspacing=0.5,labelspacing=0.25,handlelength=1.2,borderaxespad=0)
ax.set_ylabel(r"$\Delta$ SAT (K)")
ax.set_xlabel("Year")
ax.set_title("Temperature differences for C_1000C_*")


fig.tight_layout()
fig.show()
#merlinLib.saveFig(fig) # uncomment when ready to roll

