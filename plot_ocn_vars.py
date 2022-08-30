"""
Python script to plot the ocean volume average temperature and explain its behaviour.
"""

import iris.plot
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import merlinLib
import statsmodels.api as sm



plot_constraint = iris.Constraint(year=lambda cell: cell >= 2019)
refs = {}
diffs = {}
timeseries = {}
for var in ['VAOT', 'OcC_Profile']:
    rr, dd, tt = merlinLib.proc_all(var, plot_constraint=plot_constraint)
    refs[var] = rr
    diffs[var] = dd
    timeseries[var] = tt
    print(f"Processed {var}")

## now we can plot them
start_emission = 2020
"""
want 3 columns by 2 rows. 
Row 1
1) VOAT time-series from all simulations.  Linear scale
2) VOAT vs C*year to draw down. Linear scale
3) VAOT vs time after drawn down. Log scale and relative to mean for 10 years before drawdown. Nope

Row 2
1) Scatter plots of SST + error (hard to compute) for 30 year timescales relative to after draw down, 100 & 200 years latter. 
   Norm by  Carbon emitted and plotted as a fn of over shoot timescale.
2) Profile of Ocn Temp
3) Profile of C conc.


"""

fig, axes = plt.subplots(2, 2, num="ocn_vars", clear=True, figsize=[7.5, 7],gridspec_kw=dict(wspace=0.2,top=0.95,hspace=0.25,bottom=0.15))
flatten_axes = axes.flatten()
(ax_VAOT, ax_VAOT_C, ax_VAOT_post, ax_C_profile) = flatten_axes


for k, ts in diffs['VAOT'].items():
    series = merlinLib.lookup.loc[k]
    prop = merlinLib.properties(series, linewidth=1)
    name = merlinLib.gen_name(k)
    iris.plot.plot(ts.coord('year'), ts, axes=ax_VAOT, label=name, markevery=50, **prop)

ax_VAOT.set_title('Vol. Avg. Ocean T')
ax_VAOT.set_ylabel(r"$\Delta$ VOAT (K)")
ax_VAOT.set_xlabel(r"Year")
ax_VAOT.axhline(0.0, linestyle='dashed', color='black')
merlinLib.plot_emission(ax=ax_VAOT)  # Plot grey bars to show start and end of emissions.

# plot VAOT in overshoot period vs C*year post emission.
x_data=[]
vaot_data=[]
for exper, series in merlinLib.lookup.iterrows():
    ts = diffs['VAOT'][exper]  # .rolling_window('time', iris.analysis.MEAN, 11)
    x = (ts.coord('year').points - start_emission) * series.Carbon / 1000
    end_t= series.Time - 10
    x_data.append(x[10:end_t])
    vaot_data.append(ts.data[10:end_t])
    prop = merlinLib.properties(series)
    ax_VAOT_C.plot(x[0:end_t], ts.data[0:end_t], label=exper, markevery=25, **prop)

# compute regression & uncert slope from all data

x=np.hstack(x_data)
vaot = np.hstack(vaot_data)

model = sm.OLS(vaot,x)
fit=model.fit()
ce=[float(p*1000) for p in fit.conf_int().flat]

print(f"Fit is {1000*float(fit.params):3.2f} ({ce[0]:3.2f}-{ce[1]:3.2f}) mK/(Eg C Yr)")
print(fit.summary())
# decorate axis.
ax_VAOT_C.set_xlabel('Emission * time (Eg C Year)')
ax_VAOT_C.set_ylabel(r"$\Delta$ VAOT (K)")
ax_VAOT_C.set_title("Transient VAOT")
ax_VAOT_C.set_xlim(0.0)
ax_VAOT_C.set_ylim(0.0)
xr = ax_VAOT_C.get_xlim()
yr = ax_VAOT_C.get_ylim()
ax_VAOT_C.plot(xr, fit.predict(xr), linestyle='dashed', color='black')  # plot best fit

# plot norm VAOT post draw down.
ax_VAOT_post.set_yscale('symlog', linthresh=0.05)  # make VAOT a log plot

for exper, series in merlinLib.lookup.iterrows():
    prop = merlinLib.properties(series, linewidth=1)
    periods = merlinLib.named_periods(series,pnames=['-D'],constraint=True)
    start_year = merlinLib.start_emission+10+series.Time
    ts = diffs['VAOT'][exper]
    norm = ts.extract(periods['-D']).data.mean()  # 30 year mean in the D- period
    ts = ts.extract(iris.Constraint(year=lambda cell: cell >= start_year))
    y = ts / norm  # end of drawdown period.
    x = ts.coord('year').points - start_year  # years since draw down complete.
    ax_VAOT_post.plot(x, y.data, label=exper, markevery=50, **prop)
# decorate the axis
ax_VAOT_post.axhline(0.0, linestyle='dashed', color='black')

ax_VAOT_post.set_xlabel("Years post drawdown")
ax_VAOT_post.set_ylabel(f"Normalised VAOT")
ax_VAOT_post.set_title(f"Post drawdown VAOT")
ax_VAOT_post.set_yticks([-0.1, 0.1, 0.2, 0.5, 1.0])
ax_VAOT_post.set_yticklabels(ax_VAOT_post.get_yticks())
ax_VAOT_post.set_ylim([0.05, 1.2])

# plot the normalised ocean carbon just before drawdown.

for exper, series in merlinLib.lookup.iterrows():
    prop = merlinLib.properties(series, linewidth=1)
    constraints = merlinLib.named_periods(series,constraint=True,pnames=['-D'])

    profile = diffs['OcC_Profile'][exper].extract(constraints['-D']).collapsed('time', iris.analysis.MEAN)
    thick = profile.coord('depth').bounds[:, 1] - profile.coord('depth').bounds[:, 0]
    profile *= 1e6 / series.Carbon  # change per PgC
    profile /= thick # C/m
    iris.plot.plot(profile, profile.coord('depth'), axes=ax_C_profile, **prop)
ax_C_profile.axvline(0, linestyle='dashed')
ax_C_profile.set_title(r"$\Delta$ C ")
ax_C_profile.set_xlabel("Tg/(Eg m)")
ax_C_profile.set_ylabel("Depth (m)")

# done do final figure stuff.
label = merlinLib.plotLabel()  # new set of labels
for ax in flatten_axes:
    label.plot(ax)


h, l = ax_VAOT.get_legend_handles_labels()
fig.legend(h,l,loc='lower left',ncol=6,columnspacing=0.5,fontsize='small')
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)
