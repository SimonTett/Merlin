# pyton code to plot timescales -- will do fro the 1000 PG  @ 200 year seperation.
import merlinLib
import iris
import iris.plot
import numpy as np
import scipy.stats #,
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
def model_fn(xdata,A0,Ainf,timescale):
    """
    model of the form params[0]+params[1]*(1-exp(-xdata*params[2]))

    :param xadata -- xvalues
    :param *params -- parameters (being optimized)
    """
    

    y= A0 + (Ainf-A0)*(1-np.exp(-xdata/timescale))
    return y
def pwr_fit(y):
    """
    Fit a model of the form y=A_0+(A_infty-A_0)(1-exp(-t/kappa))

    :arg y -- iris timeseries.


    """
    t = y.coord('year').points
    t = t-np.min(t)
    data = y.data
    try:
        fitp,fitcov = scipy.optimize.curve_fit(model_fn,t,y.data,
                                               p0=[0,data.mean(),100],bounds=([-np.inf,-np.inf,1],[np.inf,np.inf,1000]))#,maxfev=20000)
    except RuntimeError:
        fitp = np.repeat(np.nan,3)
    
    sd = np.sqrt(np.diag(fitcov))
    return dict(A0=fitp[0],Ainf=fitp[1],timescale=fitp[2],SD0=sd[0],SDinf=sd[1],SDtimescale=sd[2])
emission_start=2020
start_up = emission_start+10
end_up = emission_start+200
start_down = emission_start+200+10
end_down = 2500
pos_constraint = iris.Constraint(year = lambda yr:  start_up < yr < end_up)
neg_constraint = iris.Constraint(year = lambda yr:  start_down < yr < end_down)
#wanted = merlinLib.lookup.query("(Carbon == 1000) & (Time==200)")
wanted = merlinLib.lookup.query("(Time==200) & (Carbon > 200) ")
# get the C stuff first -- atmosphere; veg; soil; ocean
varsC = dict()
fits = []
for var in ['AtmC','VegC','SoilC','OcC']:
    varsC[var] = merlinLib.delta(var,wanted.index)
    for exper in wanted.index: # loop over experiments

        ts = 100*varsC[var][exper]/wanted.loc[exper,'Carbon'] # per  Pg C
        fit_pos = pwr_fit(ts.extract(pos_constraint))
        fit_pos.update(name=exper,variable=var,which='pos')
        fit_neg = pwr_fit(ts.extract(neg_constraint))
        fit_neg.update(name=exper,variable=var,which='neg')
        fits.extend([pd.Series(fit_pos),pd.Series(fit_neg)])

fitsC = pd.DataFrame(fits).set_index(['name','variable','which'])

varsT = dict()
fitsT = []
for var in ['SAT','VAOT']:
    varsT[var] = merlinLib.delta(var,wanted.index)
    for exper in wanted.index: # loop over experiments

        ts = 1000*varsT[var][exper]/wanted.loc[exper,'Carbon']
        fit_pos = pwr_fit(ts.extract(pos_constraint))
        fit_pos.update(name=exper,variable=var,which='pos')
        fit_neg = pwr_fit(ts.extract(neg_constraint))
        fit_neg.update(name=exper,variable=var,which='neg')
        fitsT.extend([pd.Series(fit_pos),pd.Series(fit_neg)])

fitsT = pd.DataFrame(fitsT).set_index(['name','variable','which'])


## now to plot -

def plot_xticks(ticknames,ax):
    """:
    Plot ticknames nicely separated on an axis -- sure I can do this with the formatter!

    :arg ticknames -- list/array of ticknames to plot
    :arg ax -- axis on which to plot them

    """
    xrange = ax.get_xlim()
    pos = np.linspace(xrange[0], xrange[1], len(ticknames), endpoint=False)
    posticks = pos + (pos[1] - pos[0]) * 0.5
    ax.set_xticks(posticks)
    ax.set_xticklabels(ticknames, rotation='horizontal')
    for p in pos[1:]:
        ax.axvline(p, linestyle='dashed', color='black', linewidth=2)

def plot_stuff(posFit,negFit,ax,colors,var,title,ylabel):
    """
    Plot stuff!

    :param posFit:
    :param negFit:
    :param ax:
    :param colors:
    :param var:
    :param title:
    :param ylabel:
    :return:
    """
    nx = posFit.shape[1]
    posFit.droplevel('name').plot(kind='bar',width=-0.35,y=var,ax=ax,color=colors*nx,align='edge',yerr=varerr,legend=None)
    negFit.droplevel('name').plot(kind='bar', width=0.35, y=var, ax=ax, color=colors * nx, hatch=r'oo',align='edge',
                                  yerr=varerr,legend=None)
    ax.set_xlim(-0.5,None)
    ax.set_ylabel(ylabel)
    ax.axhline(color='black')

    ax.set_title(title)
    plot_xticks(posFit.index.get_level_values('variable').unique(),ax)
    ax.set_xlabel('')


names = {exper:series.Reference for exper,series in wanted.iterrows()}
fitsC = fitsC.rename(index=names) # rename to helpful names,
fitsT = fitsT.rename(index=names) # rename to helpful names,

colors = [merlinLib.properties(series)['color'] for name,series in wanted.iterrows()]
fig, axis = plt.subplots(nrows=3, ncols=3,clear=True,figsize=merlinLib.fsize[::-1],num='timescalesC')
titles=dict(A0='C0',Ainf=r'C$_\infty$',timescale=r'$\tau_C$')

for axx,vars_want in zip(axis[0:2],[['AtmC','OcC'],['VegC','SoilC']]):
    posFit = fitsC.loc[pd.IndexSlice[:,vars_want,'pos'],:].droplevel(['which'])
    negFit = fitsC.loc[pd.IndexSlice[:,vars_want,'neg'],:].droplevel(['which'])
    for ax,var,varerr, ylabel in zip(axx,['A0','Ainf','timescale'],
                                     ['SD0','SDinf','SDtimescale'],
                                     ['%','%','Years']):
        plot_stuff(posFit,negFit,ax,colors,var,titles[var],ylabel)

for ax in axis[0:2,0:2].flatten(): # Make sure the % of C emission is the same on both axis
    ax.set_ylim(-40,100)

for ax in axis[0:2,2].flatten(): # Make sure the timescales are 0, 150
    ax.set_ylim(0,200)

# plot T
posFit = fitsT.loc[pd.IndexSlice[:,:,'pos'],:].droplevel(['which'])
negFit = fitsT.loc[pd.IndexSlice[:,:,'neg'],:].droplevel(['which'])
titles=dict(A0='T0',Ainf=r'T$_\infty$',timescale=r'$\tau_T$')
for ax,var,varerr, ylabel in zip(axis[2].flatten(),['A0','Ainf','timescale'],
                                 ['SD0','SDinf','SDtimescale'],
                                 ['K/ExG','K/ExG','Years']):
    plot_stuff(posFit,negFit,ax,colors,var,titles[var],ylabel)

# set limits up
for ax in axis[2,0:2]:
    ax.set_ylim(-2.5,4)
axis[2,2].set_ylim(0,200) # truncate the timescales show. More than 200 years and we don't have enough data to tell...

# put the labels on.
lab = merlinLib.plotLabel()
for a in axis.flatten():
    lab.plot(ax=a)
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)
