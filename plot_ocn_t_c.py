"""
Plot the ocean  C ZM plots so can understand what is happening
test case for proplot
"""
proplot=False # set true to use proplot. Note that changes the appearance of the plot..propolot not yet ready for use
#
import iris
import iris.plot
if proplot:
    import proplot as plot
import merlinLib
import numpy as np
import matplotlib.pyplot as plt

refs = {}
diffs = {}
timeseries = {}
for var in ['ocCarbonZM','ocnT_ZM']:
    rr, dd, tt = merlinLib.proc_all(var)
    refs[var] = rr
    diffs[var] = dd
    timeseries[var] = tt

# now we can plot them
start_emission = 2020

## now to process data to give means at particular times
sub=merlinLib.lookup.query("Time == 200 & Carbon == 1000")

times = [-20,20,250]
delta=dict()
historical=dict()
control=dict()
year_start=[]
sub2 = sub.query('Reference=="Historical"')
for name,series in sub2.iterrows():
    L = (sub.Time == series.Time) & (sub.Carbon == series.Carbon) & (sub.Reference == 'Control')
    if L.sum() != 1:
        raise Exception(f"Did not find A single match. Found {L.sum()}")
    c1 = diffs['ocCarbonZM'][name]
    c2 = diffs['ocCarbonZM'][sub[L].iloc[0].name]
    diff= merlinLib.diff(c1,c2)

    delta[name]=[]
    historical[name]=[]
    control[name]=[]
    for t in times:
        yr = start_emission + t + series.Time
        tc = iris.Constraint(year=lambda cell: yr <= cell < yr + 10)
        year_start.append(yr)
        for c, l in zip([c1, c2, diff], [historical[name], control[name], delta[name]]):
            cube = c.extract(tc).collapsed('time', iris.analysis.MEAN)
            l.append(cube)

## now to plot
levels=np.arange(-1.25,1.5,0.25)
#levels=[-2,0,2]
levels=np.array([-1.5,-1,-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.75,1,1.5])*2
cmap = 'BrBG'



# make a normal pyplot figure
fig2, axes2 = plt.subplots(ncols=3,nrows=3,num='zm_plot2',figsize=[11.7,8.3],clear=True,sharex=True,sharey=True)

series=sub2.iloc[0]
exper=series.name


for indx,yr in enumerate(year_start):
    for a,c,name in zip(axes2[indx,:],
                       [historical[exper][indx],control[exper][indx],delta[exper][indx]],
                       ['Historical','Control','Difference']):
        yc = c.coord('model_level_number').points
        xc = c.coord('latitude').points
        data = c.data
        yl = [f"{int(depth):d}" for depth in c.coord('depth').points]
        title = f"{name} {yr:4d}-{yr + 10:4d} Total Pg C: {data.sum():3.0f}"
        # the proplot axis contourf/contour seem to have the wrong labels. This is a big problem...

        a.contourf(xc,yc,data,levels=levels,cmap=cmap)
        cs=a.contour(xc,yc,data,levels=levels,colors='black',linewidths=2)
        a.clabel(cs,fmt='%3.2f',fontsize='medium',inline=True,use_clabeltext=True,inline_spacing=10)
        # decoration
        a.set_yticks(yc[::2])
        a.set_yticklabels(yl[::2])
        a.set_ylim(21, 0.5)
        a.axhline(14, color='black', linestyle='dashed')
        a.xaxis.set_major_formatter(merlinLib.lat_format)
        a.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        a.set_title(title)
        a.set_xlabel('latitude')
        a.set_ylabel('depth (m)')



label=merlinLib.plotLabel()
for ax in axes2.flat:
    ax.label_outer()
    label.plot(ax)

fig2.tight_layout()
fig2.show()
merlinLib.saveFig(fig2)


# proplot version
if proplot:
    fig, axes = plot.subplots(ncols=3,nrows=3,num='zm_plot',figsize=['29.7cm','21.0cm'],tight=False)
    # problems with proplot -- clear =True fails with an error AttributeError: 'NoneType' object has no attribute 'get_text'
    for ax in axes:
        ax.cla()
    for indx, yr in enumerate(year_start):
        for a, c, name in zip(axes[indx, :],
                              [historical[exper][indx], control[exper][indx], delta[exper][indx]],
                              ['Historical', 'Control', 'Difference']):
            yc = c.coord('model_level_number').points
            xc = c.coord('latitude').points
            data = c.data
            yl = [f"{int(depth):d}" for depth in c.coord('depth').points]
            title = f"{name} {yr:4d}-{yr + 10:4d} Total Pg C: {data.sum():3.0f}"
            # the proplot axis contourf/contour seem to have the wrong labels. This is a big problem...
            a.contourf(xc, yc, data, levels=levels, cmap=cmap)
            cs = a.contour(xc, yc, data, levels=levels, colors='black', linewidths=2)
            a.clabel(cs, fmt='%3.2f', fontsize='medium', inline=True, use_clabeltext=True, inline_spacing=10)
            # decoration
            a.set_yticks(yc[::2])
            a.set_yticklabels(yl[::2])
            a.set_ylim(21, 0.5)
            a.axhline(14, color='black', linestyle='dashed')
            a.xaxis.set_major_formatter(merlinLib.lat_format)
            a.set_xticks([-90, -60, -30, 0, 30, 60, 90])
            a.format(xlabel='latitude', ylabel='depth (m)', abc=True, title=title)
    fig.show()
    merlinLib.saveFig(fig)
