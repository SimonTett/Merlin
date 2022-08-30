"""
Test PFT version of NPP
"""
import matplotlib.pyplot as plt
import merlinLib
import iris
import cf_units
import iris.quickplot as qplt
import cartopy.crs as ccrs
import numpy as np
import matplotlib.colors as mplcol

npp = merlinLib.read_file('NPP','xocyc')*merlinLib.secs_year # convert to k /year
npp.units = cf_units.Unit('1')
pft = merlinLib.read_file('PFT','xocyc')
npp_pft2 = merlinLib.read_file('NPP_PFT','xocyc') # units need fixing..
npp_pft2.units = cf_units.Unit('1')
total=npp_pft2.collapsed('pseudolevel',iris.analysis.SUM)
total2 = (npp_pft2*pft[:,0:5,:,:]).collapsed('pseudolevel',iris.analysis.SUM)
# try scaling by land-fract.
total2.data *= merlinLib.land_frac.data.squeeze()
npp.data *= merlinLib.land_frac.data.squeeze()
delta = total2-npp
npp_mean = iris.util.squeeze(npp.collapsed('time',iris.analysis.MEAN))
mn=iris.util.squeeze(delta.collapsed('time',iris.analysis.MEAN))
sd = iris.util.squeeze(delta.collapsed('time',iris.analysis.STD_DEV))
## plot
levels_mn=np.linspace(-0.2,0.2,9)
cmap_mn = 'RdYlGn'
cmap_sd='Wistia'
levels_sd = np.linspace(0,1,11)
fig,axis= plt.subplots(2,1,sharex=True,num='test_pft',figsize=merlinLib.fsize,clear=True,subplot_kw=dict(projection=ccrs.PlateCarree()))
cmap = 'Pastel2'

for ax,var,title,levels,cmap in zip(axis.flatten(),[mn,sd],
                                    ['Mean','SD'],[levels_mn,levels_sd],[cmap_mn,cmap_sd]):

    norm = mplcol.BoundaryNorm(levels,levels.size)
    #cm = iris.plot.pcolormesh(var,axes=ax,cmap=cmap,norm=norm)
    cm = iris.plot.pcolormesh(var, axes=ax, cmap=cmap)
    ax.set_title(title)
    ax.coastlines()
    fig.colorbar(cm,ax=ax,orientation='horizontal',fraction=0.1)
fig.show()