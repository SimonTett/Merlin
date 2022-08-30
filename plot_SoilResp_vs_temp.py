"""

Plot soil resp vs temperature.


"""
import merlinLib
import matplotlib.pyplot as plt
import iris
import iris.plot
import iris.quickplot
import numpy as np

import warnings
warnings.filterwarnings("ignore", message="Collapsing a non-contiguous coordinate.")
warnings.filterwarnings("ignore",message = "Collapsing spatial coordinate")

def delta_var(var,exper,smooth=None):
    series = merlinLib.lookup.loc[exper]
    ref = merlinLib.references[series.Reference]

    ref_cube = merlinLib.read_data(var, ref)
    cube = merlinLib.read_data(var, exper)

    diff = merlinLib.diff(cube, ref_cube)
    if smooth is not None:
        diff = diff.rolling_window('time',iris.analysis.MEAN,smooth)
    return diff

def delta_h_c(var,exper1,exper2,cumsum=False,smooth=None):
    cubes = []

    for exper in [exper1,exper2]:

        series = merlinLib.lookup.loc[exper]
        ref=merlinLib.references[series.Reference]

        ref_cube = merlinLib.read_data(var,ref)
        cube = merlinLib.read_data(var,exper)

        diff = merlinLib.diff(cube, ref_cube)
        cubes.append(diff)
    delta = merlinLib.diff(cubes[0],cubes[1])
    if smooth is not None:
        delta = delta.rolling_window('time',iris.analysis.MEAN,smooth)
    if cumsum:
        d2 = delta.copy()
        d2.data = np.cumsum(d2.data,axis=0)
        delta = d2
    delta.long_name='H-R diff '+var
    print("delta_h_c ",var)
    return delta

def read_resp_temp(experiment):
    """
    read the ZM SoilResp and land SAT
    :param experiment: the name of the experiment to read from
    :return:  the soilResp and the land SAT
    """
    resp = merlinLib.read_data('SoilResp_ZM', experiment)
    temp = merlinLib.read_data('Land_SAT_ZM', experiment)
    temp = temp - temp.collapsed('time',iris.analysis.MEAN)

    resp = resp/resp.collapsed('time',iris.analysis.MEAN)
    return resp,temp

def fix_cube(cube):
    """
    Fix cube
    :param cube: cube to be fixed
    :return:
    """

    result = cube.copy()
    try:
        result.data= np.ma.array(data=cube.data.data,mask=cube.data.mask,copy=True) # this seems needed...
    except AttributeError:
        pass
    result.units=''
    return result


def comp_delta(perturb_vars, ref_vars):
    """
    Compute difference between perturb and reference variables
    :param

    """

    result_diff = []
    for perturb,ref in zip(perturb_vars,ref_vars):
        delta = merlinLib.diff(perturb,ref)
        result_diff.append(delta)

    return result_diff


def part_resp_change(perturb_sim, extract=None):
    """
    Partition change in respiration into change in rate, change in soilC, from both changes and resid
    :param perturb_sim -- name of perturbation simulation

    :param extract (Optional -- default None). If not None use to extract subset.

    Fn first computes rate (assumed constant) from *mean* values then partition changes and also return resid.
    """

    perturb_resp = merlinLib.read_cube('SoilResp_ZM', perturb_sim).extract(extract)
    perturb_soil = merlinLib.read_cube('SoilC_ZM', perturb_sim).extract(extract)
    ref = merlinLib.lookup.loc[perturb_sim, 'Reference']
    ref = merlinLib.references[ref] # now have name of experiment.
    ref_resp = merlinLib.read_cube('SoilResp_ZM',ref).extract(extract)
    ref_soil = merlinLib.read_cube('SoilC_ZM',ref).extract(extract)

    perturb_rate = perturb_resp/perturb_soil
    #perturb_rate = perturb_resp / perturb_soil
    ref_rate = ref_resp/ ref_soil
    (delta_rate,delta_soil,delta_resp)= comp_delta([perturb_rate,perturb_soil,perturb_resp],[ref_rate,ref_soil,ref_resp])

    rprime_cmn = ref_soil*delta_rate
    rmn_cprime = delta_soil*ref_rate

    rprime_cprime = delta_soil*delta_rate

    resid = delta_resp- (rprime_cprime+rmn_cprime+rprime_cmn)

    result = {"Total":delta_resp,"r'_c":rprime_cmn,"r_c'":rmn_cprime,"r'_c'":rprime_cprime,"resid":resid}
    return result



def part_litter_change(perturb_sim, extract=None,var='Litter_PFT'):
    """
    Partition change in respiration into change in rate, change in soilC, from both changes and resid
    :param perturb_sim -- name of perturbation simulation

    :param extract (Optional -- default None). If not None use to extract subset.

    Fn first computes rate (assumed constant) from *mean* values then partition changes and also return resid.
    """
    pl_cons = iris.Constraint(pseudolevel = lambda p: p<=5)

   # read in the perturbation values.
    perturb_var_pft = merlinLib.read_file(var, perturb_sim).extract(extract).extract(pl_cons)/1e12 # leave in raw units
    perturb_var_pft.units='Pg/year/m^2 C'
    perturb_pft = merlinLib.read_cube('PFT', perturb_sim).extract(extract).extract(pl_cons) # read PFT -- note in m^2
    perturb_var_total = merlinLib.read_cube(var,perturb_sim).extract(extract).extract(pl_cons)
   # get in the reference values
    ref = merlinLib.lookup.loc[perturb_sim, 'Reference']
    ref = merlinLib.references[ref] # now have name of experiment.
    ref_var_pft = merlinLib.read_file(var,ref).extract(extract).extract(pl_cons)/1e12
    ref_var_pft.units = 'Pg/year/m^2 C'
    ref_pft = merlinLib.read_cube('PFT',ref).extract(extract).extract(pl_cons)
    ref_var_total = merlinLib.read_cube(var,ref).extract(extract).extract(pl_cons)
    (delta_pft,delta_var_pft,delta_total)=comp_delta([perturb_pft,perturb_var_pft,perturb_var_total],
                                                 [ref_pft,ref_var_pft,ref_var_total])

    var_prime_area_mn = delta_var_pft*ref_pft
    var_mn_area_prime = delta_pft*ref_var_pft # order matters as iris does not really handle broadcasting..
    var_prime_area_prime = delta_var_pft*delta_pft

    resid =  delta_total- (var_prime_area_prime+var_mn_area_prime+var_prime_area_mn)
    key=var[0].lower()
    result = {"Total":delta_total,key+"'_a":var_prime_area_mn,key+"_a'":var_mn_area_prime,key+"'_a'":var_prime_area_prime,"resid":resid}
    return result



ctl= merlinLib.references['Control']
hist = merlinLib.references['Historical']
histPulse=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Historical"').index[0]
ctlPulse=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Control"').index[0]

ctl_resp,ctl_temp = read_resp_temp(ctl)

hist_resp,hist_temp = read_resp_temp(hist)
min_yr = 2030
max_yr = 2080
cons = iris.Constraint(year=lambda yr: min_yr < yr < max_yr)
part_hist_1000Pulse = part_resp_change(histPulse,extract=cons)
part_ctl_1000Pulse = part_resp_change(ctlPulse,extract=cons)

part_litter_hist_1000Pulse = part_litter_change(histPulse,extract=cons)
part_litter_ctl_1000Pulse = part_litter_change(ctlPulse,extract=cons)

part_npp_hist_1000Pulse = part_litter_change(histPulse,var='NPP_PFT',extract=cons)
part_npp_ctl_1000Pulse = part_litter_change(ctlPulse,var='NPP_PFT',extract=cons)

total_litter_diff = delta_h_c('Litter_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)/5
total_litter_PFT_diff = delta_h_c('Litter_PFT_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)/5
total_litter_PFT_diff = delta_h_c('Litter_PFT_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)/5
total_resp_diff = delta_h_c('SoilResp_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)/5
total_npp_diff = delta_h_c('NPP_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)/5
## plot contributions to litter & resp.
fig,axis = plt.subplots(nrows=2,ncols=3,num='resp budget',sharex=True,sharey=True,figsize=merlinLib.fsize,clear=True)
ax=axis.flatten()
total_cs = 0.0
total_ts =0.0
total_cs_litter =0.0
total_ts_litter =0.0
ax_litter_count=1
dims = ['latitude', 'longitude']
for kresp,color in zip(part_hist_1000Pulse.keys(),['red','green','blue','orange']):
    hist_resp = part_hist_1000Pulse[kresp]/5
    ctl_resp  = part_ctl_1000Pulse[kresp]/5
    ts_v = hist_resp.collapsed(dims,iris.analysis.SUM)
    cs_hist = hist_resp.collapsed('longitude',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
    ts_c = ctl_resp.collapsed(dims,iris.analysis.SUM)
    cs_ctl  = ctl_resp.collapsed('longitude',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
    ts_diff = ts_v - ts_c
    cs_diff = cs_hist - cs_ctl
    iris.plot.plot(cs_diff, color=color, label=kresp,axes=ax[0])


    hist_v = ts_v.collapsed('time',iris.analysis.MEAN).data
    ctl_v = ts_c.collapsed('time',iris.analysis.MEAN).data
    diff = hist_v - ctl_v
    print(f"{kresp}: {hist_v: 3.1f} {ctl_v:3.1f} {diff:3.1f}")

iris.plot.plot(total_resp_diff, label='Total Resp.', axes=ax[0], color='black', linewidth=2)
ax[0].legend(ncol=2,loc='upper left')
ax[0].set_title('Soil Respiration')

for klitter,knpp,a in zip(list(part_litter_hist_1000Pulse.keys())[0:-1],part_npp_hist_1000Pulse.keys(),ax[1:]):

    hist_litter = part_litter_hist_1000Pulse[klitter]/5e12
    ctl_litter = part_litter_ctl_1000Pulse[klitter]/5e12
    ts_v = hist_litter.collapsed(dims, iris.analysis.SUM)
    cs_hist = hist_litter.collapsed('longitude', iris.analysis.SUM).collapsed('time', iris.analysis.MEAN)
    ts_c = ctl_litter.collapsed(dims, iris.analysis.SUM)
    cs_ctl = ctl_litter.collapsed('longitude', iris.analysis.SUM).collapsed('time', iris.analysis.MEAN)
    ts_diff = ts_v - ts_c
    cs_diff = cs_hist - cs_ctl
    merlinLib.plot_pft(cs_diff,ax=a,xcoord='latitude',marker=None)
    iris.plot.plot(total_litter_diff, label='Total', color='black', axes=a,linewidth=3)
    hist_v = ts_v.collapsed('time', iris.analysis.MEAN).data.sum()
    ctl_v = ts_c.collapsed('time', iris.analysis.MEAN).data.sum()
    diff = hist_v - ctl_v
    total_cs_litter += cs_diff
    total_ts_litter += ts_diff
    print(f"{klitter}: {hist_v: 3.1f} {ctl_v:3.1f} {diff:3.1f}")
    # plot the NPP terms
    hist_npp = part_npp_hist_1000Pulse[knpp] / 5e12
    ctl_npp = part_npp_ctl_1000Pulse[knpp] / 5e12

    cs_hist = hist_npp.collapsed('longitude', iris.analysis.SUM).collapsed('time', iris.analysis.MEAN)

    cs_ctl = ctl_npp.collapsed('longitude', iris.analysis.SUM).collapsed('time', iris.analysis.MEAN)
    cs_diff = cs_hist - cs_ctl
    merlinLib.plot_pft(cs_diff, ax=a, xcoord='latitude', linestyle='None',label=None)
    iris.plot.plot(total_npp_diff, color='black', marker='o',axes=a, linestyle='None',ms=10)

    a.set_title(klitter)
    iris.plot.plot(total_resp_diff, label='Total Resp.', axes=a, color='gray', linewidth=2,linestyle='dashed')



# add figure legend for PFT names. Do by grabbing legend handles from *last* plot and using those.
lines=a.get_legend_handles_labels()
legend = fig.legend(lines[0],lines[1],loc='lower left',ncol=8)
# add std axis stuff.
for a in ax:
    a.axhline(color='black',linestyle='dashed')
    a.set_ylim(-0.3,0.3)
    a.set_xlim(30,80)
    a.set_ylabel('Pg C/(yr deg)')
    a.set_xlabel('Latitude')
    a.set_xticks([30,45,60,75])
    a.xaxis.set_major_formatter(merlinLib.lat_format)
    a.margins(x=0.1)
plt.show()
## do next plot!

## let's make some pnadas dataframes then plot them
rgn_select = iris.Constraint(latitude = lambda lat:  30 <= lat <= 60)
# plot is Resp, 5 PFT + total litter, 5 PFT + total NPP.
figBudget,axis = plt.subplots(nrows=1, ncols=1, figsize=merlinLib.fsize,num='cBudget',clear=True)
width=0.07
delta_offset =0.08
last_y = 0
offset =0.0
npft = 5

## plot heatmaps
import pandas as pd
import seaborn as sns
def mn_rgn(cube):
    # no need to weight as PFT scaling includes areas.
    extract_cube = cube.extract(cons & rgn_select).collapsed('time',iris.analysis.MEAN)
    #wt = irpart_is.analysis.cartography.area_weights(extract_cube)
    result = extract_cube.collapsed(['longitude','latitude'],iris.analysis.SUM)
    return result.data

def df_pft(pft_var):
    mn_rgn_series = lambda item: pd.Series(mn_rgn(item[1]), merlinLib.pft_labels[0:5]).rename(item[0])
    df = pd.DataFrame(map(mn_rgn_series, pft_var.items()))
    # compute the total and cpts from trees and grasses.
    total = df.sum(axis=1)
    df.loc[:, 'Trees'] = df.loc[:, ['BLT', 'NLT']].sum(axis=1)
    df.loc[:, 'Grasses'] = df.loc[:, ['C3', 'C4']].sum(axis=1)
    df.loc[:, 'Total'] = total
    return df

def series_var(var):
    series = pd.Series(map(float, map(mn_rgn, var.values())), index=var.keys())
    return series

def df_delta_pft(pft_hist,pft_ctl):

    df_hist= df_pft(pft_hist)
    df_ctl = df_pft(pft_ctl)
    delta = df_hist - df_ctl

    return delta

def series_delta(hist,ctl):
    series_hist =series_var(hist)
    series_ctl = series_var(ctl)
    delta = series_hist - series_ctl
    return delta

d_litter=df_delta_pft(part_litter_hist_1000Pulse,part_litter_ctl_1000Pulse)
d_npp=df_delta_pft(part_npp_hist_1000Pulse,part_npp_ctl_1000Pulse)
d_resp = series_delta(part_hist_1000Pulse,part_ctl_1000Pulse)
cbar_kws=dict(location='top',use_gridspec=False,)
use_cbar=True
figBudget,axes = plt.subplots(2,1,num='C_budget_hm',figsize=merlinLib.fsize[::-1],clear=True,sharex=True)
for var,title,ax in zip([d_litter,d_npp],[r'$\Delta$ Litter (Pg C/yr)',r'$\Delta$ NPP (Pg C/yr)'],axes.flatten()):
    cm=sns.heatmap(var.drop('resid',axis='index'),vmin=-4.5,vmax=4.5,annot=True,fmt='2.1f',
                   square=True,center=0,ax=ax,cbar=False)
    ax.set_title(title)

    ax.set_ylabel('Decomp.')
ax.set_xlabel('PFT')
figBudget.colorbar(axes.flatten()[0].get_children()[0],ax=axes,orientation='horizontal')
#figBudget.tight_layout()
figBudget.show()

## end of heatmap plot

soil_resp =  pd.Series(map(float,map(mn_rgn,part_hist_1000Pulse.values())),part_hist_1000Pulse.keys())
litter = pd.DataFrame(map(pd.Series,map(mn_rgn,part_litter_hist_1000Pulse.values()),part_hist_1000Pulse.keys()))
colors=['red','orange','green','blue']
for respKey,litterKey,NPPKey,color in zip(part_hist_1000Pulse.keys(),
                                          part_litter_hist_1000Pulse.keys(),part_npp_hist_1000Pulse.keys(),colors):

    delta_resp = mn_rgn(part_hist_1000Pulse[respKey]-part_ctl_1000Pulse[respKey])
    delta_litter = mn_rgn(part_litter_hist_1000Pulse[litterKey]-part_litter_ctl_1000Pulse[litterKey])/1e12
    delta_npp = mn_rgn(part_npp_hist_1000Pulse[NPPKey]-part_npp_ctl_1000Pulse[NPPKey])/1e12

    tot_delta_litter = delta_litter.collapsed('pseudolevel', iris.analysis.SUM)
    tot_delta_npp = delta_npp.collapsed('pseudolevel',iris.analysis.SUM)
    litter = np.append(delta_litter.data, tot_delta_litter.data).reshape(npft+1,1)
    npp = np.append(delta_npp.data, tot_delta_npp.data).reshape(npft+1,1)
    pft= np.concatenate((litter,npp),axis=1)
    pft_x = pft.copy() # copy values
    pft_x[:,0]=np.arange(npft+1)+1
    pft_x[:, 1] =pft_x[:,0] +0.4
    pft_alpha = pft_x.copy()
    pft_alpha[:,0]=1
    pft_alpha[:,1]=0.5

    y = np.concatenate((delta_resp.data.reshape(1),pft.flatten()))
    indx = np.concatenate((np.array([0.0]).reshape(1),pft_x.flatten()))
    alpha = np.concatenate((np.array([1.0]).reshape(1),pft_alpha.flatten()))
    axis.bar(indx+offset,y,width,color=color)
    # add some annotation...
    # resp first!
    dx=0.0
    dy=0.5
    arrowprops = dict(arrowstyle=r'->')
    pos = [0,-2,-1]
    ytext = y[pos]
    xtext = indx[pos]+offset
    for text,x,y in zip([respKey,litterKey,NPPKey],xtext,ytext):
        axis.annotate(text,(x,y),(x+dx,y+np.sign(y)*dy),ha='center',
                  arrowprops=arrowprops)
    last_y += y
    offset +=delta_offset

#axis.bar(indx+offset,last_y,width,color='black')
axis.axhline(color='black',linestyle='dashed')
xlabels = ['Resp']
xlabels.extend(merlinLib.pft_labels[0:npft])
xlabels.append('Total')
axis.set_xticks(np.arange(len(xlabels)))
axis.set_xticklabels(xlabels)
figBudget.show()

##
fig,ax = plt.subplots(nrows=1,ncols=2,num='scatter',figsize=merlinLib.fsize,sharex=True,sharey=True,clear=True)
size=5
for indx in range(5,40,5,):
    ax[0].scatter(ctl_temp.data[:,indx],ctl_resp.data[:,indx]-1,s=size)
    ax[1].scatter(hist_temp.data[:, indx],hist_resp.data[:, indx] - 1,s=size)

for a,title in zip(ax,['Control','Historical']):
    a.set_title(title)
    a.set_ylabel('$\Delta$ Resp ')
    a.set_xlabel('$\Delta$ Temp')
    xlim = a.get_xlim()
    x=np.linspace(xlim[0],xlim[1])
    y=np.log(2.)*0.1*x
    a.plot(x,y,linestyle='dashed',color='black',linewidth=2)


fig.show()
PFT_c=dict()
PFT_Litter_c=dict()
litter_PFT_rate=dict()
soil_c=dict()
soil_resp=dict()
soil_rate=dict()
NPP = dict()
sfc_temp=dict()

ref_colors=['red','brown','blue','lightblue']

exper_list=[hist,histPulse,ctl,ctlPulse]
delta_Litter=delta_h_c('Litter_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)
delta_Resp = delta_h_c('SoilResp_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN)
delta_soil = delta_h_c('SoilC_ZM',histPulse,ctlPulse).extract(cons).collapsed('time',iris.analysis.MEAN) # actual change in soil carbon fluxes (dbl diff)

hist_partion = part_resp_change()



for exper in exper_list:
    # compute expected soil resp rates from ctl.
    soil_c[exper] = merlinLib.read_data('SoilC_ZM', exper).extract(cons).collapsed('time', iris.analysis.MEAN)
    soil_resp[exper] = merlinLib.read_data('SoilResp_ZM', exper).extract(cons).collapsed('time', iris.analysis.MEAN)

    soil_rate[exper] = fix_cube(soil_resp[exper]/soil_c[exper])
    sfc_temp[exper] = merlinLib.read_data('Land_SAT_ZM',exper).extract(cons).collapsed('time', iris.analysis.MEAN)
    # work out expected litter rate PER PFT.

    PFT_c[exper] = merlinLib.read_data('VegC_PFT_ZM', exper).extract(cons).collapsed('time', iris.analysis.MEAN)
    PFT_Litter_c[exper] = merlinLib.read_data('Litter_PFT_ZM', exper).extract(cons).collapsed('time', iris.analysis.MEAN)
    # "fix" -ve values (which are compensating non-cons)
    PFT_Litter_c[exper].data = PFT_Litter_c[exper].data.clip(0.0, None)

    litter_PFT_rate[exper] = fix_cube(PFT_Litter_c[exper] / PFT_c[exper])

soil_rate_diff = (soil_rate[hist]*delta_soil)
soilC_diff = (soil_c[hist]-soil_c[ctl] ) # carbon diff
# compute rate from soil_rate_diff acting on total C.

expect_delta_Resp = soil_rate[hist]*delta_soil
expect_delta_Resp2 = soil_rate[ctl]*delta_soil
fig_explain,ax = plt.subplots(nrows=3,ncols=1,num='soil',figsize=merlinLib.fsize[::-1],sharex=True,clear=True)


# pft stuff -- should  move into merlinlib...Sort of 1/2 way there already too!
colors = dict(BLT='green', NLT='darkgreen', C3='goldenrod', C4='orange', Shrub='palegreen', Urb='gray',
              Wat='blue', Soil='brown', ice='skyblue',Resid='blue',Total='grey')
linestyles = dict(BLT='solid', NLT='dashed', C3='solid', C4='dashed', Shrub='solid', Urb='dotted',
                 Wat='dotted', Soil='solid', ice='dotted')
pft_indx = {k:ind for ind,k in enumerate(colors.keys())}


for exper,colour in zip(exper_list,ref_colors):
    iris.quickplot.plot(soil_rate[exper], axes=ax[0],linewidth=2,label=exper,color=colour)
    #pass

for exper,alpha in zip([hist,ctl],[1,0.5]):
    for name in ['BLT','NLT','C3','C4','Shrub']:
        ind = pft_indx[name]
        style=dict(linewidth=2, marker='o', color=colors[name],alpha=alpha,linestyle=linestyles.get(name,'solid'))
        iris.quickplot.plot(litter_PFT_rate[exper][ind,:],label=name,axes=ax[1],**style)
ax[1].set_yscale('log')
ax[1].set_ylim(5e-2,10)
ax[1].set_title(f'Est Soil Resp. rate ')

iris.quickplot.plot(delta_Resp,label='Resp Diff',axes=ax[2])
expect_delta_Resp.units = delta_Resp.units
iris.quickplot.plot(expect_delta_Resp,label='Expect Resp Diff',axes=ax[2])
iris.quickplot.plot(expect_delta_Resp2,label='Expect Resp Diff ctl',axes=ax[2])
ax[2].legend()
for a in ax: # put some ticks
    a.set_xlabel('Latitude')
    a.set_xticks([-90,-60,-30,0,30,60,90])
    a.xaxis.set_major_formatter(merlinLib.lat_format)
    a.set_ylabel('')

fig_explain.tight_layout()
fig_explain.show()
