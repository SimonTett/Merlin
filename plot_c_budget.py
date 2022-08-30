"""
Plot contributions to Carbon budgets.

"""
import functools

import merlinLib
import matplotlib.pyplot as plt
import iris
import iris.plot
import iris.quickplot
import numpy as np
import cf_units
import functools

import warnings

warnings.filterwarnings("ignore", message="Collapsing a non-contiguous coordinate.")
warnings.filterwarnings("ignore", message="Collapsing spatial coordinate")
warnings.filterwarnings("ignore", "Using DEFAULT_SPHERICAL_EARTH_RADIUS.")


def delta_var(var, exper, smooth=None):
    series = merlinLib.lookup.loc[exper]
    ref = merlinLib.references[series.Reference]

    ref_cube = merlinLib.read_data(var, ref)
    cube = merlinLib.read_data(var, exper)

    diff = merlinLib.diff(cube, ref_cube)
    if smooth is not None:
        diff = diff.rolling_window('time', iris.analysis.MEAN, smooth)
    return diff


def delta_h_c(var, exper1, exper2, cumsum=False, smooth=None):
    cubes = []

    for exper in [exper1, exper2]:
        series = merlinLib.lookup.loc[exper]
        ref = merlinLib.references[series.Reference]

        ref_cube = merlinLib.read_data(var, ref)
        cube = merlinLib.read_data(var, exper)

        diff = merlinLib.diff(cube, ref_cube)
        cubes.append(diff)
    delta = merlinLib.diff(cubes[0], cubes[1])
    if smooth is not None:
        delta = delta.rolling_window('time', iris.analysis.MEAN, smooth)
    if cumsum:
        d2 = delta.copy()
        d2.data = np.cumsum(d2.data, axis=0)
        delta = d2
    delta.long_name = 'H-R diff ' + var
    print("delta_h_c ", var)
    return delta


def fix_cube(cube):
    """
    Fix cube
    :param cube: cube to be fixed
    :return:
    """

    result = cube.copy()
    try:
        result.data = np.ma.array(data=cube.data.data, mask=cube.data.mask, copy=True)  # this seems needed...
    except AttributeError:
        pass
    result.units = ''
    return result


def comp_delta(perturb_vars, ref_vars):
    """
    Compute difference between perturb and reference variables
    :param

    """

    result_diff = []
    for perturb, ref in zip(perturb_vars, ref_vars):
        delta = merlinLib.diff(perturb, ref)
        result_diff.append(delta)

    return result_diff


## plot heatmaps
import pandas as pd
import seaborn as sns


def mn_rgn(cube, constraint=None):
    # no need to weight as PFT scaling includes areas.
    extract_cube = cube.extract(constraint).collapsed('time', iris.analysis.MEAN)
    result = extract_cube.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    return result.data


def df_pft(pft_var, constraint=None):
    mn_rgn_series = lambda item: pd.Series(mn_rgn(item[1], constraint), merlinLib.pft_labels[0:5]).rename(item[0])
    df = pd.DataFrame(map(mn_rgn_series, pft_var.items()))
    # compute the total and cpts from trees and grasses.
    total = df.sum(axis=1)
    df.loc[:, 'Trees'] = df.loc[:, ['BLT', 'NLT']].sum(axis=1)
    df.loc[:, 'Grasses'] = df.loc[:, ['C3', 'C4']].sum(axis=1)
    df.loc[:, 'Total'] = total
    return df


def series_var(var, constraint=None):
    import functools
    fn = functools.partial(mn_rgn, constraint=constraint)
    series = pd.Series(map(float, map(fn, var.values())), index=var.keys())
    return series


"""
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
"""


def fix_negative_flux(cube):
    """
    Fix negative fluxes in a cube
    :param cube -- input cube
    returns cube with -ve values set to zero.
    """
    result = cube.copy()
    result.data = np.ma.where(cube.data < 0, 0, cube.data)
    return result


def part_resp_change(perturb_sim, extract=None):
    """
    Partition change in respiration into change in rate, change in soilC, from both changes and resid
     Important is that -ve resperation fluxes (arising from conservation moving them into resperation from litter) are set to zero
    :param perturb_sim -- name of perturbation simulation

    :param extract (Optional -- default None). If not None use to extract subset.

    Fn first computes rate (assumed constant) from *mean* values then partition changes and also return resid.
    """

    perturb_resp = fix_negative_flux(merlinLib.read_cube('SoilResp', perturb_sim).extract(extract))
    perturb_soil = merlinLib.read_cube('SoilC', perturb_sim).extract(extract)
    perturb_rate = perturb_resp / perturb_soil

    ref = merlinLib.lookup.loc[perturb_sim, 'Reference']
    ref = merlinLib.references[ref]  # now have name of experiment.
    ref_resp = fix_negative_flux(merlinLib.read_cube('SoilResp', ref).extract(extract))
    ref_soil = merlinLib.read_cube('SoilC', ref).extract(extract)
    ref_rate = (ref_resp / ref_soil)  ##.collapsed('time',iris.analysis.MEAN)
    # and compute time means of reference sim.
    ##ref_resp = ref_resp.collapsed('time',iris.analysis.MEAN)
    ##ref_soil = ref_soil.collapsed('time',iris.analysis.MEAN)

    (delta_rate, delta_soil, delta_resp) = comp_delta([perturb_rate, perturb_soil, perturb_resp],
                                                      [ref_rate, ref_soil, ref_resp])

    rprime_cmn = delta_rate * ref_soil  # alas iris does not do commutative broadcasting.,.,.
    rmn_cprime = delta_soil * ref_rate
    rprime_cprime = delta_soil * delta_rate

    resid = delta_resp - (rprime_cprime + rmn_cprime + rprime_cmn)
    r_cprime = rprime_cprime + rmn_cprime
    # want to show \alpha'C
    result = {"Total": delta_resp, "r'_c": rprime_cmn, "R_c'": r_cprime, "r_c'": rmn_cprime, "r'_c'": rprime_cprime,
              "resid": resid}
    result = pd.DataFrame(series_var(result).rename('Resp'))
    return result


def part_litter_change(perturb_sim, extract=None, var='Litter_PFT'):
    """
    Partition change in litter  into change in PFT area and change in  litter area from both changes and resid,
    Important is that -ve fluxes (arising from conservation moving them into resperation from litter) are set to zero
    :param perturb_sim -- name of perturbation simulation

    :param extract (Optional -- default None). If not None use to extract subset.

    Fn first computes rate  from *mean* values then partition changes and also return resid.
    """
    pl_cons = iris.Constraint(pseudolevel=lambda p: p <= 5)

    # read in the perturbation values.
    v2 = var.replace('_PFT', '')
    perturb_var = fix_negative_flux(merlinLib.read_cube(v2, perturb_sim).extract(extract))  # read in the full fluxes.

    perturb_var_pft = fix_negative_flux(
        merlinLib.read_file(var, perturb_sim).extract(extract & pl_cons)) / 1e12  # Convert to Pg C/year

    perturb_var_pft.units = 'Pg/year/m^2 C'
    perturb_pft = merlinLib.read_cube('PFT', perturb_sim).extract(extract & pl_cons)  # read PFT -- note in m^2
    perturb_var_total = merlinLib.read_cube(var, perturb_sim).extract(extract & pl_cons)
    perturb_resid = perturb_var - perturb_var_total.collapsed('pseudolevel', iris.analysis.SUM)

    # get in the reference values and make time means of them.
    ref = merlinLib.lookup.loc[perturb_sim, 'Reference']
    ref = merlinLib.references[ref]  # now have name of experiment.
    ref_var = fix_negative_flux(merlinLib.read_cube(v2, ref).extract(extract))  # read in the full fluxes.
    ref_var_pft = fix_negative_flux(merlinLib.read_file(var, ref).extract(extract & pl_cons)) / 1e12
    ref_var_pft.units = 'Pg/year/m^2 C'

    ref_pft = merlinLib.read_cube('PFT', ref).extract(pl_cons & extract)
    ref_var_total = merlinLib.read_cube(var, ref).extract(pl_cons & extract)
    ref_resid = ref_var - ref_var_total.collapsed('pseudolevel', iris.analysis.SUM)

    ref_var, ref_var_pft, ref_pft, ref_var_total, ref_resid = \
        [merlinLib.interp_by_coord(var, perturb_var) for var in
         (ref_var, ref_var_pft, ref_pft, ref_var_total, ref_resid)]

    (delta_pft, delta_var_pft, delta_total, delta_resid) = comp_delta(
        [perturb_pft, perturb_var_pft, perturb_var_total, perturb_resid],
        [ref_pft, ref_var_pft, ref_var_total, ref_resid])
    var_prime_area_mn = delta_var_pft * ref_pft
    var_mn_area_prime = delta_pft * ref_var_pft  # iris is non commutative for broadcasting so order matters...
    var_prime_area_prime = delta_var_pft * delta_pft

    resid = delta_total - (var_prime_area_prime + var_mn_area_prime + var_prime_area_mn)
    key = var[0].lower()
    var_area_prime = var_prime_area_prime + var_mn_area_prime
    result = {"Total": delta_total, key + "'_a": var_prime_area_mn, key.upper() + "_a'": var_area_prime,
              key + "_a'": var_mn_area_prime, key + "'_a'": var_prime_area_prime, "resid": resid}
    result = df_pft(result)  # convert to a dataframe.
    return result


ctl = merlinLib.references['Control']
hist = merlinLib.references['Historical']
histPulse = merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Historical"').iloc[0]
ctlPulse = merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Control"').iloc[0]

hist_OcC = delta_var('OcC_Profile', histPulse.name)
ctl_OcC = delta_var('OcC_Profile', ctlPulse.name)
delta_OcC = merlinLib.diff(hist_OcC, ctl_OcC)

hist_OcPrimC = delta_var('OcPrimFluxC', histPulse.name)
ctl_OcPrimC = delta_var('OcPrimFluxC', ctlPulse.name)
delta_OcPrimC = merlinLib.diff(hist_OcC, ctl_OcC)

## let's make some pandas dataframes then plot them
rename = {"N_a'": r"$a' N^P$", "n'_a": r"$a^R N'$",
          "L_a'": r"$a' L^P$", "l'_a": r"$a^R L'$",
          "r'_c": r"$\alpha' C^R$", "R_c'": r"$\alpha^P C'$"}

hperiods = merlinLib.named_periods(histPulse,constraint=True)
cperiods = merlinLib.named_periods(ctlPulse,constraint=True)
def min_mx_range(period):
    return period.start, period.stop,period.stop-period.start

@functools.cache
def land_changes(exper,constraint):
    """
    Compute changes in land fluxes; litter, npp and resp
    """
    df_Litter = part_litter_change(exper, extract=constraint)
    df_NPP = part_litter_change(exper, var='NPP_PFT', extract=constraint)
    df_resp = part_resp_change(exper, extract=constraint)
    return df_Litter,df_NPP,df_resp

#compute fluxes
delta_land =dict()
hist_land=dict()
ctl_land=dict()
for prd_name in (['+E', '+D']):

    # delta = max_yr-min_yr
    rgn_select = iris.Constraint(latitude=lambda lat: -90 <= lat <= 90)

    # plot is Resp, 5 PFT + total litter, 5 PFT + total NPP.
    hist_land[prd_name]=land_changes(histPulse.name,hperiods[prd_name])
    ctl_land[prd_name]= land_changes(ctlPulse.name,cperiods[prd_name])


    delta_land[prd_name] = [df_hist - df_ctl for (df_hist,df_ctl) in
                             zip(hist_land[prd_name],ctl_land[prd_name])]

## now to plot the values
gridspec_kw = dict(width_ratios=[1, 1, 0.25], height_ratios=[1, 1,1,1],wspace=0.175,left=0.1)
figBudget2, axes = plt.subplots(ncols=3, nrows=4, num=f'C_budget', figsize=[7.5,8],
                                    clear=True, sharex='col', gridspec_kw=gridspec_kw)
cbar_kws = dict(location='top', use_gridspec=False, )
for indx,key in enumerate(hist_land.keys()):
    ax_prd = axes[indx*2:indx*2+2,:]
    # work out fluxes to atmosphere to land, veg, soil and total
    budget_land = pd.Series(dtype='float64')
    budget_veg = pd.Series(dtype='float64')
    budget_soil = pd.Series(dtype='float64')
    budget_atmos = pd.Series(dtype='float64')

    for exper, lnd ,ax_row,range in \
            zip([histPulse.name, ctlPulse.name],
                [hist_land[key],ctl_land[key]],
                ax_prd,[[-10,25],[-10,25]]):
        (litter, npp, resp) = lnd
        name = merlinLib.gen_name(exper)
        budget_land.loc[name] = float(npp.loc['Total', 'Total'] - resp.loc['Total'])
        budget_veg.loc[name] = float(npp.loc['Total', 'Total'] - litter.loc['Total', 'Total'])
        budget_soil.loc[name] = float(litter.loc['Total', 'Total'] - resp.loc['Total'])
    # make some plots


        for df, ax in zip([npp,litter,resp],ax_row):
            dataf = df.drop('resid').iloc[0:3, :].rename(index=rename)
            cm = sns.heatmap(dataf.iloc[0:3, :], vmin=range[0], vmax=range[1],
                             annot=True,annot_kws=dict(fontsize='x-small'),
                             fmt='2.2g', center=0, ax=ax, cbar=False,)
            ax.set_yticklabels(ax.get_yticklabels(), rotation='vertical')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        ax_row[0].annotate(name+" "+key, (-0.25, 0.0), xycoords='axes fraction', rotation='vertical',va='bottom')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        #ax.set_yticklabels(ax.get_yticklabels(), rotation='vertical')
    # put land cpt budgets on top right.
    for ax_list, tname,var in zip(ax_prd.T,['Land','Veg','Soil'], [budget_land, budget_veg, budget_soil]):
        for ax, (name, value) in zip(ax_list, var.items()):
            title = f"{tname}: {value:3.1f}"
            ax.annotate(title, (1.0, 1.00), xycoords='axes fraction',ha='right')

# put titles on the top row
for ax, title in zip(axes[0,:], ['NPP', 'Litter', 'Resp']):
    ax.set_title(title,loc='center',pad=12)


# put axis labels on everything.
label = merlinLib.plotLabel()
label.plot(ax=axes)
figBudget2.tight_layout()
figBudget2.show()
merlinLib.saveFig(figBudget2)

    # end of loop over plots