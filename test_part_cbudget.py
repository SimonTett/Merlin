"""

Test code for carbon budget partitioning.

"""

import merlinLib
import matplotlib.pyplot as plt
import iris
import iris.plot
import iris.quickplot
import numpy as np
import cf_units

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

hperiods = merlinLib.named_periods(histPulse)
cperiods = merlinLib.named_periods(ctlPulse)
def min_mx_range(period):
    return period.start, period.stop,period.stop-period.start

def land_changes(exper,constraint):
    """
    Compute changes in land fluxes; litter, npp and resp
    """
    df_Litter = part_litter_change(exper, extract=constraint)
    df_NPP = part_litter_change(exper, var='NPP_PFT', extract=constraint)
    df_resp = part_resp_change(exper, extract=constraint)
    return df_Litter,df_NPP,df_resp
for prd_name in (['+E', '-D']):

    # delta = max_yr-min_yr
    rgn_select = iris.Constraint(latitude=lambda lat: -90 <= lat <= 90)
    hmin_yr,hmax_yr,htime_sep = min_mx_range(hperiods[prd_name])
    cmin_yr, cmax_yr, ctime_sep = min_mx_range(cperiods[prd_name])

    hyr_select = iris.Constraint(year=lambda yr:  hmin_yr <= yr < hmax_yr)
    cyr_select = iris.Constraint(year=lambda yr: cmin_yr <= yr < cmax_yr)

    hist_OcC_flux = (hist_OcC.extract(iris.Constraint(year=hmax_yr)) - hist_OcC.extract(
        iris.Constraint(year=hmin_yr))) / htime_sep
    ctl_OcC_flux = (ctl_OcC.extract(iris.Constraint(year=cmax_yr)) - ctl_OcC.extract(
        iris.Constraint(year=cmin_yr))) / ctime_sep
    delta_OcC_flux = hist_OcC_flux - ctl_OcC_flux

    hist_Oc_PrimC_prd = hist_OcPrimC.extract(hyr_select).collapsed('time', iris.analysis.MEAN)
    ctl_Oc_PrimC_prd = ctl_OcPrimC.extract(hyr_select).collapsed('time', iris.analysis.MEAN)

    delta_Oc_PrimC_prd = iris.analysis.maths.subtract(hist_Oc_PrimC_prd, ctl_Oc_PrimC_prd)
    # plot is Resp, 5 PFT + total litter, 5 PFT + total NPP.
    df_Litter_hist, df_NPP_hist, df_resp_hist = land_changes(histPulse.name,hyr_select)
    df_Litter_ctl, df_NPP_ctl, df_resp_ctl = land_changes(ctlPulse.name,hyr_select)


    d_Litter,d_NPP,d_resp = (df_hist - df_ctl for (df_hist,df_ctl) in
                             zip([df_Litter_hist,df_NPP_hist,df_resp_hist],[df_Litter_ctl,df_NPP_ctl,df_resp_ctl]))


    cbar_kws = dict(location='top', use_gridspec=False, )
    # work out fluxes to atmosphere to land, veg, soil and total
    budget_land = pd.Series(dtype='float64')
    budget_veg = pd.Series(dtype='float64')
    budget_soil = pd.Series(dtype='float64')
    budget_atmos = pd.Series(dtype='float64')
    for name, npp, litter, resp, ocn in zip(['Hist', 'Ctl', 'Delta'], [df_NPP_hist, df_NPP_ctl, d_NPP],
                                            [df_Litter_hist, df_Litter_ctl, d_Litter],
                                            [df_resp_hist, df_resp_ctl, d_resp],
                                            [hist_OcC_flux, ctl_OcC_flux, delta_OcC_flux]):
        budget_land.loc[name] = -float(npp.loc['Total', 'Total'] - resp.loc['Total'])
        budget_veg.loc[name] = -float(npp.loc['Total', 'Total'] - litter.loc['Total', 'Total'])
        budget_soil.loc[name] = -float(litter.loc['Total', 'Total'] - resp.loc['Total'])
        budget_atmos.loc[name] = float(
            budget_land.loc[name] - float(ocn.collapsed('depth', iris.analysis.SUM).data))  # budget to atmos
    # make some plots
    gridspec_kw = {'width_ratios': [1, 1, 0.3, 0.3], 'height_ratios': [1, 1, 1]}
    figBudget2, axes = plt.subplots(ncols=4, nrows=3, num=f'C_budget_{hmin_yr}_{hmax_yr}', figsize=merlinLib.fsize,
                                    clear=True, sharex='col', gridspec_kw=gridspec_kw)
    # aflatten = axes.flatten()
    for df, ax, title in zip([df_NPP_hist, df_Litter_hist, df_resp_hist, df_NPP_ctl, df_Litter_ctl, df_resp_ctl],
                             axes[0:2, 0:3].flatten(),
                             ['NPP Hist', 'Litter Hist', 'Resp', 'NPP Ctl', 'Litter Ctl', 'Resp']):
        dataf = df.drop('resid').iloc[0:3, :].rename(index=rename)
        cm = sns.heatmap(dataf.iloc[0:3, :], vmin=-10, vmax=25, annot=True, fmt='2.1f', center=0, ax=ax, cbar=False)
        ax.set_yticklabels(ax.get_yticklabels(), rotation='vertical')
    # plot the differences
    for df, ax, title in zip([d_NPP, d_Litter, d_resp], axes[2, 0:3].flatten(),
                             [r'$\Delta$ NPP', r'$\Delta$ Litter', r'$\Delta$ Resp']):
        dataf = df.drop('resid').iloc[0:3, :].rename(index=rename)
        cm = sns.heatmap(dataf, vmin=-4.5, vmax=4.5, annot=True, fmt='2.1f',
                         center=0, ax=ax, cbar=False)

        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        ax.set_yticklabels(ax.get_yticklabels(), rotation='vertical')

    # plot the Ocean C.
    for var, primFlux, ax in zip([hist_OcC_flux, ctl_OcC_flux, delta_OcC_flux],
                                 [hist_Oc_PrimC_prd, ctl_Oc_PrimC_prd, delta_Oc_PrimC_prd],
                                 axes[:, 3]):
        thick = var.coord('depth').bounds[:, 1] - var.coord('depth').bounds[:, 0]
        v = 1000 * var / thick  # convert to Tg/m/yr
        iris.plot.plot(v, v.coord('depth'), axes=ax)
        ax.axvline(0.0, linestyle='dashed', linewidth=2, color='black')
        ax.set_ylabel('Depth (m)')

        ax.annotate(f'{var.data.sum():3.1f} ', (0.25, 0.6), xycoords='axes fraction',
                    bbox=dict(facecolor='grey', alpha=0.2))
        ax.annotate(f'{primFlux.data:3.1f} ', (0.25, 0.3), xycoords='axes fraction',
                    bbox=dict(facecolor='green', alpha=0.2))

    # put titles on the top row
    for ax, title in zip(axes[0], ['NPP', 'Litter', 'Resp', 'Ocn']):
        ax.set_title(title)
    # label bottom ocean lab
    axes[-1][-1].set_xlabel('Tg C/(m year)')
    axes[-1][-1].set_xlim(-20, 20)  # set limit
    # put titles on the left col
    for ax, title in zip(axes[:, 0], ['Hist', 'Ctl', r'$\Delta$']):
        ax.annotate(title, (-0.12, 0.65), xycoords='axes fraction', rotation=90, fontsize='x-large')
    # put land cpt budgets on top left.
    for ax_list, var in zip(axes.T, [budget_land, budget_veg, budget_soil, budget_atmos]):
        for ax, (name, value) in zip(ax_list, var.items()):
            title = f"{value:3.1f}"
            ax.annotate(title, (-0.01, 1.00), xycoords='axes fraction', fontsize='large')
    # put axis labels on everything.
    label = merlinLib.plotLabel()
    label.plot(ax=axes)

    figBudget2.suptitle(f"C Fluxes (Pg C/year) {hmin_yr:4d} - {hmax_yr:4d}", y=0.999)

    figBudget2.tight_layout()
    merlinLib.saveFig(figBudget2)
    figBudget2.show()
    # end of loop over plots
breakpoint()
## plot time-series of Atmos, Ocean, veg, soil C
diffs = dict()
diff_mean = dict()
vars_to_plot = ['AtmC', 'OcC', 'LandC', 'Forcing', 'VAOT']
linestyles = ['solid', 'dotted', 'dashed']
expers_200 = merlinLib.lookup.query('Time == 200')


def comp_diff(var, exper, info):
    ref_exper = info.Reference
    ref_exper = merlinLib.references[ref_exper]
    ref = merlinLib.read_data(var, ref_exper)
    ts = merlinLib.read_data(var, exper)
    delta = merlinLib.diff(ts, ref)
    return delta


"""for exper, info in merlinLib.lookup.iterrows():
    min_yr = 2018 + info.Time
    max_yr = min_yr + time_seperation
    extract_time = iris.Constraint(year=lambda yr: min_yr <= yr < max_yr)
    vars = dict()
    var_mean = dict()
    for var in vars_to_plot:
        if var == 'Forcing':
            ts = merlinLib.read_data(var, exper)
        else:
            ts = comp_diff(var, exper, info)

        var_mean[var] = ts.extract(extract_time).collapsed('time', iris.analysis.MEAN) * 1000 / info.Carbon
        vars[var] = ts
    # work out the non-cons and other useful terms..
    for var in ['nonCons', 'CO2emis']:
        ts = comp_diff(var, exper, info)
        # convert units from fluxes...
        ts.data = np.cumsum(ts.data)
        ts.units = cf_units.Unit('Pg C')
        vars[var] = ts
        var_mean[var] = ts.extract(extract_time).collapsed('time', iris.analysis.MEAN) * 1000 / info.Carbon
    vars['ResidC'] = vars['AtmC'] + vars['OcC'] + vars['LandC'] - vars['CO2emis']
    var_mean['ResidC'] = vars['ResidC'].extract(extract_time).collapsed('time', iris.analysis.MEAN) * 1000 / info.Carbon
    diffs[exper] = vars
    diff_mean[exper] = var_mean"""

## generate a dataframe
series = []
for exper, info in merlinLib.lookup.iterrows():
    vars = ['AtmC', 'OcC', 'LandC', 'CO2emis', 'nonCons', 'ResidC', 'Forcing']
    data = pd.Series({var: float(diff_mean[exper][var].data) for var in vars})
    series.append(info.append(data).rename(exper))

summary_df = pd.DataFrame(series)
# and print it out (nicely formatted)
with pd.option_context("display.max_columns", 99, "display.precision", 1, 'display.width', 140):
    print(summary_df.loc[[histPulse, ctlPulse], :])

## now plot things

figDF, axis = plt.subplots(nrows=1, ncols=2, num='postDrawDown', figsize=merlinLib.fsize, clear=True)


# scatter plot the rad forcings vs CO2 & Ocn/Land C vs time. Alas pandas does not allow different markers...
def xcoord(time):
    x = np.log(time)
    return x


for e, info in summary_df.iterrows():
    jitter = np.random.uniform(0.95, 1.05)
    props = merlinLib.properties(info, linestyle=None, alpha=1, ms=12)
    axis[0].plot(info.loc['AtmC'], info.loc['Forcing'], **props)
    props = merlinLib.properties(info, linestyle=None, alpha=1, marker='d', ms=12)
    axis[1].plot(xcoord(info.loc['Time'] * 1.1 * jitter), info.loc['LandC'], **props)
    props = merlinLib.properties(info, linestyle=None, alpha=1, marker='o', ms=12)
    axis[1].plot(xcoord(info.loc['Time'] * 0.9 * jitter), info.loc['OcC'], **props)

axis[0].set_xlabel(r'$\Delta $ Atmospheric C ($10^{-3}$ Pg per Pg emitted)')
axis[0].set_ylabel(r'Forcing ($10^{-3}$ Wm$^{-2}$ per Pg emitted)')
axis[0].set_title('Forcing vs Atmospheric C')

times = [50, 100, 200]
axis[1].set_xticks(xcoord(times))
axis[1].set_xticklabels(times)
axis[1].set_xlabel('Pulse Separation (Years) ')
axis[1].set_ylabel("C (Pg per Pg Emitted)")
axis[1].set_title("Ocean/Land Carbon vs Pulse Separation")
for ax in axis:
    ax.axhline(color='black', linestyle='dashed')

figDF.show()
merlinLib.saveFig(figDF)

##
figTS, axis = plt.subplots(nrows=2, ncols=2, sharex=True, num='ts_200', clear=True, figsize=merlinLib.fsize)
for var, ax in zip(vars_to_plot, axis.flatten()):
    for exper, info in expers_200.iterrows():
        prop = merlinLib.properties(info, alpha=1, linewidth=2, marker=None)
        label = f'{info.Carbon} Pg C {info.Reference}'
        delta = diffs[exper][var]
        iris.plot.plot(delta.coord('year'), delta, axes=ax, label=label, **prop)

    ax.axhline(color='black', linestyle='dashed')
    ax.set_title(var)
    ax.set_xlabel('Year')
    merlinLib.plot_emission([2019, 2219], ax)
axis[0][0].legend()

# ax.legend()
figTS.show()
merlinLib.saveFig(figTS)

## see what time- constants are.
import statsmodels.api as sm


def exp_fn(t, *params):
    """

    fn of form C0+C_\infty (1-exp(-kt)) to support timecale estimation.

    """
    order = (len(params) - 1) // 2  # note use of interger divide as slice do not like floats..
    C0 = params[0]
    Cinf = params[1:order + 1]
    kappa = params[order + 1:]

    try:
        result = C0
        for c_val, k in zip(Cinf, kappa):
            result += c_val * (1 - np.exp(-t * k))
    except TypeError:
        result = C0 + Cinf * (1 - np.exp(-t * kappa))

    return result


def exp_decay_fn(t, *params):
    """

    fn of form C_0 exp(-kt) to support timecale estimation.

    """
    order = (len(params) - 1) // 2  # note use of interger divide as slice do not like floats..
    Cinf = params[0]
    C0 = params[1:order + 1]
    kappa = params[order + 1:]

    try:
        result = Cinf
        for c_val, k in zip(C0, kappa):
            result += c_val * np.exp(-t * k)
    except  TypeError:
        result = Cinf + C0 * np.exp(-t * kappa)
    return result


import scipy.optimize

rate_df = []
for exper, info in expers_200.iterrows():
    data = info.copy()
    values = {}
    for var, order in zip(['OcC', 'LandC'], [2, 1]):
        for min_yr, fn, name in zip([2030, 2230],
                                    [exp_fn, exp_decay_fn], ['r', 'd']):

            constraint = iris.Constraint(year=lambda yr: min_yr <= yr < (min_yr - 10 + info.Time))
            var_regress = diffs[exper][var].extract(constraint) / data.Carbon
            x = var_regress.coord('year').points - min_yr
            y = var_regress.data
            p0 = [0.1];
            p0.extend([0.1] * order);
            p0.extend([1e-3] * order)
            bnd_lower = [-1] * (order + 1);
            bnd_lower.extend([-10] * order)
            bnd_upper = [1] * (order + 1);
            bnd_upper.extend([10] * order)

            fit, cov = scipy.optimize.curve_fit(fn, x, y, p0=p0, bounds=(bnd_lower, bnd_upper))
            err = np.mean((fn(x, *fit) - y) ** 2)
            fit_sd = np.sqrt(np.diag(cov))
            # need to move things around so im terms of carbon pool size.
            ordering = np.argsort(fit[1:order + 1])
            c0 = fit[0]
            c = fit[ordering + 1]
            kappa = fit[ordering + 1 + order]
            c0_sd = fit_sd[0]
            c_sd = fit_sd[ordering + 1]
            kappa_sd = fit_sd[ordering + 1 + order]
            key = f'{var}_{name}_err'
            values[key] = err
            key = f'{var}_{name}_C0'
            values[key] = c0
            values[key + '_sd'] = c0_sd
            for indx in range(0, order):
                key = f'{var}_{name}_C{indx + 1:1d}'
                values[key] = c[indx]
                values[key + '_sd'] = c_sd[indx]
                key = f'{var}_{name}_kappa{indx + 1:1d}'
                values[key] = kappa[indx]
                values[key + '_sd'] = kappa_sd[indx]
    # end loop over variables.
    data2 = data.append(pd.Series(values).rename(exper))
    rate_df.append(data2)

rate_df = pd.DataFrame(rate_df)
## and print it out (nicely formatted)
# work out keys..
# cols_want = ['Carbon','Time','Reference']
cols_want = [col for col in rate_df.columns if not col.endswith('sd')]

# with np.printoptions(precision=1):
with pd.option_context("display.max_columns", 99, "display.precision", 0, 'display.width', 140):
    print(rate_df.loc[:, cols_want])

## now make some plots
fig, axis = plt.subplots(nrows=1, ncols=2, num='fit_C_exp', figsize=merlinLib.fsize, clear=True)
# for exper,data in rate_df.loc[[histPulse,ctlPulse],:].iterrows():
for exper, data in rate_df.iterrows():
    style = merlinLib.properties(data)
    style2 = merlinLib.properties(data, linewidth=1, linestyle='dashed', marker=None)
    print(f"{exper} {data.Carbon} {data.Reference}", end=' ')
    for ax, var, order in zip(axis, ['OcC', 'LandC'], [2, 1]):

        delta = diffs[exper][var] / data.Carbon

        for min_yr, fn, name in zip([2030, 2230], [exp_fn, exp_decay_fn], ['r', 'd']):
            constraint = iris.Constraint(year=lambda yr: min_yr <= yr < (min_yr - 10 + data.Time))
            ts = delta.extract(constraint)
            yr = ts.coord('year').points
            # work out vars
            cols = [f'{var}_{name}_C0']
            cols.extend([f'{var}_{name}_C{indx + 1:1d}' for indx in range(0, order)])
            cols.extend([f'{var}_{name}_kappa{indx + 1:1d}' for indx in range(0, order)])
            params = data.loc[cols].values.tolist()
            y = fn(yr - min_yr, *params)
            iris.plot.plot(ts.coord('year'), ts, axes=ax, markevery=20, **style)
            ax.plot(yr, y, **style2)

        ax.axhline(0, color='black', linestyle='dashed')
        ax.set_title(var)
        # print out ratio of decline,
        p0 = iris.Constraint(year=lambda yr: 2220 < yr <= 2230)
        p1 = iris.Constraint(year=lambda yr: 2300 < yr <= 2310)
        ratio = (delta.extract(p1).collapsed('time', iris.analysis.MEAN).data / delta.extract(p0).collapsed('time',
                                                                                                            iris.analysis.MEAN).data)
        print(f"{var} {ratio:3.2f}", end=' ')
    print("")

## plot the ocean carbon.
levels_sig = np.linspace(-450, 450, 10)
levels_diff = np.linspace(-180, 180, 10)
gridspec_kw = {'width_ratios': [8, 1], 'height_ratios': [1] * 3}
yrs_to_show = [2025, 2225, 2300, 2490]
line_cols = ['red', 'brown', 'purple', 'blue']
fig, axis = plt.subplots(nrows=3, ncols=2, num='OceanCarbon', figsize=merlinLib.fsize[::-1], clear=True, sharex='col',
                         gridspec_kw=gridspec_kw)
coords = ['year', 'depth']
thick = hist_OcC.coord('depth').bounds[:, 1] - hist_OcC.coord('depth').bounds[:, 0]
for ax, var, levels, title in zip(axis[:, 0], [hist_OcC, ctl_OcC, delta_OcC], [levels_sig, levels_sig, levels_diff],
                                  ['Historical', 'Control', 'Difference']):
    carbon_per_m = 1000 * var / thick
    cm = iris.plot.contourf(carbon_per_m, coords=coords, axes=ax, levels=levels)
    cs = iris.plot.contour(carbon_per_m, coords=coords, axes=ax, levels=levels, colors='black')
    ax.clabel(cs, fmt='%3.0f')
    for yr, col in zip(yrs_to_show, line_cols):
        ax.axvline(yr, linewidth=3, color=col, alpha=0.8)
    ax.set_title(f'{title} Tg C/m')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Year (CE)')

# plot ocn profiles.
for ax, var, levels, title in zip(axis[:, 1], [hist_OcC, ctl_OcC, delta_OcC], [levels_sig, levels_sig, levels_diff],
                                  ['Historical', 'Control', 'Difference']):
    carbon_per_m = 1000 * var / thick
    for yr, col in zip(yrs_to_show, line_cols):
        c = carbon_per_m.extract(iris.Constraint(year=yr))
        iris.plot.plot(c, c.coord('depth'), axes=ax, color=col)  #
    # ax.set_title(f'{title} Tg C/m')
    # ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Tg C/m')
    t_c = (c * thick / 1000).collapsed('depth', iris.analysis.SUM)
    # ax.set_title(f'{t_c.data:3.0f} Pg C')
    ax.axvline(0.0, color='black')
fig.show()
