"""

Fit simple exp model to carbon fluxes. 

"""



# need sfc ocean, deep ocean and land.
import merlinLib
import iris.plot
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import iris.pandas
import functools
##THINK I NEED STASH DIAGNOSTIC 30206 -- GBM HTN INTO OCEAN BUDGET but nwed to check/

def simple_cc_model(t, y, alpha_OLR=None, alpha_temp=None, alpha_carbon=None, beta_carbon=None, wts=None, forcing_fn=None, atmosCarbon_fn=None):
    """:
    Expression for dy/dt for simple carbon cycle model.
     to be given to scipy.integrate.solve_ivp with vectorize-True
     A
    """

    # split up the input vector

    nlayer = y.shape[0]//2
    npt = y.shape[1]
    temperature = y[0:nlayer,:]
    carbon = y[nlayer:,:]
    if len(wts) != nlayer:
        raise ValueError(f"Wts has size {len(wts)} not {nlayer}")
    if len(alpha_temp) != (nlayer-1):
        raise ValueError(f"alpha_temp should have len {nlayer-1} but has len {len(alpha_temp)}")
    if len(alpha_carbon) != nlayer:
        raise ValueError(f"alpha_carbon should have len {nlayer} but has len {len(alpha_carbon)}")
    atmos_carbon = atmosCarbon_fn(t) # should be a scaler (or maybe I can make  it a vector)
    forcing = forcing_fn(t)
    weights = wts.reshape(nlayer,1)
    secs_per_yr = 60.*60.*24*360. # UM uses 360 day year.
    dydt = np.zeros(y.shape)
    # work out temperature deriv.
    flux_in_top = np.zeros([nlayer+1,npt])
    flux_in_top[0] = (forcing-alpha_OLR*temperature[0,:])*secs_per_yr
    flux_in_top[1:-1] = (10**alpha_temp).reshape(nlayer-1,1) * (temperature[0:-1,:]-temperature[1:,:])
    flux_in_top[-1] =0 # no flux out of the bottom layer.
    # compute the temperature part of the deriv.
    dydt[0:nlayer,:] = (flux_in_top[0:-1,:]-flux_in_top[1:,:])/(weights*4.2e6) # convert J to K.
    # now lets get the carbon cpt.
    flux_in_top[0,:]= (10.0**alpha_carbon[0]) *(atmos_carbon-(carbon[0,:]-beta_carbon[0]*temperature[0,:]))
    flux_in_top[1:-1,:] = (10.0**alpha_carbon[1:]).reshape(nlayer-1,1)*(carbon[0:-1,:]-(carbon[1:,:]-beta_carbon[1:].reshape(nlayer-1,1)*temperature[1:,:]))
    flux_in_top[-1,:] =0 # no flux out of the bottom layer.
    # compute the carbon part of the deriv.
    dydt[nlayer:,:] = (flux_in_top[0:-1]-flux_in_top[1:])  # for ocn just work in Pg.

    return dydt

def unpack_params(param_values,level_names=None):
    """
    unpack parameter values into  dicts -- givens named arguments making interp easier
    See simple_cc_model for details.
    does the following:

    :param param_values: -- a nlayer*3 list or iterable.

    [0] -- alpha_OLR
    [1:nlayer-1]  -- alpha_temp  -- one for each layer but no values needed for bottom layer as diffuses to ground!
    [nlayer:nlayer*2] -- alpha_carbon    one for each layer including atmosphere but no values needed for bottom layer as diffuses to ground!
    [nlayer*2:] -- beta_carbon
    alpha_olr returned in alpha_temp[0].

    if level_names is passed which are names for ocean layers then everything will be returned as pandas Series.
    """
    nlayer = (len(param_values))//3
    if len(param_values) != 3*nlayer:
        raise ValueError(f"Len param_values = {len(param_values)} expected {3*nlayer} ")
    alpha_temp = np.array(param_values[0:nlayer])
    alpha_carbon = np.array(param_values[nlayer:nlayer*2])
    beta_carbon = np.array(param_values[nlayer*2:])

    if level_names is not None:
        # convert everything to pandas series...
        alpha_temp = pd.Series(alpha_temp,index= ['OLR']+level_names[0:-1])
        alpha_carbon = pd.Series(alpha_carbon, index = ['Atmosphere']+level_names[0:-1])
        beta_carbon =  pd.Series(beta_carbon,index=level_names)


    unpacked_args = dict(alpha_temp=alpha_temp,
         alpha_carbon=alpha_carbon,
         beta_carbon=beta_carbon)
    return unpacked_args

def construct_deriv_matrix(alpha_temp=None,alpha_carbon=None,beta_carbon=None,weights=None,index=None):
    """
    Model is dx/dt = M x -- this function computed M from various parameters
    Expect nlayer temperatures for the ocean and nlayer+1 levels for Carbon -- level[0] is the atmosphere.
    Construct matrix for dx/dt = M x from parameters   -- all of which (unless stated otherwise are pandas Series)
    :param alpha_temp -- diffusion constant for Temperature in each layer -- expect nlayer values.    Units are J/K/year
            alpha_temp.loc['OLR'] is the OLR term (OLR=alpha.loc['OLR']*temperature of top level)
    :param alpha_carbon -- diffusion coefficient for Carbon in ocean layer and in atmosphere. Units are  Pgm^-3/Pgm^-3/year. Expect nlayer values
    :param beta_carbon   -- temperature impact on Carbon concentration -- brought in via carbon concentrations. Units are   Pgm^-3/Pgm^-3/year/K
          This is an optional argument. If not specified then no interaction between temperature and carbon will be computed. Expect nlayer-1 values
    :param weights -- thickness of Ocean layers
    :param index -- index for matrix.      (pandas index NOT series) MIGHT be able to get it from data.


    """
    raise NotImplementedError
    # index is mixture of temperature & carbon & depth
    secs_per_year = 60.*60*24*360 # 360 day year for the model
    nIndex = len(index)
    nlevels = (nIndex-1)//2
    model_matrix = pd.DataFrame(np.zeros([nIndex,nIndex]),index=index,columns=index)
    # specials for surface.
    temp_names = index.loc[pd.IndexSlice[:,'Temperature']]
    C_names =  index.loc[pd.IndexSlice[:,'Carbon'] ]            # should include atmosphere.
    top = temp_names[0]
    top_t = pd.IndexSlice[top,'Temperature']
    model_matrix.loc[top_t,top_t] -= alpha_temp.loc['OLR']*secs_per_year     # OLR loss -- only for temperature.
    # diffusion part.
    for var,alpha,level_names in zip(['Temperature','Carbon'],[alpha_temp,alpha_carbon]),[temp_names,C_names]:
        # first set up the top levels
        model_matrix.loc[top_t,pd.IndexSlice[level_names[0],var]] += [-alpha[0],alpha[1]]      # sfc flux. Lose alpha[0]*sfc and gain alpha[1]*layer[1]
        for ii in range(1,len(level_names)-1):  # then go and do the interior values where we interact with levels above and below.
            values = [alpha.iloc[ii-1],-alpha.iloc[ii]-alpha.iloc[ii-1],alpha.iloc[ii] ]
            # gain alpha[ii-1]*(var[ii-1]-var[ii]) and lose (alpha[ii]*(var[ii]-var[ii+1])
            model_matrix.loc[pd.IndexSlice[level_names[ii],var],pd.IndexSlice[level_names[ii-1:ii+2],var]] += values
        # and finally do the bottom values -- zero flux out of the bottom.
        bottom_t =  pd.IndexSlice[level_names[-1],var]      # special for bottom layer-- gain from above and lose from here.
        model_matrix.loc[bottom_t,pd.IndexSlice[level_names[-2:2],var]] += [alpha.iloc[-1],-alpha.iloc[-1]]

    if beta_carbon is not None: # have interaction between temperature & carbon
        # interaction is through change in reference concentration in diffusion.
        # so modify diffusion equations
        # gain from above alpha_C[ii-1]*(C[ii-1]-(C[ii]-beta[ii]*T[ii])) and lose to below alpha_C[ii]*(C[ii]-(C[ii+1]-beta[ii+1]*temp[ii+1]))
        # change is thus +alpha_C[ii-1]*beta[ii]*T[ii]    & -alpha_C[ii]*beta[ii+1]*temp[ii+1]
        # Complication is that have atmospheric C while T has no atmosphere. (Param by beta_OLR*T)
        # so in that case alpha_carbon[1] corresponds to alpha_temp[0]
        top_c = pd.IndexSlice[C_names[0],'Carbon']
        top_t = pd.IndexSlice[temp_names[0],'Temperature']
        model_matrix.loc[top_c,top_t] += -alpha_carbon.loc[C_names[0]] * beta_carbon.loc[depth_names[0]]
        for ii in range(1,len(C_names)):
            T_indx = pd.IndexSlice[temp_names[ii-1:ii+2],'Temperature']
            values = [alpha_carbon.loc[C_names.loc[ii-1]]*beta_carbon.loc[depth_names[ii-1]],
                      -alpha_carbon.loc[C_names.loc[ii]]*beta_carbon.loc[depth_names[ii]]]
            model_matrix.loc[pd.IndexSlice[C_names[ii-1]],T_indx] += values

    model_matrix /= weights # convert E fluxes and Pg C/m^3 to J and PgC per layer.
    model_matrix.loc[pd.IndexSlice[:,'Temperature'],pd.IndexSlice[:,'Temperature']]   /=4.2e6 # scale by heat cap of m^3 of water.

    return model_matrix

def opt_fn(param_values,weights=None,forcing_fn=None,atmosCarbon_fn=None,temperature=None, carbon=None,
           temperature_scale = 1,carbon_scale=1e-7,return_values=False):
    """
    Optimise fit.

    :param param_values:
    :param kwargs:
    :return:
    """
    # extract parameters from param_values -- expect nlayer-1 temp diff params, OLR, nlayer C diff parmas (has an atmosphere), and nlayer
    nlayer = len(param_values)//3 # expect nlayer*3 params
    if temperature.shape[0] != nlayer:
        raise ValueError(f"Temperature shape {temperature.shape} unnexpected. First dim should be {nlayer}")
    if carbon.shape[0] != nlayer:
        raise ValueError(f"Carbon shape {carbon.shape} unnexpected. First dim should be {nlayer}")

    args = unpack_params(param_values)

    args.update( wts=weights,forcing_fn=forcing_fn, atmosCarbon_fn=atmosCarbon_fn)
    integrate_fn = functools.partial(simple_cc_model, **args)
    yr = np.array([v.year for v in temperature.columns.values]) # extract the years -- which is the times we want in the output.

    result = scipy.integrate.solve_ivp(integrate_fn, [yr.min(),yr.max()], np.zeros(nlayer*2), t_eval=yr,method='BDF',vectorized=True)
    if return_values:
        # wrap result up as data-array
        temp_result = temperature.copy()
        temp_result.loc[:,:] = result.y[0:nlayer,:]
        carbon_result = carbon.copy()
        carbon_result.loc[:,:] = result.y[nlayer:,:]
        return temp_result,carbon_result
    else:
        w = np.broadcast_to((weights/weights.sum()).reshape([3,1]),temperature.shape)
        temperature_resid = (result.y[0:nlayer,:]-temperature.values)*temperature_scale*w
        carbon_resid = (result.y[nlayer:, :] - carbon.values) * carbon_scale
        #carbon_resid[:,:]=0.0
        result = np.concatenate([temperature_resid,carbon_resid]).flatten()
        print(np.mean(result**2),np.mean(temperature_resid**2),np.mean(carbon_resid**2))
        return result



    
##
ctl= merlinLib.references['Control']
hist = merlinLib.references['Historical']
histPulse=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Historical"').index[0]
ctlPulse=merlinLib.lookup.query('Carbon==1000 & Time==200 & Reference=="Control"').index[0]



not_polar = iris.Constraint(latitude = lambda lat: lat < 890)
rgns = dict(#polar=iris.Constraint(latitude= lambda lat: lat >= 80),
            shallow=(iris.Constraint(depth= lambda z: z <= 200) & not_polar),
            intermediate=(iris.Constraint(depth= lambda z:  200 < z < 1000) & not_polar),
            deep=(iris.Constraint(depth= lambda z:  1000 < z) & not_polar),
            )


OcCarbon= dict()
OcTemp = dict()

for exper in [histPulse,ctlPulse]:
    OcCarbon[exper],wt=merlinLib.read_cube('OcC',exper,weight=False,return_weights=True)
    OcTemp[exper] = merlinLib.read_cube('OcT', exper,weight=False)
    ref = merlinLib.references[merlinLib.lookup.loc[exper].Reference]
    OcCarbon[ref] = merlinLib.read_cube('OcC', ref,weight=False)
    OcTemp[ref] = merlinLib.read_cube('OcT', ref,weight=False)

##
delta = []


collapse_dims=['latitude','longitude','depth']
vars = ['AtmC','OcSfcFluxC','Forcing']
wts = dict()
scale_wts = dict() # needed to convert carbon concs to total carbon.
for exper in [histPulse, ctlPulse]:
    ref = merlinLib.references[merlinLib.lookup.loc[exper].Reference]
    for name,rgn in rgns.items():
        for var_name, var,collapse in zip(['Carbon','Temperature'],[OcCarbon,OcTemp],[iris.analysis.SUM,iris.analysis.MEAN]):
            weight = wt.extract(rgn).data
            v = var[exper].extract(rgn)
            w = np.broadcast_to(weight.data,v.shape)
            ts= v.collapsed(collapse_dims,collapse,weights=w)
            v = var[ref].extract(rgn)
            w = np.broadcast_to(weight.data, v.shape)
            ref_ts = v.collapsed(collapse_dims,collapse,weights=w)
            deltaS = iris.pandas.as_series(merlinLib.diff(ts,ref_ts))
            deltaS.loc['name'] = name
            deltaS.loc['exper'] = exper
            deltaS.loc['variable'] = var_name
            # Wts for fn are mean thickness of each layer. SUm gets zero from where nothing...
            # for ease (to be done later) will just use area of earth.
            wts[name] = weight.sum()/(merlinLib.areaEarth)
            scale_wts[name] = weight.sum()
            delta.append(deltaS)

    for var in vars:
        ts = merlinLib.read_data(var,exper)
        if var is 'Forcing':
            deltaS = iris.pandas.as_series(ts)
        else:
            ref_ts = merlinLib.read_data(var,ref)
            deltaS = iris.pandas.as_series(merlinLib.diff(ts,ref_ts))
        deltaS.loc['exper'] = exper
        deltaS.loc['variable'] = var
        delta.append(deltaS)




delta = pd.DataFrame(delta).set_index(['name','exper','variable'])

# test compute_deriv_matrix

alpha_temp = pd.Series([1.0,0.5,0.5],index=['OLR']+)
## now work out fluxes... based on
fns=dict()
exper = histPulse
all_vars = vars + ['Carbon','Temperature']
weights = np.array(list(wts.values()))


arguments = dict(alpha_temp=np.repeat(11.,2),alpha_OLR=1.0,
                 alpha_carbon=np.repeat(1.,3),beta_carbon=np.repeat(0,3),wts=weights)
bounds_min =[0.5,5,5,-5,-5,-5,0,0,0]
bounds_max = [3,13,13,3,3,3,0.1,0.1,0.1]
bounds = [bounds_min,bounds_max]
for var in all_vars:
    s = delta.loc[pd.IndexSlice[:,exper,var]] # extract data
    yr = [v.year for v in s.columns.values]
    fns[var] = scipy.interpolate.interp1d(yr,s.values)
arguments.update(forcing_fn=fns['Forcing'],atmosCarbon_fn=fns['AtmC'])
temperature=delta.loc[pd.IndexSlice[:,exper,'Temperature']]
carbon=delta.loc[pd.IndexSlice[:,exper,'Carbon']]
args = dict(weights=weights,forcing_fn=fns['Forcing'],atmosCarbon_fn=fns['AtmC'],
       temperature=temperature,carbon=carbon)
# param values...
param_values =[arguments['alpha_OLR']]
for k in ['alpha_temp','alpha_carbon','beta_carbon']:
    param_values += arguments[k].tolist()
result = opt_fn(param_values,return_values=True,**args)

# try the optimization...
print("optimizing")
opt = scipy.optimize.least_squares(opt_fn,param_values,kwargs = args,verbose=1,ftol=1e-6,xtol=None,gtol=None,diff_step=1e-3,bounds=bounds)

# extract optimal params...
best_temp,best_carbon = opt_fn(opt.x,return_values=True,**args)

## and plot em.
fig,(ax_t,ax_c) = plt.subplots(nrows=2,ncols=1,sharex=True,num='fit',figsize=merlinLib.fsize,clear=True)
col=['red','green','blue']
temperature.T.plot(ax=ax_t,color=col,linestyle='dashed',label=['None']*3)
best_temp.T.plot(ax=ax_t,color=col)
ax_t.legend()

carbon.T.plot(ax=ax_c,color=col,linestyle='dashed')
best_carbon.T.plot(ax=ax_c,color=col)
ax_c.legend()


##
"""fig,ax = plt.subplots(nrows=1,ncols=1,num='ocean_carbon_ts',clear=True)
for exper,col in zip([histPulse,ctlPulse],['red','blue']):
    for name,linestyle in zip(rgns.keys(),['solid','dashed','dashdot','dotted']):
        iris.plot.plot(deltas[name][exper],color=col,linestyle=linestyle,label=exper+'_'+name,axes=ax)
ax.axhline(color='black')
ax.legend()"""

