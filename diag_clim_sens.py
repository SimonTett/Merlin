"""
Diagnose climate sensitivity of FAMOUS from difference between 2xCO2 (in radn) and fixed CO2 ctl
"""
import iris.analysis
import matplotlib.pyplot as plt
import statsmodels.api as sm
import merlinLib
import numpy as np

experiment='xovfz'
ref = merlinLib.lookup_extra.loc[experiment,'Reference']

net = merlinLib.delta('NetFlux',experiment,refName=ref).rolling_window('time',iris.analysis.MEAN,11)
temp = merlinLib.delta('SAT',experiment,refName=ref).rolling_window('time',iris.analysis.MEAN,11)
netD= net.data
tempD = temp.data
tempD = sm.add_constant(tempD)
fit = sm.RLM(netD,tempD).fit()
print(fit.summary())
fig,ax = plt.subplots(nrows=1,ncols=1,clear=True,figsize=[6,6],num='clim_sens',ms=1)
ax.scatter(net.data,temp.data)
ax.axvline(0,color='black',linestyle='dashed')
# add on best fit... from Forcing est to netFlux = 0
csense = fit.params[0]/(-fit.params[1])# climate sense
T=np.array([0,csense])
T = sm.add_constant(T)
predict_net = fit.predict(T) # predicted netflux from linear relationship
ax.plot(T[:,1],predict_net,linewidth=2,color='red')
ax.axhline(0,color='black',linestyle='dashed')
fig.tight_layout()
fig.show()
