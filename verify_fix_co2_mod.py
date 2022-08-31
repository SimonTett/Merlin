"""
Code to verify that fixed CO2 mod is working as expected
"""
import  matplotlib.pyplot as plt
import merlinLib
import iris.plot
sim = 'xovfz' # name of simulation
ref = 'Control' # reference name
fig,axes=plt.subplots(nrows=3,ncols=1,clear=True,sharex=True,figsize=[7,10],num='verify_fix_co2_mod')
# plot SAT delta
ts= merlinLib.delta('SAT',sim,refName=ref)
iris.plot.plot(ts.coord('year'),ts,axes=axes[0])
axes[0].set_title('SAT')
# then Carbon
sum=0
for var in ['AtmC','VegC','SoilC','OcC']:
    ts = merlinLib.delta(var, sim, refName=ref)
    iris.plot.plot(ts.coord('year'), ts, axes=axes[1],label=var)
    sum += ts
iris.plot.plot(sum.coord('year'), sum, axes=axes[1],label='Total')
axes[1].set_title('Carbon')
axes[1].legend(ncol=2)
axes[1].axhline(0,color='black',linestyle='dashed')
cs = merlinLib.delta('AtmT_Profile',sim, refName=ref)
ts= cs.extract(iris.Constraint(air_pressure=30))
iris.plot.plot(ts.coord('year'), ts, axes=axes[2])
#iris.plot.pcolormesh(cs,coords=['year','air_pressure'],axes=axes[2])
axes[2].set_title("30 hPa Temperature")

fig.tight_layout()
fig.show()