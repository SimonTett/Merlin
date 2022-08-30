# check emissions totals and other things about start and end.
# Ancillary is set up for 1st of July so will interpolate between those values. Has values for 2020 to 2029.
import merlinLib
import iris
start_emis=2020 # includes August to November from 2019
var='CO2emis'
for exper,series in merlinLib.lookup.iterrows():
    ts=merlinLib.read_data(var,exper)
    ref = merlinLib.references[series.Reference]
    ref_ts = merlinLib.read_data(var,ref)

    delta_ts = merlinLib.diff(ts,ref_ts)
    emit = delta_ts.extract(iris.Constraint(year=lambda cell: (start_emis-1) <= cell <= start_emis+10)).data.sum()
    drawDown = delta_ts.extract(iris.Constraint(year=lambda cell: (start_emis+series.Time-1) <= cell <= (start_emis+series.Time+10))).data.sum()
    total = delta_ts.data.sum()
    print(f"{series.Reference} Nominal: {series.Carbon:4.1f} Total Emitted: {emit:4.1f} Draw Down: {drawDown:4.1f} total: {total:4.1f}")

