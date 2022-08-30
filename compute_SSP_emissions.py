"""
Work out total emissions for SSP scenarios
"""
import pandas as pd
import pathlib
import numpy as np

file = 'iamc_db.xlsx'
C_emmis = pd.read_excel(file,index_col=0,nrows=4,usecols=[1,5,6,7,8,9,10,11,12,13,14])*12/44*1e-3 # convert to Pg C C.
#need to interpolate data
for yr in range(2015,2100):
    if yr in C_emmis.columns:
        continue # got it so skip
    C_emmis[yr]=np.nan

C_interp = C_emmis.reindex(columns=np.sort(C_emmis.columns)).interpolate(axis=1).T
