"""
test merlinLib.var_lookup
"""

import merlinLib
import copy

new_vars  ={}
for k,v in merlinLib.var_lookup.items():
    if k.endswith('ZM'): # skip zonal mean keys.
        continue
    k2=k+'_ZM'
    k3 =k+'ZM'
    if (k2 not in merlinLib.var_lookup) and (k3 not in merlinLib.var_lookup):
        collapse = v.get('collapse_dims',[])
        collapse.append('longitude')
        vv = copy.deepcopy(v) # copy the information.
        vv['collapse_dims'] = collapse # add modify the collapse info
        new_vars[k2] = vv

