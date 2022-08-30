"""
Python script data to move data to their proper place
"""

import pathlib

rootDir=pathlib.Path(r'c:\users\stett2\data\MERLIN\FAMOUS')
files=rootDir.glob('x*.nc')
new_files=[]
for f in files:
    name=f.name
    bits=name.split('-')
    new_path=rootDir/bits[0]/'agg'/bits[1]/name
    print(f"File {str(f)} -> {str(new_path)}")
    f.rename(new_path)
    new_files.append(new_path)

"""
## undo broken attempt
files = rootDir.glob('x*/op?')
for f in files:
    new_file='-'.join([f.parent.name,f.name,'UM_m02s30i292_vn405.0.nc'])
    new_file = rootDir/new_file
    print(f"{str(f)} -> {str(new_file)}")
    f.rename(new_file)
"""
