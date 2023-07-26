#!/bin/bash
# run aggregation script on annual datafiles. Assumed to be running on datastore
# run as run_agg.sh experiment 
# for example run_agg.sh xaaaa
# This will need editing to reflect where you want data to go

export HDF5_USE_FILE_LOCKING=FALSE # so netcdf reads etc work on datastore.

dirname=$1 ; shift
rootdir=/exports/csce/datastore/geos/groups/merlin/top/stett2/ # where data is found
aggdir=$rootdir/AggData # where you want data to go.
script=$aggdir/transpose_script/splitAgg.py #path to agg script
for dir in apx opx opy apy
do
    echo "Processing $aggdir $dirname $dir"
    pattern="$dirname-$dir-"
    file=${dir::-1}y_varlist_Viv.txt # varlist file
    $script -i $rootdir/famous/$dirname/nc/$dir/*.nc -p $pattern -l $aggdir/transpose_script/$file -o $dirname/agg/$dir
done



