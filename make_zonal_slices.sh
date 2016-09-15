#!/bin/bash

ROMS_FILE="/short/y99/kaa561/roms_spinup_newestbenchmark/ocean1_avg.nc"
for i in `seq 0 10 180`;
do
python -c "from sose_roms_seasonal import *; sose_roms_seasonal('$ROMS_FILE', 'temp', $i, -500, True, 'zonal_slices/temp_""$i""E.png')"
python -c "from sose_roms_seasonal import *; sose_roms_seasonal('$ROMS_FILE', 'salt', $i, -500, True, 'zonal_slices/salt_""$i""E.png')"
done
for i in `seq 170 -10 10`;
do
python -c "from sose_roms_seasonal import *; sose_roms_seasonal('$ROMS_FILE', 'temp', -$i, -500, True, 'zonal_slices/temp_""$i""W.png')"
python -c "from sose_roms_seasonal import *; sose_roms_seasonal('$ROMS_FILE', 'salt', -$i, -500, True, 'zonal_slices/salt_""$i""W.png')"
done