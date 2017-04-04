#!/bin/bash
#PBS -N era_evap
#PBS -P y99
#PBS -q normal
#PBS -l walltime=3:00:00,ncpus=1,mem=30gb
#PBS -j oe
#PBS -v YEAR,COUNT

# Call romscice_evap.py for the given year, and self-submit enough times to
# process the entire ERA-Interim file for that year (the python script only
# process 50 12-hour timesteps at once, otherwise memory overflows).

# To get it started for eg 1992, type
# qsub -v YEAR=1992,COUNT=0 era_evap.job

module unload python/2.7.3
module unload python/2.7.3-matplotlib
module load python/2.7.6
module load python/2.7.6-matplotlib

echo "YEAR = $YEAR"
echo "COUNT = $COUNT"
cd $PBS_O_WORKDIR
python -c "import romscice_evap; romscice_evap.convert_file($YEAR, $COUNT)"

if [ $COUNT -lt 1500 ]; then
    qsub -v YEAR=$YEAR,COUNT=$(($COUNT+100)) era_evap.job
fi

