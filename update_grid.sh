#!/bin/bash
#PBS -P gh8
#PBS -q express
#PBS -l walltime=2:00:00,ncpus=1,mem=16gb
#PBS -j oe

# After making changes to the grid generating file roms_grid_rtopo2.py,
# rerun this script and also regenerate the initial conditions and boundary
# conditions which depend on the water column thickness. Then submit a test
# run of ROMS-CICE-MCT with the new input files.

cd $PBS_O_WORKDIR
module unload python/2.7.3
module unload python/2.7.3-matplotlib
module load python/2.7.6
module load python/2.7.6-matplotlib

# Make new grid
python roms_grid_rtopo2.py
# Copy into model folder
cp rtopo2_circumpolar_quarterdegree.nc ../ROMS-CICE-MCT/apps/common/grid/
cd ../ROMS-CICE-MCT/python_jobs
# Make new initial conditions
python romscice_ini_woa.py
# Make new lateral boundary conditions for 1995
python -c "from romscice_nbc import convert_file; convert_file(1995)"
# Convert these to repeating annual boundary conditions
mv ../data/ecco2_cube92_lbc_1995.nc ../data/ecco2_cube92_lbc_1995_rep.nc
python romscice_nbc_rep.py
# Submit the 
cd ../apps/common/circumpolar
qsub submit.job