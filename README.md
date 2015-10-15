Python and Matlab scripts for pre- and post-processing of ROMS simulations. Read the comments at the beginning of each file to find out what they do.

Created for the circumpolar Antarctic configuration of ROMSIceShelf/ROMSSeaIce with a quarter-degree grid; should be easy to adapt to other domains.

To run roms_cice_atm_subdaily.py: open ipython or python, type "from romscice_atm_subdaily.py import *", then "convert_file(year)" for the year of ERA-Interim data you want to process. File paths will need to be hard-coded into the file.

To run any other Python script file.py: open ipython, then type "run file". You will be prompted to type in file paths and parameters at the command line.
To run Matlab script function.m: open Matlab, then type "function()". File paths and parameters will need to be hard-coded into the file.


