Python and Matlab scripts for pre- and post-processing of ROMS simulations. Read the comments at the beginning of each file to find out what they do.

Created for the circumpolar Antarctic configuration of ROMS-CICE-MCT with a quarter-degree grid; should be easy to adapt to other domains.

romscice_atm_subdaily.py is designed to be run through a self-submitting batch script; see convert_era.job for an example. File paths will need to be hard-coded into the file.

romscice_ini.py also needs file paths hard-coded, but you can run it normally (in ipython, type "run romscice_ini".

To run any other Python script file.py: open ipython, then type "run file". You will be prompted to type in file paths and parameters at the command line.

To run Matlab script function.m: open Matlab, then type "function()". File paths and parameters will need to be hard-coded into the file. NB: Use the 2008 version of Matlab on katabatic. These scripts depend on existing ROMS Matlab tools which are stored on katabatic and don't work with the newer versions of Matlab.


