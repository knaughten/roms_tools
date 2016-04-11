from os.path import exists

# Model object containing name of CMIP5 model, institution and ensemble member.
class Model:

    def __init__ (self, inst, name):

        self.name = name
        self.inst = inst
        self.ens = 'r1i1p1'
        if self.name == 'CESM1-WACCM':
            self.ens = 'r2i1p1'


    # Return the directory containing monthly output for the given variable and
    # experiment.
    # Input:
    # expt = string containing name of experiment, eg 'historical'
    # var_name = string containing name of variable, eg 'tas'
    # Output: dir = string containing directory where the model output is
    #               stored. If this model output doesn't exist, return an
    #               empty string.
    def get_directory (self, expt, var_name):

        # Figure out whether it is an atmosphere or ocean variable
        if var_name in ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']:
            realm = 'atmos'
        elif var_name in ['thetao', 'so', 'uo', 'vo', 'zos']:
            realm = 'ocean'
        else:
            print 'Unknown variable'
            # Exit early
            return ''

        # Build typical directory structure in ua6's archive on raijin
        dir = '/g/data1/ua6/drstree/CMIP5/GCM/' + self.inst + '/' + self.name + '/' + expt + '/mon/' + realm + '/' + var_name + '/' + self.ens + '/'

        # Exceptions
        if self.name == 'GFDL-ESM2M' and expt == 'rcp85' and var_name == 'clt':
            dir = '/g/data/ua6/unofficial-ESG-replica/tmp/tree/esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/NOAA-GFDL/GFDL-ESM2M/rcp85/mon/atmos/Amon/r1i1p1/v20110601/clt/'

        # Check if this directory actually exists, or if the data is missing.
        if not exists(dir):
            # Return an empty string
            dir = ''

        return dir


# Build an array of Model objects for the 39 CMIP5 models used in this project.
def build_model_list ():

    Models = []
    Models.append(Model('BCC', 'bcc-csm1-1'))
    Models.append(Model('BCC', 'bcc-csm1-1-m'))
    Models.append(Model('BNU', 'BNU-ESM'))
    Models.append(Model('CCCMA', 'CanESM2'))
    Models.append(Model('CMCC', 'CMCC-CM'))
    Models.append(Model('CMCC', 'CMCC-CMS'))
    Models.append(Model('CNRM', 'CNRM-CM5'))
    Models.append(Model('CSIRO-BOM', 'ACCESS1-0'))
    Models.append(Model('CSIRO-BOM', 'ACCESS1-3'))
    Models.append(Model('CSIRO-QCCCE', 'CSIRO-Mk3-6-0'))
    Models.append(Model('FIO', 'FIO-ESM'))
    Models.append(Model('ICHEC', 'EC-EARTH'))
    Models.append(Model('INM', 'inmcm4'))
    Models.append(Model('IPSL', 'IPSL-CM5A-LR'))
    Models.append(Model('IPSL', 'IPSL-CM5A-MR'))
    Models.append(Model('IPSL', 'IPSL-CM5B-LR'))
    Models.append(Model('LASG-CESS', 'FGOALS-g2'))
    Models.append(Model('MIROC', 'MIROC-ESM'))
    Models.append(Model('MIROC', 'MIROC-ESM-CHEM'))
    Models.append(Model('MIROC', 'MIROC5'))
    Models.append(Model('MOHC', 'HadGEM2-CC'))
    Models.append(Model('MOHC', 'HadGEM2-ES'))
    Models.append(Model('MPI-M', 'MPI-ESM-LR'))
    Models.append(Model('MPI-M', 'MPI-ESM-MR'))
    Models.append(Model('MRI', 'MRI-CGCM3'))
    Models.append(Model('NASA-GISS', 'GISS-E2-H'))
    Models.append(Model('NASA-GISS', 'GISS-E2-H-CC'))
    Models.append(Model('NASA-GISS', 'GISS-E2-R'))
    Models.append(Model('NASA-GISS', 'GISS-E2-R-CC'))
    Models.append(Model('NCAR', 'CCSM4'))
    Models.append(Model('NCC', 'NorESM1-M'))
    Models.append(Model('NCC', 'NorESM1-ME'))
    Models.append(Model('NIMR-KMA', 'HadGEM2-AO'))
    Models.append(Model('NOAA-GFDL', 'GFDL-CM3'))
    Models.append(Model('NOAA-GFDL', 'GFDL-ESM2G'))
    Models.append(Model('NOAA-GFDL', 'GFDL-ESM2M'))
    Models.append(Model('NSF-DOE-NCAR', 'CESM1-BGC'))
    Models.append(Model('NSF-DOE-NCAR', 'CESM1-CAM5'))
    Models.append(Model('NSF-DOE-NCAR', 'CESM1-WACCM'))

    return Models
