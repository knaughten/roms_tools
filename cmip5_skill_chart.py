from numpy import *

# Global variables
# List of CMIP5 variable names
var_names = ['tas', 'uas', 'vas', 'pr', 'prsn', 'huss', 'clt', 'ps', 'evspsbl', 'rsds', 'rlds', 'thetao', 'so', 'uo', 'vo']
# List of all possible season codes
season_names = ['annual', 'djf', 'mam', 'jja', 'son']


# ModelScore object containing model name and rankings for each 
# variable and season.
class ModelScore:

    def __init__ (self, name):
        self.name = name
        self.scores = zeros([len(var_names), len(season_names)])


    # Add a ranking for a specific variable and season.
    def add_score (self, var_name, season, score):
        var_index = var_names.index(var_name)
        season_index = season_names.index(season)
        self.scores[var_index, season_index] = score
    

# Summarise the results of cmip5_eraint_skill.py and cmip5_ecco2_skill.py
# by listing the model rankings in a giant table of models x variables.
# Each cell in the table will have 5 entries, representing the 5 season codes
# (annual, DJF, MAM, JJA, SON).
def cmip5_skill_chart ():    

    # Create ModelScore objects for each model
    model_scores = []
    model_scores.append(ModelScore('bcc-csm1-1'))
    model_scores.append(ModelScore('bcc-csm1-1-m'))
    model_scores.append(ModelScore('BNU-ESM'))
    model_scores.append(ModelScore('CanESM2'))
    model_scores.append(ModelScore('CMCC-CM'))
    model_scores.append(ModelScore('CMCC-CMS'))
    model_scores.append(ModelScore('CNRM-CM5'))
    model_scores.append(ModelScore('ACCESS1-0'))
    model_scores.append(ModelScore('ACCESS1-3'))
    model_scores.append(ModelScore('CSIRO-Mk3-6-0'))
    model_scores.append(ModelScore('FIO-ESM'))
    model_scores.append(ModelScore('EC-EARTH'))
    model_scores.append(ModelScore('inmcm4'))
    model_scores.append(ModelScore('IPSL-CM5A-LR'))
    model_scores.append(ModelScore('IPSL-CM5A-MR'))
    model_scores.append(ModelScore('IPSL-CM5B-LR'))
    model_scores.append(ModelScore('FGOALS-g2'))
    model_scores.append(ModelScore('MIROC-ESM'))
    model_scores.append(ModelScore('MIROC-ESM-CHEM'))
    model_scores.append(ModelScore('MIROC5'))
    model_scores.append(ModelScore('HadGEM2-CC'))
    model_scores.append(ModelScore('HadGEM2-ES'))
    model_scores.append(ModelScore('MPI-ESM-LR'))
    model_scores.append(ModelScore('MPI-ESM-MR'))
    model_scores.append(ModelScore('MRI-CGCM3'))
    model_scores.append(ModelScore('GISS-E2-H'))
    model_scores.append(ModelScore('GISS-E2-H-CC'))
    model_scores.append(ModelScore('GISS-E2-R'))
    model_scores.append(ModelScore('GISS-E2-R-CC'))
    model_scores.append(ModelScore('CCSM4'))
    model_scores.append(ModelScore('NorESM1-M'))
    model_scores.append(ModelScore('NorESM1-ME'))
    model_scores.append(ModelScore('HadGEM2-AO'))
    model_scores.append(ModelScore('GFDL-CM3'))
    model_scores.append(ModelScore('GFDL-ESM2G'))
    model_scores.append(ModelScore('GFDL-ESM2M'))
    model_scores.append(ModelScore('CESM1-BGC'))
    model_scores.append(ModelScore('CESM1-CAM5'))
    model_scores.append(ModelScore('CESM1-WACCM'))

    # Create a list of the corresponding model names (we need this because
    # it is searchable, the ModelScore array is not)
    model_names = []
    for i in range(len(model_scores)):
        model_names.append(model_scores[i].name)

    # Loop over variables and seasons
    for var in var_names:
        for season in season_names:
            # Read the rankings from skill score file
            file_name = 'cmip5_skill_scores/' + var + '_' + season + '.txt'
            f = open(file_name, 'r')
            # Models are listed in order of ranking, from best to worst
            rank = 1
            for line in f:
                # Parse the model name and find its index in the arrays
                name = line.split()[0]
                model_index = model_names.index(name)
                # Add this ranking to the ModelScore object
                model_scores[model_index].add_score(var, season, rank)
                rank += 1
            f.close()

    # Write results to a text file, tabulated nicely
    f = open('cmip5_skill_scores/summary.txt', 'w')
    f.write('{0:16}'.format('Model'))
    for var in var_names:
        f.write('{0:16}'.format(var))
    f.write('\n')
    for model in model_scores:
        f.write('{0:16}'.format(model.name))
        for i in range(size(var_names)):
            curr_scores = ''
            for j in range(size(season_names)):
                if model.scores[i,j] > 0:
                    curr_scores += str(int(model.scores[i,j]))
                    if j+1 < size(season_names):
                        curr_scores += ','
            f.write('{0:16}'.format(curr_scores))
        f.write('\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    cmip5_skill_chart()
    
        
    
                        
    
