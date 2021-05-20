# --- SCRIPT TO GENERATE INDIVIDUALIZED SURVIVAL CURVES --- #

dir_base = getwd()
# !!! SET FOLDER/FILE AS APPROPRIATE !!! #
dir_data = file.path(dir_base, 'example_run')
fn_data = file.path(dir_data, 'example_X.csv')
stopifnot(file.exists(fn_data))

source('funs_support.R')

# --- (1) LOAD IN THE DATA --- #

# (i) Load in the model coefficients


