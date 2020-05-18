import shutil

pyorbit = True
simulation_parameters = False
flat_files = False
tune_files = False
distn_gen = False

master_directory = './00_Master'
pyorbit_file = master_directory + '/pyOrbit.py'
sim_params_file = master_directory + '/simulation_parameters.py'
flat_file = master_directory + '/Flat_file.madx'
tune_file = master_directory + '/tunes.str'
distn_generator = master_directory + '/lib/pyOrbit_GenerateInitialDistribution.py'

sbs_locations = []
noSC_locations = []

sbs_locations.append('./01_01/')
sbs_locations.append('./01_02/')
sbs_locations.append('./01_03/')
sbs_locations.append('./01_04/')
sbs_locations.append('./01_05/')
sbs_locations.append('./01_06/')
sbs_locations.append('./01_07/')
sbs_locations.append('./01_08/')
sbs_locations.append('./01_09/')
sbs_locations.append('./01_10/')
sbs_locations.append('./01_11/')
sbs_locations.append('./01_12/')
sbs_locations.append('./01_13/')
sbs_locations.append('./01_14/')
sbs_locations.append('./01_15/')

sbs_locations.append('./02_01/')
sbs_locations.append('./02_02/')
sbs_locations.append('./02_03/')
sbs_locations.append('./02_04/')
sbs_locations.append('./02_05/')
sbs_locations.append('./02_06/')
sbs_locations.append('./02_07/')
sbs_locations.append('./02_08/')
sbs_locations.append('./02_09/')
sbs_locations.append('./02_10/')
sbs_locations.append('./02_11/')
sbs_locations.append('./02_12/')
sbs_locations.append('./02_13/')
sbs_locations.append('./02_14/')
sbs_locations.append('./02_15/')

sbs_locations.append('./03_01/')
sbs_locations.append('./03_02/')
sbs_locations.append('./03_03/')
sbs_locations.append('./03_04/')
sbs_locations.append('./03_05/')
sbs_locations.append('./03_06/')
sbs_locations.append('./03_07/')
sbs_locations.append('./03_08/')
sbs_locations.append('./03_09/')
sbs_locations.append('./03_10/')
sbs_locations.append('./03_11/')
sbs_locations.append('./03_12/')
sbs_locations.append('./03_13/')
sbs_locations.append('./03_14/')
sbs_locations.append('./03_15/')

if pyorbit:
	for loc in sbs_locations:
		newPath = shutil.copy(pyorbit_file, loc)
		print pyorbit_file, ' copied to ', loc

if simulation_parameters:
	for loc in sbs_locations:
		newPath = shutil.copy(sim_params_file, loc)
		print sim_params_file, ' copied to ', loc

if flat_files:
	for loc in sbs_locations:
		newPath = shutil.copy(flat_file, loc)
		print flat_file, ' copied to ', loc

if tune_files:
	for loc in sbs_locations:
		newPath = shutil.copy(tune_file, loc)
		print flat_file, ' copied to ', loc

if distn_gen:
	for loc in sbs_locations:
		loc_ = loc + 'lib/'
		newPath = shutil.copy(distn_generator, loc_)
		print distn_generator, ' copied to ', loc_
