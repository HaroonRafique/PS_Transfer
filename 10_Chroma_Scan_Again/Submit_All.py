import os

master_dir = os.getcwd()

locations = []

flag_1 = False
flag_2 = False
flag_3 = False
flag_4 = True
flag_5 = True

if flag_1:
	locations.append('/01_01')
	locations.append('/01_02')
	locations.append('/01_03')
	locations.append('/01_04')
	locations.append('/01_05')
	locations.append('/01_06')
	locations.append('/01_07')
	locations.append('/01_08')
	locations.append('/01_09')
	locations.append('/01_10')
	locations.append('/01_11')
	locations.append('/01_12')
	locations.append('/01_13')
	locations.append('/01_14')
	locations.append('/01_15')
if flag_2:
	locations.append('/02_01')
	locations.append('/02_02')
	locations.append('/02_03')
	locations.append('/02_04')
	locations.append('/02_05')
	locations.append('/02_06')
	locations.append('/02_07')
	locations.append('/02_08')
	locations.append('/02_09')
	locations.append('/02_10')
	locations.append('/02_11')
	locations.append('/02_12')
	locations.append('/02_13')
	locations.append('/02_14')
	locations.append('/02_15')
if flag_3:
	locations.append('/03_01')
	locations.append('/03_02')
	locations.append('/03_03')
	locations.append('/03_04')
	locations.append('/03_05')
	locations.append('/03_06')
	locations.append('/03_07')
	locations.append('/03_08')
	locations.append('/03_09')
	locations.append('/03_10')
	locations.append('/03_11')
	locations.append('/03_12')
	locations.append('/03_13')
	locations.append('/03_14')
	locations.append('/03_15')
if flag_4:
	locations.append('/04_01')
	locations.append('/04_02')
	locations.append('/04_03')
	locations.append('/04_04')
	locations.append('/04_05')
	locations.append('/04_06')
	locations.append('/04_07')
	locations.append('/04_08')
	locations.append('/04_09')
	locations.append('/04_10')
	locations.append('/04_11')
	locations.append('/04_12')
	locations.append('/04_13')
	locations.append('/04_14')
	locations.append('/04_15')
if flag_5:
	locations.append('/05_01')
	locations.append('/05_02')
	locations.append('/05_03')
	locations.append('/05_04')
	locations.append('/05_05')
	locations.append('/05_06')
	locations.append('/05_07')
	locations.append('/05_08')
	locations.append('/05_09')
	locations.append('/05_10')
	locations.append('/05_11')
	locations.append('/05_12')
	locations.append('/05_13')
	locations.append('/05_14')
	locations.append('/05_15')

for loc in locations:
	print '---------------------------------------------------------------------------'
	print '\t Submitting SLURM Job: PS Transfer Chroma Scan Sim ', loc[-2:]
	print '---------------------------------------------------------------------------'
	dir_ = master_dir + loc
	make_command = 'python Make_SLURM_submission_script.py'
	submit_command = 'sbatch SLURM_submission_script.sh'
	os.chdir(dir_)
	os.system(make_command)
	os.system(submit_command)
