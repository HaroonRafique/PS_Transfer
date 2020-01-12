import os

master_dir = os.getcwd()

locations = []

locations.append('/01_Lattice_NoSC_NoRF')
locations.append('/01_Op_NoSC_NoRF')
locations.append('/01_ReM_NoSC_NoRF')
locations.append('/02_Lattice_SC_RF')
locations.append('/02_Op_SC_RF')
locations.append('/02_ReM_SC_RF')

for loc in locations:
	print '--------------------------------------------------------------------------------------------'
	print '\t Submitting HPC-Batch simulation: PS Transfer Injection Oscillations'
	print '--------------------------------------------------------------------------------------------'
	dir_ = master_dir + loc
	make_command = 'python Make_SLURM_submission_script.py'
	submit_command = 'sbatch SLURM_submission_script.sh'
	os.chdir(dir_)
	os.system(make_command)
	os.system(submit_command)
