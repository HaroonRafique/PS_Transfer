import os

master_dir = os.getcwd()

locations = []

locations.append('/01_Lattice_NoSC_NoRF')
locations.append('/02_Op_NoSC_NoRF')
locations.append('/03_ReM_NoSC_NoRF')
locations.append('/04_Lattice_SC_NoRF')
locations.append('/05_Op_SC_NoRF')
locations.append('/06_ReM_SC_NoRF')
locations.append('/07_Lattice_NoSC_RF')
locations.append('/08_Op_NoSC_RF')
locations.append('/09_ReM_NoSC_RF')
locations.append('/10_Lattice_SC_RF')
locations.append('/11_Op_SC_RF')
locations.append('/12_ReM_SC_RF')

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
