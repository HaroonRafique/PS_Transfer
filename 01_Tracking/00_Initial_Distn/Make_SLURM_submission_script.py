#!/usr/bin/env python
# Python script to create a SLURM submission script for PyORBIT
# 21 March 2019 Haroon Rafique CERN BE-ABP-HSI 

import os
  
#-----------------------------------------------------------------------
#	SETTINGS
#-----------------------------------------------------------------------  
script_name = "SLURM_submission_script.sh"

# Switches
hyperthreading = False	# Enable hyperthreading
exclusive = True		# Exclusive (see SLURM documentation)
autotime = True			# 2 days for short queues, 2 weeks for long queues
autotask = True			# Automatically set nodes to maximum tasks
clean_all = True		# Clean simulation folder before running (False when resuming pickle checkpoint)

# Must be chosen

# ~ queue = 'inf-long', 'inf-short', 'batch-long', 'batch-short'
queue = 'inf-short'

n_nodes = 4

jobname = 'PS_00'

path_to_simulation = os.path.dirname(os.path.realpath(__file__)) # This directory

# Optional - have to use with correct switches
manual_time = '504:00:00' # manually set using format 'hours:minutes:seconds'
manual_tasks = 40	# manually change ntasks

# Defaults - can be changed
output_file_name = 'slurm.%N.%j.out'
error_file_name = 'slurm.%N.%j.err'
root_dir = '/hpcscratch/user/harafiqu'
simulation_file = 'pyOrbit.py'
#-----------------------------------------------------------------------
#	AUTOMATICALLY FORMAT SCRIPT
#-----------------------------------------------------------------------  

n_tasks = 0	
if autotask:
	if hyperthreading:
		if 'batch' in queue: n_tasks = 32
		elif 'inf' in queue: n_tasks = 40
		else: 
			print 'queue not recognised'
			exit(0)

	else:
		if 'batch' in queue: n_tasks = 16
		elif 'inf' in queue: n_tasks = 20
		else: 
			print 'queue not recognised'
			exit(0)
else: n_tasks = manual_tasks

time = '48:00:00'
if autotime:
	if queue == 'batch-short': time = '48:00:00'
	elif queue == 'inf-short': time = '120:00:00'
	elif queue == 'inf-long' or 'batch-long': time = '504:00:00'
	else: 
		print 'queue not recognised'
		exit(0)
else: time = manual_time

#-----------------------------------------------------------------------
#	WRITE FILE
#-----------------------------------------------------------------------  
if os.path.exists(script_name):  
	print 'SLURM submission script ' + script_name + ' already exists. Deleting'
	os.remove(script_name)

print "Creating ", script_name

f= open(script_name,"w")

f.write('#!/bin/bash')
f.write('\n#SBATCH --job-name=' + str(jobname))
f.write('\n#SBATCH --output=' + str(output_file_name))
f.write('\n#SBATCH --error=' + str(error_file_name))
f.write('\n#SBATCH --nodes=' + str(n_nodes))
f.write('\n#SBATCH --ntasks-per-node=' + str(n_tasks))
f.write('\n#SBATCH --partition=' + str(queue))
f.write('\n#SBATCH --time=' + str(time))
f.write('\n#SBATCH --mem-per-cpu=3200M')
if (exclusive): f.write('\n#SBATCH --exclusive')
if not hyperthreading: f.write('\n#SBATCH --hint=nomultithread')
f.write('\n')
f.write('\nBATCH_ROOT_DIR=' + str(root_dir))
f.write('\nRUN_DIR=' + str(path_to_simulation))
f.write('\nOrigIwd=$(pwd)')
f.write('\n')
f.write('\n# Make an output folder in the root directory to hold SLURM info file')
f.write('\ncd ${BATCH_ROOT_DIR}')
f.write('\noutput_dir="output"')
f.write('\nmkdir -p $output_dir')
f.write('\n')
f.write('\n# Fill the SLURM info file')
f.write('\nsimulation_info_file="${BATCH_ROOT_DIR}/${output_dir}/simulation_info_${SLURM_JOB_ID}.${SLURM_NODEID}.${SLURM_PROCID}.txt"')
f.write('\necho "PyOrbit path:  `readlink -f ${ORBIT_ROOT}`" >> ${simulation_info_file}')
f.write('\necho "Run path:  `readlink -f ${RUN_DIR}`" >> ${simulation_info_file}')
f.write('\necho "Submit host:  `readlink -f ${SLURM_SUBMIT_HOST}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Job name:  `readlink -f ${SLURM_JOB_NAME}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Job ID:  `readlink -f ${SLURM_JOB_ID}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Nodes allocated:  `readlink -f ${SLURM_JOB_NUM_NODES}`" >> ${simulation_info_file}')
f.write('\necho "SLURM CPUS per Node:  `readlink -f ${SLURM_CPUS_ON_NODE}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Node ID:  `readlink -f ${SLURM_NODEID}`" >> ${simulation_info_file}')
f.write('\necho "SLURM total cores for job:  `readlink -f ${SLURM_NTASKS}`" >> ${simulation_info_file}')
f.write('\necho "SLURM process ID:  `readlink -f ${SLURM_PROCID}`" >> ${simulation_info_file}')
f.write('\necho "****************************************" >> ${simulation_info_file}')
f.write('\n')
f.write('\n# Enter job directory, clean it, and setup environment -> SLURM info file')
f.write('\ncd ${RUN_DIR}')
if clean_all:f.write('\n./clean_all.sh')
f.write('\n. setup_environment.sh >> ${simulation_info_file}')
f.write('\n')
f.write('\n# Load correct MPI')
f.write('\nmodule load mpi/mvapich2/2.2')
f.write('\n')
f.write('\ntstart=$(date +%s)')
f.write('\n')
f.write('\n# Run the job')
if hyperthreading:f.write('\nsrun ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))	
else:f.write('\nsrun --hint=nomultithread ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))
f.write('\n')
f.write('\ntend=$(date +%s)')
f.write('\ndt=$(($tend - $tstart))')
f.write('\necho "total simulation time (s): " $dt >> ${simulation_info_file}')

f.close()

print 'SLURM submission script creation finished'
