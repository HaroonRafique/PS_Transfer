#!/bin/bash

if [ -e ./logfile_* ]; then
	mkdir Condor_Logfiles
	mv output/simulation_info_* Condor_Logfiles
	mv logfile_* Condor_Logfiles
	mv output_* Condor_Logfiles
fi

# slurm cleanup
if [ -e ./slurm* ]; then
    mkdir SLURM_Logfiles
	mv slurm* SLURM_Logfiles
fi

# Pickle and output files
rm -r output
rm -r bunch_output
rm -r lost
rm -r input

rm PS.seq
rm ptc_twiss
rm madx_twiss.tfs
rm PTC-PyORBIT_flat_file.flt
rm tunespread.dat

. clean_junk.sh
