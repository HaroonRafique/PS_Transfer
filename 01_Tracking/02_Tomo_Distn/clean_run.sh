#!/bin/bash

#~ if [ -e ./logfile_* ]; then
	#~ mkdir Condor_Logfiles
	#~ mv output/simulation_info_* Condor_Logfiles
	#~ mv logfile_* Condor_Logfiles
	#~ mv output_* Condor_Logfiles
#~ fi

#~ # slurm cleanup
#~ if [ -e ./slurm.* ]; then
    #~ mkdir Slurm_Logfiles
	#~ mv slurm* Slurm_Logfiles
#~ fi

# Pickle and output files
rm -r output
rm -r bunch_output
rm -r lost
rm -r input

rm PTC-PyORBIT_flat_file.flt
rm tunespread.dat

. clean_junk.sh
