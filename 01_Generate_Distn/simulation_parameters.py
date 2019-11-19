import numpy as np

parameters = {}

parameters['n_macroparticles']		= int(5E5) #int(5E4)

# Include machine (PS), tunes, lattice start position (BWS65H) for bunch output file label
parameters['tunex']					= '6218'
parameters['tuney']					= '624'
parameters['machine']				= 'PS'
parameters['lattice_start'] 		= 'BSG52'

parameters['bunch_label'] 		= parameters['machine'] + '_Lattice_Tune_' + parameters['tunex'] + '_' + parameters['tuney'] + '_' + parameters['lattice_start']
parameters['flat_file']			= '../00_Lattice_Setup/Optimised_Lattice/PTC-PyORBIT_flat_file.flt'
parameters['tomo_file']			= 'PyORBIT_Tomo_file_MD4224_HB.mat'

# PTC RF Table Parameters
harmonic_factors = [1] # this times the base harmonic defines the RF harmonics (for SPS = 4620, PS 10MHz 7, 8, or 9)
time = np.array([0,1,2])
ones = np.ones_like(time)
Ekin_GeV = 1.4*ones
RF_voltage_MV = np.array([0.0212942055190595723*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
