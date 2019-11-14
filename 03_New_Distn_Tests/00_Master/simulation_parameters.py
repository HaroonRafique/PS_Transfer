import numpy as np

parameters = {}

parameters['n_macroparticles']		= int(0.5E3) #int(50E3) # int(0.5E6)

# Include machine (PS), tunes, lattice start position (BWS65H) for bunch output file label
parameters['tunex']					= '6218'
parameters['tuney']					= '624'
parameters['machine']				= 'PS'
parameters['lattice_start'] 		= 'BSG52'

parameters['bunch_label'] 		= parameters['machine'] + '_Lattice_Tune_' + parameters['tunex'] + '_' + parameters['tuney'] + '_' + parameters['lattice_start']
parameters['flat_file']			= '../../00_Lattice_Setup/Optimised_Lattice/PTC-PyORBIT_flat_file.flt'
parameters['tomo_file']			= 'PyORBIT_Tomo_file_MD4224_HB.mat'
parameters['bunch_file']		= '../../01_Generate_Distn/Bunches/PyORBIT_Tomo_Bunch_Lattice_Twiss_Nmp_500_PS_Lattice_Tune_6218_624_BSG52.mat'
parameters['intensity']			= 72.5E+10
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

parameters['n_macroparticles']	= int(41)
parameters['gamma']				= 2.49253731343
parameters['intensity']			= 72.5E+10
parameters['bunch_length']		= 140e-9
parameters['blength']			= 140e-9
parameters['epsn_x']			= 1E-6
parameters['epsn_y']			= 1.2E-6
parameters['dpp_rms']			= 8.7e-04
parameters['LongitudinalJohoParameter'] = 1.2
parameters['LongitudinalCut'] 	= 2.4
parameters['TransverseCut']		= 5
parameters['rf_voltage']		= 0.0212942055190595723
parameters['circumference']		= 2*np.pi*100
parameters['phi_s']				= 0
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

# PS Injection 1.4 GeV
parameters['gamma'] 	= 2.49253731343
parameters['beta'] 		= np.sqrt(parameters['gamma']**2-1)/parameters['gamma']
c 						= 299792458
parameters['sig_z'] 	= (parameters['beta'] * c * parameters['blength'])/4.

parameters['turns_max'] = int(30)
parameters['turns_print'] = range(1, 101)
parameters['turns_update'] = range(1, 101)

switches = {
	'Optics':		'Matched', #'Rematched', #'Lattice',
	'CreateDistn':		True,
	'Space_Charge': 	False,
	'GridSizeX': 128,
	'GridSizeY': 128,
	'GridSizeZ': 64
}

# PTC RF Table Parameters
harmonic_factors = [1] # this times the base harmonic defines the RF harmonics (for SPS = 4620, PS 10MHz 7, 8, or 9)
time = np.array([0,1,2])
ones = np.ones_like(time)
Ekin_GeV = 1.4*ones
# ~ RF_voltage_MV = np.array([0.0212942055190595723*ones]).T # in MV
RF_voltage_MV = np.array([0.0*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
