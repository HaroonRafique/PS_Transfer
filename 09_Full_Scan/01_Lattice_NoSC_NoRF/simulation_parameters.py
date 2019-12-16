import numpy as np

parameters = {}

parameters['n_macroparticles']		= int(5E5) # int(5E4)

# Include machine (PS), tunes, lattice start position (BWS65H) for bunch output file label
parameters['tunex']					= '6218'
parameters['tuney']					= '624'
parameters['machine']				= 'PS'
parameters['lattice_start'] 		= 'BSG52'
parameters['Optics'] 		= 'Lattice' #'ReM', #'Op',

parameters['bunch_label'] 		= parameters['machine'] + '_Lattice_Tune_' + parameters['tunex'] + '_' + parameters['tuney'] + '_' + parameters['lattice_start']
parameters['flat_file']			= '../../00_Lattice_Setup/Optimised_Lattice/PTC-PyORBIT_flat_file.flt'
parameters['tomo_file']			= 'PyORBIT_Tomo_file_MD4224_HB.mat'
# ~ parameters['bunch_file']		= '../../01_Generate_Distn/Bunches/PyORBIT_Tomo_Bunch_Lattice_Twiss_Nmp_500_PS_Lattice_Tune_6218_624_BSG52.mat'
parameters['bunch_file']		= '../../01_Generate_Distn/Bunches/PyORBIT_Tomo_Bunch_Manual_Twiss_Nmp_' + str(parameters['n_macroparticles'])+'_PS_Lattice_Tune_6218_624_' + parameters['lattice_start']+'_'+parameters['Optics']+'.mat'
parameters['intensity']			= 65E+10
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

parameters['gamma']				= 2.49253731343
parameters['bunch_length']		= 140e-9
parameters['blength']			= 140e-9
parameters['epsn_x']			= 1E-6
parameters['epsn_y']			= 1E-6
parameters['dpp_rms']			= 9e-04
parameters['LongitudinalJohoParameter'] = 1.2
parameters['LongitudinalCut'] 	= 2.4
parameters['TransverseCut']		= 5
# ~ parameters['rf_voltage']		= 0.0212942055190595723
parameters['rf_voltage']		= 0.0
parameters['circumference']		= 2*np.pi*100
parameters['phi_s']				= 0
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

# PS Injection 1.4 GeV
parameters['gamma'] 	= 2.49253731343
parameters['beta'] 		= np.sqrt(parameters['gamma']**2-1)/parameters['gamma']
c 						= 299792458
parameters['sig_z'] 	= (parameters['beta'] * c * parameters['blength'])/4.

parameters['turns_max'] = int(2200)

tu1 = range(-1, parameters['turns_max'], 200)
tu2 = range(50, 100, 10) 
tu3 = range(1, 50)
tu = tu2 + tu1 + tu3 

parameters['turns_print'] = sorted(tu)
parameters['turns_update'] = sorted(tu)

switches = {
	'CreateDistn':		False,	# True will create dispersion vector
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
RF_voltage_MV = np.array([parameters['rf_voltage']*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
