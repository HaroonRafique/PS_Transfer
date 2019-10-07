import numpy as np

mismatch_factor=1.0

tomo_file = 'PyORBIT_Tomo_file_MD4224_HB.mat'
input_distn = '../../Initial_Distn/mainbunch.mat'

# PS Injection 1.4 GeV
gamma = 2.49253731343
beta = np.sqrt(gamma**2-1)/gamma
c = 299792458

# Beam Parameters from MD4224
intensity = 72.5E+10
epsn_x = 1E-6
epsn_y = 1E-6

blength = 140e-9
sig_z = (beta * c * blength)/4.
dpp_rms = 9e-04
rf_voltage = 0.0212942055190595723

# Simulation Parameters
n_macroparticles = int(2E3)
turns_max = int(20)
tu1 = range(10, turns_max, 10)
tu3 = range(1, 9)
tu = tu3 + tu1 

turns_print = sorted(tu)
turns_update = sorted(tu)
macrosize = intensity/float(n_macroparticles)

# Space Charge
grid_x = 128
grid_y = 128
grid_z = 64

# PTC RF Table Parameters
harmonic_factors = [1] # this times the base harmonic defines the RF harmonics (for SPS = 4620, PS 10MHz 7, 8, or 9)
time = np.array([0,1,2])
ones = np.ones_like(time)
Ekin_GeV = 1.4*ones
RF_voltage_MV = np.array([0.0212942055190595723*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

# Constants
circumference = 2*np.pi*100
m = 1.2					# Longitudinal Joho Parameter
TransverseCut = 5		# Used for some distributions (matched)

parameters = {
	'intensity': intensity,
	'gamma': gamma,
	'bunch_length': blength,
	'epsn_x': epsn_x,
	'epsn_y': epsn_y,
	'dpp_rms': dpp_rms,
	'tomo_file': tomo_file,
	'LongitudinalJohoParameter': m,
	'LongitudinalCut': 2.4,
	'n_macroparticles': n_macroparticles,
	'TransverseCut': TransverseCut,
	'macrosize': macrosize,
	'turns_max': turns_max,
	'turns_update': turns_update,
	'turns_print': turns_print,
	'rf_voltage': rf_voltage,
	'circumference': circumference,
	'input_distn_dir': input_distn
}

tunespread = {
	'intensity': intensity,
	'gamma': gamma,
	'sig_z': sig_z,
	'epsn_x': epsn_x,
	'epsn_y': epsn_y,
	'dpp_rms': dpp_rms
}

switches = {
	'SliceBySlice': False,
	'CreateDistn': True,
	'Mismatch': False,
	'MismatchFactor': mismatch_factor,
	'LongitudinalKick': False,
	'GridSizeX': grid_x,
	'GridSizeY': grid_y,
	'GridSizeZ': grid_z
}

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
