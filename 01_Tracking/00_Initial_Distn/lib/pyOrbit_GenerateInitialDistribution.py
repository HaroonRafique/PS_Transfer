import math
import sys
from itertools import chain
import numpy as np
import csv
import random
import orbit_mpi
import scipy.io as sio

from bunch import Bunch
from orbit.injection.joho import JohoLongitudinal
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist2D, GaussDist2D, KVDist2D
from orbit.utils.consts import mass_proton, speed_of_light, pi
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.bunch_utils import ParticleIdNumber


class LongitudinalDistributionFromTomoscope():

	def __init__(self, filename, matfile=0):
		if matfile:
			data = sio.loadmat(filename, squeeze_me=True)
			self.data = data
			self.I = data['density_array']
			self.t, self.dE = np.meshgrid(data['time_nsec'],data['energy_MeV'])		
			# ~ print self.t
			# ~ print self.dE.shape
			# ~ print self.I.shape
			
		else:
			thefile = open(filename,"r")
			lines=0
			for l in thefile:
				lines = lines + 1
				
			lines = lines - 2
			print 'grid size = ', lines, ' * ', lines, ' = ', (lines*lines)
				
			thefile.close()
			thefile = open(filename,"r")			
			
			# Need to know the size of the square numpy array from the input file
			# ~ self.I = np.zeros(shape=(lines,lines))
			self.dt_bin =  np.zeros(shape=(lines,lines))
			self.dE_bin =  np.zeros(shape=(lines,lines))
			self.dt_no_bins = 0
			self.dE_no_bins = 0
			
			for l in thefile:
				# First line: Minimum dt, maximum dt, binsize, bins
				if l == 0:
					self.dt_min = l.split()[0]
					self.dt_max = l.split()[1]
					self.dt_bin = l.split()[2]
					self.dt_no_bins = l.split()[3]
				# Second line: Minimum dE, maximum dE, binsize, bins
				elif l == 1:
					self.dE_min = l.split()[0]
					self.dE_max = l.split()[1]
					self.dE_bin = l.split()[2]
					self.dE_no_bins = l.split()[3]
				# Read grid line by line to create density array
				else:
					self.temp=[]
					# read each value and append to array
					for i in l.split():
						self.temp.append(float(i))
					# append array to density array
					self.I.append(self.temp)
					self.temp=[]
			print self.I
			
			self.dt_bins = np.arange(self.dt_no_bins)*self.dt_bin
			self.dE_bins = np.arange(self.dE_no_bins)*self.dE_bin
					
			# make a meshgrid using the min max and binsize etc
			self.t, self.dE = np.meshgrid( self.dt_bins, self.dE_bins)
			print self.t
			print self.t.keys()
			
		self.t_rand = []
		self.dE_rand = []			
		
	def getCoordinates(self, n_mp=1, noise_level=0):
		# ~ data = self.data
		I = self.I
		t = self.t
		dE = self.dE

		U_ = []
		V_ = []
		W_ = []
		t_min = np.min(t)
		t_max = np.max(t)
		dE_min = np.min(dE)
		dE_max = np.max(dE)
		I_max = np.max(I)
		while len(U_)<n_mp:
			u = np.random.uniform(t_min,t_max,n_mp)
			v = np.random.uniform(dE_min,dE_max,n_mp)
			w = np.random.uniform(0,I_max,n_mp)
			# unique creates an array of unique values in an array (for 
			# example dE). Digitize bins an array (v) according to bins
			# (np.unique(dE)) and returns an output array of indices the
			#  same size as (v).Basically this tells us which 'bin' v is
			# in, we minus 1 to give the start of the bin rather than 
			# the end. This points to the density in the 2D array 
			# generated by the tomo.
			d = I[np.digitize(v, np.unique(dE))-1, np.digitize(u, np.unique(t))-1]
			mask = np.where(w < d-noise_level*I_max)[0]
			U_.extend(u[mask])
			V_.extend(v[mask])
			W_.extend(w[mask])
			# print len(U_)
		t_rand = np.array(U_[:n_mp])
		dE_rand = np.array(V_[:n_mp])
		self.t_rand = t_rand
		self.dE_rand = dE_rand
		return t_rand, dE_rand

	def plot_Tomoscope_data(self):
		f,ax = plt.subplots(1)
		plt.scatter(self.t, self.dE, 15, self.I, lw=0)
		ax.set_xlabel('time [ns]')
		ax.set_ylabel('dE [MeV]')

	def plot_generated_distribution(self):
		f,ax = plt.subplots(1)
		plt.hist2d(self.t_rand, self.dE_rand, len(np.unique(self.t)))
		ax.set_xlabel('time [ns]')
		ax.set_ylabel('dE [MeV]')

class JohoNormalized2D:

	def __init__(self, order):
		"""
		Constructor.
		"""
		self.order = order
		self.orderinv = 1. / order
		self.emitlim = (1. + order) / 0.5
		self.__initialize()
		
	def __initialize(self):
		self.emit = self.emitlim * 2./(1. + self.order)
		self.poslength = math.sqrt(self.emitlim)
		self.momlength = math.sqrt(self.emitlim)
		self.emitrms = 0.5 * self.emitlim/(1. + self.order)

	def getCoordinates(self):
		s1 = random.random()
		s2 = random.random()
		a = math.sqrt(1 - pow(s1, self.orderinv))
		al = 2. * math.pi * s2
		u = a * math.cos(al)
		v = a * math.sin(al)
		pos = self.poslength * u
		mom = self.momlength * v		
		return (pos,mom)

class LongitudinalJohoDistributionSingleHarmonic():

	def __init__(self, Parameters_dict, Joho_order):
		"""
		Constructor.
		"""	
		self.beta = Parameters_dict['beta']
		circumference = Parameters_dict['circumference']
		self.harmonic_number = Parameters_dict['harmonic_number']
		gamma_transition = Parameters_dict['gamma_transition']
		gamma = Parameters_dict['gamma']
		bunch_length = Parameters_dict['bunch_length']
		self.phi_s = Parameters_dict['phi_s']
		self.rf_voltage = Parameters_dict['rf_voltage']
		self.LongitudinalCut = Parameters_dict['LongitudinalCut']
	
		self.rf_frequency = speed_of_light * self.beta / circumference * self.harmonic_number
		self.slip_factor = (1 / gamma_transition**2 - 1 / gamma**2)	
		self.energy = Parameters_dict['energy']
		self.sigma_dphi = self.rf_frequency * (bunch_length/2) * pi
		self.sigma_dE = np.sqrt(self.beta**2 * self.energy / self.slip_factor / self.harmonic_number * self.rf_voltage / pi * (np.cos(self.sigma_dphi + self.phi_s) - np.cos(self.phi_s) + self.sigma_dphi * np.sin(self.phi_s)))
		
		self.H_max = (self.rf_voltage/(2*pi) * (np.cos(self.LongitudinalCut * self.sigma_dphi + self.phi_s) - np.cos(self.phi_s) + (self.LongitudinalCut * self.sigma_dphi)*np.sin(self.phi_s)))
		self.distribution = JohoNormalized2D(Joho_order)
		
	def is_inside_limiting_countour(self, phi, dE):
		# ~ print '\n### Generate Initial Dsitribution 2: Entered LongitudinalJohoDistributionSingleHarmonic::is_inside_limiting_contour'
		H = abs(self.harmonic_number * self.slip_factor / (2*self.beta**2 * self.energy) * dE**2 + self.rf_voltage/(2*pi)*(np.cos(phi) - np.cos(self.phi_s) + (phi-self.phi_s)*np.sin(self.phi_s))) 
		# ~ H = abs(self.harmonic_number * self.slip_factor * dE**2 / 2  + self.rf_voltage/(2*pi * self.beta**2 * self.energy) * (np.cos(phi) - np.cos(self.phi_s) + (phi-self.phi_s)*np.sin(self.phi_s))) 
		return H <= abs(self.H_max)

	def getCoordinates(self):
		while True:
			phi, dE = self.distribution.getCoordinates()
			phi *= self.sigma_dphi
			dE *= self.sigma_dE
			# ~ print '\n### Generate Initial Dsitribution 2: Calling LongitudinalJohoDistributionSingleHarmonic::is_inside_limiting_contour'
			if self.is_inside_limiting_countour(phi, dE): 
				# ~ print '\n### Generate Initial Dsitribution 2: Returned LongitudinalJohoDistributionSingleHarmonic::is_inside_limiting_contour'
				break
		return phi, dE
	
def generate_initial_distribution_from_tomo_manual_Twiss(parameters, TwissDict, matfile=0, output_file='ParticleDistribution.in', outputFormat='pyOrbit', summary_file='ParticleDistribution_summary.txt', summary_mat_file=None):

	# Get parameters from the TwissDict dictionary
	parameters['alphax0'] = TwissDict['alpha_x']
	parameters['betax0']  = TwissDict['beta_x']
	parameters['alphay0'] = TwissDict['alpha_y']
	parameters['betay0']  = TwissDict['beta_y']
	parameters['etax0']   = TwissDict['D_x']
	parameters['etapx0']  = TwissDict['D_xp']
	parameters['etay0']   = TwissDict['D_y']
	parameters['etapy0']  = TwissDict['D_yp']
	parameters['x0']      = TwissDict['x0']
	parameters['xp0']     = TwissDict['xp0']
	parameters['y0']      = TwissDict['y0']
	parameters['yp0']     = TwissDict['yp0']
	parameters['gamma_transition'] = TwissDict['gamma_transition']
	parameters['circumference']    = TwissDict['circumference']
	parameters['length'] = TwissDict['length']
	
	# Create Twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 

	print dispersionx

	# Initialize empty particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	z = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# Instatiate the classes for longitudinal and transverse distns
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Longitudinal_distribution = LongitudinalDistributionFromTomoscope(parameters['tomo_file'], matfile)

	try: 
		noise_level = parameters['noise_level']
	except KeyError:
		noise_level = 0	
	t_rand, dE_rand = Longitudinal_distribution.getCoordinates(parameters['n_macroparticles'], noise_level) 
	z = (t_rand * 1e-9) * speed_of_light * parameters['beta'] # convert ns to s and then m
	dE = dE_rand * 1e-3 # convert from MeV to GeV
	
	# We need to convert z into phi
	h_main = np.atleast_1d(parameters['harmonic_number'])[0]
	R = parameters['circumference'] / 2 / np.pi
	phi = - z * h_main / R

	# Write the distn to a file only on one CPU
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:
		
		with open(output_file,"w") as fid:
			
			csv_writer = csv.writer(fid, delimiter=' ')
			for i in range(parameters['n_macroparticles']):
				
				# ~ (z[i], dE[i]) = Longitudinal_distribution.getCoordinates()
				# ~ z[i] = z[i] * speed_of_light * parameters['beta'] * 1e-9 # convert ns to s and then m
				# ~ dE[i] = dE[i] * 1e-3 # convert from MeV to GeV
				(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
				x[i] += closedOrbitx['x0']
				xp[i] += closedOrbitx['xp0']
				y[i] += closedOrbity['y0']
				yp[i] += closedOrbity['yp0']
				dpp = dE[i] / (parameters['energy']) / parameters['beta']**2 * 1E9
				#print '\n dpp = ', dpp
				x[i] += dpp * dispersionx['etax0']
				xp[i] += dpp * dispersionx['etapx0']
				y[i] += dpp * dispersiony['etay0']
				yp[i] += dpp * dispersiony['etapy0']
				
				# ~ if outputFormat == 'Orbit':
				x[i] *= 1000.
				xp[i] *= 1000.
				y[i] *= 1000.
				yp[i] *= 1000.
				# ~ dE[i] /= 1.e9	
						
			# ~ if outputFormat == 'Orbit':
			map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]]), range(parameters['n_macroparticles']))	
			# ~ elif outputFormat == 'pyOrbit':
				# ~ map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]]), range(parameters['n_macroparticles']))	
				
		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(parameters[key]) + '\n'), parameters)
				
		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file

	orbit_mpi.MPI_Barrier(comm)

	return output_file
	
def generate_initial_distribution_from_tomo(parameters, matfile=0, Lattice=None, output_file='ParticleDistribution.in', outputFormat='pyOrbit', summary_file='ParticleDistribution_summary.txt', summary_mat_file=None):
	
	# Get parameters from the lattice
	parameters['alphax0'] = Lattice.alphax0
	parameters['betax0']  = Lattice.betax0
	parameters['alphay0'] = Lattice.alphay0
	parameters['betay0']  = Lattice.betay0
	parameters['etax0']   = Lattice.etax0
	parameters['etapx0']  = Lattice.etapx0
	parameters['etay0']   = Lattice.etay0
	parameters['etapy0']  = Lattice.etapy0
	parameters['x0']      = Lattice.orbitx0
	parameters['xp0']     = Lattice.orbitpx0
	parameters['y0']      = Lattice.orbity0
	parameters['yp0']     = Lattice.orbitpy0
	parameters['gamma_transition'] = Lattice.gammaT
	parameters['circumference']    = Lattice.getLength()
	parameters['length'] = Lattice.getLength()/Lattice.nHarm
	
	# Create Twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 
	
	print dispersionx
	
	# Initialize empty particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	z = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# Instatiate the classes for longitudinal and transverse distns
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Longitudinal_distribution = LongitudinalDistributionFromTomoscope(parameters['tomo_file'], matfile)

	try: 
		noise_level = parameters['noise_level']
	except KeyError:
		noise_level = 0	
	t_rand, dE_rand = Longitudinal_distribution.getCoordinates(parameters['n_macroparticles'], noise_level) 
	z = (t_rand * 1e-9) * speed_of_light * parameters['beta'] # convert ns to s and then m
	dE = dE_rand * 1e-3 # convert from MeV to GeV
	
	# We need to convert z into phi
	h_main = np.atleast_1d(parameters['harmonic_number'])[0]
	R = parameters['circumference'] / 2 / np.pi
	phi = - z * h_main / R

	# Write the distn to a file only on one CPU
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:

		with open(output_file,"w") as fid:

			csv_writer = csv.writer(fid, delimiter=' ')
			for i in range(parameters['n_macroparticles']):

				# ~ (z[i], dE[i]) = Longitudinal_distribution.getCoordinates()
				# ~ z[i] = z[i] * speed_of_light * parameters['beta'] * 1e-9 # convert ns to s and then m
				# ~ dE[i] = dE[i] * 1e-3 # convert from MeV to GeV
				(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
				x[i] += closedOrbitx['x0']
				xp[i] += closedOrbitx['xp0']
				y[i] += closedOrbity['y0']
				yp[i] += closedOrbity['yp0']
				dpp = dE[i] / (parameters['energy']) / parameters['beta']**2 * 1E9
				#print '\n dpp = ', dpp
				x[i] += dpp * dispersionx['etax0']
				xp[i] += dpp * dispersionx['etapx0']
				y[i] += dpp * dispersiony['etay0']
				yp[i] += dpp * dispersiony['etapy0']

				# ~ if outputFormat == 'Orbit':
				x[i] *= 1000.
				xp[i] *= 1000.
				y[i] *= 1000.
				yp[i] *= 1000.
				# ~ dE[i] /= 1.e9	
						
			# ~ if outputFormat == 'Orbit':
			map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]]), range(parameters['n_macroparticles']))	
			# ~ elif outputFormat == 'pyOrbit':
				# ~ map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]]), range(parameters['n_macroparticles']))	

		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(parameters[key]) + '\n'), parameters)

		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file

	orbit_mpi.MPI_Barrier(comm)

	return output_file

def generate_initial_distribution_dispersion_mismatch(parameters, Lattice, mismatch = 0.1, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
	parameters['alphax0'] = Lattice.alphax0
	parameters['betax0']  = Lattice.betax0
	parameters['alphay0'] = Lattice.alphay0
	parameters['betay0']  = Lattice.betay0
	parameters['etax0']   = Lattice.etax0 * mismatch
	parameters['etapx0']  = Lattice.etapx0
	parameters['etay0']   = Lattice.etay0
	parameters['etapy0']  = Lattice.etapy0
	parameters['x0']      = Lattice.orbitx0
	parameters['xp0']     = Lattice.orbitpx0
	parameters['y0']      = Lattice.orbity0
	parameters['yp0']     = Lattice.orbitpy0
	parameters['gamma_transition'] = Lattice.gammaT
	parameters['circumference']    = Lattice.getLength()
	parameters['length'] = Lattice.getLength()/Lattice.nHarm
	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 
	
	print dispersionx

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# building the distributions
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Longitudinal_distribution = LongitudinalJohoDistributionSingleHarmonic(parameters, parameters['LongitudinalJohoParameter'])

	# only the main CPU is actually writing its distribution to a file ...
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:
		with open(output_file,"w") as fid:
			csv_writer = csv.writer(fid, delimiter=' ')


			for i in range(parameters['n_macroparticles']):
				(phi[i], dE[i]) = Longitudinal_distribution.getCoordinates()
				(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
				x[i] += closedOrbitx['x0']
				xp[i] += closedOrbitx['xp0']
				y[i] += closedOrbity['y0']
				yp[i] += closedOrbity['yp0']
				dpp = dE[i] / (parameters['energy']) / parameters['beta']**2
				x[i] += dpp * dispersionx['etax0']
				xp[i] += dpp * dispersionx['etapx0']
				y[i] += dpp * dispersiony['etay0']
				yp[i] += dpp * dispersiony['etapy0']
				
				if outputFormat == 'Orbit':
					x[i] *= 1000.
					xp[i] *= 1000.
					y[i] *= 1000.
					yp[i] *= 1000.
					dE[i] /= 1.e6
					csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
				#csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]])
		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(parameters[key]) + '\n'), parameters)
		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file

	orbit_mpi.MPI_Barrier(comm)

	return output_file

def generate_initial_distribution(parameters, Lattice,output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
	parameters['alphax0'] = Lattice.alphax0
	parameters['betax0']  = Lattice.betax0
	parameters['alphay0'] = Lattice.alphay0
	parameters['betay0']  = Lattice.betay0
	parameters['etax0']   = Lattice.etax0
	parameters['etapx0']  = Lattice.etapx0
	parameters['etay0']   = Lattice.etay0
	parameters['etapy0']  = Lattice.etapy0
	parameters['x0']      = Lattice.orbitx0
	parameters['xp0']     = Lattice.orbitpx0
	parameters['y0']      = Lattice.orbity0
	parameters['yp0']     = Lattice.orbitpy0
	parameters['gamma_transition'] = Lattice.gammaT
	parameters['circumference']    = Lattice.getLength()
	parameters['length'] = Lattice.getLength()/Lattice.nHarm
	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	# ~ dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	# ~ dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 
	
	print dispersionx

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# building the distributions
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Longitudinal_distribution = LongitudinalJohoDistributionSingleHarmonic(parameters, parameters['LongitudinalJohoParameter'])

	# only the main CPU is actually writing its distribution to a file ...
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:
		with open(output_file,"w") as fid:
			csv_writer = csv.writer(fid, delimiter=' ')


			for i in range(parameters['n_macroparticles']):
				(phi[i], dE[i]) = Longitudinal_distribution.getCoordinates()
				(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
				x[i] += closedOrbitx['x0']
				xp[i] += closedOrbitx['xp0']
				y[i] += closedOrbity['y0']
				yp[i] += closedOrbity['yp0']
				dpp = dE[i] / (parameters['energy']) / parameters['beta']**2
				x[i] += dpp * dispersionx['etax0']
				xp[i] += dpp * dispersionx['etapx0']
				y[i] += dpp * dispersiony['etay0']
				yp[i] += dpp * dispersiony['etapy0']
				
				if outputFormat == 'Orbit':
					x[i] *= 1000.
					xp[i] *= 1000.
					y[i] *= 1000.
					yp[i] *= 1000.
					dE[i] /= 1.e9
					csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
				#csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]])
		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(parameters[key]) + '\n'), parameters)
		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file

	orbit_mpi.MPI_Barrier(comm)

	return output_file

def generate_initial_distribution_FMA(parameters, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit', triangular_grid = True):

	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['beta']*parameters['etax0'], 'etapx0': parameters['beta']*parameters['etapx0']}
	dispersiony = {'etay0': parameters['beta']*parameters['etay0'], 'etapy0': parameters['beta']*parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles']); phi.fill(parameters['phi_s'])
	dE = np.zeros(parameters['n_macroparticles'])

	emittance_x = parameters['epsn_x'] / parameters['gamma'] / parameters['beta']
	emittance_y = parameters['epsn_y'] / parameters['gamma'] / parameters['beta']
	gamma_x = (1.+parameters['alphax0']**2) / parameters['betax0']
	gamma_y = (1.+parameters['alphay0']**2) / parameters['betay0']

	n_macroparticles_sqrt = np.floor(np.sqrt(parameters['n_macroparticles']))
	Jx = np.linspace(emittance_x/gamma_x/n_macroparticles_sqrt/10, emittance_x/gamma_x, n_macroparticles_sqrt)
	Jy = np.linspace(emittance_y/gamma_y, emittance_y/gamma_y/n_macroparticles_sqrt/10, n_macroparticles_sqrt)
	x, y = np.meshgrid(np.sqrt(Jx)*parameters['TransverseCut'], np.sqrt(Jy)*parameters['TransverseCut'])

	if triangular_grid:
		indcs = np.tril_indices(len(x))
		x = x[indcs]
		y = y[indcs]

	x = x.flatten()
	y = y.flatten()	

	if orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD) == 0:
		fid = open(output_file,"w")
		csv_writer = csv.writer(fid, delimiter=' ')
		for i in range(len(x)):
				
			if outputFormat == 'Orbit':
				x[i] *= 1000.
				xp[i] *= 1000.
				y[i] *= 1000.
				yp[i] *= 1000.
				dE[i] /= 1.e9		
				csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
		fid.close()

		fid = open(summary_file, 'w')
		parameter_list = ['circumference', 'rf_voltage', 'phi_s', 'harmonic_number', 'gamma_transition', 'n_macroparticles', 'energy', 'gamma', 'bunch_length', 'LongitudinalCut', 'LongitudinalJohoParameter', 'x0', 'xp0', 'betax0', 'alphax0', 'etax0', 'etapx0', 'y0', 'yp0', 'betay0', 'alphay0', 'etay0', 'etapy0', 'epsn_x', 'epsn_y', 'TransverseCut']
		for key in parameter_list:
			fid.write(key + ' = ' + str(parameters[key]) + '\n')
		fid.close()

		print '\nCreated particle distribution with ' + str(len(x)) + ' macroparticles into file: ', output_file
	
	return output_file

def generate_initial_distribution_3DGaussian(parameters, Lattice, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
	
	parameters['alphax0'] = Lattice.alphax0
	parameters['betax0']  = Lattice.betax0
	parameters['alphay0'] = Lattice.alphay0
	parameters['betay0']  = Lattice.betay0
	parameters['etax0']   = Lattice.etax0
	parameters['etapx0']  = Lattice.etapx0
	parameters['etay0']   = Lattice.etay0
	parameters['etapy0']  = Lattice.etapy0
	parameters['x0']      = Lattice.orbitx0
	parameters['xp0']     = Lattice.orbitpx0
	parameters['y0']      = Lattice.orbity0
	parameters['yp0']     = Lattice.orbitpy0
	parameters['gamma_transition'] = Lattice.gammaT
	parameters['circumference']    = Lattice.getLength()
	parameters['length'] = Lattice.getLength()/Lattice.nHarm
	
	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 
	
	print dispersionx

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	# ~ z = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# building the distributions
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	
	# We need to convert z into phi
	h_main = np.atleast_1d(parameters['harmonic_number'])[0]
	R = parameters['circumference'] / 2 / np.pi
	sig_E = (parameters['dpp_rms'] * parameters['energy'] * parameters['beta']**2)

	if orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD) == 0:
		fid = open(output_file,"w")
		csv_writer = csv.writer(fid, delimiter=' ')
		for i in range(parameters['n_macroparticles']):

			# Longitudinal distn - use 5 sigma as cut-off (manual)

			outside_limits_E = True
			while outside_limits_E:
				dE[i] = random.gauss(0., sig_E)	# Energy in eV
				if abs(dE[i]) < (5*sig_E):
					# ~ print '\n\tdE = ', dE[i]
					outside_limits_E = False

			outside_limits_z = True
			while outside_limits_z:
				z_temp = random.gauss(0., parameters['blength_rms'])
				if abs(z_temp) < (5*parameters['blength_rms']):
					# ~ print '\n\tz_temp = ', z_temp
					phi[i] = - z_temp * h_main / R 
					outside_limits_z = False

			(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
			x[i] += closedOrbitx['x0']
			xp[i] += closedOrbitx['xp0']
			y[i] += closedOrbity['y0']
			yp[i] += closedOrbity['yp0']
			dpp = dE[i] / (parameters['energy']) / parameters['beta']**2
			x[i] += dpp * dispersionx['etax0']
			xp[i] += dpp * dispersionx['etapx0']
			y[i] += dpp * dispersiony['etay0']
			yp[i] += dpp * dispersiony['etapy0']

			if outputFormat == 'Orbit':
				x[i] *= 1000.
				xp[i] *= 1000.
				y[i] *= 1000.
				yp[i] *= 1000
				dE[i] /= 1.e9
				csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
		#	else:
				# still need to convert from phi to z!!
				#csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]])
		fid.close()

		fid = open(summary_file, 'w')
		parameter_list = ['circumference', 'rf_voltage', 'phi_s', 'harmonic_number', 'gamma_transition', 'n_macroparticles', 'energy', 'gamma', 'bunch_length', 'LongitudinalCut', 'LongitudinalJohoParameter', 'x0', 'xp0', 'betax0', 'alphax0', 'etax0', 'etapx0', 'y0', 'yp0', 'betay0', 'alphay0', 'etay0', 'etapy0', 'epsn_x', 'epsn_y', 'TransverseCut']
		for key in parameter_list:
			fid.write(key + ' = ' + str(parameters[key]) + '\n')
		fid.close()

		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file

	return output_file

def generate_initial_poincare_distributionH(n_sigma, parameters, Lattice, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
        return generate_initial_poincare_distribution(n_sigma, parameters, Lattice, 1, output_file, summary_file, outputFormat)

def generate_initial_poincare_distributionV(n_sigma, parameters, Lattice, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
        return generate_initial_poincare_distribution(n_sigma, parameters, Lattice, 0, output_file, summary_file, outputFormat)

def generate_initial_poincare_distribution(n_sigma, parameters, Lattice, horizontal = 1,  output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):
        parameters['alphax0'] = Lattice.alphax0
	parameters['betax0']  = Lattice.betax0
	parameters['alphay0'] = Lattice.alphay0
	parameters['betay0']  = Lattice.betay0
	parameters['etax0']   = Lattice.etax0
	parameters['etapx0']  = Lattice.etapx0
	parameters['etay0']   = Lattice.etay0
	parameters['etapy0']  = Lattice.etapy0
	parameters['x0']      = Lattice.orbitx0
	parameters['xp0']     = Lattice.orbitpx0
	parameters['y0']      = Lattice.orbity0
	parameters['yp0']     = Lattice.orbitpy0
	parameters['gamma_transition'] = Lattice.gammaT
	parameters['circumference']    = Lattice.getLength()
	parameters['length'] = Lattice.getLength()/Lattice.nHarm
	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['beta']*parameters['etax0'], 'etapx0': parameters['beta']*parameters['etapx0']}
	dispersiony = {'etay0': parameters['beta']*parameters['etay0'], 'etapy0': parameters['beta']*parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# building the distributions
	# ~ Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Transverse_distribution = KVDist1D(twissX)
	Longitudinal_distribution = LongitudinalJohoDistributionSingleHarmonic(parameters, parameters['LongitudinalJohoParameter'])

	# only the main CPU is actually writing its distribution to a file ...
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:
		with open(output_file,"w") as fid:
			csv_writer = csv.writer(fid, delimiter=' ')


			for i in range(parameters['n_macroparticles']):
                                # RANDOM UNIFORM
                                # ~ x[i] = random.uniform(0., n_sigma) * np.sqrt(parameters['betax0'] * parameters['epsn_x'])
                                # EQUAL STEPS
                                if horizontal:
                                        # ~ print 'beta = ',parameters['beta'], ' gamma = ', parameters['gamma']
                                        x[i] = i * float(n_sigma/float(parameters['n_macroparticles'])) * np.sqrt(float(parameters['betax0']) * ( parameters['epsn_x'] / (parameters['beta'] * parameters['gamma'])))
                                        # ~ x[i] = i * float(n_sigma/float(parameters['n_macroparticles']))
                                        # ~ print 'Generate Poincare Distn: maximum sigma = ', n_sigma, ', current particle = ', i * float(n_sigma/float(parameters['n_macroparticles'])), ' sigma'
                                        # ~ print 'eps_geo_x = ', ( parameters['epsn_x'] / (parameters['beta'] * parameters['gamma']))
                                        # ~ print 'beta_x0 = ', float(parameters['betax0'])  
                                        # ~ print 'sigma_x = ', np.sqrt(float(parameters['betax0']) * ( parameters['epsn_x'] / (parameters['beta'] * parameters['gamma'])))
                                elif not horizontal:
                                        # ~ print '\nVERTICAL BUNCH: n_sigma = ',n_sigma, ', sigma = ',  (np.sqrt(parameters['betay0'] * parameters['epsn_y']))
                                        # ~ print '\ty =', i * (n_sigma/parameters['n_macroparticles']) * np.sqrt(parameters['betay0'] * parameters['epsn_y'])
                                        # ~ print '\ti = ', i, ', betay0 = ',  parameters['betay0'], ', epsn_y = ', parameters['epsn_y'], ', macroparticles = ',  parameters['n_macroparticles']
                                        # ~ print '\tsqrt(bet*eps) = ', np.sqrt(parameters['betay0'] * parameters['epsn_y'])
                                        # ~ print '\tn_sigma/macro = ', float(n_sigma/float(parameters['n_macroparticles']))
                                        y[i] = i * float(n_sigma/float(parameters['n_macroparticles'])) * np.sqrt(float(parameters['betay0']) * ( parameters['epsn_y'] / (parameters['beta'] * parameters['gamma'])))

				if outputFormat == 'Orbit':
					x[i] *= 1000.
					xp[i] *= 1000.
					y[i] *= 1000.
					yp[i] *= 1000.
					dE[i] /= 1.e9
					csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
				#csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]])
		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(parameters[key]) + '\n'), parameters)
		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file


	orbit_mpi.MPI_Barrier(comm)

	return output_file
