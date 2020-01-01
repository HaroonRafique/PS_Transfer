#-----------------------------------------------------------------------
# Class to output single particle co-ordinates, based on 
# Hannes Bartosik's (CERN BE-ABP-HSI) output dictionary.
# 24.07.2018: Created by Haroon Rafique, CERN BE-ABP-HSI 
#
# Function AddNewParticle(n) adds particle n to dictionary.
# Function Update(bunch, turn) stores data for particles at turn.
# Function PrintParticleForTurn(turn, n, filename) self-explanatory.
# Function PrintParticle(n, filename) self-explanatory.
# Function PrintAllParticles(filename) prints all particles for all
# turns for which Update() function was called.
#-----------------------------------------------------------------------

import orbit_mpi
import os
import numpy as np

class PTCLatticeFunctionsDictionary(object):
	
	def __init__(self, PTCLatticeFunctionsDictionary = None):
		if PTCLatticeFunctionsDictionary:
			print "PTCLatticeFunctionsDictionary::__init__: constructor with existing dictionary not yet implemented"
		else:
			self.twiss_dict = {}		# Top level dictionary : N : All data
			self.turn_list = [] 		# Record indices of stored turns

			print 'PTCLatticeFunctionsDictionary: Created initial twiss dictionary \'twiss_dict'

	def UpdatePTCTwiss(self, Lattice, turn, verbose=False):
		self.update_flag = 1

		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:

			# Create the turn dictionary
			self.twiss_dict[int(turn)] = {}	# Second level : N-2 : Turn
			
			# Third level: twiss
			self.twiss_dict[int(turn)]['beta_x'] = ([n.getParamsDict()['betax'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['beta_y'] = ([n.getParamsDict()['betay'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['alpha_x'] = ([n.getParamsDict()['alphax'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['alpha_y'] = ([n.getParamsDict()['alphay'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['eta_x'] = ([n.getParamsDict()['etax'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['eta_y'] = ([n.getParamsDict()['etay'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['eta_px'] = ([n.getParamsDict()['etapx'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['eta_py'] = ([n.getParamsDict()['etapy'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['orbit_x'] = ([n.getParamsDict()['orbitx'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['orbit_px'] = ([n.getParamsDict()['orbitpx'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['orbit_y'] = ([n.getParamsDict()['orbity'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['orbit_py'] = ([n.getParamsDict()['orbitpy'] for n in Lattice.getNodes()])
			self.twiss_dict[int(turn)]['s'] = np.cumsum([n.getLength() for n in Lattice.getNodes()])

		self.turn_list.append(turn)
		if verbose:
			print "PTCLatticeFunctionsDictionary::update: Added turn %i" % (turn)
				
	# Function to print 6D co-ordinates for a particle for 1 given turn
	def PrintPTCTwissForTurn(self, turn, lattice_folder='.', filename=None):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			if filename is None:
				filename = lattice_folder + '/PTC_Twiss_turn_' + str(turn) + '.dat'

			# Check that the particle exists
			if turn not in self.turn_list:
				print "PTCLatticeFunctionsDictionary::PrintPTCTwissForTurn: Turn not stored, use UpdatePTCTwiss function on this turn to store."
			else:
				# if file exists then overwrite
				if os.path.exists(filename):
					f = open(filename,"w")

				# if file doesn't exist create and add header
				else:
					f = open(filename,"w")
					f.write('# s\tbeta_x\tbeta_y\talpha_x\talpha_y\tD_x\tD_y\tD_px\tD_py\torbit_x\torbit_px\torbit_y\torbit_py')


				for i in range(0, len(self.twiss_dict[int(turn)]['s']), 1):
					f.write("\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % ( 	\
					self.twiss_dict[int(turn)]['s'][i],			\
					self.twiss_dict[int(turn)]['beta_x'][i],	\
					self.twiss_dict[int(turn)]['beta_y'][i],	\
					self.twiss_dict[int(turn)]['alpha_x'][i],	\
					self.twiss_dict[int(turn)]['alpha_y'][i],	\
					self.twiss_dict[int(turn)]['eta_x'][i],		\
					self.twiss_dict[int(turn)]['eta_y'][i],		\
					self.twiss_dict[int(turn)]['eta_px'][i],	\
					self.twiss_dict[int(turn)]['eta_py'][i],	\
					self.twiss_dict[int(turn)]['orbit_x'][i],	\
					self.twiss_dict[int(turn)]['orbit_px'][i],	\
					self.twiss_dict[int(turn)]['orbit_y'][i],	\
					self.twiss_dict[int(turn)]['orbit_py'][i]))
				f.close()

	# Function to print PTC twiss for all recorded turns
	def PrintAllPTCTwiss(self, lattice_folder='.'):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			for j in self.turn_list:
				self.PrintPTCTwissForTurn(j, lattice_folder)

	def PrintOrbitExtrema(self, lattice_folder='.'):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:

			filename = lattice_folder + '/Orbit_Extrema.dat'
			f = open(filename, "w")
			f.write('# Turn\tOrbit_Min_x\tOrbit_Max_x\tOrbit_Min_y\tOrbit_Max_y')
			for turn in self.turn_list:
				f.write("\n%i\t%f\t%f\t%f\t%f" % (turn,			\
				np.min(self.twiss_dict[int(turn)]['orbit_x']),	\
				np.max(self.twiss_dict[int(turn)]['orbit_x']),	\
				np.min(self.twiss_dict[int(turn)]['orbit_y']),	\
				np.max(self.twiss_dict[int(turn)]['orbit_y'])))
			f.close()

	def GetAverageParameter(self, name, turn): return np.mean(self.twiss_dict[int(turn)][str(name)])
	def GetMinParameter(self, name, turn): return np.min(self.twiss_dict[int(turn)][str(name)])
	def GetMaxParameter(self, name, turn): return np.max(self.twiss_dict[int(turn)][str(name)])

	def ReturnTwissDict(self): return self.twiss_dict
	
	def ReturnTurnList(self): return self.turn_list



















