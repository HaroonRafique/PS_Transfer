import math
import sys
import time
import orbit_mpi
import timeit
import numpy as np
import scipy.io as sio
from scipy.stats import moment
import os

# Use switches in simulation_parameters.py in current folder
#-------------------------------------------------------------
from simulation_parameters import switches as s
from simulation_parameters import parameters as p

# utils
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.utils.consts import mass_proton, speed_of_light, pi

# bunch
from bunch import Bunch
from bunch import BunchTwissAnalysis, BunchTuneAnalysis
from orbit.bunch_utils import ParticleIdNumber

# diagnostics
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNodeAsChild
from orbit.diagnostics import addTeapotMomentsNodeSet, addTeapotStatLatsNodeSet

# PTC lattice
from libptc_orbit import *
from ext.ptc_orbit import PTC_Lattice
from ext.ptc_orbit import PTC_Node
from ext.ptc_orbit.ptc_orbit import setBunchParamsPTC, readAccelTablePTC, readScriptPTC
from ext.ptc_orbit.ptc_orbit import updateParamsPTC, synchronousSetPTC, synchronousAfterPTC
from ext.ptc_orbit.ptc_orbit import trackBunchThroughLatticePTC, trackBunchInRangePTC
from orbit.aperture import TeapotApertureNode

# transverse space charge
from spacecharge import SpaceChargeCalcSliceBySlice2D

from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalcAnalyticGaussian
from spacecharge import InterpolatedLineDensityProfile

from lib.output_dictionary import *
from lib.pyOrbit_GenerateInitialDistribution import *
from lib.save_bunch_as_matfile import *
from lib.pyOrbit_Tunespread_Calculator import *
from lib.suppress_stdout import suppress_STDOUT
readScriptPTC_noSTDOUT = suppress_STDOUT(readScriptPTC)

# MPI stuff
#-----------------------------------------------------------------------
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
print '\n\tStart PyORBIT simulation on MPI process: ', rank

# Function to check that a file isn't empty (common PTC file bug)
def is_non_zero_file(fpath):  
	print '\n\t\t\tis_non_zero_file:: Checking file ', fpath
	print '\n\t\t\tis_non_zero_file:: File exists = ', os.path.isfile(fpath)
	print '\n\t\t\tis_non_zero_file:: Size > 3 bytes = ', os.path.getsize(fpath)
	return os.path.isfile(fpath) and os.path.getsize(fpath) > 3

# Function to check and read PTC file
def CheckAndReadPTCFile(f):
	if is_non_zero_file(f): 
		readScriptPTC_noSTDOUT(f)
	else:
		print '\n\t\t\CheckAndReadPTCFile:: ERROR: PTC file ', f, ' is empty or does not exist, exiting'
		exit(0)

# Function to open TWISS_PTC_table.OUT and return fractional tunes
def GetTunesFromPTC():
	readScriptPTC_noSTDOUT('../PTC/twiss_script.ptc')
	with open('TWISS_PTC_table.OUT') as f:
		first_line = f.readline()
		Qx = (float(first_line.split()[2]))
		Qy = (float(first_line.split()[3]))
	os.system('rm TWISS_PTC_table.OUT')
	return Qx, Qy


# Function to return second moment (mu^2) of distribution
# ~ def GetBunchMus(b, smooth=True):
	# ~ print '\n\t\t\t\t GetBunchMus called'
	# ~ print '\n\t\t\t\t GetBunchMus b.getSize() = ', b.getSize()

	# ~ # MPI stuff to run on a single node
	# ~ rank = 0
	# ~ numprocs = 1

	# ~ mpi_init = orbit_mpi.MPI_Initialized()
	# ~ comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	# ~ orbit_mpi.MPI_Barrier(comm)

	# ~ if(mpi_init):
		# ~ rank = orbit_mpi.MPI_Comm_rank(comm)
		# ~ print '\n\t\t\t\t GetBunchMus rank = ', rank 
		# ~ numprocs = orbit_mpi.MPI_Comm_size(comm)
		# ~ print '\n\t\t\t\t GetBunchMus numprocs = ', numprocs 

	# ~ nparts_arr_local = []
	# ~ for i in range(numprocs):
		# ~ nparts_arr_local.append(0)
	# ~ print '\n\t\t\t\t GetBunchMus nparts_arr_local before initialisation = ', nparts_arr_local 

	# ~ nparts_arr_local[rank] = b.getSize()
	# ~ print '\n\t\t\t\t GetBunchMus nparts_arr_local after initialisation  = ', nparts_arr_local 
	# ~ data_type = mpi_datatype.MPI_INT
	# ~ op = mpi_op.MPI_SUM

	# ~ nparts_arr = orbit_mpi.MPI_Allreduce(nparts_arr_local,data_type,op,comm)
	# ~ print '\n\t\t\t\t GetBunchMus nparts_arr  = ', nparts_arr

# ~ # Arrays to hold x and y data
	# ~ x = []
	# ~ y = []

	# ~ for i in range(b.getSize()):
		# ~ x.append(b.x(i))
		# ~ y.append(b.y(i))

# ~ # Calculate moments of the bunch
	# ~ return moment(x, 2), moment(y, 2)
	

# Function to return second moment (mu^2) of distribution
def GetBunchMus3(b, smooth=True):

	print '\n\t\t\t\t GetBunchMus called'
	#take the MPI Communicator from bunch: it could be different from MPI_COMM_WORLD
	comm = b.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	print '\n\t\t\t\t GetBunchMus rank = ', rank 
	size = orbit_mpi.MPI_Comm_size(comm)
	print '\n\t\t\t\t GetBunchMus rank = ', size 
	main_rank = 0

	print '\n\t\t\t\t GetBunchMus b.getPartAttrNames() = ', b.getPartAttrNames()

	# n_parts_arr - array of size of the number of CPUs, 
	# and have the number of macroparticles on each CPU
	n_parts_arr = [0]*size
	print '\n\t\t\t\t GetBunchMus n_parts_arr = ', n_parts_arr
	n_parts_arr[rank] = b.getSize()
	print '\n\t\t\t\t GetBunchMus n_parts_arr after rank setting = ', n_parts_arr
	n_parts_arr = orbit_mpi.MPI_Allreduce(n_parts_arr,mpi_datatype.MPI_INT,mpi_op.MPI_SUM,comm)
	print '\n\t\t\t\t GetBunchMus n_parts_arr after MPI Sum = ', n_parts_arr

	# ~ mp_array = range(n_parts_arr[rank])	#indexes of all particles
	mp_array = range(25000)	#indexes of all particles
	particles = {}
	print '\n\t\t\t\t GetBunchMus particles  = ', particles
	particles['x'] = map(b.x, mp_array)
	print '\n\t\t\t\t GetBunchMus particles[\'x\'] = ', particles['x']
	particles['y'] = map(b.y, mp_array)
	print '\n\t\t\t\t GetBunchMus particles[\'y\'] = ', particles['y']

	print '\n\t\t\t\t GetBunchMus particles.keys() = ', particles.keys()
	return 0, 0

# Function to return second moment (mu^2) of distribution
def GetBunchMus(b, smooth=True):

	print '\n\t\t\t\t GetBunchMus called'
	#take the MPI Communicator from bunch: it could be different from MPI_COMM_WORLD
	comm = b.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	print '\n\t\t\t\t GetBunchMus rank = ', rank 
	size = orbit_mpi.MPI_Comm_size(comm)
	print '\n\t\t\t\t GetBunchMus size = ', size 
	main_rank = 0

	# mu_x - 2nd moment in x
	# and have the number of macroparticles on each CPU
	mu_x = [0]*size
	mu_x_mpi = 0.
	print '\n\t\t\t\t GetBunchMus mu_x = ', mu_x

	x = []
	y = []
	for i in range(b.getSize()):
		x.append(b.x(i))
		y.append(b.y(i))

	mu_x[rank] = moment(x, 2)
	print '\n\t\t\t\t GetBunchMus mu_x after rank setting = ', mu_x

	orbit_mpi.MPI_Barrier(comm)
	# ~ mu_x = orbit_mpi.MPI_Allreduce(mu_x, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)
	orbit_mpi.MPI_Allreduce(mu_x, mu_x_mpi, 1, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)
	print '\n\t\t\t\t GetBunchMus mu_x after MPI Sum = ', mu_x

	return 0, 0

def GetBunchMus2(b, smooth=True):
	window = 40
	print '\n\t\t\t\t GetBunchMus called'
	print '\n\t\t\t\t GetBunchMus b.getSize() = ', b.getSize()

	# MPI stuff to run on a single node
	rank = 0
	numprocs = 1

	mpi_init = orbit_mpi.MPI_Initialized()
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	orbit_mpi.MPI_Barrier(comm)
	
	if(mpi_init):
		rank = orbit_mpi.MPI_Comm_rank(comm)
		print '\n\t\t\t\t GetBunchMus rank = ', rank 
		numprocs = orbit_mpi.MPI_Comm_size(comm)
		print '\n\t\t\t\t GetBunchMus numprocs = ', numprocs 

	x = []
	y = []
	for i in range(b.getSize()):
		x.append(b.x(i))
		y.append(b.y(i))

	x_mom_local = 0.
	y_mom_local = 0.
	x_mom_mpi = 0.
	y_mom_mpi = 0.
	
	# ~ x_mom_local = []
	# ~ y_mom_local = []
	# ~ x_mom_mpi = []
	# ~ y_mom_mpi = []

	# ~ for i in range(numprocs):
		# ~ x_mom_local.append(0)
		# ~ y_mom_local.append(0)
		# ~ x_mom_mpi.append(0)
		# ~ y_mom_mpi.append(0)

	x_mom_local = moment(x, 2)
	y_mom_local= moment(y, 2)

	print '\n\t\t\t\t GetBunchMus x_mom_local = ', x_mom_local
	print '\n\t\t\t\t GetBunchMus y_mom_local = ', y_mom_local
	print '\n\t\t\t\t GetBunchMus x_mom_mpi = ', x_mom_mpi
	print '\n\t\t\t\t GetBunchMus y_mom_mpi = ', y_mom_mpi

	# ~ test = orbit_mpi.MPI_Allreduce(x_mom_local, x_mom_mpi, numprocs, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)
	# ~ x_mom_mpi = orbit_mpi.MPI_Allreduce(x_mom_local, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)
	x_mom_mpi = orbit_mpi.MPI_Allreduce(x_mom_local, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)
	y_mom_mpi = orbit_mpi.MPI_Allreduce(y_mom_local, mpi_datatype.MPI_DOUBLE, mpi_op.MPI_SUM, comm)

	# ~ print '\n\t\t\t\t GetBunchMus test = ', test

	print '\n\t\t\t\t GetBunchMus x_mom_local = ', x_mom_local
	print '\n\t\t\t\t GetBunchMus y_mom_local = ', y_mom_local
	print '\n\t\t\t\t GetBunchMus x_mom_mpi = ', x_mom_mpi
	print '\n\t\t\t\t GetBunchMus y_mom_mpi = ', y_mom_mpi

# Calculate moments of the bunch
	print '\n\t\t\t\t GetBunchMus::return moments on MPI rank: ', rank
	return moment(x, 2), moment(y, 2)

# Create folder structure
#-----------------------------------------------------------------------
print '\n\t\tmkdir on MPI process: ', rank
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('input')
mpi_mkdir_p('bunch_output')
mpi_mkdir_p('output')
mpi_mkdir_p('lost')

# Dictionary for simulation status
#-----------------------------------------------------------------------
import pickle # HAVE TO CLEAN THIS FILE BEFORE RUNNING A NEW SIMULATION
status_file = 'input/simulation_status.pkl'
if not os.path.exists(status_file):
	sts = {'turn': -1}
else:
	with open(status_file) as fid:
		sts = pickle.load(fid)

# Generate PTC RF table
#-----------------------------------------------------------------------
print '\n\t\tCreate RF file on MPI process: ', rank
from lib.write_ptc_table import write_RFtable
from simulation_parameters import RFparameters as RF 
write_RFtable('input/RF_table.ptc', *[RF[k] for k in ['harmonic_factors','time','Ekin_GeV','voltage_MV','phase']])

# Initialize a Teapot-Style PTC lattice
#-----------------------------------------------------------------------
print '\n\t\tRead PTC flat file: ',p['flat_file'],' on MPI process: ', rank
PTC_File = p['flat_file']
Lattice = PTC_Lattice("PS")
Lattice.readPTC(PTC_File)

print '\n\t\tRead PTC files on MPI process: ', rank
CheckAndReadPTCFile('PTC/fringe.ptc')
CheckAndReadPTCFile('PTC/time.ptc')
CheckAndReadPTCFile('PTC/ramp_cavities.ptc')

# Create a dictionary of parameters
#-----------------------------------------------------------------------
print '\n\t\tMake paramsDict on MPI process: ', rank
paramsDict = {}
paramsDict["length"]=Lattice.getLength()/Lattice.nHarm

# Add apertures
#-----------------------------------------------------------------------
print '\n\t\tAdd apertures on MPI process: ', rank
position = 0
for node in Lattice.getNodes():
	myaperturenode = TeapotApertureNode(1, 10, 10, position)
	node.addChildNode(myaperturenode, node.ENTRANCE)
	node.addChildNode(myaperturenode, node.BODY)
	node.addChildNode(myaperturenode, node.EXIT)
	position += node.getLength()

# Import a bunch and relevant parameters for it
#-----------------------------------------------------------------------
if sts['turn'] < 0:
	print '\n\t\tCreate bunch on MPI process: ', rank
	bunch = Bunch()
	setBunchParamsPTC(bunch)

	p['harmonic_number'] = Lattice.nHarm 
	p['phi_s']           = 0
	p['gamma']           = bunch.getSyncParticle().gamma()
	p['beta']            = bunch.getSyncParticle().beta()
	p['energy']          = 1e9 * bunch.mass() * bunch.getSyncParticle().gamma()
	# ~ p['bunch_length'] = p['sig_z']/speed_of_light/bunch.getSyncParticle().beta()*4
	p['bunch_length'] = p['bunch_length']
	kin_Energy = bunch.getSyncParticle().kinEnergy()

	print '\n\t\tOutput simulation_parameters on MPI process: ', rank
	for i in p:
		print '\t', i, '\t = \t', p[i]

	print '\n\t\tLoad distribution from ', p['input_distn'] ,' on MPI process: ', rank
	path_to_distn = p['input_distn']
	bunch = bunch_from_matfile(path_to_distn)

# Add Macrosize to bunch
#-----------------------------------------------------------------------
	bunch.addPartAttr("macrosize")
	map(lambda i: bunch.partAttrValue("macrosize", i, 0, p['macrosize']), range(bunch.getSize()))
	ParticleIdNumber().addParticleIdNumbers(bunch) # Give them unique number IDs

# Dump and save as Matfile
#-----------------------------------------------------------------------
	# ~ bunch.dumpBunch("input/mainbunch_start.dat")
	print '\n\t\tSave bunch in bunch_output/mainbunch_-000001.mat on MPI process: ', rank
	saveBunchAsMatfile(bunch, "bunch_output/mainbunch_-000001")
	print '\n\t\tSave bunch in input/mainbunch.mat on MPI process: ', rank
	saveBunchAsMatfile(bunch, "input/mainbunch")
	sts['mainbunch_file'] = "input/mainbunch"

# Create empty lost bunch
#-----------------------------------------------------------------------
	lostbunch = Bunch()
	bunch.copyEmptyBunchTo(lostbunch)
	lostbunch.addPartAttr('ParticlePhaseAttributes')
	lostbunch.addPartAttr("LostParticleAttributes")	
	saveBunchAsMatfile(lostbunch, "input/lostbunch")
	sts['lostbunch_file'] = "input/lostbunch"

# Add items to pickle parameters
#-----------------------------------------------------------------------
	sts['turns_max'] = p['turns_max']
	sts['turns_update'] = p['turns_update']
	sts['turns_print'] = p['turns_print']
	sts['circumference'] = p['circumference']

bunch = bunch_from_matfile(sts['mainbunch_file'])
lostbunch = bunch_from_matfile(sts['lostbunch_file'])
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch

#############################-------------------########################
#############################	SPACE CHARGE	########################
#############################-------------------########################

# Add space charge nodes
#----------------------------------------------------
if s['SliceBySlice']:
	print '\n\t\tAdding slice-by-slice space charge nodes on MPI process: ', rank
	# Make a SC solver
	calcsbs = SpaceChargeCalcSliceBySlice2D(s['GridSizeX'], s['GridSizeY'], s['GridSizeZ'], useLongitudinalKick=s['LongitudinalKick'])
	sc_path_length_min = 1E-8
	# Add the space charge solver to the lattice as child nodes
	sc_nodes = scLatticeModifications.setSC2p5DAccNodes(Lattice, sc_path_length_min, calcsbs)
	print '\n\t\tInstalled', len(sc_nodes), 'space charge nodes ...'


# Add tune analysis child node
#-----------------------------------------------------
parentnode_number = 97
parentnode = Lattice.getNodes()[parentnode_number]
Twiss_at_parentnode_entrance = Lattice.getNodes()[parentnode_number-1].getParamsDict()
tunes = TeapotTuneAnalysisNode("tune_analysis")

tunes.assignTwiss(*[Twiss_at_parentnode_entrance[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])
addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)

# Define twiss analysis and output dictionary
#-----------------------------------------------------------------------
print '\n\t\tbunchtwissanalysis on MPI process: ', rank
bunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())

output_file = 'output/output.mat'
output = Output_dictionary()
output.addParameter('turn', lambda: turn)
output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
output.addParameter('D_x', lambda: bunchtwissanalysis.getDispersion(0))
output.addParameter('D_y', lambda: bunchtwissanalysis.getDispersion(1))
output.addParameter('mu_x', lambda: GetBunchMus(bunch)[0])
output.addParameter('mu_y', lambda: GetBunchMus(bunch)[1])
output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))
output.addParameter('beta_x', lambda: bunchtwissanalysis.getBeta(0))
output.addParameter('beta_y', lambda: bunchtwissanalysis.getBeta(1))
output.addParameter('alpha_x', lambda: bunchtwissanalysis.getAlpha(0))
output.addParameter('alpha_y', lambda: bunchtwissanalysis.getAlpha(1))
output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
output.addParameter('eff_beta_x', lambda: bunchtwissanalysis.getEffectiveBeta(0))
output.addParameter('eff_beta_y', lambda: bunchtwissanalysis.getEffectiveBeta(1))
output.addParameter('eff_epsn_x', lambda: bunchtwissanalysis.getEffectiveEmittance(0))
output.addParameter('eff_epsn_y', lambda: bunchtwissanalysis.getEffectiveEmittance(1))
output.addParameter('eff_alpha_x', lambda: bunchtwissanalysis.getEffectiveAlpha(0))
output.addParameter('eff_alpha_y', lambda: bunchtwissanalysis.getEffectiveAlpha(1))
output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())

if os.path.exists(output_file):
	output.import_from_matfile(output_file)

# Track
#-----------------------------------------------------------------------
print '\n\t\tStart tracking on MPI process: ', rank
start_time = time.time()
last_time = time.time()

#turn = -1
#bunchtwissanalysis.analyzeBunch(bunch)
# ~ output.addParameter('turn_time', lambda: time.strftime("%H:%M:%S"))
# ~ output.addParameter('turn_duration', lambda: (time.time() - last_time))
# ~ output.addParameter('cumulative_time', lambda: (time.time() - start_time))
# ~ start_time = time.time()
#output.update()
print '\n\t\tstart time = ', start_time

for turn in range(sts['turn']+1, sts['turns_max']):
	if not rank:
		print '\n\t\tTURN ', turn
		last_time = time.time()

	if turn == 0:
		output.addParameter('turn_time', lambda: time.strftime("%H:%M:%S"))
		output.addParameter('turn_duration', lambda: (time.time() - last_time))
		output.addParameter('cumulative_time', lambda: (time.time() - start_time))
		start_time = time.time()
		print '\n\t\tstart time = ', start_time

	if not rank: print '\n\t\tTURN ', turn, ' trackBunch'
	Lattice.trackBunch(bunch, paramsDict)
	if not rank: print '\n\t\tTURN ', turn, ' analyzeBunch'
	bunchtwissanalysis.analyzeBunch(bunch)  # analyze twiss and emittance

	if turn in sts['turns_update']:	sts['turn'] = turn

	if not rank: print '\n\t\tTURN ', turn, ' output.update()'
	output.update()

	if turn in sts['turns_print']:
		saveBunchAsMatfile(bunch, "input/mainbunch")
		saveBunchAsMatfile(bunch, "bunch_output/mainbunch_%s"%(str(turn).zfill(6)))
		saveBunchAsMatfile(lostbunch, "lost/lostbunch_%s"%(str(turn).zfill(6)))
		output.save_to_matfile(output_file)
		if not rank:
			with open(status_file, 'w') as fid:
				pickle.dump(sts, fid)
