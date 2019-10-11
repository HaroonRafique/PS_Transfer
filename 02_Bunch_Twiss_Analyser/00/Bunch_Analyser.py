import os
import math
import sys
import time
import orbit_mpi
import timeit
import numpy as np
import scipy.io as sio
from scipy.stats import moment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

# FUNCTION DEFINITIONS
#-----------------------------------------------------------------------

# Function to check that a file isn't empty (common PTC file bug)
def is_non_zero_file(fpath):  
	print 'Checking file ', fpath
	print 'File exists = ', os.path.isfile(fpath)
	print 'Size > 3 bytes = ', os.path.getsize(fpath)
	return os.path.isfile(fpath) and os.path.getsize(fpath) > 3

# Function to check and read PTC file
def CheckAndReadPTCFile(f):
	if is_non_zero_file(f): 
		readScriptPTC_noSTDOUT(f)
	else:
		print 'ERROR: PTC file ', f, ' is empty or does not exist, exiting'
		exit(0)

# Function to open TWISS_PTC_table.OUT and return fractional tunes
def GetTunesFromPTC():
	readScriptPTC_noSTDOUT('PTC/twiss_script.ptc')
	with open('TWISS_PTC_table.OUT') as f:
		first_line = f.readline()
		Qx = (float(first_line.split()[2]))
		Qy = (float(first_line.split()[3]))
	os.system('rm TWISS_PTC_table.OUT')
	return Qx, Qy

def GetBunchSigmas(b, smooth=True):
	window = 40

	# MPI stuff to run on a single node
	rank = 0
	numprocs = 1

	mpi_init = orbit_mpi.MPI_Initialized()
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	orbit_mpi.MPI_Barrier(comm)

	if(mpi_init):
		rank = orbit_mpi.MPI_Comm_rank(comm)
		numprocs = orbit_mpi.MPI_Comm_size(comm)

	nparts_arr_local = []
	for i in range(numprocs):
		nparts_arr_local.append(0)

	nparts_arr_local[rank] = b.getSize()
	data_type = mpi_datatype.MPI_INT
	op = mpi_op.MPI_SUM

	nparts_arr = orbit_mpi.MPI_Allreduce(nparts_arr_local,data_type,op,comm)

# Arrays to hold x and y data
	x = []
	y = []

	for i in range(b.getSize()):
		x.append(b.x(i))
		y.append(b.y(i))

# Calculate moments of the bunch
	return moment(x, 2), moment(y, 2)

# MPI stuff
#-----------------------------------------------------------------------
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
print 'Start on MPI process: ', rank

# Create folder structure
#-----------------------------------------------------------------------
# ~ print '\nmkdir on MPI process: ', rank
# ~ from lib.mpi_helpers import mpi_mkdir_p
# ~ mpi_mkdir_p('input')
# ~ mpi_mkdir_p('bunch_output')
# ~ mpi_mkdir_p('output')
# ~ mpi_mkdir_p('lost')

# load bunch from file
#-----------------------------------------------------------------------

# ~ path_to_distn = '../01_Tracking/MD4224_Nominal_WP_Tomo_Distn.mat'
path_to_distn = 'PyORBIT_Tomo_file_MD4224_HB.mat'
bunch = bunch_from_matfile(path_to_distn)

# Add Macrosize to bunch
#-----------------------------------------------------------------------
	# ~ bunch.addPartAttr("macrosize")
	# ~ map(lambda i: bunch.partAttrValue("macrosize", i, 0, p['macrosize']), range(bunch.getSize()))
	# ~ ParticleIdNumber().addParticleIdNumbers(bunch) # Give them unique number IDs

# Dump and save as Matfile
#-----------------------------------------------------------------------
	# ~ bunch.dumpBunch("input/mainbunch_start.dat")
	# ~ print 'Save bunch in bunch_output/mainbunch_-000001.mat'
	# ~ saveBunchAsMatfile(bunch, "bunch_output/mainbunch_-000001")
	# ~ print 'Save bunch in input/mainbunch.mat'
	# ~ saveBunchAsMatfile(bunch, "input/mainbunch")
	# ~ sts['mainbunch_file'] = "input/mainbunch"

# Create empty lost bunch
#-----------------------------------------------------------------------
	# ~ lostbunch = Bunch()
	# ~ bunch.copyEmptyBunchTo(lostbunch)
	# ~ lostbunch.addPartAttr('ParticlePhaseAttributes')
	# ~ lostbunch.addPartAttr("LostParticleAttributes")	
	# ~ saveBunchAsMatfile(lostbunch, "input/lostbunch")
	# ~ sts['lostbunch_file'] = "input/lostbunch"

# Add items to pickle parameters
#-----------------------------------------------------------------------
	# ~ sts['turns_max'] = p['turns_max']
	# ~ sts['turns_update'] = p['turns_update']
	# ~ sts['turns_print'] = p['turns_print']
	# ~ sts['circumference'] = p['circumference']

# ~ bunch = bunch_from_matfile(sts['mainbunch_file'])
# ~ lostbunch = bunch_from_matfile(sts['lostbunch_file'])
# ~ paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch

# Add tune analysis child node
#-----------------------------------------------------
# ~ parentnode_number = 97
# ~ parentnode = Lattice.getNodes()[parentnode_number]
# ~ Twiss_at_parentnode_entrance = Lattice.getNodes()[parentnode_number-1].getParamsDict()
# ~ tunes = TeapotTuneAnalysisNode("tune_analysis")

# ~ tunes.assignTwiss(*[Twiss_at_parentnode_entrance[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
# ~ tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])
# ~ addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)

# Define twiss analysis and output dictionary
#-----------------------------------------------------------------------
print '\nTWISS on MPI process: ', rank
bunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())

output_file = 'output/output.mat'
output = Output_dictionary()
output.addParameter('turn', lambda: turn)
output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())
output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))
output.addParameter('beta_x', lambda: bunchtwissanalysis.getBeta(0))
output.addParameter('beta_y', lambda: bunchtwissanalysis.getBeta(1))
output.addParameter('alpha_x', lambda: bunchtwissanalysis.getAlpha(0))
output.addParameter('alpha_y', lambda: bunchtwissanalysis.getAlpha(1))
output.addParameter('D_x', lambda: bunchtwissanalysis.getDispersion(0))
output.addParameter('D_y', lambda: bunchtwissanalysis.getDispersion(1))
output.addParameter('eff_beta_x', lambda: bunchtwissanalysis.getEffectiveBeta(0))
output.addParameter('eff_beta_y', lambda: bunchtwissanalysis.getEffectiveBeta(1))
output.addParameter('eff_epsn_x', lambda: bunchtwissanalysis.getEffectiveEmittance(0))
output.addParameter('eff_epsn_y', lambda: bunchtwissanalysis.getEffectiveEmittance(1))
output.addParameter('eff_alpha_x', lambda: bunchtwissanalysis.getEffectiveAlpha(0))
output.addParameter('eff_alpha_y', lambda: bunchtwissanalysis.getEffectiveAlpha(1))
output.addParameter('sigma_x', lambda: GetBunchSigmas(bunch)[0])
output.addParameter('sigma_y', lambda: GetBunchSigmas(bunch)[1])

# Track
#-----------------------------------------------------------------------
print '\nTracking on MPI process: ', rank
start_time = time.time()
last_time = time.time()

turn = -1
bunchtwissanalysis.analyzeBunch(bunch)
output.addParameter('turn_time', lambda: time.strftime("%H:%M:%S"))
output.addParameter('turn_duration', lambda: (time.time() - last_time))
output.addParameter('cumulative_time', lambda: (time.time() - start_time))
start_time = time.time()
output.update()
output.save_to_matfile(output_file)
