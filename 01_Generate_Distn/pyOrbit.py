########################################################################
#	PTC-PyORBIT simulation script that creates initial distributions   #
#	for CERN Proton Synchrotron Injection simulations of Machine       #
#	Development study MD4224.                                          #
#																	   #
#	Script creates a Gaussian, Joho, and Tomo distribution for each set#
#	of input parameters.											   #
#	Created 16.10.19: Haroon Rafique CERN BE-ABP-HSI                   #
########################################################################

import math
import sys
import time
import orbit_mpi
import timeit
import numpy as np
import scipy.io as sio
from scipy.stats import moment
import os

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
def GetBunchMus(b, smooth=True):
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

# Function to return a bunch from some input parameters
def Create_Bunch(Lattice, p=None, TwissDict=None, label=None, DistType = 'Gaussian', TwissType = 'Lattice', rank=0):

	print '\n\t\tCreate_Bunch on MPI rank ', rank
	bunch = Bunch()
	setBunchParamsPTC(bunch) # Need to check if this works here

	c = 299792458

	if label is None:
		# Include machine (PS), tunes, lattice start position (BWS65H)
		label = 'PS_Tune_607_624_BWS65H'

	if p is None:
		# If we make the parameters somewhere else they should look like this:
		p = {}
		p['n_macroparticles']	= int(0.5E6)
		p['gamma']			= 2.49253731343
		p['intensity']			= 72.5E+10
		p['bunch_length']		= 140e-9
		p['blength']			= 140e-9
		p['epsn_x']				= 1E-6
		p['epsn_y']				= 1.2E-6
		p['dpp_rms']			= 8.7e-04
		p['tomo_file']			= 'PyORBIT_Tomo_file_MD4224_HB.mat'
		p['LongitudinalJohoParameter'] = 1.2
		p['LongitudinalCut'] = 2.4
		p['TransverseCut']		= 5
		p['rf_voltage']			= 0.0212942055190595723
		p['circumference']		= 2*np.pi*100
		p['phi_s']				= 0
	
	# Have to be calculated in function
	p['DistType']			= DistType
	p['kin_Energy']			= bunch.getSyncParticle().kinEnergy()
	p['harmonic_number']	= Lattice.nHarm 
	p['gamma']				= bunch.getSyncParticle().gamma()
	p['energy']				= 1e9 * bunch.mass() * bunch.getSyncParticle().gamma()

	# Check
	p['beta']				= bunch.getSyncParticle().beta()
	beta_check 				= np.sqrt(p['gamma']**2-1)/p['gamma']
	if ((beta_check - p['beta'])/beta_check) < 0.05: print '\n\tCreate_Bunch:: bunch length check failed.\n\t\tp[\'beta\'] = ', p['beta'], ' beta_check = ', beta_check
	# ~ if beta_check is not p['beta'] :  print '\n\tCreate_Bunch:: bunch length check failed.\n\t\tp[\'beta\'] = ', p['beta'], ' beta_check = ', beta_check

	p['sig_z'] = (p['beta'] * c * p['bunch_length'])/4.
	bunch_length_check	= p['sig_z'] / c / bunch.getSyncParticle().beta()*4
	if ((bunch_length_check - p['blength'])/bunch_length_check) < 0.05: print '\n\tCreate_Bunch:: bunch length check failed.\n\t\tp[\'blength\'] = ', p['blength'], ' bunch_length_check = ', bunch_length_check
	# ~ if bunch_length_check is not p['blength']: print '\n\tCreate_Bunch:: bunch length check failed.\n\t\tp[\'blength\'] = ', p['blength'], ' bunch_length_check = ', bunch_length_check
	
	p['macrosize']			= p['intensity']/float(p['n_macroparticles'])

	p['bunch_save_name'] = 'PyORBIT_'+DistType+'_Bunch_'+TwissType+'_Twiss_Nmp_' + str(p['n_macroparticles']) + '_' + p['bunch_label']

	if TwissType is 'Lattice':
		
		if DistType is 'Gaussian':
			print '\n\tCreate_Bunch::generate_initial_distribution_3DGaussian on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution_3DGaussian(p, Lattice, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		elif DistType is 'Joho':
			print '\n\tCreate_Bunch::generate_initial_distribution on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution(p, Lattice, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		elif DistType is 'Tomo':
			print '\n\tCreate_Bunch::generate_initial_distribution_from_tomo on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution_from_tomo(p,  Lattice, 1, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		else:
			print '\n\tCreate_Bunch::Error: Distribution Type not specified. Options are \'Gaussian\', \'Joho\', and \'Tomo\'. Exiting.'
			exit(0)
			
	elif TwissType is 'Manual':
		
		if TwissDict is None:
			print '\n\tCreate_Bunch::Error: Selected Manual Twiss, but no TwissDict set. Exiting.'
			exit(0)

		if DistType is 'Gaussian':
			print '\n\tCreate_Bunch::generate_initial_distribution_3DGaussian_manual_Twiss on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution_3DGaussian_manual_Twiss(p, TwissDict, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		elif DistType is 'Joho':
			print '\n\tCreate_Bunch::generate_initial_distribution_manual_Twiss on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution_manual_Twiss(p, TwissDict, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		elif DistType is 'Tomo':
			print '\n\tCreate_Bunch::generate_initial_distribution_from_tomo_manual_Twiss on MPI process: ', rank
			Particle_distribution_file = generate_initial_distribution_from_tomo_manual_Twiss(p, TwissDict, 1, output_file=('Distributions/'+p['bunch_save_name']+'.in'), summary_file=('Distributions/'+p['bunch_save_name']+'_summary.txt'))
		else:
			print '\n\tCreate_Bunch::Error: Distribution Type not specified. Options are \'Gaussian\', \'Joho\', and \'Tomo\''
			exit(0)

	print '\n\t\tCreate_Bunch: Output bunch parameters on MPI process: ', rank
	for i in p:
		print '\t', i, '\t = \t', p[i]

	print '\n\t\tbunch_orbit_to_pyorbit on MPI process: ', rank
	bunch_orbit_to_pyorbit((Lattice.getLength()/Lattice.nHarm), p['kin_Energy'], Particle_distribution_file, bunch, p['n_macroparticles'] + 1) #read in only first N_mp particles.

	# Add Macrosize to bunch
	#-----------------------------------------------------------------------
	bunch.addPartAttr("macrosize")
	map(lambda i: bunch.partAttrValue("macrosize", i, 0, p['macrosize']), range(bunch.getSize()))
	ParticleIdNumber().addParticleIdNumbers(bunch) # Give them unique number IDs

	# Dump and save as Matfile
	#-----------------------------------------------------------------------
	print '\n\t\tSave bunch in ',p['bunch_save_name'],'.mat on MPI process: ', rank
	savename = str('Bunches/'+p['bunch_save_name'])
	saveBunchAsMatfile(bunch, savename)

	return bunch

# Function to peform bunch twiss analysis, extract bunch info, and write to output file
def Analyse_Bunch(bunch, p, rank=0):

	print '\n\t\tAnalyse_Bunch:: Analysing bunch:', p['bunch_label'] ,' of type ', p['DistType'] ,' on MPI process: ', rank
	bunchtwissanalysis = BunchTwissAnalysis() # Prepare the bunch analysis class
	get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
	get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
	get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())

	output_file = str('Distributions/'+p['bunch_save_name']+'_analysis.mat')
	output = Output_dictionary()
	output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())
	output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
	output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
	output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
	output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
	output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
	output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
	output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
	output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
	output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
	output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
	output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
	output.addParameter('beta_x', lambda: bunchtwissanalysis.getBeta(0))
	output.addParameter('beta_y', lambda: bunchtwissanalysis.getBeta(1))
	output.addParameter('alpha_x', lambda: bunchtwissanalysis.getAlpha(0))
	output.addParameter('alpha_y', lambda: bunchtwissanalysis.getAlpha(1))
	output.addParameter('D_x', lambda: bunchtwissanalysis.getDispersion(0))
	output.addParameter('D_y', lambda: bunchtwissanalysis.getDispersion(1))
	output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
	output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))
	output.addParameter('mu_x', lambda: GetBunchMus(bunch)[0])
	output.addParameter('mu_y', lambda: GetBunchMus(bunch)[1])
	output.addParameter('eff_beta_x', lambda: bunchtwissanalysis.getEffectiveBeta(0))
	output.addParameter('eff_beta_y', lambda: bunchtwissanalysis.getEffectiveBeta(1))
	output.addParameter('eff_epsn_x', lambda: bunchtwissanalysis.getEffectiveEmittance(0))
	output.addParameter('eff_epsn_y', lambda: bunchtwissanalysis.getEffectiveEmittance(1))
	output.addParameter('eff_alpha_x', lambda: bunchtwissanalysis.getEffectiveAlpha(0))
	output.addParameter('eff_alpha_y', lambda: bunchtwissanalysis.getEffectiveAlpha(1))

	bunchtwissanalysis.analyzeBunch(bunch)

	output.update()
	output.save_to_matfile(output_file)
	print '\n\t\tAnalyse_Bunch:: Analysis Complete, saved to file :', output_file ,' on MPI process: ', rank

	return output_file

def Read_Bunch_Analysis_File(filename):
	f = filename
	p = dict()
	sio.loadmat(f, mdict=p)
	print '\n\tRead_Bunch_Analysis_File::Open output data from ', filename
	return p

def Compare_Parameter(b, p, n1, n2, tol):
	# ~ print '\n\t\t Compare_Parameter:: n1 = ', n1
	# ~ print '\n\t\t Compare_Parameter:: n2 = ', n2
	# ~ print '\n\t\t Compare_Parameter:: b[n1][0][0] = ', b[n1][0][0]
	# ~ print '\n\t\t Compare_Parameter:: p[n2] = ', p[str(n2)]

	if (b[n1][0][0] - p[n2])/b[n1][0][0] < tol: 
		print '\n\t\tCompare_Parameter:: ', n1, '=', b[n1][0][0] ,'with ', n2 , '=' ,p[n2], ' exceeds tolerance of ', (tol*100), '\%'
	return
	
def Compare_Parameter2(b, p, n1, n2, tol):
	# ~ print '\n\t\t Compare_Parameter:: n1 = ', n1
	# ~ print '\n\t\t Compare_Parameter:: n2 = ', n2
	# ~ print '\n\t\t Compare_Parameter:: b[n1][0][0] = ', b[n1][0][0]
	# ~ print '\n\t\t Compare_Parameter:: p[n2] = ', p[str(n2)]
	print '\t\t\tCompare_Parameter:: \t\t', n1, '=', b[n1][0][0] ,'\twith\t\t', n2 , '=' , p[n2], '\t\tDifference =\t\t', 100*(b[n1][0][0] - p[n2])/b[n1][0][0], '%'
	return

# Function to compare parameters and analysed parameters from Analyse_Bunch
def Compare_Parameters(p, a, tolerance=0.05):
	
	# Read analysis file
	b = Read_Bunch_Analysis_File(a)
	
	print 'b[\'D_x\'][0][0]=', b['D_x'][0][0]

	# Iterate over parameters/outputs and compare
	Compare_Parameter2(b, p, 'D_x', 'etax0', tolerance)
	Compare_Parameter2(b, p, 'D_y', 'etay0', tolerance)
	Compare_Parameter2(b, p, 'beta_x', 'betax0', tolerance)
	Compare_Parameter2(b, p, 'beta_y', 'betay0', tolerance)
	Compare_Parameter2(b, p, 'alpha_x', 'alphax0', tolerance)
	Compare_Parameter2(b, p, 'alpha_y', 'alphay0', tolerance)
	Compare_Parameter2(b, p, 'n_mp', 'n_macroparticles', 0.001)
	Compare_Parameter2(b, p, 'gamma', 'gamma', 0.01)
	Compare_Parameter2(b, p, 'bunchlength', 'bunch_length', 0.01)
	Compare_Parameter2(b, p, 'epsn_x', 'epsn_x', 0.01)
	Compare_Parameter2(b, p, 'epsn_y', 'epsn_y', 0.01)

	return

########################################################################
#################			SIMULATION START			################
########################################################################
print '\n\n\tStart simulation main on MPI process: ', rank, '\n'

print '\n\t\tImport simulation parameters from input file on MPI process: ', rank
from simulation_parameters import parameters as p

for i in p:
	print '\t', i, '\t = \t', p[i]

# Create folder structure
#-----------------------------------------------------------------------
print '\n\t\tmkdir on MPI process: ', rank
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('Distributions')
mpi_mkdir_p('Bunches')

# Generate PTC RF table
#-----------------------------------------------------------------------
print '\n\t\tCreate RF file on MPI process: ', rank
from lib.write_ptc_table import write_RFtable
from simulation_parameters import RFparameters as RF 
write_RFtable('RF_table.ptc', *[RF[k] for k in ['harmonic_factors','time','Ekin_GeV','voltage_MV','phase']])

# Initialize a Teapot-Style PTC lattice
#-----------------------------------------------------------------------
print '\n\t\tRead PTC flat file: ',p['flat_file'],' on MPI process: ', rank
PTC_File = p['flat_file']
Lattice = PTC_Lattice("PS")
Lattice.readPTC(PTC_File)

print '\n\t\tRead PTC files on MPI process: ', rank
CheckAndReadPTCFile('PTC/fringe.ptc')
CheckAndReadPTCFile('PTC/time.ptc')

if p['lattice_version'] is 'Original':		CheckAndReadPTCFile('PTC/ramp_cavities.ptc')
elif p['lattice_version'] is 'Optimised':	CheckAndReadPTCFile('PTC/ramp_cavities_optimised.ptc')
else: print '\n\tERROR: p[\'lattice_version\'] = ', p['lattice_version'] ,' not recognised, options are \'Original\' and \'Optimised\''

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

########################################################################
#################			BUNCH PARAMETERS			################
########################################################################

p['gamma']				= 2.49253731343
p['intensity']			= 72.5E+10
p['bunch_length']		= 140e-9
p['blength']			= 140e-9
p['epsn_x']				= 1E-6
p['epsn_y']				= 1.2E-6
p['dpp_rms']			= 8.7e-04
p['TransverseCut']		= 5
p['rf_voltage']			= 0.0212942055190595723
p['circumference']		= 2*np.pi*100
p['phi_s']				= 0
p['bunch_length']		= 140e-9

# For Joho Distn
p['LongitudinalJohoParameter'] = 1.2
p['LongitudinalCut'] = 2.4

########################################################################
#################			TWISS DICTIONARY			################
########################################################################

TwissDict = dict()	# Dictionary used to set initial TWISS parameters
twissdict = None	# switch used to select vertical or horizontal scan data

if float(p['tunex'])/100. is not 6.21:
	twissdict = 'H'
elif float(p['tuney'])/100. is not 6.24:
	twissdict = 'V'
else: 
	print '\n\t Tune set to nominal value of ', (float(p['tunex'])/100.), ',', (float(p['tunex'])/100.)
	if p['lattice_start'][3] is 'H':
		twissdict = 'H'
	elif p['lattice_start'][3] is 'V':
		twissdict = 'V'
	else: 
		print '\n\t Error:: Start position not recognised, options are BWSH65 or BWSV64'
		exit(0)

if twissdict is 'H':
	Qx = [6.07, 6.08, 6.09, 6.1, 6.11, 6.12, 6.13, 6.14, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21]

# Don't have correct manual values for H - use PTC values (not correct - using old lattice values for new optimised lattice etc
# Fix after first mini scan
# HR 19.10.19

	beta_x = [10.0407925, 10.42665652, 10.74232159, 11.00546113, 11.22835049, 11.41974585, 11.58604376, 11.73201709, 11.86129497, 11.97668378, 12.08038684, 12.17415812, 12.25941171, 12.33730138, 12.40877927]
	beta_y = [22.97409451, 22.91670201, 22.85938091, 22.80228084, 22.74549605, 22.68908711, 22.63309333, 22.57754015, 22.52244367, 22.46781359, 22.41365514, 22.35997039, 22.30675909, 22.25401937, 22.20174814]
	# ~ beta_x_PTC = [10.0407925, 10.42665652, 10.74232159, 11.00546113, 11.22835049, 11.41974585, 11.58604376, 11.73201709, 11.86129497, 11.97668378, 12.08038684, 12.17415812, 12.25941171, 12.33730138, 12.40877927]
	# ~ beta_y_PTC = [22.97409451, 22.91670201, 22.85938091, 22.80228084, 22.74549605, 22.68908711, 22.63309333, 22.57754015, 22.52244367, 22.46781359, 22.41365514, 22.35997039, 22.30675909, 22.25401937, 22.20174814]
	# ~ beta_x = [8.394003463991467, 8.708286729183827, 9.032784558742966, 9.400138802204964, 9.81505304121718, 10.259924073989835, 10.672865484475587, 11.046993746565482, 11.359511019253874, 11.618268536499134, 11.831317088827735, 12.011757599553201, 12.165142064083902, 12.298633658877456, 12.417464319120704]
	# ~ beta_y = [24.25812972124533, 24.16839768942765, 24.073361910869952, 23.964643713301978, 23.85424051924806, 23.7272113160871, 23.606226646832162, 23.47912207724786, 23.3574957557102, 23.2453856341838, 23.138855945291812, 23.03831720890891, 22.94285747960047, 22.85183321206868, 22.763701472111915]

	D_x_PTC = [3.335109665, 3.249546421, 3.173842314, 3.106596461, 3.046447635, 2.992216686, 2.942923052, 2.897763213, 2.856081019, 2.817339676, 2.781098089, 2.74699172, 2.714717411, 2.68402141, 2.654689951]
	D_x = [3.4105513098362543, 3.283853441633427, 3.1663951421572083, 3.064909264300527, 2.9802558144715183, 2.907693019575976, 2.8456978265552255, 2.790859909624802, 2.7415926975644336, 2.6965805277242407, 2.6556961455446224, 2.617993986744784, 2.5833867377832256, 2.550996227198573, 2.5202485644118187]
	D_y = [3.894379185517653e-08, -1.6034580020930179e-06, -1.115807284249341e-06, -2.1091450414330707e-06, -4.772279146403814e-06, -7.177485201322374e-06, 7.1281520163150744e-06, -4.6800680492703475e-08, 8.644617839412422e-06, 8.666252155001297e-06, 4.960717057444899e-06, 2.4192006032738668e-05, -1.1739619500046247e-07, 7.808612298156221e-08, 2.1546832524519833e-06]

	alpha_x = [0.09494016969245705, 0.07946471956476579, 0.06621276615003152, 0.05453415926070188, 0.044041249911320454, 0.03444655377771064, 0.02646927739735752, 0.019831167246733963, 0.014717326140394312, 0.010812784493718515, 0.007887953301228174, 0.08087020614528506, 0.0038908272377804125, 0.0025206572457230746, 0.0013920706296361995]
	alpha_x_eff = [-0.0903977145803502, -0.08582703877146283, -0.07626111163734424, -0.06464854602071507, -0.05258770122560416, -0.04177356319519536, -0.03275767440563894, -0.025545783701886052, -0.019676546552782857, -0.014827600526220296, -0.010709742897535353, -0.06409810003694091, -0.004319489877500304, -0.001824179752989762, 0.00031512852672848713]
	alpha_y = [0.028911008609906703, 0.0298091467777196, 0.030729883940520928, 0.03176127978144148, 0.03279985972617165, 0.03398321534703865, 0.03512308463163282, 0.03630292440921215, 0.03743631684202831, 0.03848237473612636, 0.03945838121386767, -0.005633718017399714, 0.041238677981921885, 0.04205208697115897, 0.04283043541275512]
	alpha_y_eff = [0.02891101864671029, 0.02980914993748392, 0.030729892309368283, 0.031761289242922874, 0.03279987034323106, 0.03398322381947462, 0.035123091193636595, 0.03630293244096726, 0.037436323423418065, 0.038482376128553625, 0.039458384594416154, -0.005633686693042419, 0.04123868178535598, 0.04205209046459274, 0.04283043864772137]

	index = Qx.index(float(p['tunex'])/100)

	TwissDict['alpha_x'] 			= alpha_x[index]
	TwissDict['alpha_y'] 			= alpha_y[index]
	TwissDict['beta_x'] 			= beta_x[index]
	TwissDict['beta_y'] 			= beta_y[index]
	TwissDict['D_x'] 				= D_x[index]
	TwissDict['D_y'] 				= D_y[index]
	TwissDict['D_xp'] 				= 0.
	TwissDict['D_yp'] 				= 0.
	TwissDict['x0'] 				= 0.
	TwissDict['xp0'] 				= 0.
	TwissDict['y0'] 				= 0.
	TwissDict['yp0'] 				= 0.

elif twissdict is 'V':
	Qy = [6.1, 6.11, 6.12, 6.13, 6.14, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21, 6.22, 6.23, 6.24]

	beta_x_ptc = [11.83456792, 11.87596148, 11.91741726, 11.95887931, 12.00030736, 12.04167193, 12.08295116, 12.12412865, 12.16519202, 12.20613183, 12.24694093, 12.28761386, 12.32814649, 12.36853572, 12.40877927]
	beta_x = [11.497769168690875, 11.561219729086927, 11.621678373194436, 11.686151694872475, 11.752090514751318, 11.815069726257521, 11.881710579253923, 11.949973709644429, 12.019244092286321, 12.088000160376582, 12.155714548092652, 12.221668817192917, 12.286522995992001, 12.351792373251067, 12.41733027437319]
	beta_x_eff = [15.70727832774324, 15.746782430840286, 15.786312009995264, 15.827611157406391, 15.865523318272942, 15.900303674092575, 15.934927000400245, 15.970255546048815, 16.00642316437317, 16.044394910006403, 16.08771768319707, 16.171233397798897, 16.287840619373497, 16.365446219985902, 16.388017386590324]
	beta_y_ptc = [24.92373871, 24.42660692, 24.02924164, 23.70631133, 23.44016924, 23.21815984, 23.03097747, 22.87162883, 22.7347561, 22.61618278, 22.51260139, 22.42135386, 22.34027393, 22.26757179, 22.20174814]
	beta_y = [28.75941629371245, 27.901539601878184, 27.213290322184257, 26.55273573875031, 25.951088509809498, 25.46986237056283, 24.988794222684113, 24.545888785547756, 24.150861736726874, 23.81407818999024, 23.529804455256162, 23.280644484660698, 23.067685323867895, 22.898961752454188, 22.764817219063925]
	beta_y_eff = [28.75941363876187, 27.901537495704986, 27.21328851050792, 26.552734005497832, 25.951087908101496, 25.469861065350035, 24.988792628258622, 24.545887021172206, 24.150860438788662, 23.814077295064806, 23.529803570633504, 23.280643730316875, 23.067684825182425, 22.89896132742581, 22.764816922948334]

	alpha_x_ptc = [0.011513, 0.01069, 0.009857, 0.009017, 0.00817, 0.007316, 0.006457, 0.005594, 0.004726, 0.003854, 0.002979, 0.002101, 0.001221, 0.000338, -0.000546]
	alpha_x = [0.020149363762314573, 0.018841738892805736, 0.017628458124003044, 0.016317893976010907, 0.014987158717971134, 0.013721509703739742, 0.012370721644305303, 0.010976584559329877, 0.009552467712758283, 0.008137364733604945, 0.00674769784387647, 0.00540321374339072, 0.004084552073385008, 0.0027449656310449765, 0.001399858908938113]
	alpha_x_eff = [0.01263857994342128, 0.011787402752768992, 0.010993119054166795, 0.010138160374765592, 0.009268803518251613, 0.008440528874409346, 0.007565913699962478, 0.00665668167725413, 0.005727886071697952, 0.0047986057589267895, 0.0038721576771630377, 0.0029754572927884995, 0.00209091394149012, 0.0011996244965600173, 0.00030772719575205697]
	alpha_y_ptc = [-0.015862, -0.00566, 0.002752, 0.009781, 0.015724, 0.020797, 0.025166, 0.028959, 0.032276, 0.035195, 0.03778, 0.040084, 0.04215, 0.044012, 0.045702]
	alpha_y = [-0.0759981299041554, -0.06026690948643802, -0.04710360569302436, -0.034364222714159766, -0.0226876810922163, -0.0128778777253371, -0.0031960530794603887, 0.00573900436130687, 0.013766673300182369, 0.020700217072685663, 0.026624706918980934, 0.031779608236245774, 0.03618275788076127, 0.03980850948662125, 0.04282193036546169]
	alpha_y_eff = [-0.07599819813583442, -0.06026687716107299, -0.04710358259720198, -0.03436419121508966, -0.022687704211584163, -0.012877859092497113, -0.0031960257774123163, 0.005739023549823524, 0.013766695978317036, 0.02070022530014002, 0.026624716144330474, 0.03177961640012476, 0.036182764652471, 0.039808513074497745, 0.04282193193366175]

	D_x_ptc = [2.858624854, 2.844314013, 2.829919074, 2.815459752, 2.800950597, 2.786402598, 2.771824228, 2.757222142, 2.742601663, 2.727967112, 2.71332205, 2.698669456, 2.684011849, 2.66935139, 2.654689951]
	D_x = [2.7232597364718965, 2.711617531257042, 2.6995782883372725, 2.6866659088601867, 2.673333040270362, 2.660009904023271, 2.645587060435802, 2.6306159716576585, 2.6153104916827656, 2.5997777190545786, 2.584122745980839, 2.5682844312570325, 2.552395926698185, 2.5364870147593908, 2.5203979092578246]
	D_y_ptc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	D_y = [0.00032016495885226815, -0.0002189058945793827, 2.7764495806106146e-05, -0.00024143037602403352, 0.0009865662335166862, 0.00015608160748784814, -0.0002853837009812899, -4.574488900503671e-06, -6.670454167668475e-05, -9.192803073975801e-05, -5.7217673705390416e-05, 3.9571270250276e-05, -1.7313614618598883e-05, -2.407762236190365e-06, -2.2386320787418418e-06]

	index = Qy.index(float(p['tuney'])/100)

	TwissDict['alpha_x'] 			= alpha_x[index]
	TwissDict['alpha_y'] 			= alpha_y[index]
	TwissDict['beta_x'] 			= beta_x[index]
	TwissDict['beta_y'] 			= beta_y[index]
	TwissDict['D_x'] 				= D_x[index]
	TwissDict['D_y'] 				= D_y[index]
	TwissDict['D_xp'] 				= 0.
	TwissDict['D_yp'] 				= 0.
	TwissDict['x0'] 				= 0.
	TwissDict['xp0'] 				= 0.
	TwissDict['y0'] 				= 0.
	TwissDict['yp0'] 				= 0.


TwissDict['gamma_transition']	= Lattice.gammaT
TwissDict['circumference']		= Lattice.getLength()
TwissDict['length'] 			= Lattice.getLength()/Lattice.nHarm

#	First we create the bunch using the Create_Bunch function which 
#	creates the _summary.txt, .in, and bunch .mat file.
#	Next we analyse the bunch using the Analyse_Bunch function which
#	creates the output _analysis.mat file.
#	Finally we load output _analysis.mat file, and compare values to 
#	the input parameters.

gaussian_bunch = Create_Bunch(Lattice, p, TwissDict=None, label=p['bunch_label'], DistType = 'Gaussian', TwissType = 'Lattice', rank=rank)
gaussian_analysis = Analyse_Bunch(gaussian_bunch, p)
Compare_Parameters(p, gaussian_analysis)

joho_bunch = Create_Bunch(Lattice, p, TwissDict=None, label=p['bunch_label'], DistType = 'Joho', TwissType = 'Lattice', rank=rank)
joho_analysis = Analyse_Bunch(joho_bunch, p)
Compare_Parameters(p, joho_analysis)

tomo_bunch = Create_Bunch(Lattice, p, TwissDict=None, label=p['bunch_label'], DistType = 'Tomo', TwissType = 'Lattice', rank=rank)
tomo_analysis = Analyse_Bunch(tomo_bunch, p)
Compare_Parameters(p, tomo_analysis)

gaussian_mt_bunch = Create_Bunch(Lattice, p, TwissDict=TwissDict, label=p['bunch_label'], DistType = 'Gaussian', TwissType = 'Manual', rank=rank)
gaussian_mt_analysis = Analyse_Bunch(gaussian_mt_bunch, p)
Compare_Parameters(p, gaussian_mt_analysis)

joho_mt_bunch = Create_Bunch(Lattice, p, TwissDict=TwissDict, label=p['bunch_label'], DistType = 'Joho', TwissType = 'Manual', rank=rank)
joho_mt_analysis = Analyse_Bunch(joho_mt_bunch, p)
Compare_Parameters(p, joho_mt_analysis)

tomo_mt_bunch = Create_Bunch(Lattice, p, TwissDict=TwissDict, label=p['bunch_label'], DistType = 'Tomo', TwissType = 'Manual', rank=rank)
tomo_mt_analysis = Analyse_Bunch(tomo_mt_bunch, p)
Compare_Parameters(p, tomo_mt_analysis)

print '\n\n\tFinish simulation main on MPI process: ', rank, '\n'
