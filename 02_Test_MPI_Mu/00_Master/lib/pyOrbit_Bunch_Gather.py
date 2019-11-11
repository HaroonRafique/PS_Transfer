# 11.11.19: Created HR CERN BE-ABP-HSI
# Script to gather all bunch particles from MPI processes and calculate
# bunch moments and other quantities

import orbit_mpi
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op
import numpy as np
import scipy.io as sio
from bunch import Bunch
from scipy.stats import moment, kurtosis
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class resonance_lines(object):
	
	def __init__(self, Qx_range, Qy_range, orders, periodicity):
		
		if np.std(Qx_range):
			self.Qx_min = np.min(Qx_range)
			self.Qx_max = np.max(Qx_range)
		else:
			self.Qx_min = np.floor(Qx_range)-0.05
			self.Qx_max = np.floor(Qx_range)+1.05
		if np.std(Qy_range):
			self.Qy_min = np.min(Qy_range)
			self.Qy_max = np.max(Qy_range)
		else:
			self.Qy_min = np.floor(Qy_range)-0.05
			self.Qy_max = np.floor(Qy_range)+1.05

		self.periodicity = periodicity

		nx, ny = [], []

		for order in np.nditer(np.array(orders)):
			t = np.array(range(-order, order+1))
			nx.extend(order - np.abs(t))
			ny.extend(t)
		nx = np.array(nx)
		ny = np.array(ny)
	
		cextr = np.array([nx*np.floor(self.Qx_min)+ny*np.floor(self.Qy_min), \
						  nx*np.ceil(self.Qx_max)+ny*np.floor(self.Qy_min), \
						  nx*np.floor(self.Qx_min)+ny*np.ceil(self.Qy_max), \
						  nx*np.ceil(self.Qx_max)+ny*np.ceil(self.Qy_max)], dtype='int')
		cmin = np.min(cextr, axis=0)
		cmax = np.max(cextr, axis=0)
		res_sum = [range(cmin[i], cmax[i]+1) for i in xrange(cextr.shape[1])]
		self.resonance_list = zip(nx, ny, res_sum)
		
	def plot_resonance(self, figure_object = None):
		plt.ion()
		if figure_object:
			fig = figure_object
			plt.figure(fig.number)
		else:
			fig = plt.figure()

		Qx_min = self.Qx_min
		Qx_max = self.Qx_max
		Qy_min = self.Qy_min
		Qy_max = self.Qy_max 
		plt.xlim(Qx_min, Qx_max)
		plt.ylim(Qy_min, Qy_max)
		plt.xlabel('Qx')
		plt.ylabel('Qy')
		for resonance in self.resonance_list:
			nx = resonance[0]
			ny = resonance[1]
			# ~ print 'res ', nx, ' ', ny
			for res_sum in resonance[2]:
				if ny:
					line, = plt.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny])
				else:
					line, = plt.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max])
				if ny%2:
					plt.setp(line, linestyle='--') # for skew resonances
				if res_sum%self.periodicity:
					plt.setp(line, color='b')	# non-systematic resonances
				else:
					plt.setp(line, color='r', linewidth=2.0) # systematic resonances
		plt.draw()
		return fig
		
	def print_resonances(self):
		for resonance in self.resonance_list:
			for res_sum in resonance[2]:
				'''
				print str(resonance[0]).rjust(3), 'Qx ', ("+", "-")[resonance[1]<0], \
					  str(abs(resonance[1])).rjust(2), 'Qy = ', str(res_sum).rjust(3), \
					  '\t', ("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0]
				'''
				print '%s %s%s = %s\t%s'%(str(resonance[0]).rjust(2), ("+", "-")[resonance[1]<0], \
						str(abs(resonance[1])).rjust(2), str(res_sum).rjust(4), \
						("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0])


# Takes the bunch, gathers it from all MPI processes, and returns a 
# dictionary of output parameters.
def BunchGather(bunch, turn, p, plot_footprint=False):

	b = bunch
	verbose = True

	# take the MPI Communicator from bunch: it could be different from MPI_COMM_WORLD
	comm = b.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	size = orbit_mpi.MPI_Comm_size(comm)
	main_rank = 0

	# n_parts_arr - array of size of the number of CPUs, 
	# and have the number of macroparticles on each CPU
	n_parts_arr = [0]*size
	n_parts_arr[rank] = b.getSize()
	n_parts_arr = orbit_mpi.MPI_Allreduce(n_parts_arr,mpi_datatype.MPI_INT,mpi_op.MPI_SUM,comm)	

	if verbose:
		print 'BunchMoments:: bunch_size on MPI Rank: ', rank, ' = ', n_parts_arr[rank]
		print 'BunchMoments:: n_parts_arr on MPI Rank: ', rank, ' = ', n_parts_arr

	mp_array = range(n_parts_arr[rank])
	particles = {}
	particles['x'] = map(b.x, mp_array)
	particles['xp'] = map(b.xp, mp_array)
	particles['y'] = map(b.y, mp_array)
	particles['yp'] = map(b.yp, mp_array)
	particles['z'] = map(b.z, mp_array)
	particles['dE'] = map(b.dE, mp_array)
	phase_space_keys = particles.keys()

	for attribute in b.getPartAttrNames():
		particles[attribute] = [[] for i in range(b.getPartAttrSize(attribute))]
		for j in xrange(b.getPartAttrSize(attribute)):
			particles[attribute][j] += map(lambda i: b.partAttrValue(attribute, i, j), mp_array)

	# This is just in case. Actually, MPI_Barrier command is not necessary.
	orbit_mpi.MPI_Barrier(comm)

	for i_cpu in range(1,size):
		for key in phase_space_keys:
			if(rank == main_rank):
				# get the particle coordinates and attributes
				bunch_size_remote = orbit_mpi.MPI_Recv(mpi_datatype.MPI_INT,i_cpu,222,comm)
				if bunch_size_remote:
					particles[key] += list(np.atleast_1d(orbit_mpi.MPI_Recv(mpi_datatype.MPI_DOUBLE,i_cpu,222,comm)))
			elif(rank == i_cpu):
				# send the coordinate array if there are any particles ...
				bunch_size_local = bunch.getSize()
				orbit_mpi.MPI_Send(bunch_size_local,mpi_datatype.MPI_INT,main_rank,222,comm)
				if bunch_size_local:
					orbit_mpi.MPI_Send(particles[key],mpi_datatype.MPI_DOUBLE,main_rank,222,comm)

	for i_cpu in range(1,size):
		for attribute in b.getPartAttrNames():
			if(rank == main_rank):
				bunch_size_remote = orbit_mpi.MPI_Recv(mpi_datatype.MPI_INT,i_cpu,222,comm)
				if bunch_size_remote:
					# get the particle coordinates and attributes
					for j in xrange(b.getPartAttrSize(attribute)):
						particles[attribute][j] += list(np.atleast_1d(orbit_mpi.MPI_Recv(mpi_datatype.MPI_DOUBLE,i_cpu,222,comm)))
			elif(rank == i_cpu):
				bunch_size_local = bunch.getSize()
				orbit_mpi.MPI_Send(bunch_size_local,mpi_datatype.MPI_INT,main_rank,222,comm)
				if bunch_size_local:
					# send the coordinate array if there are any particles ...
					for j in xrange(b.getPartAttrSize(attribute)):
						orbit_mpi.MPI_Send(particles[attribute][j],mpi_datatype.MPI_DOUBLE,main_rank,222,comm)

	bunchparameters = {'classical_radius': bunch.classicalRadius(), \
					   'charge': bunch.charge(), 
					   'mass': bunch.mass(), \
					   'momentum': bunch.getSyncParticle().momentum(), \
					   'beta': bunch.getSyncParticle().beta(), \
					   'gamma': bunch.getSyncParticle().gamma(), \
					   'time': bunch.getSyncParticle().time()}

########################################################################
#                Plot tune footprint with histograms                   #
########################################################################

	if rank is 0:
		print 'Rank: ', rank
		if turn >=0:
			print 'Turn: ', turn
			if plot_footprint:
				if verbose: print 'BunchGather:: Plot tune footprint on rank', rank
			
				tunex = str(p['tunex'][0] + '.' + p['tunex'][1:])
				tuney = str(p['tuney'][0] + '.' + p['tuney'][1:])
				tunex_sav = str(p['tunex'][0] + 'p' + p['tunex'][1:])
				tuney_sav = str(p['tuney'][0] + 'p' + p['tuney'][1:])
				fontsize=15

				qx = np.array(particles['ParticlePhaseAttributes'][2])
				qy = np.array(particles['ParticlePhaseAttributes'][3])

				qx[np.where(qx>0.5)] -= 1
				qy[np.where((qy>0.6) & (qx<0.25))] -= 1

				print 'resonances'
				resonances = resonance_lines((5.75, 6.25),(5.75, 6.25),(1,2,3,4),10)
				fontsize=17

				f, ax = plt.subplots(1, figsize=(6,6))
				gridspec.GridSpec(3,3)
				#f.subplots_adjust(hspace = 0)	# Horizontal spacing between subplots
				f.subplots_adjust(wspace = 0)	# Vertical spacing between subplots

				my_cmap = plt.cm.jet
				my_cmap.set_under('w',1)

				r = resonances

				print 'title'
				title = str( tunex_sav + ' ' + tuney_sav + ' turn ' + str(turn))

				# First subplot
				print 'plot1'
				plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=1)
				plt.hist(6+qx, bins=1000, range=(r.Qx_min, r.Qx_max)) #, norm=mcolors.PowerNorm(gamma))
				plt.ylabel('Frequency')
				plt.grid(which='both')
				plt.title(title, fontsize=fontsize)

				# Main plot
				print 'plot2'
				plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
				plt.hist2d(6+qx, 6+qy, bins=1000, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]]) #, norm=mcolors.PowerNorm(gamma))
				plt.xlabel(r'Q$_x$')
				plt.ylabel(r'Q$_y$')
				
				print 'plot_resonance'
				resonances.plot_resonance(f)

				# Second subplot
				print 'plot3'
				plt.subplot2grid((3,3), (1,2), colspan=1, rowspan=2)    
				plt.hist(6+qy, bins=1000, range=(r.Qy_min, r.Qy_max), orientation=u'horizontal') #, norm=mcolors.PowerNorm(gamma))
				plt.xlabel('Frequency')
				plt.grid(which='both')

				current_axis = plt.gca()
				#current_axis.axes.get_yaxis().set_visible(False)

				ax.xaxis.label.set_size(fontsize)
				ax.yaxis.label.set_size(fontsize)
				ax.tick_params(labelsize=fontsize)

				plt.tight_layout()
				savename = str('Tune_Footprints/' + tunex_sav + '_' + tuney_sav + '_turn_' + str(turn) + '_hist.png' )
				
				print 'savefig'
				f.savefig(savename, dpi=100)
				plt.close(f)

	# Later add something to cut large amplitude particles to reduce noise for kurtosis calculation
	outputs = {
		'Mu_x' : moment(particles['x'], 2),
		'Mu_xp' : moment(particles['xp'], 2),
		'Mu_y' : moment(particles['y'], 2),
		'Mu_yp' : moment(particles['yp'], 2),
		'Mu_z' : moment(particles['z'], 2),
		'Mu_dE' : moment(particles['dE'], 2),
		'Kurtosis_x' : kurtosis(particles['x'], fisher=True, nan_policy='omit'),
		'Kurtosis_xp' : kurtosis(particles['xp'], fisher=True, nan_policy='omit'),
		'Kurtosis_y' : kurtosis(particles['y'], fisher=True, nan_policy='omit'),
		'Kurtosis_yp' : kurtosis(particles['yp'], fisher=True, nan_policy='omit'),
		'Kurtosis_z' : kurtosis(particles['z'], fisher=True, nan_policy='omit'),
		'Kurtosis_dE' : kurtosis(particles['dE'], fisher=True, nan_policy='omit')
	}


	return outputs
