import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09
from reax_utils import *

def calculate_error(dataset):
	write_reax_file(dataset) #all LAMMPS job use same reax file and same data files
	for elements,systems in dataset.by_elements.iteritems():
		commands = ('''clear
units real
atom_style full
pair_style reax/c NULL mincap 10
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary f f f
read_data	'''+systems[0].name+'''.data

neigh_modify every 1 delay 0 check no
compute atom_pe all pe/atom
compute		test_pe all reduce sum c_atom_pe
#compute atom_f all property/atom fx fy fz
#variable b equal c_atom_f[1][1]
thermo_style custom pe c_test_pe
pair_coeff * * '''+dataset.name+'''.reax Pb Cl
fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c''').splitlines()

		for line in commands:
			dataset.lmp.command(line)
		
		for s in systems:
			for i,a in enumerate(s.atoms):
				dataset.lmp.command('set atom %d x %f y %f z %f' % (i+1, a.x, a.y, a.z) )
				dataset.lmp.command('set atom %d charge 0.0' % (i+1) )
			
			dataset.lmp.command('run 1')
			
			lammps_energies_by_atom = dataset.lmp.extract_compute('atom_pe',1,1) #http://lammps.sandia.gov/doc/Section_python.html
			s.lammps_energy = sum( [lammps_energies_by_atom[i] for i in range(0,len(s.atoms)) ] )
			lammps_forces = dataset.lmp.gather_atoms('f',1,3)
			for i,a in enumerate(s.atoms):
				a.lfx, a.lfy, a.lfz = lammps_forces[(i*3)+0], lammps_forces[(i*3)+1], lammps_forces[(i*3)+2]

	#calculate energy error
	relative_energy_error, absolute_energy_error = 0.0, 0.0
	relative_force_error, absolute_force_error = 0.0, 0.0
	real_E,test_E,real_F,test_F = [],[],[],[]
	for elements,systems in dataset.by_elements.iteritems():
		baseline_energy = systems[0].lammps_energy #should be in order of increasing energy
		for s in systems:
			s.lammps_energy -= baseline_energy
			try:
				relative_energy_error += ( (s.lammps_energy-s.energy)/(s.energy+1.0) )**2
				absolute_energy_error += (s.lammps_energy-s.energy)**2
			except OverflowError: return 1e10
			#print s.name, s.energy, s.lammps_energy
			real_E.append(s.energy); test_E.append(s.lammps_energy)

			for i,a in enumerate(s.atoms):
				fx, fy, fz = a.lfx, a.lfy, a.lfz
				real_force_squared = a.fx**2 + a.fy**2 + a.fz**2
				try:
					lammps_force_squared = fx**2 + fy**2 + fz**2
					real_F.append(real_force_squared**0.5); test_F.append(lammps_force_squared**0.5)
					relative_force_error += ((fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2) / (real_force_squared + 20.0**2)
					absolute_force_error += (fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2
				except OverflowError:
					return 1e10
	
	plot_forces = False
	plot_energies = False
	if plot_forces or plot_energies: #plot energies
		import matplotlib.pyplot as plt
		if plot_energies:
			plt.plot(real_E, label='HSE06/TZVP energies')
			plt.plot(test_E, label='REAX energies')
			plt.xlabel('Runs (sorted by energy)')
			plt.ylabel('Energy (kcal/mol)')
		elif plot_forces:
			real_F,test_F = [list(x) for x in zip(*sorted(zip(real_F,test_F)))] #sort forces
			plt.plot(real_F, label='HSE06/TZVP forces')
			plt.plot(test_F, label='REAX forces')
			plt.xlabel('Atoms (sorted by force)')
			plt.ylabel('Force (kcal/mol/Angstrom)')
		plt.legend(bbox_to_anchor=(0.6, 0.9), bbox_transform=plt.gcf().transFigure)
		plt.title('%s REAX fit' % dataset.name)
		plt.show()
		exit()
	
	n_atoms = sum( [len(s.atoms) for s in dataset.systems] )
	relative_force_error = math.sqrt( relative_force_error/n_atoms )
	absolute_force_error = math.sqrt( absolute_force_error/n_atoms )
	relative_energy_error = math.sqrt( relative_energy_error/len(dataset.systems) )
	absolute_energy_error = math.sqrt( absolute_energy_error/len(dataset.systems) )

	error = relative_energy_error + relative_force_error
	
	if math.isnan(error):
		return 1e10
	else:
		return error

def pack_params(dataset):
	params, names = [], []
	
	for atom,include in zip(dataset.reax_params.atom_types, dataset.reax_includes.atom_types):
		names_list = dataset.reax_params.atom_types_names
		atom_name = atom[0][0]
		for line in range(4):
			for i,b in enumerate(include[line]):
				if b and type(b)!=str:
					params.append(atom[line][i])
					names.append(atom_name + '.' + names_list[line][i])

	for bond,include in zip(dataset.reax_params.bonds, dataset.reax_includes.bonds):
		names_list = dataset.reax_params.bond_types_names
		bondname='tbp(' + bond[0][0] + ',' + bond[0][1] + ').'
		for line in range(2):
			for i,b in enumerate(include[line]):
				if b and type(b)!=str:
					params.append(bond[line][i])
					names.append(bondname+names_list[line][i])
	
	for offdiag,include in zip(dataset.reax_params.offdiags, dataset.reax_includes.offdiags):
		names_list=dataset.reax_params.offdiags_names
		offdname='offd('+offdiag[0]+','+offdiag[1]+').'
		for i,b in enumerate(include):
			if b and type(b)!=str:
				params.append(offdiag[i])
				names.append(offdname+names_list[i])

	for thbp,include in zip(dataset.reax_params.thbps,dataset.reax_includes.thbps):
		names_list = dataset.reax_params.three_body_names
		thbpname='thbp('+thbp[0] +',' + thbp[1]+',' + thbp[2]+').'
		for i,b in enumerate(include):
			if b and type(b)!=str:
				params.append(thbp[i])
				names.append(thbpname + names_list[i])

	if len(params)!=len(names):
		print 'There are %d parameters, but %d names!' % (len(params), len(names))
		raise SystemExit
	
	return params, names

def unpack_params(params, dataset):
	p = 0
	
	for atom,include in zip(dataset.reax_params.atom_types, dataset.reax_includes.atom_types):
		for line in range(4):
			for i,b in enumerate(include[line]):
				if b and type(b)!=str:
					atom[line][i] = params[p]
					p += 1

	for bond,include in zip(dataset.reax_params.bonds, dataset.reax_includes.bonds):
		for line in range(2):
			for i,b in enumerate(include[line]):
				if b and type(b)!=str:
					bond[line][i] = params[p]
					p += 1

	for offdiag,include in zip(dataset.reax_params.offdiags, dataset.reax_includes.offdiags):
		for i,b in enumerate(include):
			if b and type(b)!=str:
				offdiag[i] = params[p]
				p += 1

	for thbp,include in zip(dataset.reax_params.thbps,dataset.reax_includes.thbps):
		for i,b in enumerate(include):
			if b and type(b)!=str:
				thbp[i] = params[p]
				p += 1

def run(run_name, other_run_names=[],restart=False):
	Cl_ = 66
	H_ = 54
	N_ = 53
	Pb_ = 111

	Pb = 907
	Cl = 838

	extra = {
		Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=0.0, vdw_r=3.0),
		Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.0, vdw_r=3.0),
	}

	dataset = utils.Struct(name=run_name, systems=[])
	filenames = []
	for root, dirs, file_list in os.walk('gaussian'):
		for ff in file_list:
			if ff.endswith('.log'):
				name = ff[:-4]
				filenames.append(name)

	for name in filenames:
				if not name.startswith('PbCl2'): continue
				if not '_vac' in name: continue
				result = g09.parse_atoms(name, check_convergence=True)
				if not result: continue #don't use the log file if not converged
				energy, atoms = result
				#if len(atoms)!=3: continue
				if any( [a!=b and utils.dist(a,b)<2.0 for a in atoms for b in atoms] ): continue
				system = utils.System(box_size=[40, 40, 40], name=name)
				total = utils.Molecule('gaussian/'+name, extra_parameters=extra, check_charges=False)
				system.energy = energy*627.509 #convert energy from Hartree to kcal/mol
				system.element_string = ' '.join( [a.element for a in total.atoms] )
				element_string_2 = ' '.join( [a.element for a in atoms] )
				if element_string_2 != system.element_string:
					#print 'Inconsistent elements in cml vs log:', name
					continue
				#print 'Read', system.element_string, 'from', name
				for i,a in enumerate(total.atoms):
					b = atoms[i]
					if a.element != b.element:
						print 'Inconsistent elements in cml vs log:', name
						exit()
					a.element=b.element
					a.x, a.y, a.z = b.x, b.y, b.z
					a.fx, a.fy, a.fz = [f*1185.8113 for f in (b.fx, b.fy, b.fz)] # convert forces from Hartree/Bohr to kcal/mol / Angstrom
				system.add(total)
				dataset.systems.append(system)
	#group systems by .element_string
	dataset.by_elements = {}
	for system in dataset.systems:
		if system.element_string not in dataset.by_elements:
			dataset.by_elements[ system.element_string ] = []
		dataset.by_elements[ system.element_string ].append(system)
	#sort molecules of same type by energy, set baseline energy as zero
	for element_string, systems in dataset.by_elements.iteritems():
		systems.sort(key=lambda s:s.energy)
		baseline_energy = systems[0].energy
		for s in systems:
			s.energy -= baseline_energy #baseline energy = 0.0

	if True: #hide screen
		dataset.lmp = lammps('',['-log','none','-screen','none'])
	else: #show screen
		dataset.lmp = lammps('',['-log','none'])
	
	os.chdir('lammps')
	for s in dataset.systems:
		files.write_lammps_data(s)

	default_input='../input.reax'

	if restart:
		best_error = 1e10
		best_run = run_name
		input_file=default_input
		for name in [run_name]+other_run_names:
			if os.path.isfile(name+'_best.reax'):
				error = float(open(name+'_best.reax').readline().split()[1])
				if error<best_error:
					best_error = error
					best_run = name
					input_file=best_run+'_best.reax'
	else:
		input_file=default_input
	commands = ('''units real
atom_style full
pair_style reax/c NULL mincap 10
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary f f f
read_data	'''+dataset.systems[0].name+'''.data

neigh_modify every 1 delay 0 check no
compute atom_pe all pe/atom
compute		test_pe all reduce sum c_atom_pe
thermo_style custom pe c_test_pe
pair_coeff * * '''+input_file+''' Pb Cl
fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c''').splitlines()

	for line in commands:
		dataset.lmp.command(line)

	dataset.reax_params = read_reax_file(input_file)
	dataset.reax_includes, bounds = read_reax_include_file('../include.reax',dataset.reax_params)

	def calculate_error_from_list(params):
		unpack_params(params, dataset)
		error = calculate_error(dataset)
		return error

	initial_params, names = pack_params(dataset)
	for i in range(len(names)):
		if names[i].endswith('gamma_w'):
			if bounds[i][0]<=.5:
				print names[i] + ' was bounded by ' +  str(bounds[i]) + 'changing to (.501,%f)' % bounds[i][1]
				bounds[i]=(.501,bounds[i][1])
			if bounds[i][1]<=.5:
				print names[i] + ' bounds were unacceptable, they have been changed to (.501,45).'
				bounds[i][1]=(bounds[i][0],45)

	if False: #print params at start
		print "Parametrization on:"
		print '      Variable         Initial         Bounds'
		for x,y,z in zip(names,initial_params,bounds):
			print '%20s' % x, '%9.4f' % y, '(%8.4f,%8.4f)' % (z[0],z[1])

	#get params from other, parallel runs
	dataset.others = [utils.Struct(name=n) for n in other_run_names]
	def check_others(dataset):
		for other in dataset.others:
			other.bond_types, other.angle_types, other.dihedral_types, = [], [], []
			if os.path.exists('../'+other.name+'_best.reax'):
				other.reax_params = [utils.Struct(reax_params=read_reax_file('../'+other.name+'_best.reax')[0]) ]
				other.reax_includes = dataset.reax_includes
				other.list, _ = pack_params(other)
			else:
				other.list = None
		dataset.how_long_since_checked_others = 0
	check_others(dataset)

	def parameter_effect(initial_params=initial_params):
		print "Testing whether included parameters have an apparent effect"
		best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
		for i,p in enumerate(initial_params):
			# print 'Testing '+ names[i]
			testparams = initial_params[:]
			if p!=0:
				testparams[i] = p*.2
			else:
				print 'PARAMETER ' + names[i] + " WAS ZERO"
				testparams[i] = -5
			test1 = utils.Struct(fun=calculate_error_from_list(testparams),x=testparams)
			if p!=0:
				testparams[i] = p*5
			else:
				testparams[i] = 5
			test2 = utils.Struct(fun=calculate_error_from_list(testparams),x=testparams)
			if (best_min.fun == test1.fun and best_min.fun == test2.fun):
				print 'Parameter ' + names[i] +' has no effect'
			else:
				print 'Parameter ' + names[i] +' has an effect' 
		raise SystemExit # Exit once all parameters are tested

	def stochastic(use_gradient=True):
		# parameter_effect(initial_params)
		if use_gradient:
			import numpy
			from scipy.optimize import minimize, fmin_l_bfgs_b
		best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
		print dataset.name, ('starting error: %.4g from ' % best_min.fun ) + input_file
		# raise SystemExit #just print starting error
		def new_param_guess(start, gauss=True):
			dataset.how_long_since_checked_others += 1
			if dataset.how_long_since_checked_others > 10:
				check_others(dataset)
			
			while True: #keep going until new params are generated
				params = []
				for p,b in zip(start, bounds):
					if gauss: #Gaussian random
						new = random.gauss(p, 0.2*(b[1]-b[0]) ) if random.random()<0.3 else p
						#reflect
						if new < b[0]:
							new += (b[0]-new)
						elif new > b[1]:
							new -= (b[1]-new)
						#just set
						if new < b[0]:
							new = b[0]
						elif new > b[1]:
							new = b[1]
					else: #uniform random
						new = b[0] + random.random()*(b[1]-b[0])
					params.append( new )
				if cmp(list(params), list(start)) == 0: #if new param list is the same as old one
					continue
				else:
					#keep different solutions apart
					if dataset.others:
						distances = []
						for other in dataset.others:
							if other.list:
								distances.append( sum([ ( (x-y)/(b[1]-b[0]) )**2 for x,y,b in zip(others.list, params, bounds)])**0.5 / len(params) )
							else:
								distances.append(1e6)
						if all( [ d>0.01*len(params) for d in distances ] ):
							return params
					else:
						return params
	
		def error_gradient(x):
			e0 = calculate_error_from_list(x)
			gradient = []
			for i in range(len(x)):
				if bounds[i]==(0.0,0.0):
					gradient.append(0.0)
					continue
				oldx = x[i]
				sign_oldx = -1 if oldx<0 else 1
				newx = oldx + 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
				if newx>bounds[i][1] or newx<bounds[i][0]:
					newx = oldx - 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
				if newx>bounds[i][1] or newx<bounds[i][0]:
					pass
					#print 'Bounds exceeded:', oldx, newx, bounds[i][0], bounds[i][1]
					#exit()
				dif = newx - oldx
				if dif != 0.0:
					x[i] = newx
					gradient.append( (calculate_error_from_list(x)-e0)/dif )
					x[i] = oldx
				else:
					print 'dif=0.0 in numerical gradient!'
					print newx, oldx, (bounds[i][0], bounds[i][1])
					exit()
			return numpy.array(gradient)

		while True:
			if use_gradient:
				params = new_param_guess(best_min.x, gauss=True)
				start_error = calculate_error_from_list(params)
				#while start_error > 2.0:
				#	params = new_param_guess(best_min.x, gauss=True)
				#	start_error = calculate_error_from_list(params)
				x, fun, stats = fmin_l_bfgs_b(calculate_error_from_list, params, fprime=error_gradient, factr=1e8)#, bounds=bounds)
				guess = utils.Struct(fun=fun,x=x)
			else: #non-gradient optimization
				params = new_param_guess(best_min.x, gauss=True)
				guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
				start_error = guess.fun
			dif_from_best = 100*sum([ ( (x-y)/(b[1]-b[0]) )**2 for x,y,b in zip(guess.x, best_min.x, bounds)])**0.5 / len(best_min.x)
			#print dataset.name, 'error', start_error, guess.fun, best_min.fun, 'Dif = %.1f%%' % dif_from_best
			if guess.fun < best_min.fun:
				best_min = guess
				unpack_params(best_min.x, dataset)
				write_reax_file(dataset,best=True,error = best_min.fun)
				print dataset.name, 'new best error = %.4g' % best_min.fun, 'Dif = %.1f%%' % dif_from_best

		return best_min

	stochastic(True)

from multiprocessing import Process, Queue
from time import sleep
def run_multiple(jobname, N):
	queue = Queue()
	procs=[]
	for i in range(N):
		procs.append(Process(target=run, args=(jobname+str(i), [jobname+str(other) for other in range(N) if other!=i],False)))
		procs[-1].start()
	while True:
		for i,p in enumerate(procs):
			if not p.is_alive():
				procs[i] = Process(target=run, args=(jobname+str(i), [jobname+str(other) for other in range(len(procs)) if other!=i],True))
				procs[i].start()
		sleep(15)

def if_multiprocessing_does_not_work():
	jobname = sys.argv[1]
	if len(sys.argv)==2:
		run(jobname)
	else:
		this_job = int(sys.argv[2])
		n_other_jobs = int(sys.argv[3])
		run(jobname+str(this_job), [jobname+str(other) for other in range(n_other_jobs) if other!=this_job])

# run('test')
run_multiple('test', 4)

