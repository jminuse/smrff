import os, copy, sys, random, cPickle, math
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

def read_tersoff_file(filename):
	# format of a single entry (one or more lines):
	#   element 1, element 2, element 3, 
	#           m, gamma, lambda3, c, d, costheta0, 
	#           n, beta, lambda2, B, R, D, lambda1, A
	tersoff_params = []
	contents = ''.join([line for line in open(filename) if not line.startswith('#')])
	lines = contents.split('\n\n') #expects entries to be separated by a blank line
	for line in lines:
		line = line.replace('\n','') #an entry can span any number of lines
		columns = line.split()
		if len(columns)==17:
			e1, e2, e3, m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A = [ (float(s) if s[-1].isdigit() else s) for s in columns]
			tersoff_params.append( utils.Struct(e1=e1, e2=e2, e3=e3, m=m, gamma=gamma, lambda3=lambda3, c=c, d=d, costheta0=costheta0, n=n, beta=beta, lambda2=lambda2, B=B, R=R, D=D, lambda1=lambda1, A=A) )
	
	return tersoff_params

def write_tersoff_file(system, best=False,error=-1):
	if best:
		f = open(system.name+'_best.tersoff', 'w')
		f.write('# Error: ' + str(error)+'\n')
	else:
		f = open(system.name+'.tersoff', 'w')

	pb_i_i = [t for t in system.tersoff_params if (t.e1,t.e2,t.e3)==('Pb','I','I')][0]
	for i,t in enumerate(system.tersoff_params):
		if (t.e1,t.e2,t.e3)==('I','Pb','Pb') and False:
			n = t
			t = pb_i_i
			f.write(('%-2s '*3+('%12.6g '*6)+'\n         '+('%12.6g '*8)+'\n\n') % (n.e1, n.e2, n.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
		else:
			f.write(('%-2s '*3+('%12.6g '*6)+'\n         '+('%12.6g '*8)+'\n\n') % (t.e1, t.e2, t.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
	f.close()

def set_lammps_parameters(system):
	write_tersoff_file(system)
	lmp.command('pair_coeff * * tersoff '+system.name+'.tersoff Pb I '+(' NULL'*(len(system.atom_types)-2)) ) #is it possible to do this with the LAMMPS set command?
	
	for t in system.atom_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			#if hasattr(t,'vdw_e'):
			#	lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
			#	lmp.command('pair_coeff %d * lj/cut/coul/cut %f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
	for t in system.bond_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
	for t in system.angle_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
	for t in system.dihedral_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

def calculate_error(system):
	for t in system.tersoff_params:
		try:
			powermint = int(t.m)
			assert t.c >= 0.0
			assert t.d >= 0.0
			assert t.n >= 0.0
			assert t.beta >= 0.0
			assert t.lambda2 >= 0.0
			assert t.B >= 0.0
			assert t.R >= 0.0
			assert t.D >= 0.0
			assert t.D <= t.R
			assert t.lambda1 >= 0.0
			assert t.A >= 0.0
			assert t.m - powermint == 0.0
			assert (powermint == 3 or powermint == 1)
			assert t.gamma >= 0.0
		except AssertionError:
			print 'Bad Tersoff parameters'
			return 1e10
	
	#run LAMMPS
	set_lammps_parameters(system)
	lmp.command('run 0')
	lammps_energies = lmp.extract_compute('atom_pe',1,1) #http://lammps.sandia.gov/doc/Section_python.html
	lammps_forces = lmp.extract_atom('f',3)
	
	#assign energies to their proper groups of atoms
	atom_count = 0
	for m in system.molecules:
		m.lammps_energy = sum( [lammps_energies[i] for i in range(atom_count, len(m.atoms)+atom_count)] )
		
		#testing
		for a in m.atoms:
			for b in m.atoms:
				if a.element!=b.element:
					r = utils.dist(a,b)
					eps = (a.type.vdw_e * b.type.vdw_e)**0.5
					sigma = (a.type.vdw_r * b.type.vdw_r)**0.5
					lj = 4*eps*( (sigma/r)**12 - (sigma/r)**6 )
					coul = 332*a.type.charge*b.type.charge/r
					R = 3.0
					D = 0.5
					lj_at_cutoff = 4*eps*( (sigma/(R+D))**12 - (sigma/(R+D))**6 )
					coul_at_cutoff = 332*a.type.charge*b.type.charge/(R+D)
					
					#if r>R+D:
					#	print r, len(m.atoms)			
					
					m.lammps_energy += (lj+coul if r>R+D else lj_at_cutoff+coul_at_cutoff)/2
		#end testing
		
		atom_count += len(m.atoms)

	#calculate energy error
	relative_energy_error, absolute_energy_error = 0.0, 0.0
	
	xx = []
	yy = []
	
	for elements,molecules in system.molecules_by_elements.iteritems():
		baseline_energy = molecules[0].lammps_energy #should be in order of increasing energy
		for m in molecules:
			m.lammps_energy -= baseline_energy
			relative_energy_error += ( (m.lammps_energy-m.energy)/(m.energy+1.0) )**2
			absolute_energy_error += (m.lammps_energy-m.energy)**2
			
			#print m.energy, m.lammps_energy
			xx += [m.energy]
			yy += [m.lammps_energy]
	
	if 0: #plot data and exit
		import matplotlib.pyplot as plt
		plt.plot(xx)
		plt.plot(yy)
		plt.show()
		exit()
	
	#calculate force error
	relative_force_error, absolute_force_error = 0.0, 0.0
	for i,a in enumerate(system.atoms):
		fx, fy, fz = lammps_forces[i][0], lammps_forces[i][1], lammps_forces[i][2]
		real_force_squared = a.fx**2 + a.fy**2 + a.fz**2
		#if real_force_squared < 20.0**2:
		relative_force_error += ((fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2) / (real_force_squared + 20.0**2)
		absolute_force_error += (fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2
	
	relative_force_error = math.sqrt( relative_force_error/len(system.atoms) )
	absolute_force_error = math.sqrt( absolute_force_error/len(system.atoms) )
	relative_energy_error = math.sqrt( relative_energy_error/len(system.molecules) )
	absolute_energy_error = math.sqrt( absolute_energy_error/len(system.molecules) )

	error = relative_energy_error# + relative_force_error
	
	#print absolute_energy_error, absolute_force_error, relative_energy_error, relative_force_error
	
	return error

def pack_params(system):
	params, bounds, names = [], [], []
	for t in system.atom_types:
		if True:
			params += [t.vdw_e, t.vdw_r]
			bounds += [(0.1,40.0), (2.0,4.0)]
			names += ['%d:vdw_e'%t.element, '%d:vdw_r'%t.element]
			if t.element==53:
				pass #defined by Pb
			else:
				params += [t.charge]
				bounds += [(0.0,2)]
				names += ['%d:charge'%t.element]
	'''
	for t in system.bond_types:
		params += [t.e, t.r]
		bounds += [(0,100), (0.9,3.0)]
	for t in system.angle_types:
		params += [t.e, t.angle]
		bounds += [(0,100), (0,180)]
	for t in system.dihedral_types:
		if len(t.e)==3:
			t.e = list(t.e)+[0.0]
		params += list(t.e)
		bounds += [(-100,100), (-50,50), (-20,20), (-10,10)]
	for t in system.tersoff_params:
		s = t.e1+t.e2+t.e3+':'
		if s=='PbII:' or s=='IPbPb:':
			names += [s+'c', s+'costheta0', s+'beta', s+'lambda2', s+'B', s+'lambda1', s+'A', s+'R', s+'D']
			params += [t.c, t.costheta0, t.beta, t.lambda2, t.B, t.lambda1, t.A, t.R, t.D]
			bounds += [(10,1e4), (-0.05,0.05), (0,1), (1,3), (1,1e6), (1,5), (1,1e6), (2.5,3.5), (0.1,0.5)]
		if s=='IPbI:':
			names += [s+'R']
			params += [t.R]
			bounds += [(5.5,10.0)]
	'''
	if len(params)!=len(bounds) or len(params)!=len(names):
		print 'There are %d parameters, but %d bounds and %d names!' % (len(params), len(bounds), len(names))
		exit()
		
	return params, bounds, names

def unpack_params(params, system):
	i = 0
	for t in system.atom_types:
		if True:
			t.vdw_e, t.vdw_r = params[i], params[i+1]
			i += 2
			if t.element==53:
				pass #defined by Pb
			else:
				t.charge = params[i]
				i += 1
	'''
	for t in system.bond_types:
		t.e, t.r = params[i], params[i+1]
		i += 2
	for t in system.angle_types:
		t.e, t.angle = params[i], params[i+1]
		i += 2
	for t in system.dihedral_types:
		t.e = tuple(params[i:i+4])
		i += 4
	for t in system.tersoff_params:
		s = t.e1+t.e2+t.e3+':'
		num_params = 0
		if s=='PbII:' or s=='IPbPb:':
			num_params = 9
			t.c, t.costheta0, t.beta, t.lambda2, t.B, t.lambda1, t.A, t.R, t.D = params[i:i+num_params]
		if s=='IPbI:':
			num_params = 1
			t.R = params[i]
		i+=num_params
	'''

	pb_type = [t for t in system.atom_types if t.element==82][0]
	i_type = [t for t in system.atom_types if t.element==53][0]
	i_type.charge = -pb_type.charge/2


I_ = 66
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838

extra = {
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=25.0, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.0, vdw_e=25.0, vdw_r=3.0),
}

system = utils.System(box_size=[1e3, 1e3, 1e3], name='test_tersoff')

for root, dirs, file_list in os.walk("gaussian"):
	count = 0
	for ff in file_list:
		if ff.endswith('.log'):
			name = ff[:-4]
	#for step in range(20):
	#		name = 'PbI2_r%d' % step
			if not name.startswith('PbI') : continue #for PbI testing
			if not name.endswith('_def2SVP'): continue
			energy, atoms = g09.parse_atoms(name, check_convergence=False)
			if len(atoms)>3: continue
			#if any([utils.dist(atoms[0], a)>3.5 for a in atoms]) and len(atoms)<6: continue
			total = utils.Molecule('gaussian/'+name, extra_parameters=extra, check_charges=False)
			total.energy = energy*627.509 #convert energy from Hartree to kcal/mol
			total.element_string = ' '.join( [a.element for a in total.atoms] )
			print total.element_string
			for i,a in enumerate(total.atoms):
				b = atoms[i]
				a.x, a.y, a.z = b.x, b.y, b.z
				a.fx, a.fy, a.fz = [f*1185.8113 for f in (b.fx, b.fy, b.fz)] # convert forces from Hartree/Bohr to kcal/mol / Angstrom
			system.add(total, count*200.0)
			count += 1

system.box_size[0] = count*400+200 #make system big enough to hold all atoms

#group molecules by .element_string
system.molecules_by_elements = {}
for m in system.molecules:
	if m.element_string not in system.molecules_by_elements:
		system.molecules_by_elements[ m.element_string ] = []
	system.molecules_by_elements[ m.element_string ].append(m)
#sort molecules of same type by energy, set baseline energy as zero
for element_string, molecules in system.molecules_by_elements.iteritems():
	molecules.sort(key=lambda m:m.energy)
	baseline_energy = molecules[0].energy
	for m in molecules:
		m.energy -= baseline_energy #baseline energy = 0.0

os.chdir('lammps')
files.write_lammps_data(system)

commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/cut 100.0 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary f f f
read_data	'''+system.name+'''.data

pair_coeff * * lj/cut/coul/cut 0.0 1.0

compute atom_pe all pe/atom
''').splitlines()
lmp = lammps('',['-log',system.name+'.log', '-screen','none'])
for line in commands:
	lmp.command(line)

system.tersoff_params = read_tersoff_file('input.tersoff')

def calculate_error_from_list(params):
	unpack_params(params, system)
	
	for t in system.atom_types + system.bond_types + system.angle_types + system.dihedral_types:
		t.written_to_lammps = False
	
	error = calculate_error(system)
	
	return error

initial_params, bounds, names = pack_params(system)

import numpy
from scipy.optimize import minimize, fmin_l_bfgs_b

def stochastic(use_gradient=True):
	best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
	print 'Error: %.4g' % best_min.fun
	
	def new_param_guess(start):
		while True: #keep going until new params are generated
			params = []
			for p,b in zip(start, bounds):
				new = random.gauss(p, abs(p)*0.5 + 0.01*(b[1]-b[0]) )
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
				params.append( new )
			if cmp(list(params), list(start)) != 0: #compare old and new param lists
				return params
	
	def error_gradient(x):
		e0 = calculate_error_from_list(x)
		gradient = []
		for i in range(len(x)):
			oldx = x[i]
			sign_oldx = -1 if oldx<0 else 1
			newx = oldx + 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
			if newx>bounds[i][1] or newx<bounds[i][0]:
				newx = oldx - 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
			if newx>bounds[i][1] or newx<bounds[i][0]:
				print oldx, newx, bounds[i][0], bounds[i][1]
				exit()
			dif = newx - oldx
			x[i] = newx
			gradient.append( (calculate_error_from_list(x)-e0)/dif )
			x[i] = oldx
		return numpy.array(gradient)
	
	while True:
		params = new_param_guess(best_min.x)
		if use_gradient:
			while calculate_error_from_list(params)>best_min.fun*3:
				params = new_param_guess(best_min.x)
			x, fun, _ = fmin_l_bfgs_b(calculate_error_from_list, params, fprime=error_gradient, bounds=bounds)
			guess = utils.Struct(fun=fun,x=x)
			print calculate_error_from_list(params), guess.fun, best_min.fun
		else:
			guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
			print guess.fun, best_min.fun
		if guess.fun < best_min.fun:
			best_min = guess
			print ''
			for n,x in zip(names,best_min.x):
				print '%-15s %-15g' % (n, x)
			unpack_params(best_min.x, system)
			write_tersoff_file(system,best=True,error=best_min.fun)
			print 'Error = %.4g' % best_min.fun
	return best_min

def try_params():
	params = []
	for p,b in zip(initial_params, bounds):
		if p < b[0]:
			p = b[0]
		elif p > b[1]:
			p = b[1]
		params.append( p )
	best_min = utils.Struct(fun=calculate_error_from_list(params),x=params)
	print 'Error: %.4g' % best_min.fun
	for param_index in range(len(best_min.x)):
		changed = []
		for step in range(10):
			params = []
			for i,p,b in zip(range(len(best_min.x)), best_min.x, bounds):
				new = p
				if i==param_index:
					new = random.gauss(p, 0.2*(b[1]-b[0]) )
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
				params.append( new )
			guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
			changed.append( guess.fun != best_min.fun )
		if any(changed):
			print names[param_index], 'has an effect'
	return best_min

#try_params()

stochastic(False)

