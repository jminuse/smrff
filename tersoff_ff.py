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

def write_tersoff_file(system):
	f = open(system.name+'.tersoff', 'w')
	pb_i_i = [t for t in system.tersoff_params if (t.e1,t.e2,t.e3)==('Pb','I','I')][0]
	for i,t in enumerate(system.tersoff_params):
		if (t.e1,t.e2,t.e3)==('I','Pb','Pb'):
			n = t
			t = pb_i_i
			f.write(('%3s '*3+('%8.8g '*6)+'\n            '+('%8.8g '*8)+'\n\n') % (n.e1, n.e2, n.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
		else:
			f.write(('%3s '*3+('%8.8g '*6)+'\n            '+('%8.8g '*8)+'\n\n') % (t.e1, t.e2, t.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
	f.close()

def set_lammps_parameters(system):
	write_tersoff_file(system)
	lmp.command('pair_coeff * * tersoff '+system.name+'.tersoff Pb I '+(' NULL'*(len(system.atom_types)-2)) ) #is it possible to do this with the LAMMPS set command?
	
	for t in system.atom_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			if hasattr(t,'vdw_e'):
				lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
				lmp.command('pair_coeff %d * lj/cut/coul/cut %f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
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
		atom_count += len(m.atoms)
	
	#calculate energy error
	relative_energy_error, absolute_energy_error = 0.0, 0.0
	for elements,molecules in system.molecules_by_elements.iteritems():
		baseline_energy = molecules[0].lammps_energy #should be in order of increasing energy
		for m in molecules:
			m.lammps_energy -= baseline_energy
			relative_energy_error += ( (m.lammps_energy-m.energy)/(m.energy+1.0) )**2
			absolute_energy_error += (m.lammps_energy-m.energy)**2
			#print m.energy, m.lammps_energy
	#exit()
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

	error = relative_energy_error + relative_force_error
	
	#print absolute_energy_error, absolute_force_error, relative_energy_error, relative_force_error
	
	return error

def pack_params(system):
	params, bounds, names = [], [], []
	for t in system.atom_types:
		#params += [t.charge, t.vdw_e, t.vdw_r] #temporarily take away Lennard-Jones
		#if t.element==53:
		#	pass #defined by Pb
		#else:
		#	params += [t.charge]
		#	bounds += [(0.0,2)]
		#	names += ['%d charge'%t.element]
		pass
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
		if t.e1=='Pb' and t.e2=='I' and t.e3=='I':
			s = t.e1+t.e2+t.e3
			#names += [s+' A', s+' B', s+' l1', s+' l2']
			#params += [ t.A, t.B, t.lambda1, t.lambda2 ]
			#bounds += [ (0,1e6), (0,1e6), (0,6), (0,3)]
			
			names = [s+'lambda3', s+'c', s+'d', s+'costheta0', s+'n', s+'beta', s+'lambda2', s+'B', s+'lambda1', s+'A', s+'D', s+'R']
			params += [t.lambda3, t.c, t.d, t.costheta0, t.beta, t.lambda2, t.B, t.lambda1, t.A, t.D, t.R]
			bounds += [  (-3,3), (0,1e6), (0,1e6), (-1,1),   (0,1),   (0,3), (0,1e6), (0,6), (0,1e6), (0, 10), (0,10 )]

		if t.e1=='I' and t.e2=='I' and t.e3=='I':
			s = t.e1+t.e2+t.e3
			names = [s+'lambda3', s+'c', s+'d', s+'costheta0', s+'n', s+'lambda2', s+'B', s+'lambda1', s+'A', s+'D', s+'R']
			params += [t.lambda3, t.c, t.d, t.costheta0, t.lambda2, t.B, t.lambda1, t.A, t.D, t.R]
			bounds += [  (-3,3), (0,1e6), (0,1e6), (-1,1),   (0,3), (0,1e6), (0,6), (0,1e6), (0, 10), (0,10 )]
		if t.e1=='I' and t.e2=='Pb' and t.e3=='I':
			s = t.e1+t.e2+t.e3
			names = [s+'lambda3', s+'c', s+'d', s+'costheta0', s+'D', s+'R']
			params += [t.lambda3, t.c, t.d, t.costheta0, t.D, t.R]
			bounds += [  (-3,3), (0,1e6), (0,1e6), (-1,1), (0, 10), (0,10 )]
		if t.e1=='I' and t.e2=='I' and t.e3=='Pb':
			s = t.e1+t.e2+t.e3
			names = [s+'lambda3', s+'c', s+'d', s+'costheta0', s+'D', s+'R']
			params += [t.lambda3, t.c, t.d, t.costheta0, t.D, t.R]
			bounds += [  (-3,3), (0,1e6), (0,1e6), (-1,1), (0, 10), (0,10 )]

	return params, bounds, names

def unpack_params(params, system):
	i = 0
	for t in system.atom_types:
		#t.charge, t.vdw_e, t.vdw_r = params[i], params[i+1], params[i+2] #temporarily take away Lennard-Jones
		#i += 3
		#if t.element==53:
		#	pass #defined by Pb
		#else:
		#	t.charge = params[i]
		#	i += 1
		pass
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
		if t.e1=='Pb' and t.e2=='I' and t.e3=='I':
			#num_params=4
			#t.A, t.B, t.lambda1, t.lambda2 = params[i:i+num_params]
			num_params=9
			t.lambda3, t.c, t.d, t.costheta0, t.beta, t.lambda2, t.B, t.lambda1, t.A = params[i:i+num_params]
			i+=num_params

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
	(H_, I_): (100.0, 2.1), 
	(N_, H_, I_): (10.0, 180.0), 
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=0.1, vdw_r=2.0, D0=5.0, alpha=1.5, r0=2.8),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=0.1, vdw_r=2.0, D0=5.0, alpha=1.5, r0=2.8),
	(Pb_, I_): (100.0, 2.9), 
	(I_, Pb_, I_): (10.0, 95.0),
	(13, 53, 54, 66): (0.0,0.0,0.0),
	(54, 53, 54, 66): (0.0,0.0,0.0),
}

system = utils.System(box_size=[1e3, 1e3, 1e3], name='test_tersoff_')

for root, dirs, file_list in os.walk("gaussian"):
	count = 0
	for ff in file_list:
		if ff.endswith('.log'):
			name = ff[:-4]
	#for step in range(20):
	#		name = 'PbI2_r%d' % step
			if not name.startswith('PbI'): continue #for PbI testing
			energy, atoms = g09.parse_atoms(name)
			total = utils.Molecule('gaussian/'+name, extra_parameters=extra, check_charges=False)
			total.energy = energy*627.509 #convert energy from Hartree to kcal/mol
			total.element_string = ' '.join( [a.element for a in total.atoms] )
			#print total.element_string
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
lmp = lammps('',['-log',system.name+'.log','-screen','none'])
for line in commands:
	lmp.command(line)

system.tersoff_params = read_tersoff_file('input.tersoff')

log = open(system.name+'.params', 'w')
def calculate_error_from_list(params):
	unpack_params(params, system)
	
	for t in system.atom_types + system.bond_types + system.angle_types + system.dihedral_types:
		t.written_to_lammps = False
	
	error = calculate_error(system)
	
	log.write( '%s %g\n' % (list(params), error))
	
	return error

initial_params, bounds, names = pack_params(system)

['PbIIlambda3', 'PbIIc', 'PbIId', 'PbIIcostheta0', 'PbIIn', 'PbIIbeta', 'PbIIlambda2', 'PbIIB', 'PbIIlambda1', 'PbIIA']
#initial_params = [-3.0, 4.658513032577969e-05, 0.7984278517358414, -1.8699527049207918e-05, 322795.143059622, 1.109701376931432, 728.467043363168, 4.234875101724356, 1163313.832162179] #error=0.7721
# initial_params = [1.8945935591480063, 1.912324201715158e-05, 2.0538627438753703, -3.301316559841533e-05, 16.319977401150283, 0.6836858011442968, 526.419534975563, 3.6302755638145436, 392278.1549681662] #Error: 0.6949
# initial_params =[0.47293110512335906, 1.3225865502970848, 0.3658076226854543, -0.38002280592624177, 1.0, 0.057386233592209873, 366.33965385242209, 4.472688308772109, 940572.01756813889, 0.92134478301331102, 6.7770354666052963, -2.0373861809573515, 0.0005194447955741172, 0.55775726382329172, -5.9713515576926485e-05, 0.43495677355222284, 4.99564519696749, -0.59630213159043055, 0.00028657196228236415, 1.116699397510942, -0.0010051240704703761, 0.46277772840479586, 8.3307641228120843, -2.1679739873084913, 0.00019456465255806936, 0.90470459269209469, -0.0024578905750086292, 1.13454496092724, 7839.8597701875433, 1.0800612906727878, 45605.191689373984, 3.045286281508103, 1.4514354004279451] #Error 1.206

initial_params = [0.28279425601204583, 1.965384767032681e-07, 1.3665418684481765e-07, -5.5171000313435071e-15, 6.052432850916392, 0.044029809020889848, 605.6088594758387, 4.8165942034667077, 3732878.2397423275, 0.0026613704450947225, 0.082704984361372857, -0.15418251337807234, 0.00048405438645680973, 0.0016544462180032736, -7.4899559919088185e-06, 2.5353410853603125e-06, 3.4918035407623299e-05, -0.015679795456745034, 2.6315550638796417e-06, 0.56067750921784132, 3.5285044070735273e-06, 3.0868693593091818e-05, 0.00045919924952706116, 0.00045624353982215882, 0.00012924810132165212, 0.0038157634148834815, -1.8390806364969089e-05, 0.0053804960347081073, 0.033099316600820453, 1.2051557308286559e-05, 12.06688882133929, 1.1456424536215685e-07, 5.1848727883446824e-06] # Error: 0.9367

import numpy
from scipy.optimize import minimize

def randomize():
	best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
	print 'Error: %.4g' % best_min.fun
	while True:
		params = []
		for p,b in zip(best_min.x, bounds):
			new = random.gauss(p, p*0.2 if p!=0.0 else 0.001)
			if new < b[0]:
				new += (b[0]-new)
			if new > b[1]:
				new -= (b[1]-new)
			params.append( new )
		# guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
		# guess = minimize(calculate_error_from_list, params, bounds=bounds) #L-BFGS-B
		guess = minimize(calculate_error_from_list, params, method='Nelder-Mead') #L-BFGS-B
		log.write('---\n')
		if guess.fun < best_min.fun:
			best_min = guess
			print list(best_min.x)
			print 'Error: %.4g' % best_min.fun
	return best_min

#best_min = minimize(calculate_error_from_list, initial_params, bounds=bounds, method='Nelder-Mead')
#best_min = minimize(calculate_error_from_list, initial_params, bounds=bounds, method='L-BFGS-B')

best_min = randomize()

print names
print list(best_min.x)
print bounds
print 'Error: %.4g' % best_min.fun

os.chdir('..')

