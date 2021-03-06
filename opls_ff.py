import os, copy, sys, random, cPickle, math
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

def set_lammps_parameters(system):
	for t in system.atom_types:
		if not t.written_to_lammps:
			t.written_to_lammps = True
			lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
			lmp.command('pair_coeff	%d *	%f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
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
	#soft constraint functions
	def softmin(x,xmin,tol=None):
		if not tol:
			tol = abs(xmin)*0.1
			if xmin==0.0:
				tol = 0.1
		if x>xmin+tol: return 0.0
		elif x>xmin: return 1e5*((xmin+tol-x)/tol)**2
		else: return 1e10
		# return xmin/(x-xmin) if x>xmin else 1e10
	
	def softmax(x,xmax,tol=None):
		if not tol:
			tol = abs(xmax)*0.1
			if xmax==0.0:
				tol = 0.1
		if x<xmax-tol: return 0.0
		elif x<xmax: return 1e5*((x-xmax-tol)/tol)**2
		else: return 1e10

	# apply soft constraints
	constraint_error = 0.0
	for t in system.atom_types:
		constraint_error += softmin(t.vdw_e,0.001)
		constraint_error += softmin(t.vdw_r,0.5)
	for t in system.bond_types:
		constraint_error += softmin(t.e,10.0)
		constraint_error += softmin(t.r,0.5)
	for t in system.angle_types:
		constraint_error += softmin(t.e,1.0)

	#Run LAMMPS
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
	relative_energy_error = 0.0
	absolute_energy_error = 0.0
	for elements,molecules in system.molecules_by_elements.iteritems():
		for m in molecules:
			m.lammps_energy -= molecules[0].lammps_energy #should be in order with minimum first
			relative_energy_error += ( (m.lammps_energy-m.energy)/(m.energy+1.0) )**2
			absolute_energy_error += (m.lammps_energy-m.energy)**2
	# energy_error = 0.0
	# for elements,molecules in system.molecules_by_elements.iteritems():
	# 	for m in molecules:
	# 		m.lammps_energy -= molecules[0].lammps_energy #should be in order with minimum first
	# 	energy_error += ( (m.lammps_energy-m.energy)/(m.energy+1.0) )**2
	
	#calculate force error
	force_error = 0.0
	for i,a in enumerate(system.atoms):
		fx, fy, fz = lammps_forces[i][0], lammps_forces[i][1], lammps_forces[i][2]
		real_force_squared = a.fx**2 + a.fy**2 + a.fz**2
		if real_force_squared < 20.0**2:
			#force_error += ((fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2) / (real_force_squared + 20)
			force_error += (fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2
	
	force_error = math.sqrt( force_error/len(system.atoms) )
	relative_energy_error = math.sqrt( relative_energy_error/len(system.molecules) )
	absolute_energy_error = math.sqrt( absolute_energy_error/len(system.molecules) )

	error = relative_energy_error + force_error + constraint_error
	
	print absolute_energy_error, force_error
	
	return error

def pack_params(system):
	params = []
	for t in system.atom_types:
		params += [t.charge, t.vdw_e, t.vdw_r]
	for t in system.bond_types:
		params += [t.e, t.r]
	for t in system.angle_types:
		params += [t.e, t.angle]
	for t in system.dihedral_types:
		if len(t.e)==3:
			t.e = list(t.e)+[0.0]
		params += list(t.e)
	return params

def unpack_params(params, system):
	i = 0
	for t in system.atom_types:
		t.charge, t.vdw_e, t.vdw_r = params[i], params[i+1], params[i+2]
		i += 3
	for t in system.bond_types:
		t.e, t.r = params[i], params[i+1]
		i += 2
	for t in system.angle_types:
		t.e, t.angle = params[i], params[i+1]
		i += 2
	for t in system.dihedral_types:
		t.e = tuple(params[i:i+4])
		i += 4

I_ = 66
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838

extra = {
	(H_, I_): (100.0, 2.1), 
	(N_, H_, I_): (10.0, 180.0), 
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=1.0, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=1.0, vdw_r=3.0),
	(Pb_, I_): (100.0, 2.9), 
	(I_, Pb_, I_): (10.0, 95.0),
	(13, 53, 54, 66): (0.0,0.0,0.0),
	(54, 53, 54, 66): (0.0,0.0,0.0),
}

system = utils.System(box_size=[1e3, 1e3, 1e3], name='test_opls_0')

for root, dirs, file_list in os.walk("gaussian"):
	count = 0
	for ff in file_list:
		if ff.endswith('.log'):
			name = ff[:-4]
			if not name.startswith('PbI2'): continue #for PbI2 testing
			energy, atoms = g09.parse_atoms(name)
			total = utils.Molecule('gaussian/'+name, extra_parameters=extra)
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
	for m in molecules:
		m.energy -= molecules[0].energy #minimum energy = first molecule = 0.0


os.chdir('lammps')
files.write_lammps_data(system)

commands = ('''units real
atom_style full
pair_style lj/cut/coul/cut 10.0 100.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary f f f
read_data	'''+system.name+'''.data

compute atom_pe all pe/atom
''').splitlines()
lmp = lammps('',['-log',system.name+'.log','-screen','none'])
for line in commands:
	lmp.command(line)

def calculate_error_from_list(params):
	unpack_params(params, system)
	
	for t in system.atom_types + system.bond_types + system.angle_types + system.dihedral_types:
		t.written_to_lammps = False
	
	error = calculate_error(system)

	return error

# from scipy.optimize import minimize
#minimize(calculate_error_from_list, pack_params(system), method='Nelder-Mead', options={'disp':True,'maxfev':1000000,'maxiter':1000000})
#minimize(calculate_error_from_list, pack_params(system), method='Powell', options={'disp':True,'maxfev':1000000}) #seems to work better, at least for small number of parameters
from scipy.optimize import fmin_powell
fmin_powell(calculate_error_from_list, pack_params(system), full_output=True, ftol=0.0)

os.chdir('..')


