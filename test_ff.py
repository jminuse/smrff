import os, copy, sys, random, cPickle, math
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0")
import utils, files, g09

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

system = utils.System(box_size=[1e3, 1e3, 1e3], name='test3')

for root, dirs, file_list in os.walk("gaussian"):
	count = 0
	for ff in file_list:
		if ff.endswith('.log'):
			name = ff[:-4]
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

system.box_size[0] = count*200+200 #make system big enough to hold all atoms

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
lmp = lammps()
for line in commands:
	lmp.command(line)

def set_lammps_parameters(mutable_parameters):
	for indices,params in mutable_parameters.iteritems():
		if type(indices)==int: #single type
			#params = utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=1.0, vdw_r=3.0)
		if len(indices)==2: #bond
			e=params[0]; r=params[1]
		elif len(indices)==3: #angle
			e=params[0]; angle=params[1]
		elif len(indices)==4: #dihedral
			e=tuple(params)

	#TODO: find LAMMPS type index for each OPLS type

	#for params in mutable_parameters:
	#	lmp.command('set type 1 charge %f' % 0.1)
	#	lmp.command('pair_coeff	1 2		%f	2.0' % 1.0)
		

def calculate_error(system):
	#set_lammps_parameters()
	lmp.command('run 0')
	lammps_energies = lmp.extract_compute('atom_pe',1,1) #http://lammps.sandia.gov/doc/Section_python.html
	lammps_forces = lmp.extract_atom('f',3)
	
	#assign energies to their proper groups of atoms
	atom_count = 0
	for m in system.molecules:
		m.lammps_energy = sum( [lammps_energies[i] for i in range(atom_count, len(m.atoms)+atom_count)] )
		atom_count += len(m.atoms)
	
	#calculate energy error
	energy_error = 0.0
	for elements,molecules in system.molecules_by_elements.iteritems():
		for m in molecules:
			m.lammps_energy -= molecules[0].lammps_energy #should be in order with minimum first
		energy_error += ( (m.lammps_energy-m.energy)/(m.energy+1.0) )**2
	
	#calculate force error
	force_error = 0.0
	for i,a in enumerate(system.atoms):
		fx, fy, fz = lammps_forces[i][0], lammps_forces[i][1], lammps_forces[i][2]
		real_force_squared = a.fx**2 + a.fy**2 + a.fz**2
		if real_force_squared < 20.0**2:
			force_error += ((fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2) / (real_force_squared + 20)
	
	error = (energy_error + force_error) / len(system.atoms)
	
	return error

print calculate_error(system)

os.chdir('..')


