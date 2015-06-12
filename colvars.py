import os, copy, sys, random, cPickle, math
sys.path.append("/fs/home/jms875/Library")
import utils, filetypes, lammps

coul_on = True
Q_params = 3.3 # lattice constant/2 for these params

def center_of_geometry(atoms):
	return utils.Struct(  x = sum(a.x for a in atoms)/len(atoms),
		y = sum(a.y for a in atoms)/len(atoms),
		z = sum(a.z for a in atoms)/len(atoms)   )

def radius_of_gyration(atoms):
		center = center_of_geometry(atoms)
		return math.sqrt( 1.0/len(atoms) * sum( [utils.dist_squared(a,center) for a in atoms] ) )


Pb = 909
I = 878
C = 80
N = 78
HC = 85
HN = 79

def write_input_header(f, run_name, atom_types, read_restart=False):
	type_indices = dict( [(t.index,i+1) for i,t in enumerate(atom_types)] )
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	pb_types = [i+1 for i,t in enumerate(atom_types) if t.element==82]
	s_types = [i+1 for i,t in enumerate(atom_types) if t.element==16]
	o_types = [i+1 for i,t in enumerate(atom_types) if t.element==8]
	c_types = [i+1 for i,t in enumerate(atom_types) if t.element==6]
	h_types = [i+1 for i,t in enumerate(atom_types) if t.element==1]
	if not read_restart:
		f.write('''units		real
atom_style	full
variable pairstyle string '''+('lj/cut/coul/long' if coul_on else 'lj/cut')+'''
pair_style hybrid/overlay lj/cut'''+('/coul/long 10.0' if coul_on else ' 10.0')+''' morse 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary p p p
read_data	'''+run_name+'''.data''')

	f.write('''
	kspace_style pppm 1.0e-3
'''+( '\n'.join([ ('#%d:\t%d\t%s' % (i+1,t.index,t.element)) for i,t in enumerate(atom_types) ]) )+'''

pair_coeff * * morse 0.0 1.0 1.0
pair_coeff * * ${pairstyle}	0.0	1.0\n''')

	#set normal LJ from OPLS
	f.write(('\n'.join(["pair_coeff %d %d ${pairstyle}	%f %f" % (atom_type_numbers[t], atom_type_numbers[t], t.vdw_e, t.vdw_r) for i,t in enumerate(atom_types) if i+1 not in pb_types+s_types])) + '\n')

	#for pair in [ (S,S), (S,OO), (S,OH), (OH,OH) ]:
	#	f.write('pair_coeff	%d	%d	${pairstyle}	0.01	3.5\n' % tuple(sorted([type_indices[i] for i in pair ])) )
	for i in c_types: #repel Pb and S from C
		for j in pb_types+s_types:
			f.write('pair_coeff	%d	%d	${pairstyle}	0.01	4.0\n' % tuple(sorted([i, j])) )
	for i in h_types: #repel Pb and S from H
		for j in pb_types+s_types:
			f.write('pair_coeff	%d	%d	${pairstyle}	0.01	3.5\n' % tuple(sorted([i, j])) )
	
	for types in params: #LJ params - overwrite any earlier ones
		if type(types)==type(1) or len(params[types])!=2: continue
		if types[0] in type_indices and types[1] in type_indices:
			f.write('pair_coeff %d %d ${pairstyle} %f %f \n' % tuple( sorted([type_indices[i] for i in types]) + params[types] ) )
	
	for types in params: #Morse params
		if type(types)==type(1) or len(params[types])!=3: continue
		if types[0] in type_indices and types[1] in type_indices:
			f.write('pair_coeff %d %d morse %f %f %f\n' % tuple( sorted([type_indices[i] for i in types]) + params[types] ) )
			f.write('pair_coeff %d %d ${pairstyle} 0.0 1.0 \n' % tuple( sorted([type_indices[i] for i in types]) ) ) #set corresponding LJ to zero

def run_job(run_name, on_queue=True):
	if on_queue:
		f = open(run_name+'.nbs', 'w')
		f.write('''#!/bin/bash
##NBS-nproc: 1
##NBS-queue: "long"

/fs/home/jms875/build/lammps/lammps-9Dec2014/src/lmp_serial -in %s.in -log %s.log &> /dev/null
''' % (run_name, run_name) )
		f.close()
		os.system("jsub "+run_name+".nbs")
	else: #run locally
		os.system('/fs/home/jms875/build/lammps/lammps-9Dec2014/src/lmp_serial_clang -in %s.in -log %s.log' % (run_name,run_name) )

def save_to_file( atoms, bonds, angles, dihedrals, atom_types, run_name ):
	#prevent loops in stored object by replacing references with interger indices
	for x in bonds+angles+dihedrals:
		x.atoms = [a.index for a in x.atoms]
	for a in atoms:
		try:
			a.bonded = [b.index for b in a.bonded]
		except: pass
		a.type = a.type.index
		a.neighbors = None
		for x in a.__dict__: #turn things with bad types into floats
			t = type(a.__dict__[x])
			if t not in [type(1), type(1.0), type(''), type(None)]:
				try:
					a.__dict__[x] = float( a.__dict__[x] )
				except:
					pass
	cPickle.dump( (atoms, bonds, angles, dihedrals, atom_types), open(run_name+'.pickle', 'wb') )
	#reset
	for x in bonds+angles+dihedrals:
		x.atoms = [atoms[ii-1] for ii in x.atoms]
	for a in atoms:
		try:
			a.bonded = [atoms[ii-1] for ii in a.bonded]
		except: pass
		a.type = [t for t in atom_types if t.index==a.type][0]

def load_from_file(run_name):
	atoms, bonds, angles, dihedrals, atom_types = cPickle.load( open('lammps/'+run_name+'.pickle', 'rb') )
	#reload references from indices
	for x in bonds+angles+dihedrals:
		x.atoms = [atoms[ii-1] for ii in x.atoms]
	for a in atoms:
		try:
			a.bonded = [atoms[ii-1] for ii in a.bonded]
		except: pass
		a.type = [t for t in atom_types if t.index==a.type][0]
	#reload cartesian coordinates from xyz file
	xyz_lines = open('lammps/'+run_name+'.xyz').readlines()[-len(atoms):]
	xyz = [ line.split()[1:] for line in xyz_lines ]
	for i in range(len(xyz)):
		atoms[i].x, atoms[i].y, atoms[i].z = [float(s) for s in xyz[i]]
	#center core atoms
	center = center_of_geometry( [a for a in atoms if a.element=='S'] )
	for a in atoms:
		a.x -= center.x; a.y -= center.y; a.z -= center.z
	return atoms, bonds, angles, dihedrals, atom_types


def make_lattice(params, system_size, prior_run_name, lc_x=1, lc_y=1, lc_z=1, random_seed=1, run_name_prefix='lattice_', on_queue=True):
	atoms, bonds, angles, dihedrals = [], [], [], []
	x, y, z = 0, 0, 0
	spacer = utils.Molecule('ma.arc')
	for xi in range(system_size):
		x, y, z = 0, 0, 0
		for yi in range(system_size):
			x, y, z = 0, 0, 0
			for zi in range(system_size):
				x=xi*lc_x
				y=yi*lc_y
				z=zi*lc_z
				unbonded_atoms = filetypes.parse_xyz('/fs/home/jar528/Desktop/nano/xyz/unbonded.xyz')
				for a in unbonded_atoms:
					a.x += x
					a.y += y
					a.z += z
				atoms = atoms + unbonded_atoms
				spacer.add_to(x, y, z, atoms, bonds, angles, dihedrals)
	
	iodine_type = [t for t in utils.Molecule.atom_types if t.element_name=='I'][-1]
	lead_type = [t for t in utils.Molecule.atom_types if t.element_name=='Pb'][-1]
	
	for a in atoms:
		if a.element=='I':
			a.type=iodine_type
		if a.element=='Pb':
			a.type=lead_type
	
	#filetypes.write_xyz('/fs/home/jar528/Desktop/nano/out', atoms)
	#exit()
	
	random_s = str(random_seed)
	random.seed(random_seed)
	run_name = run_name_prefix+str(system_size)

	# make sure atom types are the same in atom1 and atom2
	atom_types_by_type_index = dict( [(t.type.index,t.type) for t in atoms] )
	for i,a in enumerate(atoms):
		a.type = atom_types_by_type_index[ a.type.index ]
		a.index = i+1
		if coul_on and a.type.index in params:
			a.charge = params[a.type.index]
		else:
			a.charge = a.type.charge
	
	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_types.sort(key=lambda t:-t.element-t.index*0.00001)
	
	box_size = [ 3.24+max([a.x for a in atoms])-min([a.x for a in atoms]),
				 3.24+max([a.y for a in atoms])-min([a.y for a in atoms]),
				 3.24+max([a.z for a in atoms])-min([a.z for a in atoms]) ]
	
	os.chdir('lammps')
	save_to_file( atoms, bonds, angles, dihedrals, atom_types, run_name )
	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types, pair_coeffs_included=False)
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')
	
	f = open(run_name+'.in', 'w')
	write_input_header(f, run_name, atom_types, read_restart=False)
	f.write('''
dump	1 all xyz 10000 '''+run_name+'''.xyz
minimize 0.0 1.0e-8 10000 1000
neigh_modify check yes every 1 delay 0
thermo 10000
restart 100000 '''+run_name+'''.restart1 '''+run_name+'''.restart2
#fix motion all nve
#fix implicit_solvent all langevin 400.0 400.0 1000.0 '''+random_s+''' zero yes gjf yes
fix motion all npt temp 450.0 450.0 100.0 aniso 1.0 1.0 1000.0
velocity	all create 400.0 '''+random_s+''' mom yes
timestep	2.0
thermo_style	custom step temp etotal pe tpcpu
run '''+str(int(1e6))+'''
#thermo 1
#run 1000
minimize 0.0 1.0e-8 10000 1000
write_restart '''+run_name+'''.restart
''')
	f.close()
	run_job(run_name, on_queue)
	os.chdir('..')






params = {
(Pb,Pb):[5.000000,1.500000,6.484328],
(Pb,I):[5.000000,1.500000,3.232639],
(Pb,C):[5.000000,1.500000,5.690018],
(Pb,N):[5.000000,1.500000,5.330574],
(Pb,HC):[5.000000,1.500000,5.121428],
(Pb,HN):[5.000000,1.500000,4.436786],
(I,I):[5.000000,1.500000,4.436952],
(I,C):[5.000000,1.500000,4.256463],
(I,N):[5.000000,1.500000,3.787180],
(I,HC):[5.000000,1.500000,3.525462],
(I,HN):[5.000000,1.500000,2.772864],
Pb	: 0.788807,
I	:-0.547291	
}

#C	:-0.566155,
#N	:-0.634688,
#HC	: 0.279802,
#HN	: 0.400554,

utils.Molecule.set_params('oplsaa.prm')

def react_all_dots():
	for step in range(1):
		for sum_size in range(2,max_dot_size+1):
			for dot_size1 in range(1,sum_size/2+1):
				dot_size2 = sum_size-dot_size1
				react_dots(params, dot_size1, dot_size2, 'c18_'+str(dot_size1), 'c18_'+str(dot_size2), random_seed=step+1, run_name_prefix='react_%d__' % step, on_queue=True)

def plot_all():
	import matplotlib.pyplot as plt

	colors = {1000:'blue', 500:'red', 250:'green', 100:'black'}
	for seed in range(0,5):
		for h in [1000]: #[1000,500,250,100]:
			xx = []; yy = []
			for line in open('lammps/blocky0.2_s%d__5_5.pmf' % (seed)):
				a = line.split()
				if a:
					try:
						xx.append( float(a[0]) )
						yy.append( float(a[1]) )
					except ValueError: pass
			yy = [y-yy[-1] for y in yy]
			plt.plot( xx, yy, label='s%d'%seed)
		
	plt.xlabel('Distance (Angstroms)')
	plt.ylabel('Free Energy (kcal/mol)')
	#plt.xlim(6,16)
	plt.legend()
	plt.show()

best_dots = ['c18_s6_1', 'c18_s12_2', 'c18_s11_3', 'c18_s15_4', 'c18_s10_5', 'c18_s6_6', 'c18_s15_7', 'c18_s11_8', 'c18_s13_9', 'c18_s12_10']

#for seed in range(0,5):
	#for h in [1000,500,250,100]:
	#	react_dots(params, 5, 5, 'c18_5', 'c18_5', run_name_prefix='meta_s%d_h%d__' % (seed,h), random_seed=seed+1, on_queue=True, method='metadynamics', hill_frequency=h)
	
#	react_dots(params, 5, 5, 'c18_s10_5', 'c18_s10_5', run_name_prefix='blocky_heavy_s%d__' % seed, random_seed=seed+1, on_queue=True, method='metadynamics', hill_frequency=1000)
	#react_dots(params, 10, 10, 'c18_s12_10', 'c18_s12_10', run_name_prefix='long3_s%d__' % seed, random_seed=seed+1, on_queue=True, method='metadynamics', hill_frequency=500)

def make_dots_seed():
	for dot_size in range(1,11):
		dots = []
		for seed in range(16,21):
			log = open('lammps/c18_s%d_%d.log' % (seed, dot_size) ).read()
			start = log.rfind('Step Temp TotEng PotEng KinEng E_pair E_bond T/CPU')
			start = log.find('\n', start)
			start = log.find('\n', start+1)
			a = log[start:start+80].split()
			energy_min = float(a[2])
			#print 'lammps/c18_s%d_%d.log' % (seed, dot_size), energy_min
			dots.append( (energy_min, seed) )
		#print dot_size, min(dots)
		print "'c18_s%d_%d'," % (min(dots)[1], dot_size), 
		best_seed = min(dots)[1]
		#for seed in range(6,16):
		#	make_dot(params, dot_size, 'c18_s%d_%d'%(best_seed,dot_size), random_seed=seed, run_name_prefix='c18_s%d_'%seed, on_queue=True)

#plot_all()

#for dot_size in range(1,11):
#	for seed in range(19,21):
#		make_dot(params, dot_size, None, random_seed=seed, run_name_prefix='c18_s%d_'%seed, on_queue=True)

for seed in range(1,2):
	make_lattice(params, 2, None, lc_x=6.48, lc_y=6.48, lc_z=6.48, random_seed=seed, run_name_prefix='lattice_%d_'%seed, on_queue=False)

