import os, copy, sys, random
import utils, filetypes, lammps

def ligand_binding(params, n_ligands):
	utils.Molecule.set_params('oplsaa.prm')
	run_name = 'lig_'+str(n_ligands)
	
	ligand = utils.Molecule('Pb_b2.arc')
	monomer = utils.Molecule('PbS.arc')
	atoms = []
	bonds = []
	angles = []
	dihedrals = []
	
	for a in ligand.atoms + monomer.atoms:
		if a.type.element==8: a.charge = -(1.2+params[0])*0.25
		if a.type.element==16: a.charge = -params[0]
		if a.type.element==82: a.charge = params[0]
	
	box_size = [100.0]*3
	x,y,z = 0.0, 0.0, 0.0
	if n_ligands == -1:
		ligand.add_to(x, y, z, atoms, bonds, angles, dihedrals)
	else:
		monomer.add_to(x, y, z, atoms, bonds, angles, dihedrals)
		if n_ligands>0: ligand.add_to(x, y+5, z, atoms, bonds, angles, dihedrals)
		for a in ligand.atoms:
			a.y = -a.y
		if n_ligands>1: ligand.add_to(x, y-5, z, atoms, bonds, angles, dihedrals)
		for a in ligand.atoms:
			a.z,a.y = a.y,a.z
		if n_ligands>2: ligand.add_to(x, y, z-5, atoms, bonds, angles, dihedrals)
		for a in ligand.atoms:
			a.z = -a.z
		if n_ligands>3: ligand.add_to(x, y, z+5, atoms, bonds, angles, dihedrals)
		for a in ligand.atoms:
			a.z,a.x = a.x,a.z
		if n_ligands>4: ligand.add_to(x+6, y, z, atoms, bonds, angles, dihedrals)
	
	positions = []
	for line in open('lammps/lig_'+str(n_ligands)+'_.xyz'):
		columns = line.split()
		if len(columns)>3:
			index, x, y, z = columns
			positions.append( [float(s) for s in [x,y,z]] )
	
	positions = positions[-len(atoms):]
	
	for i,a in enumerate(atoms):
		#positions[i] = [ positions[i][0]-positions[0][0], positions[i][1]-positions[0][1], positions[i][2]-positions[0][2] ]
		a.x, a.y, a.z = positions[i]
	
	#for d in dihedrals:
	#	print [a.type.element for a in d.atoms], [a.type.index2 for a in d.atoms]
	
	#filetypes.write_xyz('out', atoms)
	#sys.exit(0)
	
	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_types.sort(key=lambda t:-t.element)
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	
	pb_types = [i+1 for i,t in enumerate(atom_types) if t.element==82]
	s_types = [i+1 for i,t in enumerate(atom_types) if t.element==16]
	o_types = [i+1 for i,t in enumerate(atom_types) if t.element==8]
	c_types = [i+1 for i,t in enumerate(atom_types) if t.element==6]
	h_types = [i+1 for i,t in enumerate(atom_types) if t.element==1]
	
	tersoff_element_names = {82:'Pb', 16:'S', 8:'O'}
	tersoff_coeffs_string = ' '.join([tersoff_element_names[t.element] if t.element in tersoff_element_names else 'NULL' for t in atom_types])
	
	cooh_c_types = [i+1 for i,t in enumerate(atom_types) if t.index==213]
	
	is_charged = True

	os.chdir('lammps')
	
	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types, pair_coeffs_included=False)
	#os.chdir('..')
	#return
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')


	f = open('PbSO_test.tersoff', 'w')
	PbSO = ['Pb', 'S', 'O']
	count_params = 0
	for i in range(3):
		for j in range(3):
			for k in range(3):
				m = 3.0 #can only be 1.0 or 3.0 by definition
				gamma = 1.0
				lambda3 = 0.0
				costheta0 = 0.0
				
				if i!=0 and j!=0:
					c, d, n, beta, lambda2, R, D, lambda1, A, B = params[1:11]
					A, B = 0.0, 0.0
				else:
					c, d, n, beta, lambda2, R, D, lambda1, A, B = params[count_params*10 + 1: count_params*10 + 11]
					count_params += 1
				
				f.write(('%3s '*3+('%8.8g '*6)+'\n            '+('%8.8g '*8)+'\n\n') % (PbSO[i], PbSO[j], PbSO[k], m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A))
	f.close()

	f = open(run_name+'.in', 'w')
	f.write('''units	real
atom_style	full
pair_style hybrid/overlay lj/cut/coul/long 10.0 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5
kspace_style pppm 1.0e-2

read_data	'''+run_name+'''.data

pair_coeff * * tersoff PbSO_test.tersoff '''+tersoff_coeffs_string+'\n')

	f.write('pair_coeff	*	*	lj/cut/coul/long	0.0		1.0\n')

	for i,t in enumerate(atom_types): #repel Pb and S from all but Pb, S, and O
		if (pb_types and i+1==pb_types[0]) or (s_types and i+1==s_types[0]) or (o_types and i+1==o_types[0]): continue
		if pb_types and i+1 in cooh_c_types: #don't repel COOH carbon from Pb too much
			f.write('pair_coeff	%d	%d	lj/cut/coul/long	0.01	3.0\n' % tuple(sorted([pb_types[0], i+1])) )
		elif pb_types and i+1 in c_types: #repel other C from Pb harder
			f.write('pair_coeff	%d	%d	lj/cut/coul/long	0.01	5.0\n' % tuple(sorted([pb_types[0], i+1])) )
		elif pb_types: f.write('pair_coeff	%d	%d	lj/cut/coul/long	0.01	4.0\n' % tuple(sorted([pb_types[0], i+1])) )
		if s_types:  f.write('pair_coeff	%d	%d	lj/cut/coul/long	0.01	4.0\n' % tuple(sorted([ s_types[0], i+1])) )

	f.write(('\n'.join(["pair_coeff %d %d lj/cut/coul/long %f %f" % (atom_type_numbers[t], atom_type_numbers[t], t.vdw_e, t.vdw_r) for i,t in enumerate(atom_types) if i+1 not in pb_types+s_types])) + '\n')
	
	if o_types and h_types: f.write('pair_coeff	%d	%d	lj/cut/coul/long	0.0445	3.5\n' % tuple(sorted([o_types[0], h_types[0]])) )
	
	f.write('''
velocity	all create 300.0 1337
fix motion all nve
fix implicit_solvent all langevin 300.0 10.0 100.0 1337
timestep	2.0
neigh_modify check yes every 1 delay 0
dump	1 all xyz 1000 '''+run_name+'''.xyz
minimize 0.0 1.0e-8 1000 100000

fix average all ave/time 1 1000 1000 c_thermo_pe
thermo_style	custom f_average
thermo		1000

run		20000
''')
	f.close()
	os.system('/fs/home/jms875/build/lammps/lammps2/src/lmp_mycomp -in %s.in -log %s.log' % (run_name,run_name) )
	os.chdir('..')

params = [0.2]

for i in range(3):
	for j in range(3):
		for k in range(3):
			if i!=0 and j!=0: continue
			
			c = 100390.0
			d = 16.217
			n = 0.78734
			beta = 0.0000011
			lambda2 = 1.73222
			R = 2.65
			D = 0.15
			lambda1 = 2.4799
			A, B = 2000.0, 500.0
			
			params = params + [c, d, n, beta, lambda2, R, D, lambda1, A, B]

ligand_binding(params, 5)
sys.exit(0)

for n_ligands in range(-1,6):
	print '-'*30
	print n_ligands
	print '-'*30
	ligand_binding(params, n_ligands)

goal_dE = [-42.9145925228, -30.2515354714, -32.1963192387, -19.2409360413, 0.0]
E_by_ligand = []
for n_ligands in range(-1,6):
	energies = []
	for line in open('lammps/lig_'+str(n_ligands)+'.log'):
		try:
			a = line.split()
			if len(a)==1:
				energies.append( float(a[0]) )
		except ValueError: pass
	print n_ligands
	#energies = energies[len(energies)/2:-1]
	energies = [energies[-1]]
	E_by_ligand.append( sum(energies)/len(energies) )

lone_ligand_E = E_by_ligand[0]
dE = [E_by_ligand[i]-E_by_ligand[i-1]-lone_ligand_E for i in range(2,7)]

print lone_ligand_E
print ('%9.1f '*5) % tuple(dE)
print ('%9.1f '*5) % tuple(goal_dE)

'''
Fit:
	Pb-Pb: ligand-ligand energy
	Pb-S, Pb-Pb: lattice constants
	all: ligand removal energies
	Pb-O: ligand mobility
'''

