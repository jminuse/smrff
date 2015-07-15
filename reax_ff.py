import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

class Reax_params:
	atom_types_c = '''Atom    r_s       valency   mass         r_vdw        epsilon  gamma    r_pi    valency_e    
            alpha     gamma_w   valency_boc  p_ovun5      *****    chi      eta     p_hbond    
            r_pi_pi   p_lp2     *****        b_o_131      b_o_132  b_o_133  *****   *****    
            p_ovun2   p_val3    *****        valency_val  p_val5   rcore2   ecore2  acore2'''
	bond_types_c = '''A1  A2      De_s   De_p   De_pp  p_be1  p_bo5  v13cor  p_bo6  p_ovun1       
                  p_be2  p_bo3  p_bo4  *****  p_bo1  p_bo2   ovc    *****'''
	offdiags_c='''A1  A2      D   r_vdw  alpha   r_s   r_p   r_pp   *lgcij*'''
	three_body_c='''A1  A2  A3  theta_00   p_val1   p_val2   p_coa1   p_val7   p_pen1   p_val4'''
	torsional_c='''A1  A2  A3  A4   V1   V2   V3   p_tor1   p_cot1   *****   *****'''
	hydrogen_c = '''A1  A2  A3   r0_hb   p_hb1   p_hb2   p_hb3'''
	first_line_comment="Reactive MD-force field for ..."
	
	atom_types_names=[line.split() for line in atom_types_c.splitlines()]
	bond_types_names=[line.split() for line in bond_types_c.splitlines()]
	offdiags_names=offdiags_c.split()
	three_body_names=three_body_c.split()
	torsional_names=torsional_c.split()
	hydrogen_names=hydrogen_c.split()
	gen_p_comments=[]
	
	gen_p=[]
	atom_types=[]
	bonds=[]
	offdiags=[]
	thbps=[]
	torsional=[]
	hydrogen=[]


def read_reax_file(filename):
	reax = Reax_params()
	f=open(filename)
	reax.first_line_comment = f.readline()
	
	reax.number_of_gen_params = int(re.search('\s*([0-9\.\-]+)\s+\!.*',f.readline()).group(1))
	for i in range(reax.number_of_gen_params):
		m=re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline())
		reax.gen_p.append( float(  m.group(1) ) )
		reax.gen_p_comments.append(m.group(2).strip())

	reax.number_atoms = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	f.readline(); f.readline(); f.readline()
	for i in range(reax.number_atoms):
		reax.atom_types.append([])
		col = f.readline().split()
		reax.atom_types[i].append([col[0]] + [float(s) for s in col[1:]] )

		col = f.readline().split()
		reax.atom_types[i].append([float(s) for s in col])
		
		col = f.readline().split()
		reax.atom_types[i].append([float(s) for s in col])

		col = f.readline().split()
		reax.atom_types[i].append([float(s) for s in col])

	reax.number_bonds = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	f.readline()
	for i in range(reax.number_bonds):
		reax.bonds.append([])
		col = f.readline().split()
		reax.bonds[i].append(col[:2] + [float(s) for s in col[2:]])

		col = f.readline().split()
		reax.bonds[i].append([float(s) for s in col])

	reax.number_offdiags = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_offdiags):
		col=f.readline().split()
		reax.offdiags.append(col[:2] + [float(s) for s in col[2:]])

	reax.number_threebody = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_threebody):
		col=f.readline().split()
		reax.thbps.append(col[:3] + [float(s) for s in col[3:]])
		
	reax.number_torsional = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_torsional):
		col=f.readline().split()
		reax.torsional.append(col[:4] + [float(s) for s in col[4:]])
	
	reax.number_hydrogen = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_hydrogen):
		col=f.readline().split()
		reax.hydrogen.append(col[:3] + [float(s) for s in col[3:]])
	
	return reax
def trim_spaces(num,space_count=8,i=False):
	if not i: 
		s='%.4f' % num
	else:
		s='%d' % num
	space = ' ' * (space_count - len(s))
	return space+s

def spaced_numbered_list(numlist):
	return ' '.join(map(trim_spaces, numlist))

def write_reax_file(system, best=False):
	if best:
		f = open(system.name+'_best.reax', 'w')
	else:
		f = open(system.name+'.reax', 'w')

	# For readibility
	delim = ' !'
	s = "     "
	rp=system.reax_params

	f.write(rp.first_line_comment)
	
	# Print General Parameters:
	f.write(' ' + str(rp.number_of_gen_params) + '       ! ' + 'Number of general parameters  \n')
	for i in range(rp.number_of_gen_params):
		f.write(trim_spaces(rp.gen_p[i],10) + delim + rp.gen_p_comments[i] +'  \n')

	# Print atom types:
	f.write(trim_spaces(rp.number_atoms,3,1) + '    ! ' + rp.atom_types_c + '  \n')
	for i in range(rp.number_atoms):
		s=' '*2
		if len(rp.atom_types[i][0][0]) == 2:  s = ' '
		f.write(' ' +rp.atom_types[i][0][0] + s + spaced_numbered_list(rp.atom_types[i][0][1:]) +'  \n')
		f.write(' '*4 + spaced_numbered_list(rp.atom_types[i][1]) +'  \n')
		f.write(' '*4 + spaced_numbered_list(rp.atom_types[i][2]) +'  \n')
		f.write(' '*4 + spaced_numbered_list(rp.atom_types[i][3]) +'  \n')

	# Print bond types:
	f.write(trim_spaces(rp.number_bonds,3,1) + '    ! ' + rp.bond_types_c + '  \n')
	for i in range(rp.number_bonds):
		s=' '*2
		f.write(' '*2+rp.bonds[i][0][0] + ' '*2  + rp.bonds[i][0][1] + ' ' + spaced_numbered_list(rp.bonds[i][0][2:]) +'  \n')
		f.write(' '*7 + spaced_numbered_list(rp.bonds[i][1]) + '  \n')
	
	# Print Off-Diagonal Terms:
	f.write(trim_spaces(rp.number_offdiags,3,1) + '    ! ' + rp.offdiags_c + '  \n')
	for i in range(rp.number_offdiags):
		f.write(' '*2 + rp.offdiags[i][0] + ' '*2 +  rp.offdiags[i][1] + ' ' + spaced_numbered_list(rp.offdiags[i][2:]) + '  \n')

	# Print Three-Body Parameters:
	f.write(trim_spaces(rp.number_threebody,3,1) + '    ! ' + rp.three_body_c + '  \n')
	for i in range(rp.number_threebody):
		f.write(' '*2 + rp.thbps[i][0] + ' '*2 +  rp.thbps[i][1] + ' '*2 +  rp.thbps[i][2] + ' ' + spaced_numbered_list(rp.thbps[i][3:]) + '  \n')

	# Print Torsional terms:
	if rp.number_torsional:
		f.write(trim_spaces(rp.number_torsional,3,1) + '    ! ' + rp.torsional_c + '  \n')
		for i in range(rp.number_torsional):
			f.write(' '*2 + rp.torsional[i][0] + ' '*2 +  rp.torsional[i][1] + ' '*2 +  rp.torsional[i][2] +' '*2 +  rp.torsional[i][3] + ' ' + spaced_numbered_list(rp.torsional[i][4:]) + '  \n')

	# Print Hydrogen Bonds:
	if rp.number_hydrogen:
		f.write(trim_spaces(rp.number_hydrogen,3,1) + '    ! ' + rp.hydrogen_c + '  \n')
		for i in range(rp.number_hydrogen):
			f.write(' '*2 + rp.hydrogen[i][0] + ' '*2 +  rp.hydrogen[i][1] + ' '*2 +  rp.hydrogen[i][2] + ' ' + spaced_numbered_list(rp.hydrogen[i][3:]) + '  \n')

	f.close()


def set_lammps_parameters(system):
	write_reax_file(system)
	lmp.command('unfix 1')
	lmp.command('pair_coeff * * '+system.name+'.reax Pb I '+(' NULL'*(len(system.atom_types)-2))) #is it possible to do this with the LAMMPS set command, to avoid writing the file to disk?
	lmp.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c')
	
	for t in system.atom_types:
		pass
		#if hasattr(t,'vdw_e'):
		#	lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
		#	lmp.command('pair_coeff %d * lj/cut/coul/cut %f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
	for t in system.bond_types:
		lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
	for t in system.angle_types:
		lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
	for t in system.dihedral_types:
		lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

def calculate_error(system):
	
	#run LAMMPS
	set_lammps_parameters(system)
	lmp.command('run 1')
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
	
	#print 'Error components:', absolute_energy_error, absolute_force_error, relative_energy_error, relative_force_error
	
	return error

def pack_params(system):
	params, names = [], []
	for t in system.bond_types:
		params += [t.e, t.r]
	for t in system.angle_types:
		params += [t.e, t.angle]
	for t in system.dihedral_types:
		if len(t.e)==3:
			t.e = list(t.e)+[0.0]
		params += list(t.e)
	
	for atom in system.reax_params.atom_types:
		names_list = system.reax_params.atom_types_names
		atom_name = atom[0][0]
		if atom_name == 'Pb' or atom_name == 'I':
			for i,b in enumerate(include0):
				if b:
					params.append(atom[0][i])
					names.append(atom_name + '.' + names_list[0][i])
			
			for i,b in enumerate(include1):
				if b:
					params.append(atom[1][i])
					names.append(atom_name + '.' + names_list[1][i])

			for i,b in enumerate(include2):
				if b:
					params.append(atom[2][i])
					names.append(atom_name + '.' + names_list[2][i])

			for i,b in enumerate(include3):
				if b:
					params.append(atom[3][i])
					names.append(atom_name + '.' + names_list[3][i])

	for bond in system.reax_params.bonds:
		names_list = system.reax_params.bond_types_names
		bondname='tbp(' + bond[0][0] + ',' + bond[0][1] + ').'
		
		if (bond[0][0],bond[0][1]) == ('1','2'): #only include Pb-I bonds
			for i,b in enumerate(include4):
				if b:
					params.append(bond[0][i])
					names.append(bondname+names_list[0][i])

			for i,b in enumerate(include5):
				if b:
					params.append(bond[1][i])
					names.append(bondname+names_list[1][i])

	for offdiag in system.reax_params.offdiags:
		names_list=system.reax_params.offdiags_names
		offdname='offd('+offdiag[0]+','+offdiag[1]+').'
		for i,b in enumerate(include6):
			if b:
				params.append(offdiag[i])
				names.append(offdname+names_list[i])

	for thbp in system.reax_params.thbps:
		names_list = system.reax_params.three_body_names
		thbpname='thbp('+thbp[0] +',' + thbp[1]+',' + thbp[2]+')'
		for i,b in enumerate(include7):
			if b:
				params.append(thbp[i])
				names.append(thbpname + names_list[i])

	if len(params)!=len(names):
		print 'There are %d parameters, but %d names!' % (len(params), len(names))
		raise SystemExit
	print names
	return params, names

def unpack_params(params, system):
	i = 0
	# for t in system.atom_types:
	# 	if False:
	# 		t.vdw_e, t.vdw_r = params[i], params[i+1]
	# 		i += 2
	# 		if t.element==53:
	# 			pass #defined by Pb
	# 		else:
	# 			t.charge = params[i]
	# 			i += 1
	for t in system.bond_types:
		t.e, t.r = params[i], params[i+1]
		i += 2
	for t in system.angle_types:
		t.e, t.angle = params[i], params[i+1]
		i += 2
	for t in system.dihedral_types:
		t.e = tuple(params[i:i+4])
		i += 4
	for atom_iter in range( len(system.reax_params.atom_types) ):
		if system.reax_params.atom_types[atom_iter][0][0]=='Pb' or \
		   system.reax_params.atom_types[atom_iter][0][0]=='I':
			for param_iter,b in enumerate(include0):
				if b:
					system.reax_params.atom_types[atom_iter][0][param_iter] = params[i]
					i += 1
			for param_iter,b in enumerate(include1):
				if b:
					system.reax_params.atom_types[atom_iter][1][param_iter] = params[i]
					i += 1
			for param_iter,b in enumerate(include2):
				if b:
					system.reax_params.atom_types[atom_iter][2][param_iter] = params[i]
					i += 1
			for param_iter,b in enumerate(include3):
				if b:
					system.reax_params.atom_types[atom_iter][3][param_iter] = params[i]
					i += 1
	for bond_iter in range(system.reax_params.number_bonds):
		for param_iter,b in enumerate(include4):
			if b:
				system.reax_params.bonds[bond_iter][0][param_iter] = params[i]
				i += 1
		for param_iter,b in enumerate(include5):
			if b:
				system.reax_params.bonds[bond_iter][1][param_iter] = params[i]
				i += 1

	for offdiag_iter in range(system.reax_params.number_offdiags):
		for param_iter,b in enumerate(include6):
			if b:
				system.reax_params.offdiags[offdiag_iter][param_iter] = params[i]
				i += 1

	for thbp_iter in range(system.reax_params.number_threebody):
	    for param_iter,b in enumerate(include7):
	        if b:
				system.reax_params.thbps[thbp_iter][param_iter] = params[i]
				i += 1


I_ = 66
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838

extra = {
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=0.1, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=0.0, vdw_e=0.1, vdw_r=3.0),
}

system = utils.System(box_size=[100, 100, 100], name='test_reax')

for root, dirs, file_list in os.walk("gaussian"):
	count = 0
	for ff in file_list:
		if ff.endswith('.log'):
			name = ff[:-4]
	#for step in range(20):
	#		name = 'PbI2_r%d' % step
			if not name.endswith('_def2SVP'): continue
			energy, atoms = g09.parse_atoms(name, check_convergence=False)
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
pair_style reax/c NULL
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary f f f
read_data	'''+system.name+'''.data

compute atom_pe all pe/atom
compute		test_pe all reduce sum c_atom_pe
thermo_style custom pe c_test_pe
pair_coeff * * ../input.reax Pb I '''+(' NULL'*(len(system.atom_types)-2))+'''
fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
''').splitlines()
lmp = lammps('',['-log',system.name+'.log','-screen','none'])
for line in commands:
	lmp.command(line)

system.reax_params = read_reax_file('../input.reax')

def calculate_error_from_list(params):
	unpack_params(params, system)
	
	for t in system.atom_types + system.bond_types + system.angle_types + system.dihedral_types:
		t.written_to_lammps = False
	
	error = calculate_error(system)
	
	return error

# Atom Parameters:
#          ['Atom','r_s','valency','mass','r_vdw','epsilon','gamma','r_pi','valency_e']
include0 = [     0,    1,        0,     0,      1,        1,      1,     1,          0]
#          ['alpha','gamma_w','valency_boc','p_ovun5','*****','chi','eta','p_hbond']
include1 = [ 1,      1,        0,            0,        0,      1,    1,    0]
#          ['r_pi_pi','p_lp2','*****','b_o_131','b_o_132','b_o_133','*****','*****']
include2 = [ 1,        1,      0,      1,        1,        1,        0,      0]
#          ['p_ovun2','p_val3','*****','valency_val','p_val5','rcore2','ecore2','acore2']
include3 = [ 1,        1,       0,      1,            1,       0,       0,       0]

# Bond Parameters:
#          ['A1','A2','De_s','De_p','De_pp','p_be1','p_bo5','v13cor','p_bo6','p_ovun1']
include4 = [ 0,   0,   1,     1,     1,      1,      1,      1,       1,      1]
#          ['p_be2','p_bo3','p_bo4','*****','p_bo1','p_bo2','ovc','*****']
include5 = [ 1,      1,      1,      0,      1,      1,      1,    0]

# Off-Diagonal Parameters:
#          ['A1','A2','D','r_vdw','alpha','r_s','r_p','r_pp','*lgcij*']
include6 = [ 0,   0,   1,  1,      1,      1,    1,    1,     0]

# Three-Body Parameters:
#          ['A1','A2','A3','theta_00','p_val1','p_val2','p_coa1','p_val7','p_pen1','p_val4']
include7 = [ 0,   0,   0,   1,         1,       1,       1,       1,       1,       1]

#r_s (= r_sigma), r_pi: bond order distances for sigma and pi bonding
#p_bo1-5 = exponential coefficients for bond order calculation
#p_coa parameters = pi-conjugation

initial_params, names = pack_params(system)
bounds = [ tuple(sorted([x*0.5, x*1.5])) if x!=0.0 else (-0.1,0.1) for x in initial_params ]

import numpy
from scipy.optimize import minimize, fmin_l_bfgs_b

def stochastic(use_gradient=True):
	best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
	print 'Error: %.4g' % best_min.fun
	
	def new_param_guess(start):
		while True: #keep going until new params are generated
			params = []
			for p,b in zip(start, bounds):
				new = random.gauss(p, abs(p)*0.5 + 0.01*(b[1]-b[0]) ) if random.random()<0.2 else p
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
			if bounds[i]==(0.0,0.0):
				gradient.append(0.0)
				continue
			oldx = x[i]
			sign_oldx = -1 if oldx<0 else 1
			newx = oldx + 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
			if newx>bounds[i][1] or newx<bounds[i][0]:
				newx = oldx - 0.0001*(bounds[i][1]-bounds[i][0])*sign_oldx
			if newx>bounds[i][1] or newx<bounds[i][0]:
				print oldx, newx, bounds[i][0], bounds[i][1]
				exit()
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
		params = new_param_guess(best_min.x)
		if use_gradient:
			x, fun, stats = fmin_l_bfgs_b(calculate_error_from_list, params, fprime=error_gradient, bounds=bounds)
			guess = utils.Struct(fun=fun,x=x)
			print 'Error', calculate_error_from_list(params), guess.fun, best_min.fun
		else:
			guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
			print 'Error', guess.fun, best_min.fun
		if guess.fun < best_min.fun:
			best_min = guess
			unpack_params(best_min.x, system)
			write_reax_file(system,best=True)
			print 'New best error = %.4g' % best_min.fun
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
		'''
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
		if not any(changed):
			print names[param_index], 'has no effect'
		'''
		if best_min.x[param_index]==round(best_min.x[param_index]):
			print names[param_index], 'is an integer'
	return best_min

#try_params()

stochastic(False)

