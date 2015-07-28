import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

class Reax_params:
	def __init__(self):
		self.atom_types_c = '''Atom    r_s       valency   mass         r_vdw        epsilon  gamma    r_pi    valency_e    
		        alpha     gamma_w   valency_boc  p_ovun5      *****    chi      eta     p_hbond    
		        r_pi_pi   p_lp2     *****        b_o_131      b_o_132  b_o_133  *****   *****    
		        p_ovun2   p_val3    *****        valency_val  p_val5   rcore2   ecore2  acore2'''
		self.bond_types_c = '''A1  A2      De_s   De_p   De_pp  p_be1  p_bo5  v13cor  p_bo6  p_ovun1       
		              p_be2  p_bo3  p_bo4  *****  p_bo1  p_bo2   ovc    *****'''
		self.offdiags_c='''A1  A2      D   r_vdw  alpha   r_s   r_p   r_pp   *lgcij*'''
		self.three_body_c='''A1  A2  A3  theta_00   p_val1   p_val2   p_coa1   p_val7   p_pen1   p_val4'''
		self.torsional_c='''A1  A2  A3  A4   V1   V2   V3   p_tor1   p_cot1   *****   *****'''
		self.hydrogen_c = '''A1  A2  A3   r0_hb   p_hb1   p_hb2   p_hb3'''
		self.first_line_comment="Reactive MD-force field for ..."
		self.general_parameter_comments=[]
	
		self.atom_types_names=[line.split() for line in self.atom_types_c.splitlines()]
		self.bond_types_names=[line.split() for line in self.bond_types_c.splitlines()]
		self.offdiags_names=self.offdiags_c.split()
		self.three_body_names=self.three_body_c.split()
		self.torsional_names=self.torsional_c.split()
		self.hydrogen_names=self.hydrogen_c.split()
	
		self.general_parameters=[]
		self.atom_types=[]
		self.bonds=[]
		self.offdiags=[]
		self.thbps=[]
		self.torsional=[]
		self.hydrogen=[]
	
	def __repr__(self):
		return '\n'.join([str(s) for s in [self.general_parameters, self.atom_types, self.bonds, self.offdiags, self.thbps, self.torsional, self.hydrogen]])

def read_reax_file(filename):
	reax = Reax_params()
	f=open(filename)
	reax.first_line_comment = f.readline()
	
	reax.number_of_gen_params = int(re.search('\s*([0-9\.\-]+)\s+\!.*',f.readline()).group(1))
	for i in range(reax.number_of_gen_params):
		m=re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline())
		reax.general_parameters.append( float(  m.group(1) ) )
		reax.general_parameter_comments.append(m.group(2).strip())

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

def read_reax_include_file(filename,params):
	reax = Reax_params()
	f=open(filename)
	f.readline()
	bounds=[]
	default_bound_mults=(0.5,1.5)
	default_bounds_zero_value = (-0.1,0.1)

	reax.number_of_gen_params = int(re.search('\s*([0-9\.\-]+)\s+\!.*',f.readline()).group(1))
	for i in range(reax.number_of_gen_params):
		m=re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline())
		reax.general_parameters.append( float(  m.group(1) ) )

	reax.number_atoms = int(re.search('\s*([0-9\.\-]+)\s+\!(.*)',f.readline()).group(1))
	f.readline(); f.readline(); f.readline()
	for i in range(reax.number_atoms):
		reax.atom_types.append([])
		col = f.readline().split()
		reax.atom_types[i].append([col[0]])
		for index,item in enumerate(col[1:],1):
			if item == '0':
				reax.atom_types[i][0].append(0)
			else:
			 	m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b = (params.atom_types[i][0][index] * default_bound_mults[0],params.atom_types[i][0][index] * default_bound_mults[1])
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.atom_types[i][0].append(b)

		for line in [1,2,3]:
			col = f.readline().split()
			reax.atom_types[i].append([])
			for index,item in enumerate(col):
				if item == '0':
					reax.atom_types[i][line].append(0)
				else:
					m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
					if not m:
						b = (params.atom_types[i][line][index] * default_bound_mults[0],params.atom_types[i][line][index] * default_bound_mults[1])
					else:
						b = (float(m.group(1)),float(m.group(2)))
					b = tuple(sorted(b))
					if b == (0,0):
						b = default_bounds_zero_value
					bounds.append(b)
					reax.atom_types[i][line].append(b)

	reax.number_bonds = int(re.search('\s*([0-9]+)\s+\!(.*)',f.readline()).group(1))
	f.readline()
	for i in range(reax.number_bonds):
		reax.bonds.append([])
		col = f.readline().split()
		reax.bonds[i].append(col[:2])
		for index,item in enumerate(col[2:],2):
			if item == '0':
				reax.bonds[i][0].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b = (params.bonds[i][0][index] * default_bound_mults[0],params.bonds[i][0][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.bonds[i][0].append(b)

		col = f.readline().split()
		reax.bonds[i].append([])
		for index,item in enumerate(col):
			if item == '0':
				reax.bonds[i][1].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b =  (params.bonds[i][1][index] * default_bound_mults[0],params.bonds[i][1][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.bonds[i][1].append(b)

	reax.number_offdiags = int(re.search('\s*([0-9]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_offdiags):
		col=f.readline().split()
		reax.offdiags.append(col[:2])
		for index,item in enumerate(col[2:],2):
			if item == '0':
				reax.offdiags[i].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b =  (params.offdiags[i][index] * default_bound_mults[0],params.offdiags[i][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.offdiags[i].append(b)

	reax.number_threebody = int(re.search('\s*([0-9]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_threebody):
		col=f.readline().split()
		reax.thbps.append(col[:3])
		for index,item in enumerate(col[3:],3):
			if item == '0':
				reax.thbps[i].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b =  (params.thbps[i][index] * default_bound_mults[0],params.thbps[i][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.thbps[i].append(b)
		
	reax.number_torsional = int(re.search('\s*([0-9]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_torsional):
		col=f.readline().split()
		reax.torsional.append(col[:4])
		for index,item in enumerate(col[4:],4):
			if item == '0':
				reax.torsional[i].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b =  (params.torsional[i][index] * default_bound_mults[0],params.torsional[i][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.torsional[i].append(b)
	
	reax.number_hydrogen = int(re.search('\s*([0-9]+)\s+\!(.*)',f.readline()).group(1))
	for i in range(reax.number_hydrogen):
		col=f.readline().split()
		reax.hydrogen.append(col[:3])
		for index,item in enumerate(col[3:],3):
			if item == '0':
				reax.hydrogen[i].append(0)
			else:
				m=re.search('(\-?[0-9\.]+)\,(\-?[0-9\.]+)',item)
				if not m:
					b =  (params.hydrogen[i][index] * default_bound_mults[0],params.hydrogen[i][index] * default_bound_mults[1]) 
				else:
					b = (float(m.group(1)),float(m.group(2)))
				b = tuple(sorted(b))
				if b == (0,0):
					b = default_bounds_zero_value
				bounds.append(b)
				reax.hydrogen[i].append(b)
	f.close()
	return reax, bounds

def trim_spaces(num,space_count=8,i=False):
	if not i: 
		s='%.4f' % num
	else:
		s='%d' % num
	space = ' ' * (space_count - len(s))
	return space+s

def spaced_numbered_list(numlist):
	return ' '.join(map(trim_spaces, numlist))

def write_reax_file(dataset, best=False,error=None):
	if best:
		f = open(dataset.name+'_best.reax', 'w')
	else:
		f = open(dataset.name+'.reax', 'w')

	# For readibility
	delim = ' !'
	s = "     "
	rp=dataset.reax_params
	if error:
		error_string='Error: ' + str(error) + '  '
	else:
		error_string=''
	f.write(error_string+rp.first_line_comment)
	
	# Print General Parameters:
	f.write(' ' + str(rp.number_of_gen_params) + '       ! ' + 'Number of general parameters  \n')
	for i in range(rp.number_of_gen_params):
		f.write(trim_spaces(rp.general_parameters[i],10) + delim + rp.general_parameter_comments[i] +'  \n')

	# Print atom types:
	f.write(trim_spaces(rp.number_atoms,3,1) + '    ! ' + rp.atom_types_c + '  \n')
	for i in range(rp.number_atoms):
		s=' '*2
		if len(rp.atom_types[i][0][0]) == 2:  s = ' '
		f.write(' ' +rp.atom_types[i][0][0] + s + spaced_numbered_list(rp.atom_types[i][0][1:]) +'  \n')
		f.write(' '*14 + spaced_numbered_list(rp.atom_types[i][1]) +'  \n')
		f.write(' '*14 + spaced_numbered_list(rp.atom_types[i][2]) +'  \n')
		f.write(' '*14 + spaced_numbered_list(rp.atom_types[i][3]) +'  \n')

	# Print bond types:
	f.write(trim_spaces(rp.number_bonds,3,1) + '    ! ' + rp.bond_types_c + '  \n')
	for i in range(rp.number_bonds):
		s=' '*2
		f.write(' '*2+rp.bonds[i][0][0] + ' '*2  + rp.bonds[i][0][1] + ' '*12 + spaced_numbered_list(rp.bonds[i][0][2:]) +'  \n')
		f.write(' '*20 + spaced_numbered_list(rp.bonds[i][1]) + '  \n')
	
	# Print Off-Diagonal Terms:
	f.write(trim_spaces(rp.number_offdiags,3,1) + '    ! ' + rp.offdiags_c + '  \n')
	for i in range(rp.number_offdiags):
		f.write(' '*2 + rp.offdiags[i][0] + ' '*2 +  rp.offdiags[i][1] + ' '*10 + spaced_numbered_list(rp.offdiags[i][2:]) + '  \n')

	# Print Three-Body Parameters:
	f.write(trim_spaces(rp.number_threebody,3,1) + '    ! ' + rp.three_body_c + '  \n')
	for i in range(rp.number_threebody):
		f.write(' '*2 + rp.thbps[i][0] + ' '*2 +  rp.thbps[i][1] + ' '*2 +  rp.thbps[i][2] + ' '*10 + spaced_numbered_list(rp.thbps[i][3:]) + '  \n')

	# Print Torsional terms:
	# if rp.number_torsional:
	f.write(trim_spaces(rp.number_torsional,3,1) + '    ! ' + rp.torsional_c + '  \n')
	for i in range(rp.number_torsional):
		f.write(' '*2 + rp.torsional[i][0] + ' '*2 +  rp.torsional[i][1] + ' '*2 +  rp.torsional[i][2] +' '*2 +  rp.torsional[i][3] + ' ' + spaced_numbered_list(rp.torsional[i][4:]) + '  \n')

	# Print Hydrogen Bonds:
	# if rp.number_hydrogen:
	f.write(trim_spaces(rp.number_hydrogen,3,1) + '    ! ' + rp.hydrogen_c + '  \n')
	for i in range(rp.number_hydrogen):
		f.write(' '*2 + rp.hydrogen[i][0] + ' '*2 +  rp.hydrogen[i][1] + ' '*2 +  rp.hydrogen[i][2] + ' ' + spaced_numbered_list(rp.hydrogen[i][3:]) + '  \n')

	f.close()


def calculate_error(dataset):
	write_reax_file(dataset) #all LAMMPS job use same reax file and same data files
	
	#run LAMMPS
	for s in dataset.systems:
		if True: #hide screen
			lmp = lammps('',['-log',dataset.name+'.log','-screen','none'])
		else: #show screen
			lmp = lammps('',['-log',dataset.name+'.log'])

		commands = ('''units real
atom_style full
pair_style reax/c NULL
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary f f f
read_data	'''+s.name+'''.data

compute atom_pe all pe/atom
compute		test_pe all reduce sum c_atom_pe
thermo_style custom pe c_test_pe
pair_coeff * * '''+dataset.name+'''.reax Pb Cl
fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
run 1''').splitlines()

		for line in commands:
			lmp.command(line)
		
		lammps_energies_by_atom = lmp.extract_compute('atom_pe',1,1) #http://lammps.sandia.gov/doc/Section_python.html
		s.lammps_energy = sum( [lammps_energies_by_atom[i] for i in range(0,len(s.atoms)) ] )
		lammps_forces = lmp.extract_atom('f',3)
		for i,a in enumerate(s.atoms):
			a.lfx, a.lfy, a.lfz = lammps_forces[i][0], lammps_forces[i][1], lammps_forces[i][2]
		lmp.close()
	
	#calculate energy error
	relative_energy_error, absolute_energy_error = 0.0, 0.0
	relative_force_error, absolute_force_error = 0.0, 0.0
	for elements,systems in dataset.by_elements.iteritems():
		baseline_energy = systems[0].lammps_energy #should be in order of increasing energy
		for s in systems:
			s.lammps_energy -= baseline_energy
			try:
				relative_energy_error += ( (s.lammps_energy-s.energy)/(s.energy+1.0) )**2
				absolute_energy_error += (s.lammps_energy-s.energy)**2
			except OverflowError: return 1e10
			#print s.energy, s.lammps_energy

			for i,a in enumerate(s.atoms):
				fx, fy, fz = a.lfx, a.lfy, a.lfz
				real_force_squared = a.fx**2 + a.fy**2 + a.fz**2
				try:
					relative_force_error += ((fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2) / (real_force_squared + 20.0**2)
					absolute_force_error += (fx-a.fx)**2 + (fy-a.fy)**2 + (fz-a.fz)**2
				except OverflowError:
					return 1e10
	
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

def run(run_name, other_run_names=[]):
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
	#filenames = filenames[:5]
	#random.seed(10)
	#random.shuffle(filenames)
	for name in filenames:
				if not name.startswith('PbCl2'): continue
				if not name.endswith('_def2SVP'): continue
				energy, atoms = g09.parse_atoms(name, check_convergence=True)
				if len(atoms)!=3: continue
				system = utils.System(box_size=[40, 40, 40], name=name)
				total = utils.Molecule('gaussian/'+name, extra_parameters=extra, check_charges=False)
				system.energy = energy*627.509 #convert energy from Hartree to kcal/mol
				system.element_string = ' '.join( [a.element for a in total.atoms] )
				print 'Read', system.element_string, 'from', name
				for i,a in enumerate(total.atoms):
					b = atoms[i]
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

	os.chdir('lammps')
	for s in dataset.systems:
		files.write_lammps_data(s)

	dataset.reax_params = read_reax_file('../input.reax')
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
	
	#optimize
	import numpy
	from scipy.optimize import minimize, fmin_l_bfgs_b

	def stochastic(use_gradient=True):
		best_min = utils.Struct(fun=calculate_error_from_list(initial_params),x=initial_params)
		print dataset.name, 'starting error: %.4g' % best_min.fun
		#exit() #just print starting error
		def new_param_guess(start, gauss=True):
			dataset.how_long_since_checked_others += 1
			if dataset.how_long_since_checked_others > 10:
				check_others(dataset)
			
			while True: #keep going until new params are generated
				params = []
				for p,b in zip(start, bounds):
					if gauss: #Gaussian random
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
								distances.append( sum([ ( (x-y)/(b[1]-b[0]) )**2 for x,y,b in zip(others.list, params, bounds)])**0.5 )
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
					print 'Bounds exceeded:', oldx, newx, bounds[i][0], bounds[i][1]
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
				params = new_param_guess(best_min.x, gauss=False)
				start_error = calculate_error_from_list(params)
				#while start_error > 2.0:
				#	params = new_param_guess(best_min.x, gauss=True)
				#	start_error = calculate_error_from_list(params)
				x, fun, stats = fmin_l_bfgs_b(calculate_error_from_list, params, fprime=error_gradient, factr=1e8, bounds=bounds)
				guess = utils.Struct(fun=fun,x=x)
				print dataset.name, 'error', start_error, guess.fun, best_min.fun
			else: #non-gradient optimization
				params = new_param_guess(best_min.x, gauss=False)
				guess = utils.Struct(fun=calculate_error_from_list(params),x=params)
				#print dataset.name, 'error', guess.fun, best_min.fun
			if guess.fun < best_min.fun:
				best_min = guess
				unpack_params(best_min.x, dataset)
				write_reax_file(dataset,best=True,error = best_min.fun)
				print dataset.name, 'new best error = %.4g' % best_min.fun

		return best_min

	stochastic(False)

from multiprocessing import Process, Queue
def run_multiple(jobname, N):
	queue = Queue()
	for i in range(N):
		p = Process(target=run, args=(jobname+str(i), [jobname+str(other) for other in range(N) if other!=i]))
		p.start()

run('test')
#run_multiple('test', 8)

