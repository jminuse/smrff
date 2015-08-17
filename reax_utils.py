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
	f.write(error_string+' Reactive MD_force field') 
	# Too long of a first line results in errors (ex: WARNING: number of globals in ffield file is 0!)
	
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

