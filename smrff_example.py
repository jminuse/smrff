import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09
from reax_utils import *
from utils import elements_by_atomic_number

def run(run_name):
	print 'starting ' + run_name
	Cl_ = 66
	H_ = 54
	N_ = 53
	Pb_ = 111

	Pb = 907
	Cl = 838

	extra = {
		Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=2.0, vdw_e=0.1, vdw_r=2.3),
		Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-1.0, vdw_e=0.1, vdw_r=2.3),
	}
	factor = 1
	system = utils.System(box_size=(20* factor, 20* factor, 20* factor), name=run_name)
	
	# acet = utils.Molecule('cml/acetone')
	# dmf = utils.Molecule('cml/DMF')
	# mai = utils.Molecule('cml/MAI')
	# macl = utils.Molecule('cml/MACl', extra_parameters=extra, check_charges=False)
	# pbcl2 = utils.Molecule('cml/PbCl2', extra_parameters=extra)
	# tiled_pbcl2_mai = utils.Molecule('cml/tiledCl.cml', extra_parameters=extra)
	tiled_pbcl2_macl = utils.Molecule('cml/tiledCl_MACl.cml', extra_parameters=extra,check_charges=False)
	# pbcl_clump = utils.Molecule('cml/PbCl_clump.cml', extra_parameters=extra,check_charges=False)
	
	# pure_densities = (5.85, 0.944, 0.944)
	# ratio = (1,1,10)
	# density = sum([d*r for d,r in zip(pure_densities,ratio)])/sum(ratio)
	# print 'Estimated density =', density
	# files.packmol(system, (pbcl2, mai, dmf), ratio, density)
	# system.add(tiled_pbcl2_mai)
	system.add(tiled_pbcl2_macl)
	# system.add(pbcl_clump)

	# ## --TRYING THE SYSTEM WITH ONLY Pb & Cl
	# atoms=system.atoms
	# system.atoms=[]
	# for i in range(len(atoms)):
	# 	if atoms[i].type.element_name=='Pb' or atoms[i].type.element_name=='Cl':
	# 		system.atoms.append(atoms[i])

	for a in system.atoms:
		a.x*=factor
		a.y*=factor
		a.z*=factor
	# files.write_xyz(system.atoms)
	# exit()
	
	if False: #hide screen
		lmp = lammps('',['-log','none','-screen','none'])
	else: #show screen
		lmp = lammps('',['-log',run_name+'.log'])
	
	os.chdir('lammps')
	files.write_lammps_data(system)
	
	input_file='../input.reax'
	# input_file='testcl_new2_best.reax'
	# input_file='../new_cl_input.reax'
	# input_file='testcl1_best.reax'
	# input_file='../test_smrff.reax'
	# input_file = '../Nov_4_best_mod.reax'
	input_file = '../Nov_9_best.reax'
	# input_file='forces_correction_24_included_1_best.reax'
	# input_file='Nov_2_PbCl_24_minimization2_best.reax'

	[lmp.command(line) for line in ('''units real
atom_style full
pair_style hybrid lj/cut/coul/dsf 0.05 10.0 10.0 reax/c NULL safezone 5.0 mincap 1000
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5 

boundary p p p
read_data	'''+run_name+'''.data''').splitlines()]

	for t1 in system.atom_types:
		for t2 in system.atom_types:
			if t1.lammps_type<=t2.lammps_type and (t1.lammps_type,t2.lammps_type) not in [(1,1),(1,2),(2,2)]:
			# if t1.lammps_type<=t2.lammps_type and (t1.lammps_type,t2.lammps_type) not in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
				vdw_e = (t1.vdw_e*t2.vdw_e)**0.5
				vdw_r = (t1.vdw_r*t2.vdw_r)**0.5
				pair_elements = [utils.elements_by_atomic_number[i] for i in sorted([t1.element,t2.element])]
				if pair_elements==['O','Pb']:
					vdw_e = 10.0
					vdw_r = 2.7
				if pair_elements==['O','Cl']:
					vdw_e = 0.01
					vdw_r = 3.2
				#if pair_elements==['O','O']: #modifies original OPLS parameters
				#	vdw_e = 0.01
				#	vdw_r = 4.5
				lmp.command('pair_coeff %d %d lj/cut/coul/dsf %f %f' % (t1.lammps_type, t2.lammps_type, vdw_e, vdw_r))

	# Forming a string in order of the atoms, element names for REAX atoms, NULL otherwise
	a_types = []
	for t in system.atom_types:
		if t.element == 82 or t.element == 17 or t.index == 233:
			a_types+=[elements_by_atomic_number[t.element]]
		else:
			a_types+=['NULL']
		print "Charge on " + t.element_name + " is " + str(t.charge)
	a_types = ' '.join(a_types)
	print '''pair_coeff * * reax/c ''' + input_file + ' ' + a_types

	# Printing the charges of the system types
	charges = dict()
	for a in system.atoms:
		if not charges.has_key(a.type.element_name):
			charges[a.type.element_name]=0
		charges[a.type.element_name]+=a.type.charge;
	print charges



	[lmp.command(line) for line in ('''
pair_coeff * * reax/c ''' + input_file + ' ' + a_types + '''
# pair_coeff * * reax/c '''+input_file+''' Pb Cl '''+(' NULL'*(len(system.atom_types)-2))+'''

group qeq_atoms type 1 2 3
fix 1 qeq_atoms qeq/reax 1 0.0 10.0 1.0e-6 reax/c''').splitlines()]

	commands = ('''dump 1 all xyz 100 '''+run_name+'''.xyz
thermo 100

compute reax all pair reax/c

fix av all ave/time 1 100 100 c_thermo_pe
group pbs type 1
group cls type 2
group hs type 6
group other type 3


compute pbs_pe pbs pe/atom
compute cls_pe cls pe/atom
compute hs_pe hs pe/atom
compute other_pe other pe/atom

compute charges all property/atom q
compute sum_charges all reduce sum c_charges

compute	sum_pbs_pe pbs reduce sum c_pbs_pe
compute	sum_cls_pe cls reduce sum c_cls_pe
compute	sum_hs_pe hs reduce sum c_hs_pe
compute	sum_other_pe other reduce sum c_other_pe

thermo_style custom step temp f_av density tpcpu c_reax[1] c_reax[5] ecoul c_sum_pbs_pe c_sum_cls_pe c_sum_hs_pe c_sum_other_pe emol epair c_sum_charges
minimize 0.0 1.0e-8 1000 100000
velocity all create 300.0 1337 rot yes dist gaussian
fix motion all npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0
#fix motion all nvt temp 600.0 300.0 100.0
timestep 0.5
run 10000
write_restart '''+run_name+'''.restart''').splitlines()

	for line in commands:
		lmp.command(line)
	print run_name + ' has finished.'
run('Jan_27_testing')