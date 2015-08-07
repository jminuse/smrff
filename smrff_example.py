import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09
from reax_utils import *

def run(run_name):
	Cl_ = 66
	H_ = 54
	N_ = 53
	Pb_ = 111

	Pb = 907
	Cl = 838

	extra = {
		Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=0.01, vdw_r=4.0),
		Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.01, vdw_r=4.0),
	}

	system = utils.System(box_size=(30, 30, 30), name=run_name)
	
	acet = utils.Molecule('cml/acetone')
	dmf = utils.Molecule('cml/DMF')
	macl = utils.Molecule('cml/MACl', extra_parameters=extra, check_charges=False)
	pbcl2 = utils.Molecule('cml/PbCl2', extra_parameters=extra)
	
	#for m in [dmf, macl, pbcl2]:
	#	print [a.type.element for a in m.atoms]
	#exit()
	
	pure_densities = (5.85, 0.944, 0.944)
	ratio = (1,1,10)
	density = sum([d*r for d,r in zip(pure_densities,ratio)])/sum(ratio)
	print 'Estimated density =', density
	
	#density = 0.2
	#ratio = (10,1,1)
	#files.packmol(system, (pbcl2, macl, dmf), ratio, density)
	files.packmol(system, (dmf,), (1,), 0.8)
	
	system.add(pbcl2)
	
	if False: #hide screen
		lmp = lammps('',['-log','none','-screen','none'])
	else: #show screen
		lmp = lammps('',['-log',run_name+'.log'])
	
	os.chdir('lammps')
	files.write_lammps_data(system)
	
	input_file='../input.reax'
	
	[lmp.command(line) for line in ('''units real
atom_style full
pair_style hybrid lj/cut/coul/dsf 0.05 10.0 10.0 reax/c NULL
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data	'''+run_name+'''.data''').splitlines()]

	for t1 in system.atom_types:
		for t2 in system.atom_types:
			if t1.lammps_type<=t2.lammps_type and (t1.lammps_type,t2.lammps_type) not in [(1,1),(1,2),(2,2)]:
				vdw_e = (t1.vdw_e*t2.vdw_e)**0.5
				vdw_r = (t1.vdw_r*t2.vdw_r)**0.5
				pair_elements = [utils.elements_by_atomic_number[i] for i in sorted([t1.element,t2.element])]
				print pair_elements
				if pair_elements==['O','Pb']:
					vdw_e = 10.0
					vdw_r = 2.7
				if pair_elements==['O','Cl']:
					vdw_e = 0.01
					vdw_r = 4.2
				#if pair_elements==['O','O']: #modifies original OPLS parameters
				#	vdw_e = 0.01
				#	vdw_r = 4.5
				lmp.command('pair_coeff %d %d lj/cut/coul/dsf %f %f' % (t1.lammps_type, t2.lammps_type, vdw_e, vdw_r))
	
	[lmp.command(line) for line in ('''
pair_coeff * * reax/c '''+input_file+''' Pb Cl'''+(' NULL'*(len(system.atom_types)-2))+'''
group qeq_atoms type 1 2
fix 1 qeq_atoms qeq/reax 1 0.0 10.0 1.0e-6 reax/c''').splitlines()]

	commands = ('''dump 1 all xyz 100 '''+run_name+'''.xyz
thermo 100

compute reax all pair reax/c

fix av all ave/time 1 100 100 c_thermo_pe

thermo_style custom step temp f_av density tpcpu c_reax[1] c_reax[5]
minimize 0.0 1.0e-8 1000 100000
velocity all create 300.0 1337 rot yes dist gaussian
#fix motion all npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0
fix motion all nvt temp 300.0 300.0 100.0
timestep 2.0
run 1000
write_restart '''+run_name+'''.restart''').splitlines()

	for line in commands:
		lmp.command(line)

run('test_cl')
