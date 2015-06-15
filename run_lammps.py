import os, copy, sys, random, cPickle, math
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
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

system = utils.System(box_size=(30., 30., 30.), name='test_tersoff')

pbi2 = utils.Molecule('cml/PbI2', extra_parameters=extra)
pbi2.bonds = []
pbi2.angles = []

system.add(pbi2, 0, 2)

os.chdir('lammps')
files.write_lammps_data(system)

commands = ('''units real
atom_style full
pair_style tersoff

read_data	'''+system.name+'''.data

pair_coeff * * tersoff SiO.tersoff O Si

dump 1 all xyz 1000 '''+system.name+'''.xyz

min_style fire
minimize 0.0 1.0e-8 1000 100000
''').splitlines()
lmp = lammps()
for line in commands:
	lmp.command(line)

xyz = files.read_xyz(system.name+'.xyz')[-1] # take last frame from LAMMPS job
for i,a in enumerate(system.atoms):
	a.x, a.y, a.z = xyz[i].x, xyz[i].y, xyz[i].z

