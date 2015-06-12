import os, copy, sys, random, cPickle, math
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0")
import utils, files, g09

system = utils.System(box_size=(30., 30., 30.), name='test_propane')

propane = utils.Molecule('gaussian/propane')

system.add(propane, 0, 2)
system.add(propane, 0, -2)

os.chdir('lammps')
files.write_lammps_data(system)

commands = ('''units real
atom_style full
pair_style lj/cut/coul/cut 10.0 100.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

read_data	'''+system.name+'''.data

dump 1 all xyz 1000 '''+system.name+'''.xyz

min_style fire
minimize 0.0 1.0e-8 1000 100000

#fix motion all nve
#fix solvent all langevin 200.0 100.0 100.0 1337
#timestep 2.0
#run 10000
''').splitlines()
lmp = lammps()
for line in commands:
	lmp.command(line)

xyz = files.read_xyz(system.name+'.xyz')[-1] # take last frame from LAMMPS job
for i,a in enumerate(system.atoms):
	a.x, a.y, a.z = xyz[i].x, xyz[i].y, xyz[i].z

os.chdir('..')

files.write_cml(system.atoms, system.bonds, 'gaussian/'+system.name)
g09.job(system.atoms, 'M062X/cc-pVDZ', 'batch', system.name, 'SP Force')


