import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

os.chdir('lammps')

commands = ('''units		real

atom_style	full
read_data	tatb.data

pair_style      reax/c NULL
pair_coeff      * * ../input.reax Pb I

fix     1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

run		1
''').splitlines()
lmp = lammps('',['-log', 'tatb.log'])
for line in commands:
	lmp.command(line)


