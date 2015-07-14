import os, sys, random, math, re
from lammps import lammps
sys.path.append("/fs/home/jms875/Library/2.0/tools")
import utils, files, g09

os.chdir('lammps')

commands = ('''units		real

atom_style	charge
read_data	data.tatb

pair_style      reax/c NULL
pair_coeff      * * tatb.reax C H O N

compute reax all pair reax/c

variable eb  	 equal c_reax[1]
variable ea  	 equal c_reax[2]
variable elp 	 equal c_reax[3]
variable emol 	 equal c_reax[4]
variable ev 	 equal c_reax[5]
variable epen 	 equal c_reax[6]
variable ecoa 	 equal c_reax[7]
variable ehb 	 equal c_reax[8]
variable et 	 equal c_reax[9]
variable eco 	 equal c_reax[10]
variable ew 	 equal c_reax[11]
variable ep 	 equal c_reax[12]
variable efi 	 equal c_reax[13]
variable eqeq 	 equal c_reax[14]

neighbor	2.5 bin
neigh_modify	delay 0 every 5 check no

fix		1 all nve
fix     2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

thermo		5
thermo_style 	custom step temp epair etotal press v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq

timestep	0.0625

run		25
''').splitlines()
lmp = lammps('',['-log', 'tatb.log'])
for line in commands:
	lmp.command(line)


