import math, copy, sys, random, re, os, cPickle, shutil
sys.path.append("/fs/home/jms875/Library/2.0/tools") 
import g09, files, utils

jobs = []
for i in range(5,10):
	atoms = g09.atoms('PbI2')
	atoms = [atoms[1], atoms[0], atoms[2]]
	for a in atoms:
		a.x += 0.2*(2-random.random())
		a.y += 0.2*(2-random.random())
		a.z += 0.2*(2-random.random())
	jobs.append( g09.job('PbI2_%d' % i, 'HSEH1PBE/LANL2DZ Force SCRF(Solvent=Water)', atoms, queue=None, force=True) )
	shutil.copyfile('gaussian/PbI2.cml', 'gaussian/PbI2_%d.cml' % i)

for j in jobs: j.wait()

