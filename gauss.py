import math, copy, sys, random, re, os, cPickle, shutil
sys.path.append("/fs/home/jms875/Library/2.0/tools") 
import g09, files, utils

for i in range(0,20):
	atoms = g09.atoms('PbI+')
	for a in atoms:
		a.x += 0.5*(2-random.random())
		a.y += 0.5*(2-random.random())
		a.z += 0.5*(2-random.random())
	g09.job('PbI+_%d' % i, 'HSEH1PBE/LANL2DZ Force SCRF(Solvent=Water)', atoms, queue=None, force=True, charge_and_multiplicity='1,1').wait()
	shutil.copyfile('gaussian/PbI+.cml', 'gaussian/PbI+_%d.cml' % i)

#for j in jobs: j.wait()

