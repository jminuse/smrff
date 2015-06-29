import math, copy, sys, random, re, os, cPickle, shutil
sys.path.append("/fs/home/jms875/Library/2.0/tools") 
import g09, files, utils

for i in range(0,20):
	atoms = [utils.Atom('Pb',0.,0.,0.), utils.Atom('I',2.,0.,0.)]
	#for a in atoms:
	#	a.x += 0.5*(2-random.random())
	#	a.y += 0.5*(2-random.random())
	#	a.z += 0.5*(2-random.random())
	atoms[1].x += i*0.25
	g09.job('PbI+_r%d' % i, 'HSEH1PBE/LANL2DZ Force SCRF(Solvent=Water)', atoms, queue=None, force=True, charge_and_multiplicity='1,1').wait()
	shutil.copyfile('gaussian/PbI+.cml', 'gaussian/PbI+_r%d.cml' % i)

#for j in jobs: j.wait()

