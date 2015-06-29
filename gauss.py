import math, copy, sys, random, re, os, cPickle, shutil
sys.path.append("/fs/home/jms875/Library/2.0/tools") 
import g09, files, utils

for i in range(0,2):
	atoms = [utils.Atom('Pb',0.,0.,0.), utils.Atom('I',0.,2.7,0.), utils.Atom('I',2.2,0.,0.)]
	#for a in atoms:
	#	a.x += 0.5*(2-random.random())
	#	a.y += 0.5*(2-random.random())
	#	a.z += 0.5*(2-random.random())
	atoms[-1].x += i*0.2
	g09.job('PbI2_r%d' % i, 'HSEH1PBE/LANL2DZ Force SCRF(Solvent=Water)', atoms, queue=None, force=True).wait()
	shutil.copyfile('gaussian/PbI2.cml', 'gaussian/PbI2_r%d.cml' % i)

#for j in jobs: j.wait()

