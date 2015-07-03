import math, copy, sys, random, re, os, cPickle, shutil
sys.path.append("/fs/home/jms875/Library/2.0/tools") 
import g09, files, utils

def new_jobs():
	for i in range(0,20):
		atoms = [utils.Atom('Pb',0.,0.,0.), utils.Atom('I',0.,2.7,0.), utils.Atom('I',2.7,0.,0.), utils.Atom('I',-2.7,0.,0.)]
		for a in atoms:
			a.x += 0.5*(2-random.random())
			a.y += 0.5*(2-random.random())
			a.z += 0.5*(2-random.random())
		#atoms[-1].x += i*0.2
		g09.job('PbI3_%d' % i, 'HSEH1PBE/LANL2DZ Force SCRF(Solvent=Water)', atoms, queue=None,charge_and_multiplicity='-1,1').wait()
		shutil.copyfile('gaussian/PbI3.cml', 'gaussian/PbI3_%d.cml' % i)

	#for j in jobs: j.wait()

def restart_jobs():
	for root, dirs, file_list in os.walk("gaussian"):
		for ff in file_list:
			if ff.endswith('.log'):
				name = ff[:-4]
				if not name.startswith('PbI'): continue
				if name.endswith('_def2SVP'): continue
				if os.path.exists('gaussian/%s_def2SVP.log' % name): continue
				print name, '...',
				g09.job('%s_def2SVP' % name, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water) Geom=AllCheck', queue=None, force=True, previous=name).wait()
				print 'Done'
				shutil.copyfile('gaussian/%s.cml' % name, 'gaussian/%s_def2SVP.cml' % name)

def new_jobs_2():
	for i in range(20,40):
		atoms = [utils.Atom('Pb',0.,0.,0.), utils.Atom('I',0.,2.7,0.)]
		atoms[1].x += 0.3*(2-random.random())
		g09.job('PbI+_%d_def2SVP' % i, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue=None, charge_and_multiplicity='1,1').wait()
		shutil.copyfile('gaussian/PbI+.cml', 'gaussian/PbI+_%d_def2SVP.cml' % i)

def from_md():
	frames = files.read_xyz('out.xyz')
	#g09.job('PbI2_PbI2_def2SVP', 'HSEH1PBE/Def2SVP Opt SCRF(Solvent=Water)', frames[1], queue=None, force=True).wait()
	shutil.copyfile('gaussian/PbI2_PbI2.cml', 'gaussian/PbI2_PbI2_def2SVP.cml')
	for i,atoms in enumerate(frames[2:]):
		#print 'Running %d' % i
		#g09.job('PbI2_PbI2_%d_def2SVP' % i, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue=None).wait()
		shutil.copyfile('gaussian/PbI2_PbI2.cml', 'gaussian/PbI2_PbI2_%d_def2SVP.cml' % i)

from_md()

