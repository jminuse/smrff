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

def make_cl():
	for root, dirs, file_list in os.walk("gaussian"):
		for ff in file_list:
			if ff.endswith('.log'):
				name = ff[:-4]
				if not name.startswith('PbI2'): continue
				if not name.endswith('_def2SVP'): continue
				new_name = name.replace('I', 'Cl')
				atoms = g09.atoms(name)
				for a in atoms:
					if a.element=='I':
						a.element = 'Cl'
				if not os.path.exists('gaussian/'+new_name+'.log'):
					g09.job(new_name, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue=None).wait()
				#shutil.copyfile('gaussian/'+name+'.cml', 'gaussian/'+new_name+'.cml')
				print new_name, 'done'
	
def change_cml():
	for root, dirs, file_list in os.walk("gaussian"):
		for ff in file_list:
			if ff.endswith('.cml') and 'Cl2' in ff:
				try:
					new_contents = open('gaussian/'+ff.replace('Cl', 'I')).read().replace('I', 'Cl')
					open('gaussian/'+ff, 'w').write(new_contents)
				except:
					pass
def new_cl_jobs():
	for i in range(1,20):
		atoms = g09.atoms('PbCl2_opt_def2SVP')
		atoms[2].x += 5.0*random.random()
		atoms[2].y += 1.0*(2-random.random())
		atoms[2].z += 1.0*(2-random.random())
		name = 'PbCl2_p%d_def2SVP' % i
		g09.job(name, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue=None).wait()
		shutil.copyfile('gaussian/PbCl2_opt_def2SVP.cml', 'gaussian/'+name+'.cml')

def cl_no_solv():
	for root, dirs, file_list in os.walk("gaussian"):
		for ff in file_list:
			if ff.endswith('.cml') and 'Cl2' in ff and 'PbCl2_PbCl2' not in ff:
				old_job = ff[:-4]
				name = old_job.replace('def2SVP', 'vac')
				if os.path.exists('gaussian/%s.log' % name):
					print name
					continue
				print '%20s' % old_job, name
				g09.job(name, 'HSEH1PBE/Def2TZVP Force Geom=Check Guess=Read', previous=old_job, queue=None).wait()
				shutil.copyfile('gaussian/PbCl2_opt_def2SVP.cml', 'gaussian/'+name+'.cml')

def closer_cl_jobs():
	for i in range(10):
		atoms = g09.atoms('PbCl2_opt_vac')
		for a in [atoms[0], atoms[2]]:
			a.x += 0.5*(2-random.random())
			a.y += 0.5*(2-random.random())
			a.z += 0.5*(2-random.random())
		name = 'PbCl2_w%d_vac' % i
		g09.job(name, 'HSEH1PBE/Def2TZVP Force Guess=Read', atoms, previous='PbCl2_opt_vac', queue=None).wait()
		shutil.copyfile('gaussian/PbCl2_opt_def2SVP.cml', 'gaussian/'+name+'.cml')
		#g09.job('PbCl2_opt_vac', 'HSEH1PBE/Def2TZVP Opt Guess=Read Geom=Check', previous='PbCl2_opt_def2SVP', queue=None).wait()

def rotate_cl_jobs():
	#frames = []
	for degrees in range(80,190,10):
		atoms = [utils.Atom('Pb',0.,0.,0.), utils.Atom('Cl',2.459898,0.,0.), utils.Atom('Cl',2.459898,0.,0.)  ]
		t = degrees*math.pi/180
		R = [[math.cos(t), -math.sin(t), 0.],
			 [math.sin(t),  math.cos(t), 0.],
			 [0., 0., 1.]]
		atoms[2].x, atoms[2].y, atoms[2].z = utils.matvec(R, (atoms[2].x, atoms[2].y, atoms[2].z) )
		center = copy.deepcopy(atoms[1])
		for a in atoms:
			a.x -= center.x
			a.y -= center.y
			a.z -= center.z
		
		name = 'PbCl2_rot%d_vac' % degrees
		print name
		#frames.append(atoms)
		g09.job(name, 'HSEH1PBE/Def2TZVP Force Guess=Read', atoms, previous='PbCl2_opt_vac', queue=None).wait()
		shutil.copyfile('gaussian/PbCl2_opt_def2SVP.cml', 'gaussian/'+name+'.cml')
	#files.write_xyz(frames)

def PbCl2_2_vac():
	#g09.job('PbCl2_2_opt_vac', 'HSEH1PBE/Def2TZVP Opt Guess=Read Geom=Check', previous='PbCl2_PbCl2_def2SVP', queue=None).wait()
	for i in range(10,18):
		atoms = g09.atoms('PbCl2_2_opt_vac')
		for a in atoms:
			if a.element != 'Pb':
				a.x += 0.25*(2-random.random())
				a.y += 0.25*(2-random.random())
				a.z += 0.25*(2-random.random())
		name = 'PbCl2_2_w%d_vac' % i
		g09.job(name, 'HSEH1PBE/Def2TZVP Force Guess=Read', atoms, previous='PbCl2_2_opt_vac', queue=None)
		shutil.copyfile('gaussian/PbCl2_2_opt_vac.cml', 'gaussian/'+name+'.cml')
	
PbCl2_2_vac()

