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

Cl_ = 66
H_ = 54
N_ = 53
Pb_ = 111
Pb = 907
Cl = 838
extra = {
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.0, vdw_e=0.0, vdw_r=3.0),
	Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.0, vdw_r=3.0),
}
def PbCl2_acetone():
	#g09.job('PbCl2_4acetone', 'HSEH1PBE/Def2SVP Opt', atoms=files.read_xyz('out.xyz'), queue=None, force=True)
	acet = utils.Molecule('cml/acetone')
	pbcl2 = utils.Molecule('cml/PbCl2', extra_parameters=extra)
	system = utils.System()
	system.add(pbcl2)
	for i in range(4): system.add(acet)
	files.write_cml(system.atoms, system.bonds)


def wiggle_PbCl2_acetone():
	g09.job('PbCl2_4acetone_tzvp', 'HSEH1PBE/Def2TZVP Opt Guess=Read Geom=(Check,NewDefinition)', previous='PbCl2_4acetone', queue=None)
	shutil.copyfile('gaussian/PbCl2_4acetone.cml', 'gaussian/PbCl2_4acetone_tzvp.cml')
	for i in range(20):
		atoms = g09.atoms('PbCl2_4acetone')
		for a in atoms:
			if a.element != 'Pb':
				a.x += 0.1*(2-random.random())
				a.y += 0.1*(2-random.random())
				a.z += 0.1*(2-random.random())
		name = 'PbCl2_4acetone_w%d' % i
		g09.job(name, 'HSEH1PBE/Def2TZVP Force Guess=Read', atoms, previous='PbCl2_4acetone', queue=None).wait()
		shutil.copyfile('gaussian/PbCl2_4acetone.cml', 'gaussian/'+name+'.cml')

def tzvp_blaire_vac():
	for name in ['5_acetone_5', '5_acetone_5_PbCl2', 'pbi2_6methc_9', 'pbcl2_6methc_9', '4ACN_s_4', '4ACN_s_4_pbcl2', '5buty22111.2_2', '5buty22111.2_2_PbCl2']:
		shutil.copyfile('/fs/home/bas348/perovskites/gaussian/'+name+'.chk', 'gaussian/'+name+'.chk')
		g09.job(name+'_svp', 'HSEH1PBE/Def2SVP Opt Guess=Read Geom=(Check,NewDefinition)', previous=name, queue='batch')

'''
Running:
['5_acetone_5', '5_acetone_5_PbCl2', 'pbi2_6methc_9', 'pbcl2_6methc_9', '4ACN_s_4', '4ACN_s_4_pbcl2', '5buty22111.2_2', '5buty22111.2_2_PbCl2']+'svp'
PbCl2_4acetone_tzvp
Done:
PbCl2_4acetone_w%d' % range(10)
'''

from multiprocessing import Process

def run_cl2_pair(i):
	frame = files.read_xyz('pbcl2.xyz')
	frame[0].y -= i*0.4
	name = 'PbCl2_c_%d_vac' % i
	g09.job(name, 'HSEH1PBE/Def2TZVP Force SCF=(IntRep,QC)', frame, queue=None, force=True)

def cl2_pair():
	for i in range(5,10):
		Process(target=run_cl2_pair, args=(i,)).start()

def PbCl2_2():
	atoms=files.read_xyz('gaussian/PbCl2_2.xyz')
	g09.job('PbCl2t_2_def2SVP', 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue=None,force=True).wait()
	shutil.copyfile('gaussian/PbCl2_2.cml', 'gaussian/PbCl2t_2_def2SVP.cml')




def run_PbCl2_reax_min(i,atoms):
	g09.job('PbCl2_reax_min_vac_%d' % i, 'HSEH1PBE/Def2TZVP Force Guess=Read', atoms, previous='PbCl2_opt_vac', queue='batch',force=True).wait()
	shutil.copyfile('gaussian/PbCl2_reax_min.cml', 'gaussian/PbCl2_reax_min_vac_%d.cml' % i)

def PbCl2_reax_min():
	frames = files.read_xyz('lammps/testing.xyz')
	#g09.job('PbI2_PbI2_def2SVP', 'HSEH1PBE/Def2SVP Opt SCRF(Solvent=Water)', frames[1], queue=None, force=True).wait()
	shutil.copyfile('gaussian/PbCl2_rot80_vac.cml', 'gaussian/PbCl2_reax_min.cml')
	for i,atoms in enumerate(frames[2:]):
		procs = Process(target=run_PbCl2_reax_min, args=(i,atoms))
		procs.start()
		print 'Running %d' % i
		# g09.job('PbCl2_reax_min_%d' % i, 'HSEH1PBE/Def2SVP Force SCRF(Solvent=Water)', atoms, queue='batch',force=True).wait()
		# shutil.copyfile('gaussian/PbCl2_reax_min.cml', 'gaussian/PbCl2_reax_min_%d.cml' % i)


PbCl2_reax_min()
