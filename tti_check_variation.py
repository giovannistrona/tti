from os import listdir,chdir
import subprocess
from random import randrange,sample
from shutil import copy,copytree,rmtree
from os import mkdir
from time import sleep

sleep(randrange(300))

fff = listdir('./')
ff = []
for i in fff:
	try:
		ff.append(int(i))
	except:
		pass




check_ok = 'yes'
sc = 0
t = str(sample(ff,1)[0])
cdir = t+'-'+str(sc)+'_control'
if cdir not in fff:
	check_ok = 'no'




while (cdir in fff) and sc<10:
	for f in sample(ff,len(ff)):
		t = str(f)
		cdir = t+'-'+str(sc)+'_control'
		if cdir not in fff:
			check_ok = 'no'
			break
	sc+=1


while check_ok == 'no':
	par_dir = './'+cdir+'/parasites/'
	invdir_p = cdir.replace('_control','_invasion_p')
	mkdir(cdir)
	mkdir(invdir_p)
	mkdir(par_dir)
	y1 = randrange(150,200)*1000 #used later, random time for tt invasion
	r_seed = str(randrange(100000000))
	for cfd in [cdir,invdir_p]:
		copy('./'+t+'/environment.cfg','./'+cfd+'/environment.cfg')
		cfg_i = open('./'+t+'/avida.cfg','r')
		cfg_o = open('./'+cfd+'/avida.cfg','w')
		for i in cfg_i:
			if i[:11]=='RANDOM_SEED':
				cfg_o.write('RANDOM_SEED '+r_seed+'\n')
			else:
				cfg_o.write(i)
		for cf in ['avida','instset-transsmt.cfg','host-smt.org','parasite-smt.org']:
			copy(cf,'./'+cfd+'/'+cf)
	cfg_i.close()
	cfg_o.close()
	#WSIZE
	cfg_file = [i.split(' ') for i in open('./'+cdir+'/avida.cfg','r')]
	x = [int(i[1]) for i in cfg_file if len(i)>0 and i[0]=='WORLD_X'][0]
	y = [int(i[1]) for i in cfg_file if len(i)>0 and i[0]=='WORLD_Y'][0]
	wsize = x*y
	###GENERATE EVENTS FILES
	#reference scenario
	out = open('./'+cdir+'/events.cfg','w')
	out.write('u begin Inject host-smt.org\n')
	out.write('u 2000 InjectParasite parasite-smt.org ABB 0 400\n')
	out.write('u 0:1000:end PrintHostPhenotypeData\n')
	out.write('u 0:1000:end PrintHostTasksData\n')
	out.write('u 0:1000:end PrintParasitePhenotypeData\n')
	out.write('u 0:1000:end PrintParasiteTasksData\n')
	out.write('u 0:1000:end SavePopulation\n')
	out.write('u '+str(y1)+':100:'+str(y1+10000)+' SavePopulation\n')
	out.write('u 250000 Exit') ###modified to 500k
	out.close()
	###RUN CONTROL
	chdir('./'+cdir)
	subprocess.call(['./avida'])
	chdir('..')
	#####
	check_ok = 'yes'
	try:
		check_control = [i.split() for i in open('./'+cdir+'/data/detail-250000.spop','r')][25:]
	except:
		check_ok = 'no'
		rmtree('./'+cdir+'/data')
	if check_ok == 'yes':
		###invasion p
		out = open('./'+invdir_p+'/events.cfg','w')
		out.write('u begin Inject host-smt.org\n')
		out.write('u 2000 InjectParasite parasite-smt.org ABB 0 400\n')
		out.write('u '+str(y1)+' InjectRandom 1\n')
		out.write('u '+str(y1)+':1000:end PrintHostPhenotypeData\n')
		out.write('u '+str(y1)+':1000:end PrintHostTasksData\n')
		out.write('u '+str(y1)+':1000:end PrintParasitePhenotypeData\n')
		out.write('u '+str(y1)+':1000:end PrintParasiteTasksData\n')
		out.write('u '+str(y1)+':100:'+str(y1+10000)+' SavePopulation\n')
		out.write('u '+str(y1)+':1000:end SavePopulation\n')
		out.write('u 250000 Exit')
		out.close()
		####
		for cfd in [invdir_p]:
			copytree(par_dir,'./'+cfd+'/parasites')
		for cfd in [invdir_p]:
			chdir('./'+cfd)
			subprocess.call(['./avida'])
			chdir('..')
		try:
			check_control = [i.split() for i in open('./'+invdir_p+'/data/detail-250000.spop','r')][25:]
		except:
			check_ok = 'no'
			rmtree('./'+invdir_p+'/data')


