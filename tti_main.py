from os import listdir,chdir
import subprocess
from GenomeManipulator import GenomeManipulator
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

t = str(sample(ff,1)[0])
sc = 0
cdir = t+'-'+str(sc)+'_control'
check_ok = 'yes'
while (cdir in fff) and sc<11:
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
#	out.write('u 0:1000:end DumpHostTaskGrid\n')
#	out.write('u 0:1000:end DumpParasiteTaskGrid\n')
#	out.write('u 0:1000:end DumpParasiteVirulenceGrid\n')
#	out.write('u 0:1000:end DumpFitnessGrid\n')
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
		pop_n_s = [i for i in listdir('./'+cdir+'/data') if i[-4:]=='spop' if 50000<=int(i[7:-5])<150000]
		h,p = [],[]
		while (h==[] and p==[]):
			y0 = sample(pop_n_s,1)[0]
			out = open('./'+cdir+'/invader_age.txt','w')
			out.write(y0)
			out.close()
			pop_f = [i.split(' ') for i in open('./'+cdir+'/data/'+y0,'r')][25:]
			ccc = set([])
			for i in pop_f:
				if float(i[4])>0:
					ccc|=set(i[17].split(','))
			ccc = list(ccc)
			c_d = dict([i[::-1] for i in enumerate(ccc)])
			occ = [[] for i in ccc]
			for i in pop_f:
				if float(i[4])>0:
					for j in i[17].split(','):
						x = c_d[j]
						occ[x].append([i[2],i[16]])
			occ_p = [i for i in occ if len(i)==2]
			if occ_p!=[]:
				rhp = sample(occ_p,1)[0]
				p = [i for i in rhp if i[0]!='(none)']
				h = [i for i in rhp if i[0]=='(none)']
				p = p[0]
				h = h[0]
		out = open('./'+cdir+'/invader_seq.txt','w')
		out.write(','.join(h)+'\n')
		out.write(','.join(p)+'\n')
		out.close()
		g=GenomeManipulator('./instset-transsmt.cfg')
		p_reg,p_seq = p[0],p[1]
		dna=g.sequence_to_genome(p_seq)
		outp=open(par_dir+'p_invader.org','w')
		for k in dna:
			outp.write(k+'\n')
		outp.close()
		h_reg,h_seq = h[0],h[1]
		dna=g.sequence_to_genome(h_seq)
		outh=open(par_dir+'h_invader.org','w')
		for k in dna:
			outh.write(k+'\n')
		outh.close()
		###invasion p
		out = open('./'+invdir_p+'/events.cfg','w')
		out.write('u begin Inject host-smt.org\n')
		out.write('u 2000 InjectParasite parasite-smt.org ABB 0 400\n')
		for j in sample(range(wsize),int(round(wsize*0.05))):
			out.write('u '+str(y1)+' InjectParasite ./parasites/p_invader.org '+p_reg+' '+str(j)+' '+str(j+1)+'\n')
#		out.write('u '+str(y1)+':1000:end DumpHostTaskGrid\n')
#		out.write('u '+str(y1)+':1000:end DumpParasiteTaskGrid\n')
#		out.write('u '+str(y1)+':1000:end DumpParasiteVirulenceGrid\n')
#		out.write('u '+str(y1)+':1000:end DumpFitnessGrid\n')
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


