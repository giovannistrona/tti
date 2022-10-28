import os,subprocess
import datetime
import shutil

def modification_date(filename):
	t = os.path.getmtime(filename)
	return datetime.datetime.fromtimestamp(t)


#status = [i.split() for i in str(subprocess.check_output(['squeue', '-u','stronagi'])).split('\\n')]
#for i in status:
#	if ('tti_format' in i and 'RUNNING' not in i) or i[-2:]==['launch','fai']:
#		proc = subprocess.Popen(['scancel '+i[0]],shell=True)
#



#####get environmental variables
fff = os.listdir('./')
ff = []
for i in fff:
	try:
		ff.append(int(i))
	except:
		pass



env_dict = dict([])
for i in ff:
	env_file = open('./'+str(i)+'/environment.cfg','r')
	res = []
	for line in env_file:
		row = line.split()
		if len(row)>1 and row[0]=='RESOURCE':
			inflow = float(row[1].split(':')[1].split('=')[1])
			outflow = float(row[1].split(':')[2].split('=')[1])
			res.append([inflow,outflow])
	env_file.close()
	cfg_file = open('./'+str(i)+'/avida.cfg','r')
	x,y = 'na','na'
	for line in cfg_file:
		row = line.split()
		if len(row)>1 and row[0]=='WORLD_X':
			x = int(row[1])
		if len(row)>1 and row[0]=='WORLD_Y':
			y = int(row[1])
	env_dict[str(i)] = ','.join(map(str,[x*y,len(res),sum([j[0] for j in res]),
							 sum([j[1] for j in res]),
							 sum([j[0]*j[1] for j in res]),
							 sum([j[0]*j[1] for j in res])/float(len(res))
							 ]))





fff = [i for i in os.listdir('./') if 'tti_summary_' in i and 'complete' not in i]
col_names = 'size,res_n,inflow,outflow,eff_res,eff_res_mean,rep,i_age,i_t,scenario,step,sp,tot,sh,equ,sp_h,tot_h,sh_h,equ_h,sp_p,tot_p,sh_p,equ_p,sp_tti,tot_tti,sh_tti,equ_tti,sp_tti_p,tot_tti_p,sh_tti_p,equ_tti_p,sp_tti_h,tot_tti_h,sh_tti_h,equ_tti_h,sp_nat,tot_nat,sh_nat,equ_nat,sp_nat_p,tot_nat_p,sh_nat_p,equ_nat_p,sp_nat_h,tot_nat_h,sh_nat_h,equ_nat_h,h_n_net,p_n_net,net_l,net_fill,p_mean_deg,h_mean_deg,p_min_deg,h_min_deg,p_max_deg,h_max_deg,p_sd_deg,h_sd_deg,tot_h_inf_by_i,i_h_mean_deg,i_p_mean_deg,i_h_min_deg,i_p_min_deg,i_h_max_deg,i_p_max_deg,i_h_std_deg,i_p_std_deg,tot_h_inf_by_i0,i_h0_mean_deg,i_p0_mean_deg,i_h0_min_deg,i_p0_min_deg,i_h0_max_deg,i_p0_max_deg,i_h0_std_deg,i_p0_std_deg,p_i_lin_div,p_i_lin_ab,inv_p_sp_ind_n,phylo_div,phylo_depth,turn,ext_rate,turn0,turn_comp\n'
out = open('tti_summary_complete.csv','w')
out.write(col_names)
to_del,done_sum = set([]),set([])
for t in fff:
	inp = open(t,'r')
	row = 'null'
	for i in inp:
		row = i.split(',')
	inp.close()
	if row!='null' and len(row)==80 and row[3] == 'invasion_p' and row[4] == '250000': ##if you add fields to the format file, modify the 80
		env = env_dict[t.split('_')[2].split('-')[0]]
		inp = open(t,'r')
		for row in inp:
			out.write(env+','+row)
		inp.close()
		done_sum.add(t.split('_')[-1][:-4])
	else:
		to_del.add(t.split('_')[-1])


out.close()





fff = [i for i in os.listdir('./') if 'tti_pers_' in i and 'complete' not in i]
col_names = 'rep,scenario,sp_id,vir_free,pers,abund\n'
out = open('tti_pers_complete.csv','w')
out.write(col_names)
done_pers = set([])
for t in fff:
	inp = open(t,'r')
	row = 'null'
	for i in inp:
		row = i.split(',')
	inp.close()
	if row!='null' and len(row)==6:# and row[1] == 'invasion_p':
		inp = open(t,'r')
		for row in inp:
			out.write(row)
		inp.close()
		done_pers.add(t.split('_')[-1][:-4])
	else:
		to_del.add(t.split('_')[-1])


out.close()

#for t in to_del:
#	os.remove('tti_summary_'+t)
#	os.remove('tti_pers_'+t)



#fff = [i for i in os.listdir('./') if 'core' in i]
#for t in fff:
#	os.remove(t)
#
#
#out = open('done.csv','w')
#for f in done_sum & done_pers:
#	out.write(f+'\n')
#
#
#out.close()
#
#
##for i in range(100):
##		proc = subprocess.Popen(['sbatch format_results.job'],shell=True)
#
#
#
#fff = [i for i in os.listdir('./') if i[-8:] == '_control']
#to_del = []
#for ff in fff:
#	ff_inv = ff.replace('control','invasion_p')
#	now = datetime.datetime.now()
#	d0 = os.listdir('./'+ff)
#	if 'core' in [i[:4] for i in d0] or ('data' not in d0 and (datetime.datetime.now()-modification_date(ff)).total_seconds()>3600):
#		to_del.append(ff.split('_')[0])
#	elif 'data' in d0:
#		d1 = [int(i.split('-')[1].split('.')[0]) for i in os.listdir('./'+ff+'/data/') if i[-5:]=='.spop']
#		if 250000 not in d1:
#			mod = modification_date('./'+ff+'/data')
#			if (datetime.datetime.now()-mod).total_seconds()>3600:
#				to_del.append(ff.split('_')[0])
#		else:
#			d0_inv = os.listdir('./'+ff_inv)
#			if 'core' in [i[:4] for i in d0_inv] or ('data' not in d0_inv and (datetime.datetime.now()-modification_date(ff_inv)).total_seconds()>3600):
#				to_del.append(ff.split('_')[0])
#			elif 'data' in d0_inv:
#				d2 = [int(i.split('-')[1].split('.')[0]) for i in os.listdir('./'+ff_inv+'/data/') if i[-5:]=='.spop']
#				if 250000 not in d2:
#					mod = modification_date('./'+ff_inv+'/data')
#					if (datetime.datetime.now()-mod).total_seconds()>(3600*48):
#						to_del.append(ff.split('_')[0])
#
#
#
#for f in to_del:
#	try:
#		#proc = subprocess.Popen(['sbatch tti_3.job'],shell=True)
#		shutil.rmtree(f+'_control')
#		shutil.rmtree(f+'_invasion_p')
#	except:
#		pass
