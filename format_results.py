import csv
import os
from numpy import mean,std,array
from math import log
from collections import Counter
from time import sleep
from random import randrange,sample,random
import datetime
global ab_tre

ab_tre = 0

def modification_date(filename):
	t = os.path.getmtime(filename)
	return datetime.datetime.fromtimestamp(t)


def get_ext_rate(pop_0,pop_1):
	pop_h0 = [i for i in pop_0 if i[1].startswith('div')]
	alive_0 = set([i[0] for i in pop_h0 if int(i[4])>ab_tre])
	pop_h1 = [i for i in pop_1 if i[1].startswith('div')]
	alive_1 = [i for i in pop_h1 if int(i[4])>ab_tre]
	d = dict([])
	for i in pop_h1:
		if i[3]!='(none)':
			d[i[0]] = i[3]
		else:
			d[i[0]] = i[0]
	ddd = set([])
	for i in alive_1:
		par = i[3]
		if par!='(none)':
			while d.get(par,'no') not in ['no',par]:
				par = d[par]
				ddd.add(par)
	if len(alive_0)>0:
		ext_rate = len(alive_0-ddd)/float(len(alive_0))
	else:
		ext_rate = 'na'
	return ext_rate


def turnover(com1,com2):
	com1 = set([i[0] for i in com1 if int(i[4])>ab_tre])
	com2 = set([i[0] for i in com2 if int(i[4])>ab_tre])
	a = com1&com2
	b = com1-a
	c = com2-a
	a,b,c = map(float,[len(a),len(b),len(c)])
	return 1-a/(a+b+c)


def turnover_w(pop_0,pop_1):
	c0 = dict([[i[0],int(i[4])] for i in pop_0])
	c1 = dict([[i[0],int(i[4])] for i in pop_1])
	all_spp = sorted(list(set(c0.keys()) | set(c1.keys())))
	p0 = array([c0.get(i,0.0) for i in all_spp])
	p1 = array([c1.get(i,0.0) for i in all_spp])
	den = float(sum((p0)**2)+sum((p1)**2)-sum(p0*p1))
	if den!=0:
		return sum((p0-p1)**2)/den
	else:
		return 'na'


def div_met(x):
	sp = len(x)
	tot = float(sum(x))
	if sp>1:
		shannon = -sum([(i/tot)*log(i/tot) for i in x])
		equ = shannon/float(log(sp))
	else:
		shannon = 0.0
		equ = 1.0
	return (sp,tot,shannon,equ)


def read_dat(dat_f):
	return [i.split() for i in open(dat_f,'r') if i[0] not in ['#','\n']]


def get_pd(pop_):
	pop_h = [i for i in pop_ if i[1].startswith('div')]
	alive = [i for i in pop_h if int(i[4])>ab_tre]
	b_dict = dict([[i[0],int(i[11])] for i in pop_h])
	d = dict([])
	for i in pop_h:
		if i[3]!='(none)':
			d[i[0]] = i[3]
		else:
			d[i[0]] = i[0]
	ddd,p_depth = [],[]
	for i in alive:
		p_depth.append(int(i[13]))
		b0 = int(i[11])
		par = i[3]
		row = []
		if par!='(none)':
			while d.get(par,'no') not in ['no',par]:
				par = d[par]
				b1 = b_dict[par]
				row.append(b0-b1)
		ddd.append(sum(row))
	return (sum(ddd),sum(p_depth))


#####
tot_gen_n = 250000
exp = ['_control','_invasion_p']

for rep in range(10):
	sleep(random()*randrange(300)) ##minimize chances of two processes overwriting the same file
	fff = [i for i in os.listdir('./') if '_control' in i]
	fff_ok = []
	for f in fff:
		f_inv = f.replace('_control','_invasion_p')
		check = os.listdir('./'+f_inv)
		if 'data' in check:
			check_1 = os.listdir('./'+f_inv+'/data/')
			ff = [int(i.split('.')[0].split('-')[1]) for i in check_1 if i[-5:]=='.spop']
			if ff!=[] and max(ff)>(tot_gen_n-1000):
				fff_ok.append(f.split('_')[0])
	for rand_att in range(randrange(1,60)): ##minimize chances of two processes overwriting the same file
		try:
			done = set([i[0] for i in csv.reader(open('done.csv','r'))])
		except:
			done = set([])
		sleep(random())
	to_do = set(fff_ok)-done
	if len(to_do)>0:
		t = sample(to_do,1)[0]
		done_file = open('done.csv','a')
		done_file.write(t+'\n')
		done_file.close()
		out = open('tti_summary_'+t+'.csv','w')
		out_pers =open('tti_pers_'+t+'.csv','w')
		dirs = [t+e for e in exp]
		age = list(open('./'+t+'_control/invader_age.txt','r'))[0].split('-')[1].split('.')[0]
		h_seq,p_seq = list(open('./'+t+'_control/invader_seq.txt','r'))
		h_seq,p_seq = h_seq.split(',')[1].replace('\n',''),p_seq.split(',')[1].replace('\n','')
		inv_time = list(open('./'+t+'_invasion_p/events.cfg','r'))[2].split(' ')[1]
		steps_con = sorted([int(i.split('-')[1].split('.')[0]) for i in os.listdir('./'+t+'_control/data/') if i[-5:]=='.spop'])
		steps_inv = sorted([int(i.split('-')[1].split('.')[0]) for i in os.listdir('./'+t+'_invasion_p/data/') if i[-5:]=='.spop'])
		#make pre-dictionary to track the invader lineage pre-invasion
		pop_0 = [i.split() for i in open('./'+t+'_control/data/detail-'+str(age)+'.spop','r')][25:]
		id_i,ab_i = 'no',0
		for i in pop_0:
			if i[16]==p_seq:
				id_i = i[0]
				ab_i = int(i[4])
				break
		if ab_i == 0:
			ab_i = 1
		d_i_lin = dict()
		for i in pop_0:
			if i[3]!='(none)':
				d_i_lin[i[0]] = i[3]
			else:
				d_i_lin[i[0]] = i[0]
		d_i_lin[id_i]='p_i'
		d_i_lin['p_i']='p_i'
		for i in pop_0:
			if i[0] != id_i and i[3]!='(none)':
				par = i[3]
				while d_i_lin.get(par,'no') not in ['no',par]:
					par = d_i_lin[par]
				d_i_lin[i[0]] = par
		#make pre-dictionary to track the invader lineage post-invasion
		pop_inv = [i.split() for i in open('./'+t+'_invasion_p/data/detail-'+str(inv_time)+'.spop','r')][25:]
		ids = [int(i[0]) for i in pop_inv if int(i[4])>ab_tre]
		d = dict([])
		for i in pop_inv:
			if i[3]!='(none)':
				d[i[0]] = i[3]
			else:
				if int(i[0])>=min(ids):
					if i[2]!='(none)':
						d[i[0]] = 'tti_p'
					else:
						d[i[0]] = 'tti_h'
				else:
					d[i[0]] = i[0]
		for i in pop_inv:
			if d.get(i[3],'no') != 'no':
				par = d[i[3]]
				while d.get(par,'no') not in ['no',par]:
					par = d[par]
				d[i[0]] = par
		for ddd in range(2): #only control + invasion_p
			if ddd == 0:
				steps_ = steps_con[:]
			else:
				steps_ = steps_inv[:]
			#####new code persistency
			pop_inv_pers = [i.split() for i in open(dirs[ddd]+'/data/detail-'+str(inv_time)+'.spop','r')][25:]
			ids = [i[0] for i in pop_inv_pers]
			id_inv = 'no'
			if ddd>0:
				for i in pop_inv_pers:
					if i[16]==p_seq:
						id_inv = i[0] ##identify the id of the invader
						break
			d_pers = dict([]) ##make dictionary of lineages at the time of the invasion
			for i in pop_inv_pers:
				if i[3]!='(none)':
					d_pers[i[0]] = i[3]
				else:
					d_pers[i[0]] = i[0]
			for i in pop_inv_pers:
				if d_pers.get(i[3],'no') != 'no':
					par = d_pers[i[3]]
					while d_pers.get(par,'no') not in ['no',par]:
						par = d_pers[par]
					d_pers[i[0]] = par
			for id_ in ids:
				d_pers[id_] = id_+'_o'
			d_pers[id_inv] = id_inv+'_invader'
			ab_dict = dict([])
			for i in pop_inv_pers:
				ind_n = int(i[4])
				if ind_n>0:
					sp = d_pers[i[0]].split('_')[0]
					ab_dict[sp] = ab_dict.get(sp,0)+ind_n
			pers = []
			for step in steps_[steps_.index(int(inv_time)):steps_.index(int(inv_time))+101]:
				pop_f_p = [i.split() for i in open(dirs[ddd]+'/data/detail-'+str(step)+'.spop','r')][25:]
				for i in pop_f_p:
					if d_pers.get(i[0],'no') == 'no':
						par = d_pers.get(i[3],'no')
						if par != 'no':
							d_pers[i[0]] = par
						else:
							d_pers[i[0]] = i[3]
				for i in pop_f_p:
					if d_pers.get(i[0],'no') != 'no':
						par = d_pers[i[0]]
					else:
						par = d_pers[i[3]]
					while d_pers.get(par,'no') not in ['no',par]:
						par = d_pers[par]
					d_pers[i[0]] = par
				sp_pres = set([])
				for i in pop_f_p:
					if int(i[4])>ab_tre:
						vf = '_vir'
						if i[2]=='(none)':
							vf = '_free'#
						sp = i[0].split('_')[0]
						sp_s = d_pers[i[0]].split('_')[0]
						if id_inv in [sp,sp_s]:
							sp_pres.add(id_inv+'_invader')
						elif sp in ids:
							sp_pres.add(sp+vf)
						elif sp_s in ids:
							sp_pres.add(sp_s+vf)
						else:
							pass
				pers+=list(sp_pres)
			cp = Counter(pers)
			for i in cp.keys():
				out_pers.write(','.join(map(str,[t,exp[ddd][1:]]+i.split('_')+[cp[i],ab_dict[i.split('_')[0]]]))+'\n')
			####end new code persistency
			pop_inv_ok = [i for i in pop_inv if int(i[4])>ab_tre and i[2]=='(none)'][25:]
			for step in steps_:
				pop_f = [i.split() for i in open(dirs[ddd]+'/data/detail-'+str(step)+'.spop','r')][25:]
				phylo_div,phylo_depth = get_pd(pop_f)
				pop_f_h = [i for i in pop_f if int(i[4])>ab_tre and i[2]=='(none)']
				if step>min(steps_):
					pop_pre = [i.split() for i in open(dirs[ddd]+'/data/detail-'+str(steps_[steps_.index(step)-1])+'.spop','r')][25:]
					pop_pre = [i for i in pop_pre if int(i[4])>ab_tre and i[2]=='(none)']
					turn = turnover_w(pop_pre,pop_f_h)
					ext_rate = get_ext_rate(pop_pre,pop_f_h)
				elif ddd==1 and step==min(steps_):
					pop_pre = [i.split() for i in open(dirs[0]+'/data/detail-'+str(steps_con[steps_con.index(step)-1])+'.spop','r')][25:]
					pop_pre = [i for i in pop_pre if int(i[4])>ab_tre and i[2]=='(none)']
					turn = turnover_w(pop_pre,pop_f_h)
					ext_rate = get_ext_rate(pop_pre,pop_f_h)
				else:
					turn = 0
					ext_rate = 0
				if step>=int(inv_time):
					turn0 = turnover_w(pop_inv_ok,pop_f_h)
					if ddd==1:
						pop_control = [i.split() for i in open(dirs[0]+'/data/detail-'+str(step)+'.spop','r')][25:]
						pop_control = [i for i in pop_pre if int(i[4])>ab_tre and i[2]=='(none)']
						turn_comp = turnover_w(pop_control,pop_f_h)
				else:
					turn0 = 0
					turn_comp = 0
				inv_p_sp_ind_n = 0
				div = [int(i[4]) for i in pop_f if int(i[4])>ab_tre]
				div_p = [int(i[4]) for i in pop_f if (int(i[4])>ab_tre and i[2]!='(none)')]
				div_h = [int(i[4]) for i in pop_f if (int(i[4])>ab_tre and i[2]=='(none)')]
				for i in pop_f:
					if step>=int(age):
						if i[0]==id_i:
							inv_p_sp_ind_n = int(i[4])
					if d_i_lin.get(i[0],'no') == 'no':
						par = d_i_lin.get(i[3],'no')
						if par != 'no':
							d_i_lin[i[0]] = par
						else:
							d_i_lin[i[0]] = i[3]
				for i in pop_f:
					if i[0] != id_i and i[3]!='(none)':
						par = i[3]
						while d_i_lin.get(par,'no') not in ['no',par]:
							par = d_i_lin[par]
						d_i_lin[i[0]] = par
				p_i_count = [int(i[4]) for i in pop_f if d_i_lin[i[0]]=='p_i' and int(i[4])>ab_tre]
				p_i_lin_div = len(p_i_count)
				p_i_lin_ab = sum(p_i_count)
				if ddd>0:
					for i in pop_f:
						if d.get(i[0],'no') == 'no':
							par = d.get(i[3],'no')
							if par != 'no':
								d[i[0]] = par
							else:
								d[i[0]] = i[3]
					for i in pop_f:
						if d.get(i[0],'no') != 'no':
							par = d[i[0]]
						else:
							par = d[i[3]]
						while d.get(par,'no') not in ['no',par]:
							par = d[par]
						d[i[0]] = par
					tti_p_ids = [i for i in range(len(pop_f)) if (int(pop_f[i][4])>ab_tre and pop_f[i][2]!='(none)' and d[pop_f[i][0]]=='tti_p')]
					tti_h_ids = [i for i in range(len(pop_f)) if (int(pop_f[i][4])>ab_tre and pop_f[i][2]=='(none)' and d[pop_f[i][0]]=='tti_h')]
					div_tti_p = [int(pop_f[i][4]) for i in tti_p_ids]
					div_tti_h = [int(pop_f[i][4]) for i in tti_h_ids]
					div_tti = div_tti_p+div_tti_h
					div_nat_p = [int(i[4]) for i in pop_f if (int(i[4])>ab_tre and i[2]!='(none)' and d[i[0]] not in ['tti_h','tti_p'])]
					div_nat_h = [int(i[4]) for i in pop_f if (int(i[4])>ab_tre and i[2]=='(none)' and d[i[0]] not in ['tti_h','tti_p'])]
					div_nat = div_nat_p+div_nat_h
					h_pos,p_pos = [],[]
					p_sc,h_sc = 0,0
					for i in range(len(pop_f)):
						if int(pop_f[i][4])>ab_tre and pop_f[i][2]!='(none)':
							if d[pop_f[i][0]] == 'tti_p':
								p_pos.append(['p_i_'+str(p_sc)]+pop_f[i][17].split(','))
							elif d_i_lin[pop_f[i][0]]=='p_i':
								p_pos.append(['p_s_'+str(p_sc)]+pop_f[i][17].split(','))
							else:
								p_pos.append(['p_'+str(p_sc)]+pop_f[i][17].split(','))
							p_sc+=1
						if int(pop_f[i][4])>ab_tre and pop_f[i][2]=='(none)':
							h_pos.append(['h_'+str(h_sc)]+pop_f[i][17].split(','))
							h_sc+=1
				else:
					div_tti = []
					div_tti_p = []
					div_tti_h = []
					div_nat = div
					div_nat_p = div_p
					div_nat_h = div_h
					h_pos,p_pos = [],[]
					p_sc,h_sc = 0,0
					for i in range(len(pop_f)):
						if int(pop_f[i][4])>ab_tre and pop_f[i][2]!='(none)':
							if d_i_lin[pop_f[i][0]]=='p_i':
								p_pos.append(['p_s_'+str(p_sc)]+pop_f[i][17].split(','))
							else:
								p_pos.append(['p_'+str(p_sc)]+pop_f[i][17].split(','))
							p_sc+=1
						if int(pop_f[i][4])>ab_tre and pop_f[i][2]=='(none)':
							h_pos.append(['h_'+str(h_sc)]+pop_f[i][17].split(','))
							h_sc+=1
				net = dict([])
				for i in h_pos:
					for j in p_pos:
						ov = len(set(i[1:])&set(j[1:]))
						if ov>0:
							net[(i[0],j[0])] = ov+net.get((i[0],j[0]),0)
				net_ok = []
				for i in net.keys():
					if net[i]>0:
						net_ok.append(list(i))
				if len(net_ok)>0:
					hhh_net = [i[0] for i in net_ok]
					ppp_net = [i[1] for i in net_ok]
					h_n_net = len(set(hhh_net))
					p_n_net = len(set(ppp_net))
					net_l = len(net_ok)
					net_fill = net_l/float(h_n_net*p_n_net)
					co_h_net = Counter(hhh_net)
					co_p_net = Counter(ppp_net)
					deg_h_net = list(co_h_net.values())
					deg_p_net = list(co_p_net.values())
					p_mean_deg = mean(deg_p_net)
					h_mean_deg = mean(deg_h_net)
					p_min_deg = min(deg_p_net)
					h_min_deg = min(deg_h_net)
					p_max_deg = max(deg_p_net)
					h_max_deg = max(deg_h_net)
					p_sd_deg = std(deg_p_net)
					h_sd_deg = std(deg_h_net)
					i_h_id = set([i for i in hhh_net if i[2]=='i'])
					i_p_id = set([i for i in ppp_net if i[2]=='i'])
					tot_h_inf_by_i = len(set([i[0] for i in net_ok if i[1] in i_p_id]))
					i_h_deg = [co_h_net[i] for i in i_h_id]
					i_p_deg = [co_p_net[i] for i in i_p_id]
					i_h_mean_deg = mean(i_h_deg) if i_h_deg!=[] else 0
					i_p_mean_deg = mean(i_p_deg) if i_p_deg!=[] else 0
					i_h_min_deg = min(i_h_deg) if i_h_deg!=[] else 0
					i_p_min_deg = min(i_p_deg) if i_p_deg!=[] else 0
					i_h_max_deg = max(i_h_deg) if i_h_deg!=[] else 0
					i_p_max_deg = max(i_p_deg) if i_p_deg!=[] else 0
					i_h_std_deg = std(i_h_deg) if i_h_deg!=[] else 0
					i_p_std_deg = std(i_p_deg) if i_p_deg!=[] else 0
					i_h0_id = set([i for i in hhh_net if i[2]=='s'])
					i_p0_id = set([i for i in ppp_net if i[2]=='s'])
					tot_h_inf_by_i0 = len(set([i[0] for i in net_ok if i[1] in i_p0_id]))
					i_h0_deg = [co_h_net[i] for i in i_h0_id]
					i_p0_deg = [co_p_net[i] for i in i_p0_id]
					i_h0_mean_deg = mean(i_h0_deg) if i_h0_deg!=[] else 0
					i_p0_mean_deg = mean(i_p0_deg) if i_p0_deg!=[] else 0
					i_h0_min_deg = min(i_h0_deg) if i_h0_deg!=[] else 0
					i_p0_min_deg = min(i_p0_deg) if i_p0_deg!=[] else 0
					i_h0_max_deg = max(i_h0_deg) if i_h0_deg!=[] else 0
					i_p0_max_deg = max(i_p0_deg) if i_p0_deg!=[] else 0
					i_h0_std_deg = std(i_h0_deg) if i_h0_deg!=[] else 0
					i_p0_std_deg = std(i_p0_deg) if i_p0_deg!=[] else 0
				else:
					h_n_net = 0
					p_n_net = 0
					net_l = 0
					net_fill = 0
					p_mean_deg = 0
					h_mean_deg = 0
					p_min_deg = 0
					h_min_deg = 0
					p_max_deg = 0
					h_max_deg = 0
					p_sd_deg = 0
					h_sd_deg = 0
					tot_h_inf_by_i = 0
					i_h_mean_deg = 0
					i_p_mean_deg = 0
					i_h_min_deg = 0
					i_p_min_deg = 0
					i_h_max_deg = 0
					i_p_max_deg = 0
					i_h_std_deg = 0
					i_p_std_deg = 0
					tot_h_inf_by_i0 = 0
					i_h0_mean_deg = 0
					i_p0_mean_deg = 0
					i_h0_min_deg = 0
					i_p0_min_deg = 0
					i_h0_max_deg = 0
					i_p0_max_deg = 0
					i_h0_std_deg = 0
					i_p0_std_deg = 0
				sp,tot,sh,equ = div_met(div)
				sp_h,tot_h,sh_h,equ_h = div_met(div_h)
				sp_p,tot_p,sh_p,equ_p = div_met(div_p)
				sp_tti,tot_tti,sh_tti,equ_tti = div_met(div_tti)
				sp_tti_p,tot_tti_p,sh_tti_p,equ_tti_p = div_met(div_tti_p)
				sp_tti_h,tot_tti_h,sh_tti_h,equ_tti_h = div_met(div_tti_h)
				sp_nat,tot_nat,sh_nat,equ_nat = div_met(div_nat)
				sp_nat_p,tot_nat_p,sh_nat_p,equ_nat_p = div_met(div_nat_p)
				sp_nat_h,tot_nat_h,sh_nat_h,equ_nat_h = div_met(div_nat_h)
				res = [step,sp,tot,sh,equ,sp_h,tot_h,sh_h,equ_h,
						sp_p,tot_p,sh_p,equ_p,sp_tti,tot_tti,sh_tti,equ_tti,sp_tti_p,
						tot_tti_p,sh_tti_p,equ_tti_p,sp_tti_h,tot_tti_h,sh_tti_h,equ_tti_h,
						sp_nat,tot_nat,sh_nat,equ_nat,sp_nat_p,tot_nat_p,sh_nat_p,equ_nat_p,
						sp_nat_h,tot_nat_h,sh_nat_h,equ_nat_h,
						h_n_net,p_n_net,net_l,net_fill,p_mean_deg,h_mean_deg,
						p_min_deg,h_min_deg,p_max_deg,h_max_deg,p_sd_deg,h_sd_deg,
						tot_h_inf_by_i,i_h_mean_deg,i_p_mean_deg,i_h_min_deg,
						i_p_min_deg,i_h_max_deg,i_p_max_deg,i_h_std_deg,i_p_std_deg,
						tot_h_inf_by_i0,i_h0_mean_deg,i_p0_mean_deg,i_h0_min_deg,i_p0_min_deg,
						i_h0_max_deg,i_p0_max_deg,i_h0_std_deg,i_p0_std_deg,
						p_i_lin_div,p_i_lin_ab,inv_p_sp_ind_n,phylo_div,
						phylo_depth,turn,ext_rate,turn0,turn_comp]
				wl = out.write(','.join(map(str,[t,age,inv_time,exp[ddd][1:]]+res))+'\n')
		out.close()
		out_pers.close()

