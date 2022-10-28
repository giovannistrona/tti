import os
fff = [i for i in os.listdir('./') if 'tti_summary_' in i]
head = 'rep,i_t,scenario,step,sp,tot,sh,equ,sp_h,tot_h,sh_h,equ_h,sp_p,tot_p,sh_p,equ_p\n'
out = open('summary_simple.csv','w')
out.write(head)
sc = 0
for ff in fff:
	inp = open(ff,'r')
	for row in inp:
		out.write(row)
	sc+=1
	print (sc)


out.close()