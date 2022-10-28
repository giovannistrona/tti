###Analyses in the same order as presented in the draft
library(vioplot)
library(viridis)
library(pracma)
library(corrplot)
library(randomForest)
library(Rmpfr)

###Function to plot two series in one plot with confidence intervals

plot_ci_line<-function(x_,y_,ylab,color_n=1,
                       color_id=1,
                       xlab='time',add='no',
                       x0='na',x1='na',y0='na',y1='na'){
  xy<-aggregate(y_~x_,FUN='mean',na.rm=TRUE, na.action=NULL)
  xy_max<-aggregate(y_~x_,FUN='max',na.rm=TRUE, na.action=NULL)[,2]
  xy_min<-aggregate(y_~x_,FUN='min',na.rm=TRUE, na.action=NULL)[,2]
  x<-xy[,1]
  y<-xy[,2]
  xy_sd<-aggregate(y_~x_,FUN='sd',na.rm=TRUE, na.action=NULL)[,2]
  xy_n<-aggregate(y_~x_,FUN='length')[,2]
  ci<-qnorm(0.99)*xy_sd/sqrt(xy_n)
  y_min<-min(y-ci)
  y_max<-max(y+ci)
  x_min<-min(x)
  x_max<-max(x)
  if (y0!='na'){y_min<-y0}
  if (x0!='na'){x_min<-x0}
  if (y1!='na'){y_max<-y1}
  if (x1!='na'){x_max<-x1}
  if (add=='no'){
    plot(x,y,type='n',las=1,xlab=xlab,ylab=ylab,xlim=c(x_min,x_max),ylim=c(y_min,y_max))}
  polygon(c(rev(x), x), c(rev(y-ci),(y+ci)),
          col = plasma(color_n,alpha=0.5)[color_id], border = NA)
  lines(x,y,lwd=1.5,col=plasma(color_n)[color_id])
  print (paste('mean =',round(y[length(y)],2),'sd =',round(xy_sd[length(y)],2),'min = ',round(xy_min[length(y)],2),'max = ',round(xy_max[length(y)],2)))
}


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  #dev.new(width=1.75, height=5)
  plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,xlim=c(0,5))
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,1,y+1/scale, col=lut[i], border=NA)
  }
}

#######################################################################################
##############1 - Compare invader vs native species persistence
#######################################################################################

sum_data<-read.csv('./results/tti_summary_complete.csv',header=T)
reps_id<-unique(sum_data$rep)
####identify simulations where parasites went extinct; those will not be considered
rep_p_ext<-c()
sc<-0
for (rep_id in reps_id){
  sub_a<-sum_data[sum_data$rep==rep_id,]
  inv_time<-sub_a$i_t[1]
  sub_a<-sub_a[sub_a$step>2000,]
  sub_a<-sub_a[sub_a$step<inv_time,]
  if (min(sub_a$tot_p)==0){rep_p_ext<-c(rep_p_ext,rep_id)}
  sc<-sc+1
  print (c(sc,length(rep_p_ext)))
}


#######Overall persistence of invaders across all simulations
pers_data<-read.csv('./results/tti_pers_complete.csv',header=T)
pers_data<-pers_data[which(!(pers_data$rep%in%rep_p_ext)),]

all_vir<-pers_data[pers_data$vir_free=='vir',]
all_vir_c<-all_vir[all_vir$scenario=='control',]
all_vir_i<-all_vir[all_vir$scenario=='invasion_p',]

100*(sum(all_vir_i$pers>=50)/nrow(all_vir_i))
hist(all_vir_i$pers)
summary(all_vir_i$pers)
sd(all_vir_i$pers)
sd(all_vir_c$pers)


###null expectation 
pers_data_c<-pers_data[pers_data$scenario=='control',]
vir_nat<-pers_data_c[pers_data_c$vir_free=='vir',]$pers

pers_data_i<-pers_data[pers_data$scenario=='invasion_p',]
vir_tt<-pers_data_i[pers_data_i$vir_free=='invader',]$pers


k<-sum(vir_tt>=50)
n<-length(vir_tt)
p<-sum(vir_nat>50)/length(vir_nat)

k = mpfr(as.character(k), 128)
p = mpfr(as.character(p), 128)
n = mpfr(as.character(n), 128)

####probability of having as many highly persistent pathogens as observed
p00 = chooseMpfr(n,k)*(p**k)*(1-p)**(n-k)

sink('various_results.txt',append = F)
print (paste('probability of having as many highly persistent pathogens as observed:',p00))
sink()


#######Comparisons within simulations
reps<-unique(pers_data$rep)
res<-c()
all_pers_change<-c()
for (replicate in reps){
  a<-pers_data[pers_data$rep==replicate,]
  pers_i<-a[which(a[,4]=='invader'),]
  vir<-a[a[,4]=='vir',]
  vir<-vir[vir[,2]=='control',]
  pers_comp_1<-sum(vir[,5]<=as.numeric(pers_i[5]))/nrow(vir) #% species less persistent than the invader
  vir<-vir[vir[,6]>=as.numeric(pers_i[6]),]
  mean_pers_c<-mean(vir[,5])
  sd_pers_c<-sd(vir[,5])
  z_pers<-(pers_i[5]-mean_pers_c)/sd_pers_c
  pers_comp_2<-sum(vir[,5]<=as.numeric(pers_i[5]))/nrow(vir) #% species less persistent than the invader (as much abundant as the invader)
  a<-a[a[,4]=='free',]
  a1<-a[a[,2]=='control',]
  a2<-a[a[,2]=='invasion_p',]
  comm<-intersect(a1[,3],a2[,3])
  perc_comm1<-length(comm)/nrow(a1) #species both in control and invasion
  perc_comm2<-length(comm)/nrow(a2) #same
  a1<-a1[which(a1[,3]%in%comm),]
  a1<-a1[order(a1[,3]),]
  a2<-a2[order(a2[,3]%in%comm),]
  a2<-a2[order(a2[,3]),]
  pers_r<-(a2[,5]/a1[,5])###species by species comparison persistence invasion/control
  pers_change_<-(a2[,5]-a1[,5])/100
  pers_change<-mean(abs(pers_change_))
  all_pers_change<-rbind(all_pers_change,cbind(rep(as.numeric(as.factor(replicate)), nrow(a)),pers_change_))
  cor_pers<-cor.test(a2[,5],a1[,5],method='spearman')##persistence correlation
  red_n<-100*sum(pers_r<1)/length(pers_r)#percentage of species with reduced persistence
  print (replicate)
  res<-rbind(res,c(red_n,pers_comp_1,pers_comp_2,mean_pers_c,sd_pers_c,as.numeric(pers_i[5]),as.numeric(z_pers),mean(pers_r),cor_pers$estimate,pers_change))
  }


sink('various_results.txt',append = F)
pers_red<-100*(sum(all_pers_change[,2]<0)/nrow(all_pers_change))
print (paste('percentage of species with reduction in persistence in invasion vs control:',pers_red))
#37.29

pers_incr<-100*(sum(all_pers_change[,2]>0)/nrow(all_pers_change))
print (paste('percentage of species with increase in persistence in invasion vs control:',pers_incr))
#37.56

pers_eq<-100*(sum(all_pers_change[,2]==0)/nrow(all_pers_change))
print (paste('percentage of species with equal persistence in invasion vs control:',pers_eq))
#25.15

pers_less_mean<-100*mean(res[,3],na.rm=T)
pers_less_sd<-100*sd(res[,3],na.rm=T)
print (paste('mean % of native viruses less persistent than the invader:',pers_less_mean,'+/-',pers_less_sd))
#33.6 +/- 39


ok_pers_comp<-res[which(is.finite(res[,3])),]

pers_50<-100*sum(ok_pers_comp[,3]>0.5)/nrow(ok_pers_comp) 
print (paste('% of simulations where invader persisted more than 50% of other viruses:',pers_50))
#29.8%

pers_90<-100*sum(ok_pers_comp[,3]>0.9)/nrow(ok_pers_comp) #% of simulations where invader persisted more than 90% of other viruses
print (paste('% of simulations where invader persisted more than 90% of other viruses:',pers_90))
sink()


pdf('./plots/comparison_persistence_inv_vs_native.pdf',height=4,width=4)
hist(res[,3],xlab='fraction of native viruses\n less persistent than the invader',
     las=1,cex.axis=1.2,cex.lab=1.2,main='')
dev.off()


#######################################################################################
##############2 - Explore invaders' persistence and effects on native communities
#######################################################################################

#limit to simulations where parasites where present 
a<-sum_data[which(!(sum_data$rep%in%rep_p_ext)),]

a$step<-a$step-a$i_t
a<-a[a$step>=0,]
a<-a[is.finite(a$step),]
steps<-unique(a$step)
reps<-unique(a$rep)

###plot the fraction of simulation reaching a certain step
sce<-'invasion_p'
n_rep<-length(unique(a$rep))

sink('various_results.txt',append = T)
print (paste('number of simulations where pathogen did not go extinct:', n_rep,'/1100'))
sink()

pers<-data.frame()
b<-a[a$scenario==sce,]
for (step in steps){
  bb<-b[b$step==step,]
  pers_p<-sum(bb$sp_tti_p>0)/n_rep
  pers<-rbind(pers,data.frame('step'=step,'pers_p'=pers_p))
}


pdf('./plots/invader_persistence.pdf',width=16,height=4)
par(mfrow=c(1,4))
par(cex.axis=1.2, cex.lab=1.2, cex.main=1.2,las=1)
plot(pers$step,100*pers$pers_p,type='l',xlim=c(0,50000),ylim=c(0,100),xlab='time to/from invasion',ylab='% rep invaders persisted',main='a',cex.main=2)

#then plot invaders' abundance across all the different invasion experiments 
plot(b$step,log(b$tot_tti_p+1),type='n',xlab='time to/from invasion',ylab="log(invaders' abundance)",xlim=c(0,50000),main='b',cex.main=2)
for (rep_id in unique(a$rep)){
  bb<-b[b$rep==rep_id,]
  pers_col<-viridis(11,alpha=0.7)[12-round(10*max(bb$step[bb$sp_tti_p>0])/50000)+1]
  lines(bb$step,log(bb$tot_tti_p+1),col=pers_col,lwd=1.5)
}

plot(b$step,b$sp_tti_p,type='n',xlim=c(0,50000),xlab='time to/from invasion',ylab="invaders' genetic diversity",main='c',cex.main=2)
for (rep_id in unique(a$rep)){
  bb<-b[b$rep==rep_id,]
  pers_col<-viridis(11,alpha=0.7)[12-round(10*max(bb$step[bb$sp_tti_p>0])/50000)+1]
  lines(bb$step,bb$sp_tti_p,col=pers_col,lwd=1.5)
}


color.bar(viridis(100), min=0,max=50000,title='steps',nticks=5)
dev.off()


#####
pers_vals<-c()
for (rep in reps){
  bb<-b[b$rep==rep,]
  pers_vals<-c(pers_vals,max(bb$step[bb$sp_tti_p>0]))
  }
pers_vals[!is.finite(pers_vals)]<-0

tre<-quantile(pers_vals,0.975)

sink('various_results.txt',append = T)
print (paste('percentage of simulations where the invader persisted more than',tre,'steps:', 100*sum(pers_vals>=tre)/length(unique(b$rep))))

sink()
#3.07

#filter for the experiments where invaders persisted
reps<-unique(a$rep)
reps_ok<-c()
for (sce in c('control','invasion_p')){
  b<-a[a$scenario==sce,]
  for (rep in reps){
    bb<-b[b$rep==rep,]
    pers_p<-max(bb$step[bb$sp_tti_p>0])
    if (pers_p>=tre){
      reps_ok<-c(reps_ok,rep)
    }
  }
}

reps_ok<-unique(reps_ok) #those are the replicates where invaders persisted for more than 50k steps
rep_ok_n<-length(reps_ok)


###plot invaders' prevalence in the invaded communities
a<-sum_data[which(!(sum_data$rep%in%rep_p_ext)),]
a$step<-a$step-a$i_t
a<-a[a$step>=0,]
a<-a[a$step<=50000,]
inv_all<-a[a$scenario=='invasion_p',]


prev_div<-100*inv_all$sp_tti_p/inv_all$sp_p
prev_ab<-100*inv_all$tot_tti_p/inv_all$tot_p

ok<-which(prev_div>0)
steps<-inv_all$step[ok]

prev_div<-prev_div[ok]
prev_ab<-prev_ab[ok]


pdf('./plots/fig2_invader_prevalence_succ_sims_mean.pdf',width=5,height=5)
par(mfrow=c(1,1))
plot_ci_line(steps,prev_div,xlab='steps',ylab='average prevalence of persistent invaders (%)',y0=0,y1=80) #only successful invasions
#plot_ci_line(steps,prev_ab,xlab='steps',ylab='average abundance of persistent invaders',y0=0,y1=75)
dev.off()


sink('various_results.txt',append = T)
print (paste('average prevalence of survival invaders at step 50k:',round(mean(prev_div[which(inv_succ$step==50000)]),1),'S.D.:',
       round(sd(prev_div[which(inv_succ$step==50000)]),1)))
sink()



###Explore diversity loss
####plot mean loss of diversity and abundance for 2500 steps after the invasion
a<-sum_data[which(!(sum_data$rep%in%rep_p_ext)),]
a$step<-a$step-a$i_t
a<-a[a$step>=0,]
a<-a[a$step<=2500,]
inv_all<-a[a$scenario=='invasion_p',]
###simulations where invaders persisted more than 50k
inv_succ<-inv_all[which(inv_all$rep %in% reps_ok),]
#all simulations
inv_no_succ<-inv_all

control_all<-a[a$scenario=='control',]
control_succ<-control_all[which(control_all$rep %in% reps_ok),]
control_no_succ<-control_all


pdf('./plots/diversity_loss_.pdf',width=5,height=5)
par(mfrow=c(1,1))
par(cex.axis=1.2, cex.lab=1.2, cex.main=1.2,las=1)
sink('various_results.txt',append = T)
div_loss_succ<-100*(control_succ$sp_h-inv_succ$sp_h)/control_succ$sp_h
div_loss_no_succ<-100*(control_no_succ$sp_h-inv_no_succ$sp_h)/control_no_succ$sp_h
print ('diversity loss invasions >50k:')
plot_ci_line(inv_succ$step,div_loss_succ,ylab='diversity loss')
print ('diversity loss invasions <50k:')
plot_ci_line(inv_no_succ$step,div_loss_no_succ,add='y',color_n = 2,color_id = 2)
sink()
dev.off()


pdf('./plots/PD_and_abundance_loss_.pdf',width=8,height=4)
par(mfrow=c(1,2))
par(cex.axis=1.2, cex.lab=1.2, cex.main=1.2,las=1)
sink('various_results.txt',append = T)

ab_loss_succ<-100*(control_succ$tot_h-inv_succ$tot_h)/control_succ$tot_h
ab_loss_no_succ<-100*(control_no_succ$tot_h-inv_no_succ$tot_h)/control_no_succ$tot_h
print ('abundance loss invasions >50k:')
plot_ci_line(inv_succ$step,ab_loss_succ,ylab='abundance loss')
print ('abundance loss invasions <50k:')
plot_ci_line(inv_no_succ$step,ab_loss_no_succ,add='y',color_n = 2,color_id = 2)

pd_loss_succ<-100*(control_succ$phylo_div-inv_succ$phylo_div)/control_succ$phylo_div
pd_loss_no_succ<-100*(control_no_succ$phylo_div-inv_no_succ$phylo_div)/control_no_succ$phylo_div
print ('phylo div loss invasions >50k:')
plot_ci_line(inv_succ$step,pd_loss_succ,ylab='phylo div loss')
print ('phylo div loss invasions <50k:')
plot_ci_line(inv_no_succ$step,pd_loss_no_succ,add='y',color_n = 2,color_id = 2)
sink()
dev.off()


###compute relative loss of abundance, phylo_div and H for comparison
div_loss_all<-100*(control_all$sp_h-inv_all$sp_h)/control_all$sp_h
pd_loss_all<-100*(control_all$phylo_div-inv_all$phylo_div)/control_all$phylo_div
ab_loss_all<-100*(control_all$tot_h-inv_all$tot_h)/control_all$tot_h
sh_loss_all<-100*(control_all$sh_h-inv_all$sh_h)/control_all$sh_h

div_comp_data<-cbind(div_loss_all,pd_loss_all,ab_loss_all,sh_loss_all)

cor(div_comp_data)**2

pdf('./plots/comparison_diversity_measures_loss.pdf',width=10,height=4)
par(mfrow=c(1,3))
plot(div_loss_all,pd_loss_all,xlab='species richness',ylab='PD',cex.lab=1.2,cex.axis=1.2,las=1,pch=16)
plot(div_loss_all,ab_loss_all,xlab='species richness',ylab='abundance',cex.lab=1.2,cex.axis=1.2,las=1,pch=16)
plot(div_loss_all,sh_loss_all,xlab='species richness',ylab='Shannon',cex.lab=1.2,cex.axis=1.2,las=1,pch=16)
#plot(div_loss_all,equ_loss_all,xlab='species richness',ylab='equitability',cex.lab=1.2,cex.axis=1.2,las=1,pch=16)

dev.off()


#####DIVERSITY LOSS QUANTIFICATION
###Compare AUCs of diversity in control vs. invasion simulations to assess potential diversity loss due to invasions

###This part includes also the comparisons of CONTROL VS CONTROL, made to verify that the results cannot be explained by 
#expected variability between replicates of the same simulation (with different random seed)

a<-sum_data[which(!(sum_data$rep%in%rep_p_ext)),]
reps_a<-unique(a$rep)

cc<-read.csv('./results/summary_simple.csv',header=T)
#reps_cc<-unique(cc$rep)

# ####identify simulations where parasites went extinct; those will not be considered
# rep_p_ext_cc<-c()
# sc<-0
# for (rep_id in reps_cc){
#   sub_cc<-cc[cc$rep==rep_id,]
#   inv_time<-sub_cc$i_t[1]
#   sub_cc<-sub_cc[sub_cc$step>2000,]
#   sub_cc<-sub_cc[sub_cc$step<inv_time,]
#   if (min(sub_cc$tot_p)==0){rep_p_ext_cc<-c(rep_p_ext_cc,rep_id)}
#   sc<-sc+1
#   print (c(sc,length(rep_p_ext_cc)))
# }
#cc<-cc[which(!(cc$rep%in%rep_p_ext_cc)),]

reps_cc<-sample(unique(cc$rep),length(reps_a))
div_loss_cc<-c()
for (rep in reps_cc){
  data<-cc[cc$rep==rep,]
  i_time<-data$i_t
  con<-data[data$scenario=='control',]
  con<-con[con$step>=con$i_t,]
  inv<-data[data$scenario=='invasion_p',]
  con<-con[1:25,]
  inv<-inv[1:25,]
  x<-con$step
  y1<-con$sp_h
  y2<-inv$sp_h
  a1<-trapz(1:25,y1)
  a2<-trapz(1:25,y2)
  l_auc<-100*(a1-a2)/a1
  l_mean<-mean(100*(y1-y2)/y1)
  l_max<-max(100*(y1-y2)/y1)
  l_min<-min(100*(y1-y2)/y1)
  div_loss_cc<-rbind(div_loss_cc,c(l_auc,l_mean,l_max,l_min,rep))
}




div_loss_ic<-c()
for (rep in reps_a){
  data<-a[a$rep==rep,]
  i_time<-data$i_t
  con<-data[data$scenario=='control',]
  con<-con[con$step>=con$i_t,]
  inv<-data[data$scenario=='invasion_p',]
  con<-con[1:25,]
  inv<-inv[1:25,]
  x<-con$step
  y1<-con$sp_h
  y2<-inv$sp_h
  a1<-trapz(1:25,y1)
  a2<-trapz(1:25,y2)
  l<-100*(a1-a2)/a1
  l_auc<-100*(a1-a2)/a1
  l_mean<-mean(100*(y1-y2)/y1)
  l_max<-max(100*(y1-y2)/y1)
  l_min<-min(100*(y1-y2)/y1)
  div_loss_ic<-rbind(div_loss_ic,c(l_auc,l_mean,l_max,l_min,rep))
}



###proportion leading to high losses or gains
###invasion vs control
aucs<-as.numeric(div_loss_ic[,1])
#note that since we are computing change as  (control-invasion)/control, positive values indicate losses of diversity, while negative values are gains 
sink('various_results.txt',append = T)
print (paste('% of losses>10%',100*(sum(aucs>10)/nrow(div_loss_ic))))
print (paste('% of gains>10%',100*(sum(aucs<(-10))/nrow(div_loss_ic))))
print (paste('tot % of gains or losses >10%',100*(sum(aucs<(-10))/nrow(div_loss_ic))+
               100*(sum(aucs>10)/nrow(div_loss_ic))))
print (paste('average loss',mean(aucs)))
print (paste('max loss',max(aucs)))
print (paste('max gain',min(aucs)))
print (paste('sd loss',sd(aucs)))
print (paste('% sims with loss',100*(sum(aucs>0)/nrow(div_loss_ic))))
print (paste('% sims with gains',100*(sum(aucs<=0)/nrow(div_loss_ic))))


###only succ
hpers_inv<-aucs[which(div_loss_ic[,5]%in% reps_ok)]
lowpers_inv<-aucs[which(!(div_loss_ic[,5]%in% reps_ok))]
print (paste('only high persistent invaders - % sims with loss',100*(sum(hpers_inv>0)/length(hpers_inv))))
print (paste('only high persistent invaders - % sims with gains',100*(sum(hpers_inv<=0)/length(hpers_inv))))

print (paste('only high persistent invaders - mean losses',mean(hpers_inv)))
print (paste('only high persistent invaders - sd losses',sd(hpers_inv)))
print (paste('only high persistent invaders - max losses',max(hpers_inv)))
print (paste('only high persistent invaders - max gains',min(hpers_inv)))



###control VS control
cc_auc<-as.numeric(div_loss_cc[,1])
ic_auc<-as.numeric(div_loss_ic[,1])
range_diff<-(max(ic_auc)-min(ic_auc))/(max(cc_auc)-min(cc_auc))

print (paste('control vs control - mean losses',mean(cc_auc)))
print (paste('control vs control - sd losses',sd(cc_auc)))
print (paste('control vs control - max losses',max(cc_auc)))
print (paste('control vs control - max gains',min(cc_auc)))

sink()


pdf('./plots/cc_vs_ic.pdf',width=5,height=6)
par(mfrow=c(1,1))
vioplot(-cc_auc,-ic_auc,names=c('control vs\n control','control vs\n invasion'),las=1,cex.axis=1.2,cex.lab=1.2,
        ylab='relative change (%)\n<--losses - gains-->')
dev.off()




##############################################################################################
##############3 - Explore invaders' features increasing the chances of successful invasion?
##############################################################################################

#limit to simulations where parasites where present 
a<-sum_data[which(!(sum_data$rep%in%rep_p_ext)),]

a<-a[is.finite(a$step),]
reps<-unique(a$rep)

b<-a[a$scenario=='invasion_p',]
a<-a[a$scenario=='control',]

pers_data<-c()
for (rep in reps){
 a1<-a[a$rep==rep,]
 b1<-b[b$rep==rep,]
 i_time<-ceiling(b1$i_t[1]/1000)*1000
 wsize<-a1$size[1]
 res_n<-a1$res_n[1]
 eff_res<-a1$eff_res[1]
 source_time<-b1$i_age[1]
 x_<-b1$step[b1$step>=i_time & b1$step<=(i_time+2500)]
 y1<-a1[a1$step %in% x_,]$sp_h
 y2<-b1[b1$step %in% x_,]$sp_h
 auc_con<-trapz(seq(0,2500,100),y1)
 auc_inv<-trapz(seq(0,2500,100),y2)
 auc_loss<-100*(auc_inv-auc_con)/auc_con
 if (sum(b1$sp_tti_p>0)>0){
 pers_p<-max(b1$step[b1$sp_tti_p>0])-i_time} else {pers_p<-0}#persistence of the invader
 pers_p_con<-sum(a1$p_i_lin_ab[which(a1$step<i_time)])/sum(a1$p_i_lin_ab[which(a1$step<i_time)]>0)#mean abundance of the invader in the steps it was present;
 pers_p_l<-sum(a1$p_i_lin_ab[which(a1$step<i_time)]>0)# n steps invader was present pre-invasion;
 pers_p_tot<-sum(a1$p_i_lin_ab[which(a1$step<i_time)])#tot ind n invader pre invasion
 inv_p_div0<-a1[a1$step==source_time,]$p_i_lin_ab ###note that this is the abundance of the invader in the source community!
 ###compute time diff
 last_step <- max(a1[which(a1$inv_p_sp_ind_n>0),]$step)
 t_diff<-i_time-last_step
 
 pers_data<-rbind(pers_data,data.frame('sce'='invasion_p',
                                       'rep'=rep,
                                       'wsize' = wsize,
                                       'eff_res' = eff_res,
                                       'res_n' = res_n,
                                       'auc_loss'=auc_loss,
                                       'pers_p'=pers_p,
                                       'pers_p_con'=pers_p_con,
                                       'pers_p_l'=pers_p_l,
                                       'pers_p_tot'=pers_p_tot,
                                       'tdiff'=t_diff,#(b1$i_t[1]-b1$i_age[1]),
                                       'age'=b1$i_age[1],
                                       'i_t'=b1$i_t[1],
                                       'inv_p_div0'=inv_p_div0,
                                       'net_l_i'=a1[a1$step == i_time-1000,]$net_l,
                                       'net_fill_i'=a1[a1$step == i_time-1000,]$net_fill,
                                       'mean_deg_p_i'=a1[a1$step == source_time,]$i_p0_mean_deg, ###this is the mean degree of the target invader+its descendants in the source community
                                       'sh_h0'=a1[a1$step == i_time-1000,]$sh_h,
                                       'sh_p0'=a1[a1$step == i_time-1000,]$sh_p,
                                       'div_h0'=a1[a1$step == i_time-1000,]$sp_h,
                                       'div_p0'=a1[a1$step == i_time-1000,]$sp_p,
                                       'n_h0'=a1[a1$step == i_time-1000,]$tot_h,
                                       'n_p0'=a1[a1$step == i_time-1000,]$tot_p,
                                       'p_spec_i_time'=b1$i_p_mean_deg[b1$step == b1$i_t[1]][1]
 ))
}

###Add phylogenetic distance and phenotype information
phyl<-read.csv('./results/phylo_dist.csv',header=T)
phyl<-phyl[which(!(phyl$rep%in%rep_p_ext)),]
phyl<-phyl[which(phyl$rep_id%in%pers_data$rep),]
phyl<-phyl[match(as.character(pers_data$rep),as.character(phyl$rep_id)),]
all.equal(as.character(phyl$rep_id),as.character(pers_data$rep))

pers_data<-cbind(pers_data,phyl)
pers_data$mean_dist<-as.numeric(as.character(pers_data$mean_dist))
pers_data$max_dist<-as.numeric(as.character(pers_data$max_dist))
pers_data$min_dist<-as.numeric(as.character(pers_data$min_dist))
pers_data$std_dist<-as.numeric(as.character(pers_data$std_dist))

tasks<-read.csv('./results/inv_tasks.csv',header=F)
tasks<-tasks[which(!(tasks[,1]%in%rep_p_ext)),]
tasks<-tasks[which(tasks[,1]%in%pers_data$rep),]
tasks<-tasks[match(as.character(pers_data$rep),as.character(tasks[,1])),]
all.equal(as.character(tasks[,1]),as.character(pers_data$rep))
pers_data<-cbind(pers_data,data.frame('task_n'=tasks[,2],'task_comp'=tasks[,3]))


###############Explore loss and change vs. invaders' persistence

pdf('./plots/pers_vs_loss_boxplots.pdf',height=4,width=9)
par(mfrow=c(1,2))
rho<-round(cor.test(pers_data$auc_loss,pers_data$pers_p,method='spearman')$estimate,2)
p<-round(cor.test(pers_data$auc_loss,pers_data$pers_p,method='spearman')$p.value,3)
if (p==0){p<-'<0.001'} else (p<-paste('p =',p))
boxplot(pers_data$auc_loss~round(log(pers_data$pers_p+1)),outline=T,
        xlab='log(invader persistence)',ylab='change in species richness (%)',main=paste('rs =',rho,p))
abline(h=0,lty=2)
rho<-round(cor.test(abs(pers_data$auc_loss),pers_data$pers_p,method='spearman')$estimate,2)
p<-round(cor.test(abs(pers_data$auc_loss),pers_data$pers_p,method='spearman')$p.value,3)
if (p==0){p<-'p<0.001'} else (p<-paste('p =',p))

boxplot(abs(pers_data$auc_loss)~round(log(pers_data$pers_p+1)),outline=T,
        xlab='log(invader persistence)',ylab='absolute change in species richness (%)',main=paste('rs =',rho,p))
abline(h=0,lty=2)

dev.off()






#####################################################
#######Pairwise comparisons
####################################################

x<-cbind(pers_data$div_h0,
         pers_data$div_p0,
         pers_data$n_h0,
         pers_data$n_p0,
         pers_data$n_h0/pers_data$wsize,
         pers_data$n_p0/pers_data$wsize,
         pers_data$pers_p_l,
         pers_data$pers_p_con,
         pers_data$pers_p_tot,
         pers_data$inv_p_div0,
         pers_data$mean_deg_p_i,
         pers_data$task_n,
         pers_data$task_comp,
         pers_data$mean_dist,
         pers_data$tdiff,
         pers_data$age,
         pers_data$i_t
)

dim(x)
colnames(x)<-c('free-living richness',
               'pathogen richness',
               'free-living abundance',
               'pathogen abundance',
               'free-living density',
               'pathogen density',
               'invader persistence (pre-invasion)',#n steps invader was present pre-invasion;
               'invader mean abudance (pre-invasion)',#mean abundance of the invader in the steps it was present;
               'invader total abundance (pre-invasion)',#tot ind n invader pre invasion
               'invader abundance in source community',#abundance of the invader in the source community
               'invader generalism',
               'task diversity',
               'max task complexity',
               'average phylogenetic distance',
               'time difference',
               'age of invaded community',
               'invasion time'
          )


y<-cbind(pers_data$pers_p,pers_data$auc_loss,abs(pers_data$auc_loss))
colnames(y)<-c('invader_persistence','diversity change','abs(diversity change)')

cor_vals<-c()
for (i in 1:ncol(x)){
  row<-c()
  for (j in 1:ncol(y)){
  row<-c(row,round(cor.test(x[,i],y[,j],method='spearman')$estimate,2),
         round(cor.test(x[,i],y[,j],method='spearman')$p.value,4))
      }
  cor_vals<-rbind(cor_vals,row)
  }

row.names(cor_vals)<-colnames(x)
colnames(cor_vals)<-c("invader_persistence rs","invader_persistence p","div_change rs","div_change p","abs_div_change rs","abs_div_change p")

write.table(cor_vals,'./plots/pairwise_table.csv',quote=F,sep=',',row.names=T,col.names=T)

xy<-cbind(x,y)
pdf('./plots/correlations.pdf')
corrplot(cor(xy,method='spearman',use="complete.obs"),is.corr=T,method = 'square')
dev.off()


#####################################
###check high lossess and gains compare to persistence
h_losses<-which(y[,2]<(-5))
h_gains<-which(y[,2]>=5)
length(h_losses)
length(h_gains)

pdf('./plots/hloss_hgain_vars_comparison_main.pdf',width = 4.5,height=4.5)
par(mfrow=c(1,1))
hlv<-y[h_losses,1]
hgv<-x[h_gains,1]
p<-round(wilcox.test(hgv, hlv, alternative = "two.sided")$p.value,3)
boxplot(hlv,hgv,names=c('>5% loss','>5% gain'),
        ylab='',main=paste(colnames(y)[1],'p = ',p),las=1,cex.lab=1.2,cex.axis=1.2)

# var<-which(colnames(x)=='n_p0')
# hlv<-x[h_losses,var]
# hgv<-x[h_gains,var]
# p<-round(wilcox.test(hgv, hlv, alternative = "two.sided")$p.value,3)
# boxplot(hlv,hgv,names=c('>10% loss','>10% gain'),
#         ylab='',main=paste(colnames(x)[var],' p =',p),las=1,cex.lab=1.2,cex.axis=1.2)
dev.off()


pdf('./plots/hloss_hgain_vars_comparison_SI.pdf',width = 10,height=10)
par(mfrow=c(4,4))
for (var in 1:ncol(x)){
  if (colnames(x)[var]!='n_p0'){
    hlv<-x[h_losses,var]
    hgv<-x[h_gains,var]
    p<-round(wilcox.test(hgv, hlv, alternative = "two.sided")$p.value,3)
    boxplot(hlv,hgv,names=c('>5% loss','>5% gain'),
            ylab='',main=paste(colnames(x)[var],' p =',p),las=1,cex.lab=1.2,cex.axis=1.2,outline=F)
  }
  }
dev.off()






# ####Random Forest
library(randomForest)
citation("randomForest")
dens_h<-pers_data$n_h0/pers_data$wsize
dens_p<-pers_data$n_p0/pers_data$wsize
pers_data<-cbind(pers_data,'dens_h'=dens_h,'dens_p'=dens_p)


sink('random_forest_results.txt')
pdf('random_forest_figure.pdf',width=21,height=7)
par(mfrow=c(1,3))
rf<-randomForest(pers_p~
                   div_h0+
                   div_p0+
                   n_h0+
                   n_p0+
                   dens_h+
                   dens_p+
                   pers_p_l+
                   pers_p_con+
                   pers_p_tot+
                   inv_p_div0+
                   mean_deg_p_i+
                   task_n+
                   task_comp+
                   mean_dist+
                   tdiff+
                   age+
                   i_t,ntree=1000,localImp=T,data=pers_data,na.action = na.roughfix)#,data=train)


imp<-importance(rf)[,1]
names(imp)<-c('free-living richness',
               'pathogen richness',
               'free-living abundance',
               'pathogen abundance',
               'free-living density',
               'pathogen density',
               'invader persistence (pre-invasion)',#n steps invader was present pre-invasion;
               'invader mean abudance (pre-invasion)',#mean abundance of the invader in the steps it was present;
               'invader total abundance (pre-invasion)',#tot ind n invader pre invasion
               'invader abundance in source community',#abundance of the invader in the source community
               'invader generalism',
               'task diversity',
               'max task complexity',
               'average phylogenetic distance',
               'time difference',
               'age of invaded community',
               'invasion time'
)


barplot(rev(sort(imp[1:10])),horiz=T,las=1,xlab='% MSE Increase',main="invaders' persistence")
print(rf)

####signed change
rf<-randomForest(auc_loss~
                   pers_p+
                   div_h0+
                   div_p0+
                   n_h0+
                   n_p0+
                   dens_h+
                   dens_p+
                   pers_p_l+
                   pers_p_con+
                   pers_p_tot+
                   inv_p_div0+
                   mean_deg_p_i+
                   task_n+
                   task_comp+
                   mean_dist+
                   tdiff+
                   age+
                   i_t,ntree=1000,localImp=T,data=pers_data,na.action = na.roughfix)#,data=train)


imp<-importance(rf)[,1]
names(imp)<-c('invader persistence (post-invasion)',
              'free-living richness',
              'pathogen richness',
              'free-living abundance',
              'pathogen abundance',
              'free-living density',
              'pathogen density',
              'invader persistence (pre-invasion)',#n steps invader was present pre-invasion;
              'invader mean abudance (pre-invasion)',#mean abundance of the invader in the steps it was present;
              'invader total abundance (pre-invasion)',#tot ind n invader pre invasion
              'invader abundance in source community',#abundance of the invader in the source community
              'invader generalism',
              'task diversity',
              'max task complexity',
              'average phylogenetic distance',
              'time difference',
              'age of invaded community',
              'invasion time'
)


barplot(rev(sort(imp[1:10])),horiz=T,las=1,xlab='% MSE Increase',main='signed diversity change')
print(rf)


####

rf<-randomForest(abs(auc_loss)~
                   pers_p+
                   div_h0+
                   div_p0+
                   n_h0+
                   n_p0+
                   dens_h+
                   dens_p+
                   pers_p_l+
                   pers_p_con+
                   pers_p_tot+
                   inv_p_div0+
                   mean_deg_p_i+
                   task_n+
                   task_comp+
                   mean_dist+
                   tdiff+
                   age+
                   i_t,ntree=1000,localImp=T,data=pers_data,na.action = na.roughfix)#,data=train)


imp<-importance(rf)[,1]
names(imp)<-c('invader persistence (post-invasion)',
              'free-living richness',
              'pathogen richness',
              'free-living abundance',
              'pathogen abundance',
              'free-living density',
              'pathogen density',
              'invader persistence (pre-invasion)',#n steps invader was present pre-invasion;
              'invader mean abudance (pre-invasion)',#mean abundance of the invader in the steps it was present;
              'invader total abundance (pre-invasion)',#tot ind n invader pre invasion
              'invader abundance in source community',#abundance of the invader in the source community
              'invader generalism',
              'task diversity',
              'max task complexity',
              'average phylogenetic distance',
              'time difference',
              'age of invaded community',
              'invasion time'
)


barplot(rev(sort(imp[1:10])),horiz=T,las=1,xlab='% MSE Increase',main = 'absolute diversity change')
print(rf)

dev.off()
sink()


