#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 02:46:38 2019

@author: saeedeh
"""
#######################################################################
# avg OR

new_xx= new_x[(new_x>=-2000) & (new_x<=4500)]
arr2=ar2[:,(new_x>=-2000) & (new_x<=4500)]

new_xx= new_x
arr2=ar2

zero_ind=np.where(new_xx==0)#min= -1900
min_ind=np.where(new_xx==-1300)#
t=np.nanmean(arr2,axis=0)

t_avg=t#-t[zero_ind]
plt.plot(new_xx, t_avg, color='k')
plt.xlabel('Time (ms)', weight='bold')
plt.ylabel('RR (ms)', weight='bold')
max_ind=t.argmax() #max=1000
plt.tight_layout()

plt.xticks([-1300,0])
plt.plot(new_xx[max_ind], t_avg[max_ind],'or')
plt.plot(new_xx[zero_ind], t_avg[zero_ind],'or')
plt.plot(new_xx[min_ind], t_avg[min_ind],'or')
#plt.axvline(linewidth=1, color='r',ls='--')
plt.ylim([770,805])

##
# Avg by subject
#subplots
range0=-400; range1=3200
new_xx= new_x[(new_x>=range0) & (new_x<=range1)]
arr2=ar2[:,(new_x>=range0) & (new_x<=range1)]
zero_ind=np.where(new_xx==0)#min= -1900
t=np.nanmean(arr2,axis=0)
t_avg=t-t[zero_ind]
for idx,subj in enumerate(subj_ids):
    sel_inds=all_trials.id==subj
    t=np.nanmean(arr2[sel_inds,:], axis=0)
    t=t-t[zero_ind]
    plt.subplot(5,9,idx+1)
    plt.plot(new_xx, t, color='b')
    plt.plot(new_xx, t_avg,ls='--', color='k')
    #plt.ylim([-23, 38])
    plt.xticks([0],[''])
    plt.yticks([],[])
    #plt.savefig('/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/Analysis/images/OR-avg-subj/'+str(subj)+'.png')
   # plt.close('all')
# all in one
range0=-400; range1=3200
new_xx= new_x[(new_x>=range0) & (new_x<=range1)]
arr2=ar2[:,(new_x>=range0) & (new_x<=range1)]
zero_ind=np.where(new_xx==0)#min= -1900
t=np.nanmean(arr2,axis=0)
t_avg=t-t[zero_ind]
for idx,subj in enumerate(subj_ids):
    sel_inds=all_trials.id==subj
    t=np.nanmean(arr2[sel_inds,:], axis=0)
    t=t-t[zero_ind]
    plt.plot(new_xx, t, color='gray', alpha=0.2, linewidth=2)
    #plt.ylim([-23, 48])
#plt.xticks([0],[''])
#plt.yticks([],[])
arr3=np.copy(arr2)
for idx in range(arr2.shape[0]):
    arr3[idx,:]=arr2[idx,:]-arr2[idx, zero_ind]
    
plt.plot(new_xx, t_avg,ls='--', color='k')
plt.xlabel('Time')
plt.ylabel('RR')

range0=-1300; range1=3500
new_xx= new_x[(new_x>=range0) & (new_x<=range1)]
arr2=ar2[:,(new_x>=range0) & (new_x<=range1)]
new_xx2=np.tile(new_xx, (arr2.shape[0], 1))
g=sns.lineplot(new_xx2.flatten(), arr2.flatten(), color="black")
g.set(xticks=[-1300, 0, 1400, 3500])
#plt.ylim([-20, 60])
    

############################################################################
###   


###### RR signal for the 7 tone durations
colors=['red', 'orange', 'yellow', 'm', 'c', 'b', 'green']
durations=np.unique(all_trials.duration)
for i in range(7):
    inds= (all_trials.duration==durations[i])
    t=subj_mean_mean(arr2, inds)
    plt.plot(new_xx, t, color=colors[i])
plt.legend(['1','2','3','4','5','6','7'])


###########################################
######## making inds #########################



# 1)contraction, consistent, dilation relative to 134 ms
inds1=all_trials.TP_bias2==-1 #contraction
#inds2=all_trials.TP_bias2==0 #consistent
inds2= all_trials.TP_bias2==1#dilation

# 2) consistent, distortion relative to BP (exclude near BP)
inds1=all_trials.TP_acc2==1 #consistent
inds2= (all_trials.TP_acc2==0)

# 3) objective short, long (exclude 134 ms)
inds1= all_trials.duration<134
inds2=all_trials.duration>134

# 4) short vs long min_rr
inds1= (all_trials.high_min_rr==True)
inds2=(all_trials.high_min_rr==False)

# 5) short vs long min_post_avg_rr
inds1= (all_trials.high_min_post_avg_rr==True)
inds2=(all_trials.high_min_post_avg_rr==False)

 #6) first vs last block
inds1= (all_trials.trial_num==1)
inds2= (all_trials.trial_num==2)
inds3= (all_trials.trial_num==3)
inds4= (all_trials.trial_num==4)
#################################################



new_LB= -1500; new_UB=0
new_LB= -1800; new_UB=4900
new_LB= -1300; new_UB=4000

new_xx= new_x[(new_x>=new_LB) & (new_x<=new_UB)]
arr2=ar2[:,(new_x>=new_LB) & (new_x<=new_UB)]

t1=subj_mean_mean(arr2, inds1) #long stim, long resp
t2=subj_mean_mean(arr2, inds2) #short stim, long resp
#t3=subj_mean_mean(arr2, inds3) # long stim, short resp
# t4=subj_mean_mean(arr2, inds4) # short stim, short resp

fig, ax=plt.subplots(figsize=(4*0.9,3*0.9))
plt.locator_params(axis='y', nbins=4)
plt.locator_params(axis='x', nbins=4)
plt.plot(new_xx, t1, color='k',ls='-', alpha=0.5, linewidth=2)
plt.plot(new_xx, t2, color='k',ls='--', linewidth=2)
#plt.plot(new_xx, t3, color='k', ls='-', linewidth=2)
# plt.plot(new_xx, t4, color='red', ls=':')
#plt.ylim([-15,30])
plt.ylim([770,820])
plt.xlim([new_LB,new_UB])
plt.xlabel('time',weight='bold'); plt.ylabel('RR',weight='bold')
plt.tight_layout()
plt.tick_params(axis = "x", bottom = True, top=False, length=3)
#plt.legend(['inds1','inds2'])
save_img_dir='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/Analysis/images/RR-plots2'
img_fname=ylab+'-'+type_str+'.png'
plt.savefig(os.path.join(save_img_dir, img_fname))

## 4 plots
t1=subj_mean_mean(arr2, inds1) 
t2=subj_mean_mean(arr2, inds2) 
t3=subj_mean_mean(arr2, inds3)
t4=subj_mean_mean(arr2, inds4)
plt.plot(new_xx, t1, color='k',ls='-', alpha=1, linewidth=2)
plt.plot(new_xx, t2, color='k',ls='-',alpha=0.7,  linewidth=2)
plt.plot(new_xx, t3, color='k', ls='-', linewidth=2, alpha=0.5)
plt.plot(new_xx, t4, color='k', ls='-', linewidth=2, alpha=0.2)


#with seaborn
new_xx1=np.tile(new_xx, (np.sum(inds1), 1))
new_xx2=np.tile(new_xx, (np.sum(inds2), 1))
df=pd.dataFrame()
g=sns.lineplot(new_xx1.flatten(), arr2[inds1].flatten(), color="black", ci=68)

g=sns.lineplot(new_xx2.flatten(), arr2[inds2].flatten(), color="black", ci=68)


leg=plt.legend(['objective short', 'objective long'])
leg=plt.legend(['contraction', 'dilation'])
leg=plt.legend(['consistent', 'inconsistent'])
leg=plt.legend(['high pre-RR', 'low pre-RR'])
leg=plt.legend(['high pre-post-avg-RR', 'low pre-post-avg-RR'])

leg.get_frame().set_linewidth(0.0)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['left'].set_visible(False)

