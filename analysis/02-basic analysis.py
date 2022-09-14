#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 19:41:21 2019

@author: saeedeh
"""
all_trials.min_rr.isna().groupby(all_trials.id).mean()
all_trials.near_threshold.groupby(all_trials.id).mean()


#


################
####### My last

#behavioral
print(ro.r('summary(lmer(  min_to_max_OR_fixed~(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  min_to_max_OR_flex~RT_real+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  min_to_max_OR_flex~duration+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  RT_real~ duration+(1|id),  data=rdf))'))  

#bias
print(ro.r('summary(lmer(  TP_bias2~min_to_max_OR_flex+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  TP_bias2~max_time+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  max_time~onset_to_max_OR_flex+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  TP_bias2 ~ max_time+onset_to_max_OR_flex+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  TP_bias2 ~ min_rr_fixed+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  min_to_max_OR_flex~TP_bias2 + fixation_time+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  max_time ~ fixation_time+(1|id),  data=rdf))'))  

print(ro.r('summary(lmer(  acc~min_to_max_OR_flex+(1|id),  data=rdf))'))  
print(ro.r('summary(lmer(  TP_acc2~RT_real+(1|id),  data=rdf))'))  


#consistency
print(ro.r('summary(glmer(  TP_acc2~minRR_controlledCOR+(1|id), family="binomial"  , data=rdf))'))  
print(ro.r('summary(glmer(  TP_acc2~max_time+(1|id), family="binomial"  , data=rdf))'))  
print(ro.r('summary(glmer(  TP_acc2~min_to_max_OR_flex+(1|id), family="binomial"  , data=rdf))'))  
print(ro.r('summary(glmer(  TP_acc2~post_rr+(1|id), family="binomial"  , data=rdf))'))  
print(ro.r('summary(glmer(  TP_acc2~min_rr_fixed+max_time+RT_real+(1|id), family="binomial"  , data=rdf))')) 
print(ro.r('summary(glmer(  TP_acc2~RT_real+(1|id), family="binomial"  , data=rdf))')) 

print(ro.r('summary(lmer(  min_post_avg_rr ~ RT_real+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  RT_real~min_rr_fixed+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  max_time~RT_real+(1|id) , data=rdf))')) 

print(ro.r('summary(glmer(  TP_acc2~min_post_avg_rr2 +(1|id), family="binomial"  , data=rdf))')) 
print(ro.r('summary(glmer(  TP_acc2~max_time +(1|id), family="binomial"  , data=rdf))')) 


print(ro.r('summary(lmer( TP_acc2~ max_time+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  min_post_avg_rr2~max_time+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  min_rr~max_time+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  min_rr_fixed~min_to_max_OR_flex+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  TP_bias2~min_to_max_OR_flex*min_post_avg_rr+(1|id) , data=rdf))')) 

#higher heart rate means more orienting and for shorter!!!

# high-min-rr
print(ro.r('summary(lmer(  TP_acc2~high_min_post_avg_rr+(1|id) , data=rdf))')) 

### COR and minRR decorrelated
print(ro.r('summary(lmer(  long_resp~minRR_controlledCOR+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  min_rr_fixed~CORslope_controlledMinRR+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  TP_bias2 ~ onset_to_max_OR_flex+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  long_resp~duration + COR_onset_controlledMinRR+(1|id) , data=rdf))')) 
#this is worse than not controlling. we'd better not control

### block analysis
print(ro.r('summary(lmer(  mean_rr~mean_COR+(1|id) , data=block_rdf))')) 
print(ro.r('summary(lmer(  TP_bias2~mean_COR+(1|id) , data=block_rdf))')) 

### onset to RT COR
print(ro.r('summary(lmer(  TP_bias~onset_to_RT_COR+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  TP_bias2~COR_onset_slope+(1|id) , data=rdf))')) 

### predicting long from minRR!
print(ro.r('summary(lmer(  TP_bias2~min_to_max_OR_flex+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(acc~TP_bias+(1|id), data=rdf))')) 

print(ro.r('summary(lmer(  long_resp~signed_minRR_thr+ duration+correct_long_thr+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  long_resp~signed_minRR+(1|id)+(1|duration) , data=rdf))')) 
ro.r('rdf$RT_times_slope <- rdf$COR_onset_slope*rdf$RT')
print(ro.r('summary(lmer(  long_resp~ duration+COR_onset_slope+(1|id) , data=rdf))')) 

print(ro.r('summary(lmer(  RT_real~ min_to_max_OR_flex+(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  RT_real~ min_rr_fixed+(1|id) , data=rdf))')) 

## set by step

#H1:
# =============================================================================
# Larger RR: 
# -less accuracy
# -larger RT
# -more bias towards long
# Larger COR: 
# -bias towards long
# =============================================================================

print(ro.r('summary(lmer(  long_resp~ duration+min_rr_fixed +(1|id) , data=rdf))')) 
print(ro.r('summary(lmer(  long_resp~ duration+min_rr_fixed*duration+onset_to_max_OR_flex +(1|id) , data=rdf))')) 
  
#####################################################
### Bisection plot
df=all_trials
res=df.groupby([levels_name,'id'])[resp_name].mean();
res=res.reset_index()    
g=sns.lineplot(data=res, x="duration", y="long_resp", marker='o', color="black")
g.set(xticks=np.unique(res.duration))
plt.xlabel('stimulus duration (ms)')
plt.ylabel('p(long)')
plt.tight_layout()
    
#generic
mname='onset_to_max_OR_flex'
mname_label='COR'

mname='RT_real'
mname_label='RT'

mname='min_rr_fixed'
mname_label='pre-RR'

mname='max_time'
mname_label='peak time'

mname='onset_to_RT_COR'
mname_label='COR until RT'

avg_m= all_trials.groupby(['id','duration'])[mname].median()

avg_m_UB= all_trials.groupby(['id'])[mname].quantile(0.50)
avg_m_LB= all_trials.groupby(['id'])[mname].quantile(0.50)
all_trials['high_m']=False
all_trials['low_m']=False

for idx, meas in enumerate(all_trials[mname]):
    #if OR>avg_OR_UB[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
    if meas>avg_m_UB[all_trials.id.iloc[idx]]:
    #if OR>0:
    #if OR>avg_avg_OR:
        all_trials.high_m.iloc[idx]=True
    #if OR<avg_OR_LB[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
    if meas<avg_m_LB[all_trials.id.iloc[idx]]:
    #if OR>0:
    #if OR>avg_avg_OR:
        all_trials.low_m.iloc[idx]=True
plot_bisection(all_trials[all_trials.high_m],'r')  
plot_bisection(all_trials[all_trials.low_m],'b')  
plt.legend(['high'+mname_label,'low'+mname_label])


#1) OR
plt.subplot(1,2,1)
avg_OR= all_trials.groupby(['id','duration'])['onset_to_max_OR_fixed'].median()
avg_OR_UB= all_trials.groupby(['id','duration'])['onset_to_max_OR_fixed'].quantile(0.85)
avg_OR_LB= all_trials.groupby(['id','duration'])['onset_to_max_OR_fixed'].quantile(0.15)

avg_OR_UB= all_trials.groupby(['id'])['onset_to_max_OR_fixed'].quantile(0.75)
avg_OR_LB= all_trials.groupby(['id'])['onset_to_max_OR_fixed'].quantile(0.25)

avg_avg_OR=all_trials.onset_to_max_OR_fixed.mean()
all_trials['high_OR']=False
all_trials['low_OR']=False

for idx, OR in enumerate(all_trials.onset_to_max_OR_fixed):
    #if OR>avg_OR_UB[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
    if OR>avg_OR_UB[all_trials.id.iloc[idx]]:
    #if OR>0:
    #if OR>avg_avg_OR:
        all_trials.high_OR.iloc[idx]=True
    #if OR<avg_OR_LB[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
    if OR<avg_OR_LB[all_trials.id.iloc[idx]]:
    #if OR>0:
    #if OR>avg_avg_OR:
        all_trials.low_OR.iloc[idx]=True
plot_bisection(all_trials[all_trials.high_OR],'r')  
plot_bisection(all_trials[all_trials.low_OR],'b')  
plt.legend(['high COR','low COR'])
##
#2) HR
plt.subplot(1,2,2)
plt.figure(figsize=(3,2.5))
avg_min_rr= all_trials.groupby('id')['min_rr_fixed'].median()
all_trials['high_min_rr']=False
for idx, rr in enumerate(all_trials.min_rr_fixed):
    if rr>avg_min_rr[all_trials.id.iloc[idx]]: all_trials.high_min_rr.iloc[idx]=True
plot_bisection(all_trials[all_trials.high_min_rr],'b')  
plot_bisection(all_trials[~all_trials.high_min_rr],'r')
plt.legend(['large','small'], title="pre-trial RR")  
plt.tight_layout()

#3) RT
avg_RT= all_trials.groupby(['id', 'duration'])['RT_real'].median()
all_trials['high_RT']=False
for idx, rt in enumerate(all_trials.RT_real):
    #if rt>avg_RT[all_trials.id.iloc[idx]]: 
    if rt>avg_RT[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
        all_trials.high_RT.iloc[idx]=True
plt.figure(figsize=(3,2.5))
plot_bisection(all_trials[all_trials.high_RT],'b')  
plot_bisection(all_trials[~all_trials.high_RT],'r')
plt.legend(['slow','fast'], title='RT')  
plt.tight_layout()


#4) max-time
avg_MT= all_trials.groupby(['id', 'duration'])['max_time'].median()
all_trials['high_MT']=False
for idx, MT in enumerate(all_trials.max_time):
    #if rt>avg_RT[all_trials.id.iloc[idx]]: 
    if MT>avg_MT[all_trials.id.iloc[idx]][all_trials.duration.iloc[idx]]:
        all_trials.high_MT.iloc[idx]=True
plt.figure(figsize=(3,2.5))
plot_bisection(all_trials[all_trials.high_MT],'b')  
plot_bisection(all_trials[~all_trials.high_MT],'r')
plt.legend(['late','early'], title='COR Peak time')  
plt.tight_layout()


nan_ttest_ind(all_trials.onset_to_max_OR_fixed[inds2], all_trials.onset_to_max_OR_fixed[np.logical_or(inds4, inds3, inds1)])


