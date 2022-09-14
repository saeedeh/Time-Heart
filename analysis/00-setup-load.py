#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 12:29:23 2019

@author: saeedeh
"""
import pandas as pd
from peakdet import Physio, operations, load_rtpeaks
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
import seaborn as sns

all_trials=pd.read_pickle('/home/saeedeh/Desktop/from lab PC/tp heart/all_trials.pkl')
all_subj=pd.read_pickle('/home/saeedeh/Desktop/from lab PC/tp heart/all_subjects.pkl')

#base_data_dir='/mnt/Aclab/Studies/38_TimeHeart'
base_data_dir='/home/saeedeh/Desktop/aclab-server/Studies/38_TimeHeart'

# Average performance plots!

##
ITI=3.5

data_dir=base_data_dir+'/raw data'

fname=base_data_dir+'/analysis/all_subjects.pkl';
all_subj=pd.read_pickle(fname)

fname=base_data_dir+'/analysis/all_trials.pkl';
all_trials=pd.read_pickle(fname)

fname=base_data_dir+'/analysis/df_blocks.pkl';
df_blocks=pd.read_pickle(fname)

if all_trials.RT_real.mean()<10:
    all_trials.RT_real = all_trials.RT_real*1000
if all_trials.RT.mean()<10:
    all_trials.RT = all_trials.RT*1000
all_trials.RT[all_trials.RT<88]=np.nan

if all_subj.thresh.mean()<10:
    all_subj.thresh= all_subj.thresh*1000

if all_trials.duration.mean()<10:
    all_trials.duration *=1000
#
#####################
# Filtering
#####################
all_trials['trial_num']=all_trials['.thisN']
all_trials['block']= (all_trials['.thisN'])//14  #14 trials in each block
all_trials['trial_in_block']= (all_trials['.thisN']-all_trials['block']*14) 
all_trials['block']= all_trials['block']+1
all_trials['trial_in_block']=all_trials['trial_in_block']+1

#filter_inds=np.array( (all_trials['.thisN']<5) | (all_trials['trial_in_block']==1 ))

#filter trials with inaccurate duration (Where do they come from? It's 7% of trials!!!!!)-> that's trial=14!
all_trials['trial_time']=np.nan
all_trials.trial_time.iloc[:-1]=np.array(all_trials.physioTime.iloc[1:])-np.array(all_trials.physioTime.iloc[:-1])
invalid_inds_RT=(np.array(all_trials.trial_time<3500) | np.array(all_trials.trial_time>9000)) 
all_trials.trial_time[np.array(invalid_inds_RT)]=np.nan
all_trials['ITI']=all_trials.trial_time-all_trials.RT #This is ITI+fixation of next trials!
t=all_trials.ITI-3500
all_trials['fixation_time']=np.nan
all_trials.fixation_time.iloc[1:]=np.array(t.iloc[0:-1])
all_trials.fixation_time[all_trials.fixation_time>2500]=np.nan
all_trials.fixation_time[all_trials.fixation_time<-500]= np.nan
#plt.hist(all_trials.fixation_time[all_trials.systole==1],bins=50)
#plt.hist(all_trials.fixation_time[all_trials.systole==0],bins=50)


all_trials.fixation_time[all_trials.systole==1].mean()
all_trials.fixation_time[all_trials.systole==0].mean()

all_trials.fixation_time[all_trials.fixation_time<0]=0
all_trials.fixation_time[all_trials.fixation_time>1300]=np.nan
all_trials=all_trials[~all_trials.fixation_time.isna()]
#### Filter
filter_inds=np.array( (all_trials['trial_in_block']==14 )) | np.array(all_trials.trial_in_block==14)
filter_inds=np.array( (all_trials['.thisN']<3) ) | filter_inds

nan16=np.zeros(16); nan16[:]=np.nan
inds=np.where(filter_inds==True)[0]
for ind in inds:
    all_trials.peaks_times.iloc[ind][all_trials.peaks_times.iloc[ind]-all_trials.physioTime.iloc[ind]>-100]=np.nan
######

subj_ids=all_trials.id.unique()

### ar2 = RR interpolated signal
ar=all_trials.peaks_times-all_trials.physioTime
ar=np.stack(ar)
new_x=np.arange(-2500, 5000, 100)
ar2=np.empty([ar.shape[0],len(new_x)])
naned_ar=np.empty(ar.shape[0])
for i in range(all_trials.shape[0]):
    subj_id=all_trials.id.iloc[i]
    subj_idx=all_subj[all_subj.id==subj_id].index[0]
    subj_rrs=np.array(all_subj.task_rrs[all_subj.id==subj_id])[0]
    t=ar[i,:]
    rrs=np.diff(t)
    #min_rr=np.min([np.nanquantile(subj_rrs.tolist(), 0.025),545])
    #max_rr=np.max([np.nanquantile(subj_rrs.tolist(), 0.975),1100])
    #q1=all_subj.rr_q1[idx]; q3==all_subj.rr_q3[idx]
    #min_rr=q1-3*(q3-q1); max_rr=q3+3*(q3-q1)
    #m=all_subj.rr_mean[idx]; s=all_subj.rr_sd[idx]
    #min_rr=m-4*s; max_rr=m+4*s
    #min_rr = 400; max_rr=1500
    min_rr = 0; max_rr=1314 #m+-4sd for all data combined
    
    naned_ar[i]=np.sum(np.logical_or(rrs<min_rr, rrs>max_rr))
    rrs[np.logical_or(rrs<min_rr, rrs>max_rr)]=np.nan
    ar2[i,:]=np.interp(new_x, t[1:], rrs)
all_trials['min_rr_fixed']=ar2[:,new_x==-1300] #!!!!!!!!!!!!!!!!!!!!
all_trials['max_rr_fixed']=ar2[:,new_x==1400]
all_trials['post_rr_fixed']=ar2[:,new_x==3600]

all_trials['min_rr_flex']=np.nan
all_trials['min_rr_time']=np.nan

for i in range(all_trials.shape[0]):
    fix_time=(all_trials.fixation_time.iloc[i])
    pre_time=-np.round((fix_time+700)/100)*100
    all_trials.min_rr_time.iloc[i]=pre_time
    all_trials.min_rr_flex.iloc[i]=ar2[i, new_x==pre_time]

### onset to RT COR
all_trials['onset_to_RT_COR']=np.nan
all_trials['RT_RR']=np.nan
for idx in range(all_trials.shape[0]):
     RT=all_trials.RT.iloc[idx]
     if np.isnan(RT):
         continue
     all_trials.RT_RR.iloc[idx] = ar2[idx, new_x==round((RT)/100)*100]



max_or_range=(new_x>200) & (new_x<2600)
onset_time= (new_x<=0) & (new_x>=0)
min_or_range=(new_x>=-1300) & (new_x<=0)
min_time=(new_x>=-1300) & (new_x<=0)

all_trials['max_rr_range']=np.nanmean(ar2[:,max_or_range], axis=1)
all_trials['onset_rr']=np.nanmean(ar2[:,onset_time], axis=1)
all_trials['min_rr_range']=np.nanmean(ar2[:,min_or_range], axis=1)
all_trials.onset_to_RT_COR= all_trials.RT_RR-all_trials.onset_rr    


all_trials['onset_rr']=ar2[:,new_x==0]
all_trials['min_to_onset_OR_fixed']=all_trials.onset_rr-all_trials.min_rr_fixed;
all_trials['onset_to_max_OR_fixed']=all_trials.max_rr_fixed-all_trials.onset_rr;
all_trials['min_to_max_OR_fixed']=all_trials.max_rr_fixed-all_trials.min_rr_fixed;
all_trials['onset_to_max_OR_range']=all_trials.max_rr_range-all_trials.onset_rr;
all_trials['min_to_onset_OR_range']=all_trials.min_rr_range-all_trials.min_rr_fixed;
#all_trials['min_to_max_OR_range']=all_trials.max_rr_range-all_trials.min_rr_fixed;
#ar2=ar2[:,(new_x>=-1500) & (new_x<=4000)]
#new_x= new_x[(new_x>=-1500) & (new_x<=4000)]

ar2_bpm=60000/ar2
all_trials['min_bpm_fixed']=ar2_bpm[:,new_x==-1300]
all_trials['max_bpm_fixed']=ar2_bpm[:,new_x==1400]
all_trials['onset_bpm']=ar2_bpm[:,new_x==0]
all_trials['min_to_onset_OR_fixed_bpm']=all_trials.onset_bpm-all_trials.min_bpm_fixed;
all_trials['onset_to_max_OR_fixed_bpm']=all_trials.max_bpm_fixed-all_trials.onset_bpm;
all_trials['min_to_max_OR_fixed_bpm']=all_trials.max_bpm_fixed-all_trials.min_bpm_fixed;



################
# Acc, bias
all_trials['near_threshold']=np.nan
all_trials['above_threshold']=np.nan
all_trials['below_threshold']=np.nan
all_trials['duration_rel']=np.nan

all_trials['long_bias']=np.nan
all_trials['short_bias']=np.nan
all_trials['TP_bias']=np.nan
all_trials['subj_thresh']=np.nan

thresh_limit=10
all_trials['acc']=np.nan #even near threshold trials are evaluated
all_trials['correct_long_thr']=np.nan #correct response is long based on threshold
for idx,subj_id in enumerate(subj_ids):
    sel_inds=all_trials.id==subj_id
    all_trials.subj_thresh[sel_inds]=all_subj.thresh[idx]
    correct_long=all_trials.duration[sel_inds]>(all_subj.thresh[idx])
    all_trials.correct_long_thr[sel_inds]=correct_long.copy()
    all_trials.acc[sel_inds]= (correct_long==all_trials.long_resp[sel_inds])
    all_trials.TP_bias[sel_inds]= all_subj.thresh[idx] - all_trials.duration[sel_inds]
    all_trials.TP_bias[sel_inds]=np.sign(all_trials.TP_bias[sel_inds])
    all_trials.TP_bias[sel_inds & all_trials.acc==1]=0
    all_trials.short_bias[sel_inds]= (correct_long & ~all_trials.long_resp[sel_inds] )
    all_trials.long_bias[sel_inds]= (~correct_long & all_trials.long_resp[sel_inds])
    all_trials.near_threshold[sel_inds]=(np.abs(all_trials.duration-all_subj.thresh[idx])< thresh_limit)[sel_inds]
    all_trials.below_threshold[sel_inds]=((all_subj.thresh[idx]-all_trials.duration)> thresh_limit)[sel_inds]
    all_trials.above_threshold[sel_inds]=((all_trials.duration-all_subj.thresh[idx])>thresh_limit)[sel_inds]
    all_trials.duration_rel[sel_inds]=all_trials.duration[sel_inds]-all_subj.thresh[idx]
    
all_trials.acc=all_trials.acc.astype(float)
all_trials.long_bias=all_trials.long_bias.astype(float)
all_trials.short_bias=all_trials.short_bias.astype(float)
all_trials.near_threshold=all_trials.near_threshold.astype(float)
all_trials.below_threshold=all_trials.below_threshold.astype(float)
all_trials.above_threshold=all_trials.above_threshold.astype(float)
all_trials.correct_long_thr=all_trials.correct_long_thr.astype(float)
###
# Easiness
all_trials['difficulty']=np.nan
for idx, subj in enumerate(subj_ids):
    sel_inds=all_trials.id==subj
    t=np.abs(all_trials.duration-all_subj.thresh[idx])
    all_trials.difficulty[sel_inds]=t[sel_inds]

all_trials['easiness']=np.abs(all_trials.duration-134)
###
# OR based on subbj
subj_min_time=np.empty(len(subj_ids)); subj_min_time[:]=np.nan
subj_max_time=np.empty(len(subj_ids)); subj_max_time[:]=np.nan
min_inds=np.where((new_x>=-2000) & (new_x<=0))[0]
#<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
max_inds=np.where((new_x>=900) & (new_x<=1900))[0]
#<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
post_inds=np.where((new_x>=3000) & (new_x<=4000))[0]
#<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
min_range=new_x[min_inds]
max_range=new_x[max_inds]
post_range=new_x[post_inds]
for idx,subj in enumerate(subj_ids):
   sel_inds=all_trials.id==subj
   t=np.nanmean(ar2[sel_inds,:], axis=0)  
   if np.isnan(t[min_inds]).mean()>0.5:
       subj_min_time[idx]=np.nan
   else:
       min_ind=np.nanargmin(t[min_inds])
       subj_min_time[idx]=min_range[min_ind]
   if np.isnan(t[max_inds]).mean()>0.5:
        subj_max_time[idx]=np.nan
   else:
        max_ind=np.nanargmax(t[max_inds])
        subj_max_time[idx]=max_range[max_ind]

all_subj['minRR_time']=subj_min_time
all_subj['maxRR_time']=subj_max_time

all_trials['min_rr']=np.nan 
all_trials['min_time']=np.nan 
all_trials['max_rr']=np.nan 
all_trials['max_time']=np.nan 
all_trials['post_rr']=np.nan 
all_trials['post_time']=np.nan 

for idx in range(all_trials.shape[0]):
    #sel_inds=all_trials.id==subj
    #all_trials.min_rr[sel_inds]=ar2[sel_inds, new_x== subj_min_time[idx]]
    #all_trials.max_rr[sel_inds]=ar2[sel_inds, new_x== subj_max_time[idx]]
    t=ar2[idx,:]
    if np.isnan(t[min_inds]).mean()>0.5: #or np.nanargmin(t[min_inds])==0 or np.nanargmin(t[min_inds])==len(min_inds)-1:
        min_time=np.nan; min_rr=np.nan;
    else:
        min_ind=np.nanargmin(t[min_inds])
        min_rr=np.min(t[min_inds])
        min_time=min_range[min_ind]
    all_trials.min_time.iloc[idx]=min_time
    all_trials.min_rr.iloc[idx]=min_rr
    
    if np.isnan(t[max_inds]).mean()>0.5:# or np.nanargmax(t[max_inds])==0 or np.nanargmax(t[max_inds])==(len(max_inds)-1):
        max_time=np.nan; max_rr=np.nan;
    else:
        max_ind=np.nanargmax(t[max_inds])
        max_rr=np.max(t[max_inds])
        max_time=max_range[max_ind]
    all_trials.max_time.iloc[idx]=max_time
    all_trials.max_rr.iloc[idx]=max_rr
    
    
    if np.isnan(t[post_inds]).mean()>0.5:# or np.nanargmax(t[max_inds])==0 or np.nanargmax(t[max_inds])==(len(max_inds)-1):
        post_time=np.nan; post_rr=np.nan;
    else:
        post_ind=np.nanargmin(t[post_inds])
        post_rr=np.min(t[post_inds])
        post_time=post_range[post_ind]
    all_trials.post_time.iloc[idx]=post_time
    all_trials.post_rr.iloc[idx]=post_rr
    
all_trials['onset_to_max_OR_flex']=all_trials.max_rr-all_trials.onset_rr
all_trials['min_to_max_OR_flex']=all_trials.max_rr-all_trials.min_rr_fixed

#types are defined based on middle point (134 ms) rather than BP


correct_resp=np.sign(all_trials.duration-134)
all_trials['stim_type']=correct_resp
all_trials['TP_bias2']=-correct_resp
all_trials.TP_bias2[all_trials.acc==1]=0
all_trials.TP_bias2[all_trials.duration==134]=np.nan#all_trials.long_resp*2-1

#TP_acc2 is defined based on the middle point. The middle point itself has NaN accuracy
all_trials['TP_acc2']=np.copy(all_trials.acc)#np.sign(all_trials.duration-0.134)*(all_trials.long_resp-0.5)*2

all_trials.TP_acc2[(all_trials.near_threshold==1)]=np.nan

all_trials['error_cont']=(1-all_trials.TP_acc2)*np.abs(all_trials.duration-134)
all_trials['min_post_avg_rr']= (all_trials.min_rr_fixed + all_trials.post_rr)/2
all_trials['min_post_avg_rr2']= (all_trials.min_rr_flex + all_trials.post_rr)/2

all_subj['mean_onset_to_max_OR_flex']= all_trials.groupby('id')['onset_to_max_OR_flex'].mean().values
all_subj['mean_min_to_max_OR_flex']= all_trials.groupby('id')['min_to_max_OR_flex'].mean().values
all_subj['median_min_rr_fixed']= all_trials.groupby('id')['min_rr_fixed'].median().values
all_subj['mean_min_rr_fixed']= all_trials.groupby('id')['min_rr_fixed'].mean().values
all_subj['mean_max_time']= all_trials.groupby('id')['max_time'].mean().values
all_subj['median_min_post_avg_rr']= all_trials.groupby('id')['min_post_avg_rr'].median().values

## high vs low within subject
all_trials['high_min_rr']=np.nan
all_trials['high_min_post_avg_rr']=np.nan
for idx,subj_id in enumerate(subj_ids):    
    # high vs low min_rr_fixed 
    sel_inds=all_trials.id==subj_id
    thr=all_subj.median_min_rr_fixed[idx]
    all_trials.high_min_rr[sel_inds]= all_trials.min_rr_fixed[sel_inds]>thr
    
    thr=all_subj.median_min_post_avg_rr[idx]
    all_trials.high_min_post_avg_rr[sel_inds]= all_trials.min_post_avg_rr[sel_inds]>thr

all_trials.high_min_rr=all_trials.high_min_rr.astype(float)
all_trials.high_min_post_avg_rr=all_trials.high_min_post_avg_rr.astype(float)

#### decorrelating min_RR and COR!
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
all_trials['COR_controlledMinRR']=np.nan
all_trials['COR_onset_controlledMinRR']=np.nan #from onset to max
all_trials['RTCOR_controlledMinRR']=np.nan
all_trials['minRR_controlledCOR']=np.nan
all_trials['PCA1']=np.nan
all_trials['PCA2']=np.nan

all_subj['coef_COR_minRR']=np.nan
all_subj['coef_COR_onset_minRR']=np.nan
for idx,subj_id in enumerate(subj_ids):        
    sel_inds=all_trials.id==subj_id
    x=np.array(all_trials.min_rr_fixed[sel_inds])
    y=np.array(all_trials.onset_to_RT_COR[sel_inds]);
    nan_inds=np.isnan(x) | np.isnan(y)
    X=x[~nan_inds]; Y=y[~nan_inds]
    #regressing out
    m, b = np.polyfit(X, Y, 1)
    resid=y-m*x-b
    all_trials.COR_controlledMinRR[sel_inds]=resid
    all_subj.coef_COR_minRR.iloc[idx]=m

    m, b = np.polyfit(Y, X, 1)
    resid=x-m*y-b
    all_trials.minRR_controlledCOR[sel_inds]=resid
    
    #COR‌ from onset
    y=np.array(all_trials.onset_to_max_OR_flex[sel_inds]);
    x=np.array(all_trials.min_rr_fixed[sel_inds])
    nan_inds=np.isnan(x) | np.isnan(y)
    X=x[~nan_inds]; Y=y[~nan_inds]
    m, b = np.polyfit(X, Y, 1)
    resid=y-m*x-b
    all_trials.COR_onset_controlledMinRR[sel_inds]=resid
    all_subj.coef_COR_onset_minRR.iloc[idx]=m

    #PCA
    XX=np.column_stack((X,Y))
    pca = PCA(n_components=2)   
    pca.fit(XX)
    res=pca.transform(XX)
    res1=np.copy(x); res1[~nan_inds]=res[:,0]; res1[nan_inds]=np.nan
    res2=np.copy(x); res2[~nan_inds]=res[:,1]; res2[nan_inds]=np.nan
    all_trials.PCA1[sel_inds]=res1
    all_trials.PCA2[sel_inds]=res2
    
#
### slope
all_trials['COR_slope']=all_trials.min_to_max_OR_flex/all_trials.max_time

all_trials['COR_onset_slope']=all_trials.onset_to_max_OR_flex/all_trials.max_time

all_trials['CORslope_controlledMinRR']=all_trials.COR_controlledMinRR/all_trials.max_time

all_trials['CORonsetSlope_controlledMinRR']=all_trials.COR_onset_controlledMinRR/all_trials.max_time

all_trials['COR_onsetRT_slope']=all_trials.onset_to_RT_COR/all_trials.RT

### signed minRR
all_trials['signed_minRR_thr']=(all_trials.correct_long_thr-0.5)*all_trials.min_rr_fixed

all_trials['signed_minRR']=np.sign(all_trials.duration-134)*all_trials.min_rr_fixed

all_trials['duration_sign']=np.sign(all_trials.duration-134);
all_trials['dur_times_minRR']=(all_trials.duration-134)*all_trials.min_rr_fixed/1000
all_trials['dur_times_minRR_thresh']=(all_trials.duration-all_trials.subj_thresh)*all_trials.min_rr_fixed/1000

min_rr_stan=(all_trials.min_rr_fixed-all_trials.min_rr_fixed.mean())/all_trials.min_rr_fixed.std()

all_trials['signed_minRR_stan']=all_trials.duration_sign*(min_rr_stan-min_rr_stan.min())

all_trials['dur_times_minRR_stan']=(all_trials.duration-134)*(min_rr_stan-min_rr_stan.min())


### import to R
import os
os.environ['R_HOME'] = '/usr/lib/R'

import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
ro.r("rm(list = ls(all.names=TRUE))")
base = importr('base')
lme4 = importr('lme4')
lmerTest=importr('lmerTest')
#bayesR=importr('bayestestR')
#emmeans=importr('emmeans')
rdf=pandas2ri.py2rpy_pandasdataframe(
        all_trials.loc[:,['id','systole', 'duration','long_resp','acc', 
                          'trial_num','block', 'easiness','difficulty','RT_real','RT',
                          'duration_rel','long_bias','short_bias','TP_bias','pre_systole',
                          'min_rr_fixed','max_rr_fixed', 'onset_rr', 'min_to_onset_OR_fixed', 
                          'min_to_max_OR_fixed', 'onset_to_max_OR_fixed',
                          'min_rr','max_rr','min_time','max_time','near_threshold',
                          'below_threshold', 'above_threshold','onset_to_max_OR_flex',
                          'onset_to_max_OR_range','min_rr_range','max_rr_range', 'min_rr_flex',
                          'min_to_onset_OR_range', 'min_bpm_fixed','max_bpm_fixed','onset_to_RT_COR',
                          'onset_bpm','min_to_onset_OR_fixed_bpm','onset_to_max_OR_fixed_bpm',
                          'min_to_max_OR_fixed_bpm','error_cont','fixation_time', 'TP_bias2',
                          'TP_acc2' , 'min_to_max_OR_flex', 'post_rr_fixed', 'min_post_avg_rr','min_post_avg_rr2',
                          'post_rr', 'post_time', 'high_min_rr','high_min_post_avg_rr','COR_controlledMinRR','minRR_controlledCOR',
                          'PCA1','PCA2',
                          'CORslope_controlledMinRR','COR_slope','COR_onset_controlledMinRR','COR_onset_slope' ,'COR_onsetRT_slope','CORonsetSlope_controlledMinRR','signed_minRR_thr', 'signed_minRR','dur_times_minRR','correct_long_thr','duration_sign', 'signed_minRR_stan', 'dur_times_minRR_stan','dur_times_minRR_thresh','stim_type']])

ro.globalenv['rdf'] = rdf

cnames=df_blocks.columns.tolist()
block_rdf=pandas2ri.py2rpy_pandasdataframe(df_blocks.loc[:,cnames])
ro.globalenv['block_rdf'] = block_rdf

#sjplot=importr('sjPlot')

#tidyr=importr('tidyr')
#ro.r("save(rdf,file='/mnt/Aclab/Studies/38_TimeHeart/analysis/all_trials')")

#TP_bias2 is -1, 0, 1 reflecting contractio, dilation and consistent. However, these trial 

all_trials_for_matlab=all_trials.loc[:,['id', 'duration','long_resp','block','min_rr_fixed','RT_real','min_to_max_OR_flex']]

all_trials_for_matlab.to_csv(base_data_dir+'/analysis/all_trials_for_matlab.csv', index=False)

