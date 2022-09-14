#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 08:08:20 2019

Initial setup of the data
processes trial data and ECG features and saves them in 2 dataframes:

1)all_trials: each row contains info of one trial for one subject
2)all_subjects: each row contains info of overall data for one subject

Running this file takes a while.
Results are saved as a pickle file and then loaded and further processed at "00-setup-load.py"



@author: saeedeh
"""

import pandas as pd
from peakdet import Physio, operations, load_rtpeaks
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
# Average performance plots!
subj_ids=list(range(200,258))
subj_ids.remove(202) #data not saved
subj_ids.remove(203) #peaks detected wrong
subj_ids.remove(214) #below chance acc

#subj_ids.remove(217) #only 8 blocks of data saved
subj_ids.remove(220) #Saeedeh recording=false!!!!
subj_ids.remove(223) #Saeedeh recording=false!!!!
subj_ids.remove(228) #Saeedeh recording=false!!!!
subj_ids.remove(232) #Saeedeh recording=false!!!!
subj_ids.remove(233) 
subj_ids.remove(237) 
subj_ids.remove(246)
subj_ids.remove(247)
subj_ids.remove(248)
subj_ids.remove(256)


##
ITI=3.5
base_data_dir='/mnt/Aclab/Studies/38_TimeHeart'
data_dir=base_data_dir+'/raw data'
all_trials=[];
all_subj=pd.DataFrame();
for subj_id in subj_ids:
    print(str(subj_id))
    subj_info=pd.DataFrame();
    subj=str(subj_id)
    task_fname=data_dir+'/'+subj+'/taskSubj_'+subj+'_run_1.csv'
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'

    task=pd.read_csv(task_fname);
    sys_inds= (task.systole==1)
    task['long_resp']=task.response=='l';
    trials=task[['systole', 'duration', 'long_resp', 'RT']]
    trials['id']=subj_id
     ## estimate fixation onset
    rt1=task.RT[0:-1]
    onTime1=task.physioTime[0:-1]
    onTime2=task.physioTime[1:]
    fix_physio2=np.array(onTime1)+(np.array(rt1)+ITI)*1000
    fix_dur=np.array(onTime2)-fix_physio2
    invalid_fix=fix_dur>1200
    fix_physio2[invalid_fix]=np.nan
    fix_physio2=np.insert(fix_physio2, 0, np.nan)
    task['fixation_physio']=fix_physio2
    ##
    time_stamps=task.physioTime;
    RR_lists_stimulus=get_RR_preCurPost_n(time_stamps, ecg_fname,7)
    RR_lists_stimulus.name='rr_list_stim'
    trials=pd.concat([trials, RR_lists_stimulus], axis=1)
    time_stamps=task.fixation_physio;
    RR_lists_fixation=get_RR_preCurPost_n(time_stamps, ecg_fname,7)
    RR_lists_fixation.name='rr_list_fix'
    trials=pd.concat([trials, RR_lists_fixation], axis=1)
    
    ##
    time_stamps=task.physioTime;
    physio_course_list=get_ecg_course_around_timestamps(time_stamps, ecg_fname,1, 4)
    trials=pd.concat([trials, physio_course_list], axis=1)
    #
    trials['acc']=np.nan
    trials.acc[trials.duration>trials.duration.median()]= trials.long_resp[trials.duration>trials.duration.median()]
    trials.acc[trials.duration<trials.duration.median()]= ~trials.long_resp[trials.duration<trials.duration.median()]
    trials.acc=pd.to_numeric(trials.acc)
    acc=np.nanmean(list(trials.acc.values))
    all_subj.loc[subj_id,'acc']=acc
    trials['mean_acc']=acc
    
    rrs=np.stack(np.array(trials.rr_list_stim)) #converted to 2d array
    trials['pre_rr']=rrs[:,0]
    trials['post_rr']=rrs[:,2]
    trials['cur_rr']=rrs[:,1]
    trials['rr_change1']=rrs[:,1]-rrs[:,0]
    trials['rr_change2']=rrs[:,2]-rrs[:,0]
    trials['orienting1']=trials.rr_change1>trials.rr_change1.median()
    trials['orienting2']=trials.rr_change2>trials.rr_change2.median()
    mean_rr_change1=np.nanmean(trials['rr_change1'])
    mean_rr_change2=np.nanmean(trials['rr_change2'])
    trials['mean_rr_change1']=mean_rr_change1
    trials['mean_rr_change2']=mean_rr_change2
    all_subj.loc[subj_id,'mean_rr_change1']=mean_rr_change1
    all_subj.loc[subj_id,'mean_rr_change2']=mean_rr_change2
    
    trials['trial_number']=list(range(len(trials)))
    if len(all_trials)==0:
        all_trials=trials
    else:
        all_trials= all_trials.append(trials)
#remove long RT's    
all_trials.RT[all_trials.RT>3]=np.nan
all_trials['RT_real']=all_trials.RT-all_trials.duration
all_trials.RT_real[all_trials.RT_real<0.2]=np.nan


all_trials['easiness']=np.abs(all_trials.duration-all_trials.duration.median())


#filtering nans
nan_inds=np.isnan(all_trials.pre_rr) | np.isnan(all_trials.cur_rr)| np.isnan(all_trials.post_rr) 
all_trials_filtered= all_trials[~nan_inds]

for subj_id in subj_ids:
    trials=all_trials[all_trials.id==subj_id]
    trials_filtered=all_trials_filtered[all_trials_filtered.id==subj_id]
    # high/low preRR 
    mean_rr=np.nanmean(trials.pre_rr)
    all_trials.loc[all_trials.id==subj_id, 'high_pre_rr']=trials.pre_rr>mean_rr
    
    mean_rt=np.nanmean(trials.RT)
    all_subj.loc[subj_id,'mean_RT']=mean_rt
    #orienting strength
    [t,p]=stats.ttest_rel(trials_filtered.pre_rr, trials_filtered.post_rr)
    all_subj.loc[subj_id,'orienting_tVal']=t
    all_trials.loc[all_trials.id==subj_id, 'orienting_tVal']=t
    #mean RT
    trials=all_trials[all_trials.id==subj_id]
    mean_rt=np.nanmean(trials.RT_real)
    all_trials.loc[all_trials.id==subj_id, 'mean_RT_real']=mean_rt
    all_subj.loc[subj_id,'mean_RT']=mean_rt
    #mean number of trials RR increased
    inds=trials.rr_change2>0
    inds[trials.rr_change2.isna()]=np.nan
    mean_count_rr_inc=inds.mean()
    all_trials.loc[all_trials.id==subj_id, 'mean_count_orienting']=mean_count_rr_inc
    all_subj.loc[subj_id,'mean_count_orienting']=mean_count_rr_inc

#filtering nans
nan_inds=np.isnan(all_trials.pre_rr) | np.isnan(all_trials.cur_rr)| np.isnan(all_trials.post_rr) 
all_trials_filtered= all_trials[~nan_inds]

all_task_records=pd.DataFrame()
for subj_id in np.unique(all_trials.id):
    subj=str(subj_id)
    task_fname=data_dir+'/'+subj+'/taskSubj_'+subj+'_run_1.csv'
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
    task=pd.read_csv(task_fname);
    all_task_records=all_task_records.append(task)
###################################
# Estimate Peaks-times
n_pre=6; n_after=10
all_trials=all_trials.reset_index(drop=True)

all_trials['peaks_times']=np.nan
for subj_id in np.unique(all_trials.id):
    subj=str(subj_id)
    subj
    sel_inds=all_trials.id==subj_id
    time_stamps=all_trials.physioTime.loc[sel_inds]
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
    res=get_peak_Times_n(time_stamps, ecg_fname, n_pre, n_after)
    res=res.reset_index(drop=True)
    all_trials.peaks_times[all_trials[sel_inds].index]=res
tmp=np.empty(n_pre+n_after); tmp[:]=np.nan
len(tmp)
len(all_trials.peaks_times[3000])
ar=np.array(all_trials.peaks_times)
t=[type(x) for x in ar]
float_inds=[i for i,x in enumerate(t) if x==float]
list_inds=[i for i,x in enumerate(t) if x==list]
all_trials['onset_beatInd']=n_pre

# mean and std of HR 
all_subj['id']=np.unique(all_trials.id)
all_subj['mean_task_rr']=np.nan
all_subj['sd_task_rr']=np.nan
all_subj['task_rrs']=np.nan
for idx,subj_id in enumerate(np.unique(all_trials.id)):
    subj=str(subj_id)
    subj
    sel_inds=all_trials.id==subj_id
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
    t=get_ecg(ecg_fname)
    t=np.diff(t['peak_times'])
    t[np.array(t<200) | np.array(t>2000)]=np.nan
    all_subj.task_rrs.iloc[idx]=t.astype(object)
    all_subj.mean_task_rr.iloc[idx]=t.mean()
    all_subj.sd_task_rr.iloc[idx]=t.std()
#pre-systole
all_trials['pre_systole']=np.nan #was the pre trial sys?
for idx,subj_id in enumerate(subj_ids):
   trials_sys=all_trials.systole.loc[inds]
   all_trials.pre_systole[inds[1:]]=trials_sys[:-1]

df_blocks=pd.DataFrame()
## block HR
for idx,subj_id in enumerate(np.unique(all_trials.id)):
    sel_inds=all_trials.id==subj_id
    subj=str(subj_id)
    print(subj)
    task_fname=data_dir+'/'+subj+'/taskSubj_'+subj+'_run_1.csv'
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
    task=pd.read_csv(task_fname);
    ecg_df=pd.read_csv(ecg_fname)
    bl_num= np.round(task.shape[0]/14).astype(int)
    blocks=pd.DataFrame(range(bl_num), columns=['block_num'])
    blocks['id']=subj_id
    mean_hr=np.zeros(bl_num)
    mean_hr[0:]=np.nan
    for bl in range(bl_num):
        start_tr=bl*14
        end_tr=(bl+1)*14-1
        start_physio=task.physioTime[start_tr]-4000
        end_physio=task.physioTime[end_tr]+1000
        t=get_ecg(ecg_df, start_time=start_physio, end_time=end_physio)
        tdiff=np.diff(t['peak_times'])
        tdiff[np.array(tdiff<200) | np.array(tdiff>2000)]=np.nan
        mean_hr[bl]=np.nanmean(tdiff)
    blocks['mean_rr']=mean_hr
    if len(df_blocks)==0:
        df_blocks=blocks
    else:
        df_blocks= df_blocks.append(blocks)
   

# extract outlider threshold for RR
all_subj['rr_q1']=np.nan
all_subj['rr_q3']=np.nan
all_subj['rr_mean']=np.nan
all_subj['rr_sd']=np.nan

for idx,subj_id in enumerate(np.unique(all_trials.id)):
    subj=str(subj_id)
    subj
    sel_inds=all_trials.id==subj_id
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
    t=get_ecg(ecg_fname)
    t=np.diff(t['peak_times'])
    q1=np.nanquantile(t, 0.25)
    q3=np.nanquantile(t, 0.75)
    iqr=q3-q1
    all_subj.rr_q1[idx]=q1
    all_subj.rr_q3[idx]=q3
    all_subj.rr_sd[idx]= np.nanstd(t)
    all_subj.rr_mean[idx]= np.nanmean(t)
    
        

tt=all_trials.groupby(["id", "block"])
df_blocks['TP_acc2']=tt["TP_acc2"].mean().values
df_blocks['TP_bias2']=tt["TP_bias2"].mean().values
df_blocks['mean_COR']=tt["min_to_max_OR_flex"].mean().values
df_blocks['mean_RT']=tt["RT_real"].mean().values


import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
base = importr('base')
lme4 = importr('lme4')
lme4t = importr('lmerTest')
rdf=pandas2ri.py2rpy_pandasdataframe(
        all_trials.loc[:,['id','systole', 'duration','long_resp','pre_rr','high_pre_rr', 'cur_rr', 'post_rr','acc', 'rr_change1', 'rr_change2', 'orienting1', 'orienting2', 'easiness', 'RT', 'RT_real','high_pre_rr']])

ro.globalenv['rdf'] = rdf


all_thresholds=np.array([0.1385,	0.1494,	0.1386,	0.1298,	0.1225,	0.1315,	0.1169,	0.1226,	0.1434,	0.1312,	0.1460,	0.1258,	0.1267,	0.1561,	0.1250,	0.1365,	0.1322,	0.1466,	0.1387,	0.1160,	0.1269,	0.1384,	0.1275,	0.1250,	0.1364,	0.1304,	0.1139,	0.1407,	0.1362,	0.1400,	0.1234,	0.1309,	0.1152,	0.1226,	0.1479,	0.1358,	0.1220,	0.1374,	0.1368,	0.1472,	0.1309,	0.1149,	0.1375,	0.1517,	0.1487])

all_slopes=np.array([60.2,	93.0,	85.4,	128.0,	105.9,	133.2,	98.3,	175.8,	91.6,	84.5,	66.3,	96.1,	119.4,	42.3,	84.6,	124.0,	117.8,	83.4,	102.4,	182.8,	127.0,	110.1,	135.1,	112.5,	86.4,	101.4,	106.9,	130.6,	99.9,	103.6,	108.9,	121.2,	95.4,	149.2,	108.4,	131.5,	116.6,	74.8,	93.6,	122.5,	115.1,	126.5,	100.8,	66.1,	77.9])

all_subj['slope']=all_slopes
all_subj['thresh']=all_thresholds

all_trials=all_trials.reset_index(drop=True)
all_task_records=all_task_records.reset_index(drop=True)
all_subj=all_subj.reset_index(drop=True)
all_trials=pd.concat([all_trials, all_task_records], axis=1)
all_trials = all_trials.loc[:,~all_trials.columns.duplicated()]

fname=base_data_dir+'/analysis/all_subjects.pkl';
all_subj.to_pickle(fname)

fname=base_data_dir+'/analysis/all_trials.pkl';
all_trials.to_pickle(fname)

fname=base_data_dir+'/analysis/df_blocks.pkl';
df_blocks.to_pickle(fname)
