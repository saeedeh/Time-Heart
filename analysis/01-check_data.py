#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 09:14:07 2019

@author: saeedeh
"""

import pandas as pd
from peakdet import Physio, operations, load_rtpeaks
import matplotlib.pyplot as plt
import numpy as np
data_dir='/mnt/Aclab/Studies/38_TimeHeart/raw data'

subj='244'
task_fname=data_dir+'/'+subj+'/taskSubj_'+subj+'_run_1.csv'
physio_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
peaks_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_peaks.csv'

task=pd.read_csv(task_fname);
bio_d=pd.read_csv(physio_fname)
peaks=pd.read_csv(peaks_fname)

sys_inds= (task.systole==1)
sys_physio=task.physioTime[sys_inds]
dias_physio=task.physioTime[np.logical_not(sys_inds)]
task['long_resp']=task.response=='l';
#if task.key_order_reversed[0]==True:
#    task['long_resp']=task.response=='s';
#else:
#    task['long_resp']=task.response=='l';
## Check peak finding

physio_inds1= np.isin(bio_d.time,sys_physio)
physio_inds2= np.isin(bio_d.time,dias_physio)
plt.plot(bio_d.time, bio_d.channel9)
plt.plot(bio_d.time[physio_inds1], bio_d.channel9[physio_inds1],'or')
plt.plot(bio_d.time[physio_inds2], bio_d.channel9[physio_inds2],'og')

rtp_peaks_x= peaks.time[peaks.peak==1]
rtp_peaks_y= peaks.amplitude[peaks.peak==1]
plt.plot(rtp_peaks_x, rtp_peaks_y,'*')

## Check bisection systole/diastole

res_sys=task[sys_inds].groupby('duration')['long_resp'].mean();
res_dias=task[~sys_inds].groupby('duration')['long_resp'].mean();
fig, ax = plt.subplots()
ax.plot(res_dias,'-b',label='diastole')
ax.plot(res_sys,'-r',label='systole')
legend=ax.legend()
plt.xlabel('t')
plt.ylabel('p(long)')

# bisection without physio
res=task.groupby('duration')['long_resp'].mean();
fig, ax = plt.subplots()
ax.plot(res,'b')
plt.xlabel('t')
plt.ylabel('p(long)')

# Check operations plot_physio
subj='236'

physio_fname=data_dir+'/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
bio_d=pd.read_csv(physio_fname)

sel_inds=all_trials.id==int(subj)
all_trials.min_rr[sel_inds].isna().mean()
ecg=Physio(bio_d.channel9,fs=50.)
ecg=operations.peakfind_physio(ecg, thresh=0.4, dist=300)
operations.plot_physio(ecg)
plt.hist(np.diff(ecg.peaks))

plt.hist(np.diff(peaks.time[peaks.peak==1]))

all_trials.pre_rr[sel_inds]