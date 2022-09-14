#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 12:32:11 2019

@author: saeedeh
"""
import pandas as pd
from peakdet import Physio, operations, load_rtpeaks
import matplotlib.pyplot as plt
import numpy as np
subj='007'
task_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/007/taskSubj_007_run_1_2019_Sep_15_1242.csv'
physio_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/007/007_log-run1_MP150_data.csv'
peaks_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/007/007_log-run1_MP150_peaks.csv'

subj='001'
task_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/001/taskSubj_001_run_1_2019_Sep_11_2220.csv'
physio_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/001/001_log-run1_MP150_data.csv'
peaks_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/001/001_log-run1_MP150_peaks.csv'

subj='100'
task_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/'+subj+'/taskSubj_'+subj+'_run_1.csv'
physio_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/'+subj+'/'+subj+'_log-run1_MP150_data.csv'
peaks_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/'+subj+'/'+subj+'_log-run1_MP150_peaks.csv'

task=pd.read_csv(task_fname);
bio_d=pd.read_csv(physio_fname)
peaks=pd.read_csv(peaks_fname)
#ecg = load_rtpeaks(physio_fname, fs=1000., channel=9)
#ax1=operations.plot_physio(ecg)
sys_inds= (task.systole==1)
sys_physio=task.physioTime[sys_inds]
dias_physio=task.physioTime[np.logical_not(sys_inds)]
physio_inds1= np.isin(bio_d.time,sys_physio)
physio_inds2= np.isin(bio_d.time,dias_physio)
plt.plot(bio_d.time, bio_d.channel9)
plt.plot(bio_d.time[physio_inds1], bio_d.channel9[physio_inds1],'og')
plt.plot(bio_d.time[physio_inds2], bio_d.channel9[physio_inds2],'or')

rtp_peaks_x= peaks.time[peaks.peak==1]
rtp_peaks_y= peaks.amplitude[peaks.peak==1]
plt.plot(rtp_peaks_x, rtp_peaks_y,'*')


# plot bisection
task['long_resp']=task.response=='l';
res_sys=task[sys_inds].groupby('duration')['long_resp'].mean();
res_dias=task[~sys_inds].groupby('duration')['long_resp'].mean();
fig, ax = plt.subplots()
ax.plot(res_dias,'*b',label='diastole')
ax.plot(res_sys,'*r',label='systole')
legend=ax.legend()
plt.xlabel('t')
plt.ylabel('p(long)')


#bisection without physio
subj='100'
task_fname='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/data/0_210Trial/taskSubj_000_run_2.csv'
task=pd.read_csv(task_fname);
task['long_resp']=task.response=='l';
res=task.groupby('duration')['long_resp'].mean();
fig, ax = plt.subplots()
ax.plot(res,'b')
plt.xlabel('t')
plt.ylabel('p(long)')