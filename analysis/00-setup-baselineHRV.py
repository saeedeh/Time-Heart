#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 16:17:48 2020

@author: saeedeh
"""
import glob
#############################################
# inspection
############################################
# check data
for subj_id in subj_ids:
    print(str(subj_id))
    subj=str(subj_id)
    bl_fname=data_dir+'/'+subj+'/'+subj+'_main_baseline*'
    fname=glob.glob(bl_fname)
    if(len(fname)==0):
        print('----no data:'+str(subj))
        continue
    if(len(fname)>1):
        print('***more than one file: '+str(subj))
        
bl_no_data=[209, 213,239,242,251, 257]

#RHRV preperation
rhrv = importr('RHRV')
# call an R function on a Pandas DataFrame
ro.r(('get_hrv_time_analysis<- function(time_beats){'+
      'hrv.data<-CreateHRVData();'+
      'hrv.data <- SetVerbose(hrv.data, FALSE);'+
      'hrv.data$Beat$Time<-time_beats;'+
      'hrv.data <- BuildNIHR(hrv.data);'+
      'hrv.data <- FilterNIHR(hrv.data);'+
      'hrv.data <- InterpolateNIHR(hrv.data);'+
      'hrv.data <- CreateTimeAnalysis(hrv.data);'+
      'res<- hrv.data$TimeAnalysis[[1]];'+
      'res$beats_num<-dim(hrv.data$Beat)[1];'+
      'return(res);'+
      '}'))
rhrv_time_analysis=ro.r('get_hrv_time_analysis')
rhrv_time_analysis=ro.r('get_hrv_time_analysis')

# read data and get ecg
TOO_SMALL_RR_THRESHOLD=0.45
TOO_LARGE_RR_THRESHOLD=1.3
for idx,subj_id in enumerate(subj_ids[29:]):
    print(str(subj_id))
    subj=str(subj_id)
    bl_fname=data_dir+'/'+subj+'/'+subj+'_main_baseline*'
    fname=glob.glob(bl_fname)
    if(len(fname)==0):
        continue
    ecg = get_ecg(fname[0])  
    t=ecg['HRV']
    beats=t.data.peaks/1000
    beats_list= ro.FloatVector(beats.tolist())
    res=rhrv_time_analysis(beats_list)
    res=dict(zip(res.names, list(res)))
    if(np.mean(t.rrint<TOO_SMALL_RR_THRESHOLD)>0.2 or np.mean(t.rrint>TOO_LARGE_RR_THRESHOLD)>0.2):
        print('weired subject: '+ str(subj_id))
        break
        #insepcetion
        plt.hist(t.rrint, bins=50)
        file_name=fname[0]
        ecg = load_rtpeaks(file_name, fs=500., channel=9)
        ecg=operations.peakfind_physio(ecg, thresh=0.4, dist=300)
        operations.plot_physio(ecg)
      #     res=[np.nan, np.nan, np.nan]
      # else:
      #     res=[res['SDNN'][0], res['pNN50'][0], res['rMSSD'][0]]

############################################################
###### MAIN setup
rhrv = importr('RHRV')
# call an R function on a Pandas DataFrame
ro.r(('get_hrv_time_analysis<- function(time_beats){'+
      'hrv.data<-CreateHRVData();'+
      'hrv.data <- SetVerbose(hrv.data, FALSE);'+
      'hrv.data$Beat$Time<-time_beats;'+
      'hrv.data <- BuildNIHR(hrv.data);'+
      'hrv.data <- FilterNIHR(hrv.data);'+
      'hrv.data <- InterpolateNIHR(hrv.data);'+
      'hrv.data <- CreateTimeAnalysis(hrv.data);'+
      'res<- hrv.data$TimeAnalysis[[1]];'+
      'res$beats_num<-dim(hrv.data$Beat)[1];'+
      'return(res);'+
      '}'))
rhrv_time_analysis=ro.r('get_hrv_time_analysis')

# read data and get ecg
measure_names=['SDNN','pNN50','rMSSD']
all_subj['mean_bl_rr']=np.nan
for m in measure_names: all_subj[m]=np.nan
for idx,subj_id in enumerate(subj_ids):
    print(str(subj_id))
    subj=str(subj_id)
    bl_fname=data_dir+'/'+subj+'/'+subj+'_main_baseline*'
    fname=glob.glob(bl_fname)
    if(len(fname)==0):
        continue
    ecg = get_ecg(fname[0])  
    t=ecg['HRV']
    all_subj.mean_bl_rr.iloc[idx]=t.avgnn
    beats=t.data.peaks/1000
    beats_list= ro.FloatVector(beats.tolist())
    res=rhrv_time_analysis(beats_list)
    res2=rhrv_
    res=dict(zip(res.names, list(res)))
    for m in measure_names:
        all_subj[m].iloc[idx]=res[m][0]


############################################################
###### Frequency  analysis
rhrv = importr('RHRV')
rhrv = importr('pracma')

# call an R function on a Pandas DataFrame
ro.r(('get_hrv_freq_analysis<- function(time_beats){'+
      'hrv.data<-CreateHRVData();'+
      'hrv.data <- SetVerbose(hrv.data, FALSE);'+
      'hrv.data$Beat$Time<-time_beats;'+
      'hrv.data <- BuildNIHR(hrv.data);'+
      'hrv.data <- FilterNIHR(hrv.data);'+
      'hrv.data <- InterpolateNIHR(hrv.data);'+
      'hrv.data <- CreateFreqAnalysis(hrv.data);'+
      'hrv.data <- CalculatePSD(hrv.data, indexFreqAnalysis = 1, method = "lomb", doPlot = F);'+
      't<-hrv.data$FreqAnalysis[[1]]$periodogram;'+
      'hf_inds<-t$freq>0.15 & t$freq<0.4;'+
      'hf_mean<- trapz(t$freq[hf_inds] ,t$spec[hf_inds] );'+
      'return(hf_mean);'+
      '}'))
rhrv_freq_analysis=ro.r('get_hrv_freq_analysis')

# read data and get ecg
all_subj['HF']=np.nan
for idx,subj_id in enumerate(subj_ids):
    print(str(subj_id))
    subj=str(subj_id)
    bl_fname=data_dir+'/'+subj+'/'+subj+'_main_baseline*'
    fname=glob.glob(bl_fname)
    if(len(fname)==0):
        continue
    ecg = get_ecg(fname[0])  
    t=ecg['HRV']
    all_subj.mean_bl_rr.iloc[idx]=t.avgnn
    beats=t.data.peaks/1000
    beats_list= ro.FloatVector(beats.tolist())
    hf=rhrv_freq_analysis(beats_list)[0]
    all_subj['HF'].iloc[idx]=hf
nan_pearsonr(all_subj.HF, all_subj.slope)
nan_pearsonr(all_subj.rMSSD, all_subj.slope)