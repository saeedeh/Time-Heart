#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 07:23:53 2019

@author: saeedeh
"""
from peakdet import Physio, operations, load_rtpeaks, HRV

def get_ecg(fname_or_dataframe, start_time=-1, end_time=-1):
#get ECG results
    #******** Don't forget to check the sample rate!
    if(type(fname_or_dataframe)==str):    
        d=pd.read_csv(fname_or_dataframe)
    else:
        d=fname_or_dataframe
    if len(d)==0:
        return np.nan
    if(start_time==-1 and end_time==-1):
        ind1=0; ind2=len(d)
    else:
        [ind1, ind2]=d.time.searchsorted([start_time, end_time])
    ecg=Physio(d.channel9[ind1:ind2],fs=500.)
    ecg=operations.peakfind_physio(ecg, thresh=0.4, dist=300)
    hrv=HRV(ecg)
    peak_times=d.time[ecg.peaks] #peak timestamps
    trough_times=d.time[ecg.troughs]
    #peak_times=peak_times.reset_index(drop=True)
    #trough_times=trough_times.reset_index(drop=True)
    res={'HRV':hrv, 'peak_times': peak_times, 'trough_times':trough_times}
    return res;

def get_ecg_course_around_timestamps(time_stamps, ecg_fname, len_before, len_after):
    #len_before len_after are in seconds, 
    # **** len_before must be positive
    #get physio of time_stamps given the ecg file name
    d=pd.read_csv(ecg_fname)
    physio_vals=d['channel9']
    physio_times=d['time']
    len_before=len_before*1000; len_after=len_after*1000 #sampling rate is 500 points per second
    n_points= len_after+len_before;
    res=pd.Series()
    for ind,time_stamp in enumerate(time_stamps):
        t_ind=physio_times.searchsorted(time_stamp, side='left')
        if(t_ind-len_before<0 or t_ind+len_after>=len(physio_times)):
            tmp=np.empty(n_points)
            tmp[0:n_points]=np.nan
            res=res.append(pd.Series([tmp]))
            continue
        vals=np.array(physio_vals[t_ind-len_before: t_ind+len_after])
        res=res.append( pd.Series( [vals] ) )
    res=res.reset_index(drop=True);
    res.name='physio_course'
    return res
def get_RR_preCurPost(time_stamps, ecg_fname):
    #get RR of time_stamps given the ecg file name
    ecg = get_ecg(ecg_fname)
    if(pd.isnull(ecg)):
        return pd.Series(np.nan, index=range(len(time_stamps)))
    p_times=ecg['peak_times']
    cur_res=pd.Series()
    pre_res=pd.Series()
    post_res=pd.Series()
    for ind,time_stamp in enumerate(time_stamps):
        p_ind=p_times.searchsorted(time_stamp, side='left')
        if(p_ind-2==0 or p_ind+2>=len(p_times)):
            cur_res=cur_res.append(pd.Series(np.nan))
            pre_res=pre_res.append(pd.Series(np.nan))
            post_res=post_res.append(pd.Series(np.nan))
            continue
        i1_pre=p_times.searchsorted(time_stamp)-2
        i2_pre=i1_pre+1
        i1_cur=i1_pre+1
        i2_cur=i1_cur+1
        i1_post=i1_cur+1;
        i2_post=i1_post+1
        rr_pre=p_times.iloc[i2_pre]-p_times.iloc[i1_pre]
        rr_cur=p_times.iloc[i2_cur]-p_times.iloc[i1_cur]
        rr_post=p_times.iloc[i2_post]-p_times.iloc[i1_post]
#        if( (p2-p1)>2000 or (p2-p1)<300):
#            cur_res=cur_res.append(pd.Series(np.nan))
#            continue
        cur_res=cur_res.append( pd.Series( (rr_cur) ) )
        pre_res=pre_res.append( pd.Series( (rr_pre) ) )
        post_res=post_res.append( pd.Series( (rr_post) ) )
    pre_res=pre_res.reset_index(drop=True);cur_res=cur_res.reset_index(drop=True);post_res=post_res.reset_index(drop=True)    
    res=pd.concat([pre_res,cur_res,post_res], axis=1, ignore_index=True)
    res.columns=['pre_rr', 'cur_rr', 'post_rr']
    return res



def get_RR_preCurPost_n(time_stamps, ecg_fname, n_rrs):
    #get RR of time_stamps given the ecg file name)
    ecg = get_ecg(ecg_fname)
    if(pd.isnull(ecg)):
        return pd.Series(np.nan, index=range(len(time_stamps)))
    p_times=ecg['peak_times']
    res=pd.Series()
    for ind,time_stamp in enumerate(time_stamps):
        p_ind=p_times.searchsorted(time_stamp, side='left')
        if(p_ind-2==0 or p_ind+n_rrs>len(p_times)):
            tmp=np.empty(n_rrs)
            tmp[0:n_rrs]=np.nan
            res=res.append(pd.Series([tmp]))
            continue
        i1=np.zeros(n_rrs)
        i2=np.zeros(n_rrs)
        rrs=np.zeros(n_rrs)
        for i in range(n_rrs):  
            i1[i]=p_times.searchsorted(time_stamp)-2+i
            i2[i]=i1[i]+1
            rrs[i]=p_times.iloc[int(i2[i])]-p_times.iloc[int(i1[i])]
        res=res.append( pd.Series( [rrs] ) )
    res=res.reset_index(drop=True);
    res.name='rr_list'
    return res

def get_peak_Times_n(time_stamps, ecg_fname, n_pre, n_after):
    #get time_stamp of peaks, from previous beat of the time_stamp to n_peaks after it (including current beat)
    # x=(p_ind-2)--------x=(p_ind-1)----[time_stamp]-----x=(p_ind)--------x--------x---------x(p_ind+n_peaks)
    # res is a series of np.arrays
    n_peaks=n_after
    ecg = get_ecg(ecg_fname)
    time_stamps=time_stamps.reset_index(drop=True)
    if(pd.isnull(ecg)):
        tmp=np.empty([len(time_stamps), n_peaks+n_pre])
        tmp[:]=np.nan
        return pd.Series(tmp.tolist(), index=range(len(time_stamps)))
    p_times=ecg['peak_times']
    p_times=p_times.reset_index(drop=True)
    res=pd.Series()
    for ind,time_stamp in enumerate(time_stamps):
        p_ind=p_times.searchsorted(time_stamp, side='left')
        if(p_ind-n_pre<0):
            tmp=np.empty(n_pre-p_ind); tmp[:]=np.nan
            ar=np.array(p_times.iloc[0:p_ind+n_peaks])
            ar=np.insert(ar, 0, tmp)
        if (p_ind+n_peaks>len(p_times)):
            ar=np.array(p_times.iloc[p_ind-n_pre:len(p_times)])
            tmp=np.empty(n_peaks+p_ind-len(p_times))
            tmp[:]=np.nan
            ar=np.insert(ar, len(ar), tmp)   
        else:
            ar=np.array(p_times.iloc[p_ind-n_pre:p_ind+n_peaks])
        res=res.append( pd.Series([ar]))
    res=res.reset_index(drop=True);
    res.name='peak_list'
    return res

