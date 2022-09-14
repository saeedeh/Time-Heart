#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 09:56:11 2019

@author: saeedeh
"""

def plot_bisection(df, color='b',levels_name='duration', resp_name='long_resp'):
    #levels_name: column name of level
    #resp_name: column name of response
    res=df.groupby(levels_name)[resp_name].mean();
    res_sd=df.groupby(levels_name)[resp_name].std();
    res_count=df.groupby(levels_name)[resp_name].count();
    res_se=res_sd/np.sqrt(res_count)
    #fig, ax = plt.subplots()

    plt.errorbar(x=res.keys(), y=res.values, yerr=res_se.values, fmt=color, alpha=0.7)
    #plt.plot(res.keys(), res.values, color=color,ls='--', marker='o', markersize=3)
    plt.xlabel('t')
    plt.ylabel('p(long)')

def plot_bisection2(df, color='b',levels_name='duration', resp_name='long_resp'):
    #levels_name: column name of level
    #resp_name: column name of response
    #first avg by subj, then avg all
    res=df.groupby([levels_name,'id'])[resp_name].mean();
    res=res.groupby(levels_name).mean()

def plot_bisection3(df, color='b',levels_name='duration', resp_name='long_resp'):
    #plot for each subj
    res=df.groupby(['id',levels_name])[resp_name].mean();
    
    ## plot subj by subj
    for idx, ind in enumerate(subj_ids):
         plt.plot(res[ind].keys(), res[ind].values, color='silver',ls='-', markersize=1)#,fillstyle='none', markersize=4)
    
    
    ## plot mean
    res=df.groupby(levels_name)[resp_name].mean();
    res_sd=df.groupby(levels_name)[resp_name].std();
    res_count=df.groupby(levels_name)[resp_name].count();
    res_se=res_sd/np.sqrt(res_count)
    #fig, ax = plt.subplots()

    res=res.groupby(levels_name).mean()     
    plt.plot(res.keys(), res.values, color="black",ls='-', marker='o', markersize=2)#,fillstyle='none', markersize=4)
    plt.xlabel('time (ms)')
    plt.ylabel('p(long)')
    plt.xticks(res.keys())
    plt.tight_layout()
    plot_2D_ar_errorbar(df,'black')

def estimate_BP(df, levels_name='duration', resp_name='long_resp'):
     res=df.groupby(levels_name)[resp_name].mean();
     
def plot_2D_ar_errorbar(ar,color):
    #each column is one xtick on the x axis, its mean and SEM are plotted as errorbar
    sz=ar.shape
    X=np.zeros(sz[1])
    Y=np.zeros(sz[1])
    sem=np.zeros(sz[1])
    for col in range(sz[1]):
        X[col]=col-1
        Y[col]=np.nanmean(ar[:,col])
        sem[col]=np.nanstd(ar[:,col])/np.sqrt(sz[0])
    plt.errorbar(x=X, y=Y,yerr=sem, fmt=color)
    
def plot_2D_ar_errorbar2(ar,color):
    #ar is a list of arrays. Each array will be one xtick, its mean and SEM are plotted as errorbar
    n=len(ar)
    X=np.zeros(n)
    Y=np.zeros(n)
    sem=np.zeros(n)
    for col in range(n):
        X[col]=col-1
        Y[col]=np.nanmean(ar[col])
        sem[col]=np.nanstd(ar[col])/np.sqrt(len(ar[col]))
    plt.errorbar(x=X, y=Y,yerr=sem, fmt=color)

from scipy import stats
def nan_pearsonr(ar1, ar2):
    nan_inds=np.isnan(ar1) | np.isnan(ar2)
    [r,p]=stats.pearsonr(ar1[~nan_inds], ar2[~nan_inds])
    return r,p

def nan_ttest_rel(ar1, ar2):
    nan_inds=np.isnan(ar1) | np.isnan(ar2)
    n= len(ar1[~nan_inds])
    [t,p]=stats.ttest_rel(ar1[~nan_inds], ar2[~nan_inds])
    return t,p,n-1

def nan_ttest_1samp(ar1, val=0):
    nan_inds=np.isnan(ar1)
    n= len(ar1[~nan_inds])
    [t,p]=stats.ttest_1samp(ar1[~nan_inds],val)
    return t,p,n-1

def nan_ttest_ind(ar1, ar2):
    nan_inds1=np.isnan(ar1)
    nan_inds2=np.isnan(ar2)
    
    [t,p]=stats.ttest_ind(ar1[~nan_inds1],ar2[~nan_inds2])
    return t,p

def subj_mean_mean(ar2, inds):
    if(len(ar2.shape))==2:
        res=np.empty([len(subj_ids), ar2.shape[1]])
        for idx,subj_id in enumerate(subj_ids):
          #  rr0=ar2[(all_trials.id==subj_id) & inds, new_xx==0]
          #  ars=[x-rr0[i] for i, x in enumerate(ar2[(all_trials.id==subj_id) & inds,:])]
            ars=ar2[(all_trials.id==subj_id) & inds]
            res[idx,:]=np.nanmean( ars,axis=0)
            #res[idx,:]=np.nanmean( ar2[(all_trials.id==subj_id) & inds,:],axis=0) #-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        #res=ar2[inds,:]
        return np.nanmean(res, axis=0)
    else:
        res=np.empty([len(subj_ids)])
        for idx,subj_id in enumerate(subj_ids):
            res[idx]=np.nanmean( ar2[(all_trials.id==subj_id) & inds])#-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        #res=ar2[inds,:]
        return np.nanmean(res)

def subj_mean_mean0(ar2, inds):
    if(len(ar2.shape))==2:
        res=np.empty([len(subj_ids), ar2.shape[1]])
        for idx,subj_id in enumerate(subj_ids):
            rr0=ar2[(all_trials.id==subj_id) & inds, new_xx==0]
            ars=[x-rr0[i] for i, x in enumerate(ar2[(all_trials.id==subj_id) & inds,:])]
            res[idx,:]=np.nanmean( ars,axis=0)
            #res[idx,:]=np.nanmean( ar2[(all_trials.id==subj_id) & inds,:],axis=0) #-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        #res=ar2[inds,:]
        return np.nanmean(res, axis=0)
    else:
        res=np.empty([len(subj_ids)])
        for idx,subj_id in enumerate(subj_ids):
            res[idx]=np.nanmean( ar2[(all_trials.id==subj_id) & inds])#-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        #res=ar2[inds,:]
        return np.nanmean(res)

def subj_mean_mean_se(ar, inds):
        res=np.empty([len(subj_ids)])
        for idx,subj_id in enumerate(subj_ids):
            res[idx]=np.nanmean( ar[(all_trials.id==subj_id) & inds])#-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        #res=ar2[inds,:]
        res2=ar[inds]
        return [np.nanmean(res), np.nanstd(res2)/np.sqrt(len(res2))]

def subj_mean_se(ar, inds):
        #res=np.empty([len(subj_ids)])
        #for idx,subj_id in enumerate(subj_ids):
       #     res[idx]=np.nanmean( ar[(all_trials.id==subj_id) & inds])#-np.nanmean( ar2[(all_trials.id==subj_id),:],axis=0)
        res=ar[inds]
        return [np.nanmean(res), np.nanstd(res)/np.sqrt(len(res))]

def balance_subj_inds(inds1, inds2, all_trials):
    subj_in_inds1=all_trials.id.isin(np.unique(all_trials.id[inds1]))
    subj_in_inds2=all_trials.id.isin(np.unique(all_trials.id[inds2]))
    subj_in_inds= subj_in_inds1 & subj_in_inds2
    #subj_in_inds=all_trials.id>0
    inds1=inds1&subj_in_inds; inds2=inds2&subj_in_inds
    return inds1, inds2

def fitting_line(x, y):
    nan_inds=np.isnan(x) | np.isnan(y)
    x=x[~nan_inds]; y=y[~nan_inds]
    range_LB_UB=np.array([np.min(x), np.max(x)])
    m, b = np.polyfit(x, y, 1)
    xx=range_LB_UB
    yy=m*range_LB_UB + b
    plt.plot(xx,yy, alpha=0.7, linestyle="-", color="black")