#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:56:20 2020

@author: saeedeh
"""
import math
all_subj['interocept_err']=np.nan
black_list=[213,226]
base_data_dir='/mnt/Aclab/Studies/38_TimeHeart'
data_dir=base_data_dir+'/raw data'
invalid_count=0;
for idx,subj_id in enumerate(subj_ids):
    
    if (np.isin(subj_id, black_list)):
        continue
    print(str(subj_id))
    subj_info=pd.DataFrame();
    subj=str(subj_id)
    task_fname=data_dir+'/'+subj+'/'+subj+'_Heartbeat.csv'
    ecg_fname=data_dir+'/'+subj+'/'+subj+'_heartBeatRecording-run1_MP150_data.csv'
    task=pd.read_csv(task_fname);
    
    ecg = get_ecg(ecg_fname, task.Count25_start_physio[0], task.Count25_stop_physio[0])
    beats25=ecg['peak_times'].shape[0]
    ecg = get_ecg(ecg_fname, task.Count35_start_physio[1], task.Count35_stop_physio[1])
    beats35=ecg['peak_times'].shape[0]
    ecg = get_ecg(ecg_fname, task.Count45_start_physio[2], task.Count45_stop_physio[2])
    beats45=ecg['peak_times'].shape[0]
    lens=np.array([25,35,45])
    beats=np.array([beats25,beats35,beats45])/lens*60
    invalid_beat= ((beats < 45) | (beats > 120) )
    beats25= np.nan if invalid_beat[0] else beats25
    beats35= np.nan if invalid_beat[0] else beats35
    beats45= np.nan if invalid_beat[0] else beats45
    invalid_count=invalid_count+ np.isnan(beats25)+np.isnan(beats35)+np.isnan(beats45)
    resp25=task['SliderBar_25.response'][0]
    resp35=task['SliderBar_35.response'][1]
    resp45=task['SliderBar_45.response'][2]
    
    err25= np.abs(beats25-resp25)/resp25
    err35= np.abs(beats35-resp35)/resp35
    err45= np.abs(beats45-resp45)/resp45
    all_subj.interocept_err.iloc[idx] = np.abs(np.nanmean([err25,err35,err45]))

all_subj['interocept_acc']=pd.Series([-math.log(x) for x in np.array(all_subj.interocept_err)])  
nan_pearsonr(all_subj.interocept_acc, all_subj.mean_onset_rr-all_subj.mean_task_rr) # Higher interoception is correlated with larger orienting reflex!

y=all_subj.interocept_err
x=all_subj.mean_onset_to_max_OR
idx = np.isfinite(x) & np.isfinite(y)

plt.plot(x,y , '.', markersize=12, markeredgecolor='k')
#plt.ylim([-0.2,1.9])
#plt.xlim([-8,38])

m, b = np.polyfit(x[idx], y[idx], 1)
plt.plot(x, m*x + b)



### Estimate COR importance for each subject
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import Imputer
all_subj['COR_importance']=np.nan    
all_subj['min_rr_importance']=np.nan    
all_subj['mean_COR_dilation']=np.nan    
all_subj['mean_COR_contraction']=np.nan    
all_subj['mean_COR_shortResp']=np.nan    
all_subj['mean_COR_longResp']=np.nan    
model1 = LinearRegression()
model2 = LogisticRegression()
imputer = Imputer()
for subj in subj_ids:
   print(subj)
   inds = (all_trials.id==subj)
   
   ## COR importance
   y=all_trials.TP_bias2[inds];
   x1=all_trials.onset_to_max_OR_flex[inds]
   x2=all_trials.duration[inds]
   X=pd.DataFrame({'COR':x1})
   
   #y_imputed = imputer.fit_transform(y)
   X_imputed = imputer.fit_transform(X)
   res = model1.fit(X_imputed,y)
   ind = all_subj.loc[all_subj.id == subj].index[0]
   all_subj.COR_importance.loc[ind]=res.coef_[0]
   
   #### min_rr_importance
   y=all_trials.TP_acc2[inds];
   x1=all_trials.min_rr_fixed[inds]
   X=pd.DataFrame({'min_rr':x1})
   keep_inds= ~np.isnan(y)
   X_imputed = imputer.fit_transform(X)
   res = model1.fit(X_imputed[keep_inds],y[keep_inds])
   ind = all_subj.loc[all_subj.id == subj].index[0]
   all_subj.min_rr_importance.loc[ind]=res.coef_[0]
   
      ###
   dilate_inds= (all_trials.acc==0) & (all_trials.long_resp==1) & inds
   contract_inds= (all_trials.acc==0) & (all_trials.long_resp==0) & inds
   
   long_inds= (all_trials.long_resp==1) & inds
   short_inds= (all_trials.long_resp==0) & inds
   
   all_subj.mean_COR_dilation.loc[ind]= np.nanmean(all_trials.onset_to_max_OR_flex[dilate_inds])
   all_subj.mean_COR_contraction.loc[ind]= np.nanmean(all_trials.onset_to_max_OR_flex[contract_inds])

   all_subj.mean_COR_longResp.loc[ind]= np.nanmean(all_trials.onset_to_max_OR_flex[long_inds])
   all_subj.mean_COR_shortResp.loc[ind]= np.nanmean(all_trials.onset_to_max_OR_flex[short_inds])


nan_pearsonr(all_subj.thresh, all_subj.mean_min_to_max_OR_flex)
nan_pearsonr(all_subj.interocept_acc, all_subj.slope)
nan_pearsonr(all_subj.interocept_acc, all_subj.mean_min_to_max_OR_flex) 
nan_pearsonr(all_subj.interocept_acc, all_subj.mean_min_rr_fixed) 
nan_pearsonr(all_subj.interocept_acc, all_subj.mean_max_time) 
nan_pearsonr(all_subj.thresh, all_subj.mean_RT)
 

plt.plot(all_subj.interocept_acc, all_subj.mean_onset_to_max_OR_flex, 'ok')
plt.xlabel('Heartbeat-Counting Performance')
plt.ylabel('Average MCOR in Temporal Task')
fitting_line(all_subj.interocept_acc, all_subj.mean_onset_to_max_OR_flex)
plt.tight_layout()
