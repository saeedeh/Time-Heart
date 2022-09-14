#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 09:36:54 2021

@author: saeedeh

divide duration by the minimum duration
stim is duration relative to thr
"""

def normalize_col(df, col):
    df.loc[:,col]=(df.loc[:,col])/np.nanstd(df.loc[:,col])

from pathlib import Path
cwd = str(Path.cwd())

val_inds=~(all_trials.RT.isna()| all_trials.min_rr_flex.isna() | all_trials.min_to_max_OR_flex.isna() |all_trials.min_rr_fixed.isna())
trials_sel = all_trials.loc[val_inds,['id','RT','long_resp', 'stim_type',
                            'duration','duration_rel','min_rr_fixed', 'min_rr_flex','max_time','min_to_max_OR_flex',
                            'COR_controlledMinRR','CORslope_controlledMinRR','correct_long_thr','easiness', 'COR_slope','COR_onset_slope','COR_onsetRT_slope', 'CORonsetSlope_controlledMinRR', 'onset_to_max_OR_flex','dur_times_minRR','dur_times_minRR_thresh', 'signed_minRR','duration_sign','signed_minRR_thr','onset_to_RT_COR','subj_thresh']]

trials_sel['stim']=trials_sel.duration_rel/88 #duration 
normalize_col(trials_sel,'min_rr_flex')
normalize_col(trials_sel,'min_rr_fixed')
normalize_col(trials_sel,'onset_to_max_OR_flex')
normalize_col(trials_sel,'min_to_max_OR_flex')
normalize_col(trials_sel,'dur_times_minRR')
normalize_col(trials_sel,'dur_times_minRR_thresh')
normalize_col(trials_sel,'onset_to_max_OR_flex')
normalize_col(trials_sel,'dur_times_minRR')
normalize_col(trials_sel,'dur_times_minRR_thresh')
normalize_col(trials_sel,'signed_minRR')
normalize_col(trials_sel,'signed_minRR_thr')


import seaborn as sns

if trials_sel.RT.mean()>10:
    trials_sel.RT=trials_sel.RT/1000

df1=trials_sel.rename({'id':'subj_idx','RT':'rt', 'long_resp':'response'}, axis=1) 
df1.to_csv(cwd+'/04-data/trials-resp_addedTHR.csv', index=False) 
# I added subj_thr to the fields after fitting the model.
# Just for being cautions I am saving the result in a new fild (with THR added to the end) instead of overwriting the old one


     
          
        