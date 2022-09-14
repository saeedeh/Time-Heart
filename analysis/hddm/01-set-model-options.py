#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 14:04:07 2021

@author: saeedeh
"""

import hddm
from pathlib import Path
import matplotlib.pyplot as plt


cwd = str(Path.cwd())


data_resp = hddm.load_csv(cwd+'/04-data/trials-resp.csv')

data = hddm.load_csv(cwd+'/0-data/all_trials_for_hddm.csv')



N=2000; nburn=1000



#model_no=4.12
#next: 2.3, 4.3 again with rr_thr
have_bias=True
################
# Response modeling
################

elif model_no==1.3: 
    #### 1) base 
    mdl =  hddm.HDDMRegressor(data_resp, ['v ~ stim', 't ~ stim'],bias=have_bias,
                               p_outlier=0.02)
    mdl_name='Base1.3'


elif model_no==4.9:
    #### 3)
    mdl = hddm.HDDMRegressor(data_resp, ['v ~ stim  + dur_times_minRR_thresh+ CORonsetSlope_controlledMinRR',
                                         't ~ stim + min_rr_fixed',
                                         'z ~ min_rr_fixed'], bias=True,p_outlier=0.02)
    mdl_name='heart4.9-both'

    
####################
mname='N'+str(N)+'_'+mdl_name;
if have_bias==False:
    mname="noBias_"+mname
trace_name=cwd+'/04-data/TRACE-'+mname+'.db'
out_name=cwd+'/04-data/'+mname

