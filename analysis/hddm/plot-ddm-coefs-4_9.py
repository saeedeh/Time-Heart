#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 11:26:18 2022

@author: saeedeh
"""
from pathlib import Path
import hddm
import numpy as np
import matplotlib.pyplot as plt
import os
cwd = str(Path.cwd())

fig_size=(3,2)

data_resp = hddm.load_csv(cwd+'/04-data/trials-resp.csv')
mname='N2000_heart4.9-both'

fname=cwd+'/04-data/'+mname
mdl=hddm.load(fname)
data_all=data_resp.copy()

import seaborn as sns

## DIC
mdl.dic
## plot coefficients
t=mdl.get_group_nodes()

save_dir=cwd+'/images/4.9'

###########################################################################
###########################drift rate
###########################################################################
################################################
##linear
################################################
i=0
for i in range(2):
    fig, ax = plt.subplots(figsize=fig_size)
    color='black'
    shade_color="indianred"
    data_resp=data_all
    
    min_rr = data_resp.min_rr_fixed.mean()
    stim = data_resp.stim.mean()
    corSlope = data_resp.CORslope_controlledMinRR.mean()
    dur_times_minRR_thresh = data_resp.dur_times_minRR_thresh.mean();
    
    v_i=mdl.nodes_db.node['v_Intercept'].trace().mean()
    z_rr=mdl.nodes_db.node['z_min_rr_fixed'].trace().mean()
    v_cor=mdl.nodes_db.node['v_CORonsetSlope_controlledMinRR'].trace().mean()
    v_interaction=mdl.nodes_db.node['v_dur_times_minRR_thresh'].trace().mean()
    v_stim=mdl.nodes_db.node['v_stim'].trace().mean()
    
    ###########
    if i==0:
        x = data_resp.stim
        param_str='v_stim'
        param=v_stim
        x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
        stim=x2
       # ytick=np.linspace(-4,4,5)
    elif i==1:
        x = data_resp.CORslope_controlledMinRR
        param_str='v_CORonsetSlope_controlledMinRR'
        param=v_cor
        x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
        corSlope=x2
      #  ytick=np.linspace(-0.2,0.2,5);
    #
    
    
    ###########
    
    fname=param_str;
    
    v_par_lb=t.loc[param_str, '2.5q']
    v_par_ub=t.loc[param_str, '97.5q']
    y= v_i  + stim*v_stim + corSlope*v_cor + dur_times_minRR_thresh*v_interaction
    y=y-y.mean(); # this is done because the interaction term is messing up with the intercept
        
    y_lb= y + x2*(v_par_lb - param)
    y_ub= y + x2*(v_par_ub - param)
    
    c=(v_par_ub-v_par_lb)/2
    sigma=x.std()
    ci2=np.sqrt(c*c*sigma*sigma+c*c*(x2-np.mean(x))*(x2-np.mean(x)))
    
    ax.plot(x2, y, "-", color=color, linewidth=1.5, alpha=0.7)
    ax.fill_between(x2, y-ci2, y+ci2, color=shade_color,alpha=0.3, edgecolor="")
    #ax.fill_between(x2, y_lb, y_ub, color=color,alpha=0.2, edgecolor="")
    #ax.plot(x2, y_lb, "--", color="0.5", label="95% Prediction Limits")
    #ax.plot(x2, y_ub, "--", color="0.5", label="95% Prediction Limits")
    ax.set_xticks([x.mean()-x.std(),x.mean(), x.mean()+x.std()])
    #ax.set_yticks(ytick)
    #plt.axhline(y=0, color='grey',alpha=0.4, linestyle='dotted')
    plt.savefig(os.path.join(save_dir,fname+'.png'))

###############################################
## Interatction
##############################################
fig, ax = plt.subplots(figsize=fig_size) 
fname='interaction'

for i in range(2):
#choose either one of these
    if i==0:
        data_resp=data_all[data_all.duration_sign<0]
        lbl='short'; color='--k'
    else:
        data_resp=data_all[data_all.duration_sign>0]
        lbl='long'; color='-k'
    ##
    ##### Reading model data
    
    min_rr = data_resp.min_rr_fixed.mean()
    stim = data_resp.stim.mean()
    corSlope = data_resp.CORslope_controlledMinRR.mean()
    dur_times_minRR_thresh = data_resp.dur_times_minRR_thresh.mean();
    
    v_i=mdl.nodes_db.node['v_Intercept'].trace().mean()
    z_rr=mdl.nodes_db.node['z_min_rr_fixed'].trace().mean()
    v_cor=mdl.nodes_db.node['v_CORonsetSlope_controlledMinRR'].trace().mean()
    v_interaction=mdl.nodes_db.node['v_dur_times_minRR_thresh'].trace().mean()
    v_stim=mdl.nodes_db.node['v_stim'].trace().mean()
    
    ##### fitting the model
    
    xx = data_resp.duration_rel
    x =data_resp.min_rr_fixed
    standard_coef=(data_resp.duration_rel*data_resp.min_rr_fixed).std()
    xx=xx/standard_coef
    
    v_par_lb=t.loc['v_dur_times_minRR_thresh', '2.5q']
    v_par_ub=t.loc['v_dur_times_minRR_thresh', '97.5q']
    
    param=v_interaction
    x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
    dur_times_minRR_thresh=xx
    y= v_i  + stim*v_stim + corSlope*v_cor + v_interaction*xx.mean()*x2
    y_lb= y + x2*xx.mean()*(v_par_lb - param)
    y_ub= y + x2*xx.mean()*(v_par_ub - param)
    
    c=(v_par_ub - v_par_lb)/2*xx.mean()
    sigma=x.std()
    ci2=np.sqrt(c*c*sigma*sigma+c*c*(x2-np.mean(x))*(x2-np.mean(x)))

    ax.plot(x2, y, color, linewidth=1.5, alpha=0.7)
    ax.fill_between(x2, y-ci2, y+ci2, color=shade_color,alpha=0.3, edgecolor="")
    ax.set_xticks([x.mean()-x.std(),x.mean(), x.mean()+x.std()])
plt.savefig(os.path.join(save_dir,fname+'.png'))
###########################################################################
###########################non-decision time
###########################################################################  

## 't ~ stim + min_rr_fixed'

for i in range(2):
    color='black'
    data_resp=data_all
    
    min_rr = data_resp.min_rr_fixed.mean()
    stim = data_resp.stim.mean()
    corSlope = data_resp.CORslope_controlledMinRR.mean()
    
    t_i=mdl.nodes_db.node['t_Intercept'].trace().mean()
    t_rr=mdl.nodes_db.node['t_min_rr_fixed'].trace().mean()
    t_stim=mdl.nodes_db.node['t_stim'].trace().mean()
    
    ##
    ###########
    if i==0:
        x = data_resp.min_rr_fixed
        param_str='t_min_rr_fixed'
        param=t_rr
        x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
        min_rr=x2
    #
    elif i==1:
        x = data_resp.stim
        param_str='t_stim'
        param=t_stim
        x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
        stim=x2
    #
    
    
    ###########        
    fig, ax = plt.subplots(figsize=fig_size)
    
    fname=param_str;
    
    v_par_lb=t.loc[param_str, '2.5q']
    v_par_ub=t.loc[param_str, '97.5q']
    y= t_i +  min_rr*t_rr + stim*t_stim
    y_lb= y + x2*(v_par_lb - param)
    y_ub= y + x2*(v_par_ub - param)
    
    c=(v_par_ub-v_par_lb)/2
    sigma=x.std()
    ci2=np.sqrt(c*c*sigma*sigma+c*c*(x2-np.mean(x))*(x2-np.mean(x)))
    
    ax.plot(x2, y, "-", color=color, linewidth=1.5, alpha=0.7)
    
    ax.fill_between(x2, y-ci2, y+ci2, color=shade_color,alpha=0.3, edgecolor="")
    #ax.fill_between(x2, y_lb, y_ub, color=color,alpha=0.2, edgecolor="")
    #ax.plot(x2, y_lb, "--", color="0.5", label="95% Prediction Limits")
    #ax.plot(x2, y_ub, "--", color="0.5", label="95% Prediction Limits")
    ax.set_xticks([x.mean()-x.std(),x.mean(), x.mean()+x.std()])
    ax.set_yticks([0.52, 0.54, 0.56])
#plt.axhline(y=0, color='grey',alpha=0.4, linestyle='dotted')
    plt.savefig(os.path.join(save_dir,fname+'.png'))
#############################################
######## z (initial bias)
#############################################
fig, ax = plt.subplots(figsize=fig_size)

z_i=mdl.nodes_db.node['z_Intercept'].trace().mean()
z_rr=mdl.nodes_db.node['z_min_rr_fixed'].trace().mean()
min_rr = data_resp.min_rr_fixed.mean()
x = data_resp.min_rr_fixed
param_str='z_min_rr_fixed'
param=z_rr
fname=param_str;

x2 = np.linspace(x.quantile(0.025), x.quantile(0.975), 100)
min_rr=x2
v_par_lb=t.loc[param_str, '2.5q']
v_par_ub=t.loc[param_str, '97.5q']
y= z_i +  min_rr*z_rr

c=(v_par_ub-v_par_lb)/2
sigma=x.std()
ci2=np.sqrt(c*c*sigma*sigma+c*c*(x2-np.mean(x))*(x2-np.mean(x)))

ax.plot(x2, y, "-", color=color, linewidth=1.5, alpha=0.7)

ax.fill_between(x2, y-ci2, y+ci2, color=shade_color,alpha=0.3, edgecolor="")
#ax.fill_between(x2, y_lb, y_ub, color=color,alpha=0.2, edgecolor="")
#ax.plot(x2, y_lb, "--", color="0.5", label="95% Prediction Limits")
#ax.plot(x2, y_ub, "--", color="0.5", label="95% Prediction Limits")
ax.set_xticks([x.mean()-x.std(),x.mean(), x.mean()+x.std()])
#plt.axhline(y=0, color='grey',alpha=0.4, linestyle='dotted')
plt.savefig(os.path.join(save_dir,fname+'.png'))

