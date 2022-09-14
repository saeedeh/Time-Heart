#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 12:06:55 2022

@author: saeedeh
"""

mdl.dic
v1=mdl.nodes_db.node[m4].trace()
(v1>0).mean()
 
m=mdl.get_stochastics()
np.unique(m.knode_name)
m[np.array(m.knode_name=='t_min_rr_fixed')]
#model prediction check: http://ski.clps.brown.edu/hddm_docs/tutorial_post_pred.html
sim_data=hddm.utils.post_pred_gen(mdl, samples=50)
ppc_compare = hddm.utils.post_pred_stats(trials_resp, sim_data)
print ppc_compare


## plot coefficients

v1=mdl.nodes_db.node['v_signed_minRR_thr'].trace()
(v1<0).mean()
v2=mdl.nodes_db.node['v_CORonsetSlope_controlledMinRR'].trace()
(v2>0).mean()


x0=0
v=v1
fig, ax = plt.subplots(figsize=[4,3])
ax=sns.distplot(v,color='black', hist=False,kde_kws = {'shade': False, 'linewidth': 1.5})
kde_x, kde_y = ax.lines[0].get_data()
p1 = plt.axvline(x=x0,linestyle='--',color='grey')
ax.fill_between(kde_x, kde_y, where=(kde_x<x0) , interpolate=True, color='#ffc7cf')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


sns.distplot(v2, hist=False,kde_kws = {'shade': True, 'linewidth': 3},color='orange')


