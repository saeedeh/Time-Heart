#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:52:01 2021

@author: saeedeh
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 09:35:27 2021

@author: saeedeh
"""

models=[1.3 4.9]
N=len(models)
for model_no in models:
    exec(open("set-model-options-4.py").read())
    print("model = "+ str(model_no))
    import time
    t1=time.time()
    mdl.find_starting_values()
    mdl.sample(N, burn=nburn, dbname=trace_name, db='pickle')
#   mdl.print_stats()
    mdl.save(out_name)
    t2=time.time()
    tdif=(t2-t1)/60/60
    import MyEmail
    t=mdl.get_group_nodes()
    t=t.loc[:,['knode_name','mean','std','2.5q','97.5q','mc err']]
    msg="Model "+ str(model_no)+ " is done!\nrun time = "+str(tdif)+" hours \nOutput file name is: "+ mname+"\nN = "+ str(N)+"\nnburn = "+ str(nburn)+ "\nDIC="+ str(mdl.dic)+"\n\n"+t.to_string()
#    MyEmail.sendMail(msg)


#### investigate model outputs
na='N2000_Heart2.3-driftAcc'
na='TRACE-N2000_heart4.9-both'

fname=cwd+'/04-data/'+na
mdl=hddm.load(fname)
mdl=hddm.load(out_name)
t=mdl.get_group_nodes()
t=mdl.get_non_observeds()

mdl.dic
t.knode_name
mname='v';
mname='v_Intercept'
mname='z'
vs=t[t.knode_name==mname]
plt.plot(vs['mean'])

t.loc['v_CORonsetSlope_controlledMinRR']
mdl.plot_posteriors(mname)
mdl.plot_posteriors('v_min_rr_fixed')
mdl.plot_posteriors('v_CORslope_controlledMinRR')
mdl.plot_posteriors('v_COR_onsetRT_slope')
mdl.plot_posteriors('z_onset_to_RT_COR')
mdl.plot_posteriors('v_Intercept')
mdl.plot_posteriors('v_duration')
mdl.plot_posteriors('v_min_rr_fixed')
mdl.plot_posteriors('v_dur_times_minRR_thresh')
mdl.plot_posteriors('v_signed_minRR_thr')
mdl.plot_posteriors('v_stim')
mdl.plot_posteriors('v_CORonsetSlope_controlledMinRR')

tt=mdl.get_group_traces()
(tt.v_CORonsetSlope_controlledMinRR>0).mean()

v1=mdl.nodes_db.node['v_signed_minRR_thr'].trace()
(v1<0).mean()
v2=mdl.nodes_db.node['v_CORonsetSlope_controlledMinRR'].trace()
(v2>0).mean()


vd=mdl.nodes_db.node['v_duration'].trace()
