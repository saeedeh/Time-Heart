#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 17:41:49 2021

@author: saeedeh
"""

##################################
## bar plots!
##################################

m=1 #preRR
m=2 #cor time
m=3 #RT
m=4 #cor magnitude
m=5 #onset_to_RT_COR

tar=1 #consistency
tar=2 #bias
cmap = plt.get_cmap('PiYG', 11)    # 11 discrete colors
m=2; tar=1
for m in range(1, 5):
    for tar in range(1,3):
        R=1
        fig_R=0.6
        if m==1:
            measure='min_rr_fixed'
            yrange=[750, 810]
            ylab='preRR'
        elif m==2:
            measure='max_time'
            yrange=[1220, 1550]
            ylab='OR peak time'
        elif m==3:
            measure = 'RT_real'
            ylab='RT'
            yrange=[790,1400]
        elif m==4:
            measure='min_to_max_OR_flex'
            ylab='COR'
            yrange=[30,55]

        
        ####  acc & bias bar plots
        #acc
        if tar==1:
            type_str='acc';labels=['consistant','distortion']
            inds1=all_trials.TP_acc2==1 #consistent
            inds2= (all_trials.TP_acc2==0)
            col1='blue';col2='red'
        elif tar==2:
        #bias
            type_str='bias'; labels=['contraction','dilation']
            inds1=all_trials.TP_bias2==-1 #contraction
            inds2= all_trials.TP_bias2==1#dilation
            col1='yellow';col2='green'
        col1='k';col2='k'
        Y=all_trials[measure]
        [t1,s1]=subj_mean_mean_se(Y, inds1)
        [t2,s2]=subj_mean_mean_se(Y, inds2)
        #plt.style.use('fivethirtyeight')
        plt.style.use('default')
        #plt.rcParams['axes.facecolor'] = 'white'
        fig, ax=plt.subplots(figsize=(4*fig_R,3*fig_R))
        plt.bar(range(2),[t1*R, t2*R ], width=0.3, color=[col1,col2])
        plt.errorbar(range(2),[t1*R, t2*R],yerr= [s1*R, s2*R ], ls='None', elinewidth=3,color='indianred')
        plt.xticks(range(2), labels, size='small')
        #plt.ylabel(ylab)
        #ax.spines['right'].set_visible(False)
        #ax.spines['top'].set_visible(False)
        plt.tight_layout()
        plt.xlim([-0.5,1.5])
        plt.ylim(yrange)
        plt.locator_params(axis='y', nbins=4)
        
        save_img_dir='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/Analysis/images/barplots2'
        img_fname=ylab+'-'+type_str+'.png'
        plt.savefig(os.path.join(save_img_dir, img_fname))

## RT vs COR & min_RR
for i in range(3):
    if i==0:
        mname='min_to_max_OR_flex'
        yrange=[30,55]
        ylab='COR'
    elif i==1:
        mname='min_rr_fixed'
        yrange=[750, 800]
        ylab="preRR"
    elif i==2:
        mname='max_time'
        yrange=[1220, 1550]
        ylab="latency"

    fig, ax=plt.subplots(figsize=(3,2))
    fitting_line( all_trials.RT_real, all_trials[mname])
    plt.ylim(yrange)
    plt.xlabel('RT (ms)')
    plt.ylabel(ylab)
    plt.tight_layout()
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    
    save_img_dir='/home/saeedeh/Dropbox/Cornell Research/Codes/TP-heart/Analysis/images/barplots2'
    img_fname=ylab+'-RT'+'.png'
    plt.savefig(os.path.join(save_img_dir, img_fname))