base_data_dir='/mnt/Aclab/Studies/38_TimeHeart'
base_data_dir='/home/saeedeh/Desktop/aclab-server/Studies/38_TimeHeart'

fname=[base_data_dir,'/analysis/all_trials.csv'];
opts = detectImportOptions(fname);
all_trials=readtable(fname,opts);
fname=[base_data_dir,'/analysis/all_subjects.csv'];
opts = detectImportOptions(fname);
all_subjects=readtable(fname,opts);


sys_inds=all_trials.systole==1;
dias_inds=all_trials.systole==0;
all_trials.long_resp= strcmp(all_trials.long_resp,'True');
res_sys=varfun(@sum, all_trials(sys_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration')
res_dias=varfun(@sum, all_trials(~sys_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration')

StimLevels=[res_sys.duration';res_dias.duration']; 
NumPos=[res_sys.sum_long_resp';res_dias.sum_long_resp'];
OutOfNum=[res_sys.GroupCount';res_dias.GroupCount'];



% fitting parameters
PF = @PAL_Logistic;
options = PAL_minimize('options');   %PAL_minimize search options
options.TolFun = 1e-12;     %Increase desired precision on LL
options.TolX = 1e-12;       %Increase desired precision on parameters
options.MaxFunEvals = 5000; %Allow more function evals
options.MaxIter = 5000;     %Allow more iterations
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 0.06];        %Range on lapse rates. Will go ignored here
                            %since lapse rate is not a free parameter                            
params_prior=[0.134, 100,0,0.02; 
    0.134, 100,0,0.02];



%% Fit to each subject
subj_ids=unique(all_trials.id);
for ind=1:length(subj_ids)
   subj_id=subj_ids(ind)
   subj_inds=all_trials.id==subj_id;

       %all
   subj_all=varfun(@sum, all_trials(subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
  
StimLevels=[subj_all.duration']; 
 NumPos=[subj_all.sum_long_resp'];
OutOfNum=[subj_all.GroupCount'];
[paramsF LL exitflag output] = PAL_PFML_Fit(StimLevels, NumPos, ...
    OutOfNum, params_prior(1,:), [1, 1, 0, 0], PF, 'SearchOptions',options);
    allParams(ind,:)=paramsF;
 
    %sys/dias
   subj_sys=varfun(@sum, all_trials(sys_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
   subj_dias=varfun(@sum, all_trials(~sys_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');

StimLevels=[subj_sys.duration';subj_dias.duration']; 
 NumPos=[subj_sys.sum_long_resp';subj_dias.sum_long_resp'];
OutOfNum=[subj_sys.GroupCount';subj_dias.GroupCount'];
[paramsF LL exitflag output trash numParams2T2S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params_prior, PF, 'thresholds','unconstrained','slopes',...
    'unconstrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);
    allParams_sys(ind,:)=paramsF(1,:);
    allParams_dias(ind,:)=paramsF(2,:);
   
    % OR
       highOR_inds=strcmp(all_trials.orienting2,'True');
   subj_highOR=varfun(@sum, all_trials(highOR_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
   subj_lowOR=varfun(@sum, all_trials(~highOR_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
    StimLevels=[subj_highOR.duration';subj_lowOR.duration']; 
 NumPos=[subj_highOR.sum_long_resp';subj_lowOR.sum_long_resp'];
OutOfNum=[subj_highOR.GroupCount';subj_lowOR.GroupCount'];
[paramsF LL exitflag output trash numParams2T2S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params_prior, PF, 'thresholds','unconstrained','slopes',...
    'unconstrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);
    allParams_highOR(ind,:)=paramsF(1,:);
    allParams_lowOR(ind,:)=paramsF(2,:);
    
    
    % pre_rr
    subj_median_rr= nanmedian(all_trials.pre_rr(subj_inds));
     highPreRR_inds= all_trials.pre_rr>=subj_median_rr;
     lowPreRR_inds= all_trials.pre_rr<subj_median_rr;
   subj_highRR=varfun(@sum, all_trials(highPreRR_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
   subj_lowRR=varfun(@sum, all_trials(lowPreRR_inds & subj_inds,:),'InputVariables','long_resp',...
       'GroupingVariables','duration');
    StimLevels=[subj_highRR.duration';subj_lowRR.duration']; 
 NumPos=[subj_highRR.sum_long_resp';subj_lowRR.sum_long_resp'];
OutOfNum=[subj_highRR.GroupCount';subj_lowRR.GroupCount'];
[paramsF LL exitflag output trash numParams2T2S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params_prior, PF, 'thresholds','unconstrained','slopes',...
    'unconstrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);
    allParams_highRR(ind,:)=paramsF(1,:);
    allParams_lowRR(ind,:)=paramsF(2,:);
end

all_subjects.thresh_all=allParams(:,1);
all_subjects.thresh_sys=allParams_sys(:,1);
all_subjects.thresh_dias=allParams_dias(:,1);
all_subjects.slope_all=allParams(:,2);
all_subjects.slope_sys=allParams_sys(:,2);
all_subjects.slope_dias=allParams_dias(:,2);
all_subjects.thresh_sys_dias=allParams_sys(:,1)-allParams_dias(:,1);
all_subjects.slope_sys_dias=allParams_sys(:,2)-allParams_dias(:,2);

plot(all_subjects.thresh_all, all_subjects.slope_all,'o', 'MarkerFaceColor', 'b')
xlabel('threshold (bias)'); ylabel('slope (accuracy)')
[h,p,c,st]=ttest(allParams_sys(:,2)-allParams_dias(:,2))
[h,p]=ttest(allParams_highOR(:,2)-allParams_lowOR(:,2))
[h,p]=ttest(allParams_highRR(:,1)-allParams_lowRR(:,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% between subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r,p]=corrcoef(all_subjects.orienting_tVal, allParams(:,1))
[r,p]=corrcoef(all_subjects.orienting_tVal, allParams(:,2))
[r,p]=corrcoef( allParams(:,1), allParams(:,2))
plot( allParams(:,1), allParams(:,2),'o') %bias to choose short means poorer accuracy

[r,p]=corrcoef(all_subjects.orienting_tVal, all_subjects.mean_RT)
%people who orient more are faster in responding

[r,p]=corrcoef(all_subjects.orienting_tVal, allParams_sys(:,1)-allParams_dias(:,1))
%people who orient more have a bias to perceive systole as longer

[r,p]=corrcoef(all_subjects.mean_RT, allParams_highRR(:,2)-allParams_lowRR(:,2))

[r,p]=corrcoef(allParams(:,), all_subjects.mean_RT)


plot(all_subjects.mean_count_orienting, allParams(:,1),'o')
median_rr_change=median(all_subjects.mean_rr_change2)
plotData(all_trials(all_trials.mean_rr_change2>median_rr_change,:))
hold on;
plotData(all_trials(all_trials.mean_rr_change2<=median_rr_change,:))
hold off


fitlm(all_subjects,'orienting_tVal~thresh_sys_dias')
[r,p]=corrcoef(all_subjects.orienting_tVal, allParams_sys(:,1)-allParams_dias(:,1))
