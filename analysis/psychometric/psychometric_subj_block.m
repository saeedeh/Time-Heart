
base_data_dir='/home/saeedeh/Desktop/aclab-server/Studies/38_TimeHeart'

fname=[base_data_dir,'/analysis/all_trials_for_matlab.csv'];
opts = detectImportOptions(fname);
all_trials=readtable(fname,opts);
all_trials.long_resp= strcmp(all_trials.long_resp,'True');
all_trials.duration= all_trials.duration/1000;

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


% Fit to each subject and block
all_blocks_slope=[];
all_blocks_thresh=[];
all_blocks_mean_pre_rr=[]
subj_ids=unique(all_trials.id);
for ind=1:length(subj_ids)
    for block=1:16
        subj_id=subj_ids(ind)
        subj_id
        block
        inds= (all_trials.id==subj_id) & (all_trials.block==block);% | all_trials.block==block+1 | all_trials.block==block+2);
        if(sum(inds)==0)
            continue;
        end
        trials=varfun(@sum, all_trials(inds,:),'InputVariables','long_resp',...
               'GroupingVariables','duration');  
        StimLevels=[trials.duration']; 
        NumPos=[trials.sum_long_resp'];
        OutOfNum=[trials.GroupCount'];
        [paramsF LL exitflag output] = PAL_PFML_Fit(StimLevels, NumPos, ...
        OutOfNum, params_prior(1,:), [1, 1, 0, 0], PF, 'SearchOptions',options);
        thresh=paramsF(1);
        slope= paramsF(2); 
        all_blocks_slope=[all_blocks_slope, slope];
        all_blocks_thresh=[all_blocks_thresh, thresh];
        all_blocks_mean_pre_rr=[all_blocks_mean_pre_rr, nanmean(all_trials.min_rr_fixed(inds))];
    end
end

all_blocks_slope(all_blocks_slope>250)=170;
valid_inds=~isnan(all_blocks_slope)
[r,p]=corrcoef(all_blocks_slope(valid_inds), all_blocks_mean_pre_rr(valid_inds))
