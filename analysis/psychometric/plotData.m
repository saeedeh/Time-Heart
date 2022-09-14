function [] = plotData(all_trials, subj_id)
if nargin==1
    trials=all_trials;
else
    trials=all_trials(all_trials.id==subj_id, :);
end
    res=varfun(@mean, trials,'InputVariables','long_resp',...
       'GroupingVariables','duration');
    plot(res.duration, res.mean_long_resp,'ko')

    
end

