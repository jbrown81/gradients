function [all_betas all_tstats all_ps]=coupling_parameters(grad_slopes,grad_slope_deltas,grad_slope_2deltas,demean)
% coupling_parameters.m
% Jesse Brown
% 01/2022
% jesse.brown@ucsf.edu

%   [all_betas]=coupling_parameters(grad_slopes,grad_slope_deltas,grad_slope_2deltas)
%   returns the coupling parameters based on a linear model of a given
%   dimension's second derivative (grad_slope_2deltas) as a function of all
%   gradients' timeseries (grad_slopes) and first derivative (grad_slope_deltas)

n_comps=size(grad_slopes,2);
all_betas=[];
all_tstats=[];
all_ps=[];
for i=1:n_comps
    cur_comp=i;
    cur_other_comps=setdiff(1:n_comps,cur_comp);
    cur_X=[];
    cur_y=grad_slope_2deltas(:,cur_comp);
    for j=1:n_comps
        cur_comp_2=j;
        cur_X=[cur_X grad_slopes(:,cur_comp_2) grad_slope_deltas(:,cur_comp_2)];
    end
    if demean
        cur_y=cur_y-mean(cur_y);
        cur_X=cur_X-mean(cur_X);
    end
    cur_mdl=fitlm(cur_X,cur_y);
    %B=robustfit(cur_X,cur_y);
    cur_betas=cur_mdl.Coefficients.Estimate;
    %cur_betas=B;
    cur_tstats=cur_mdl.Coefficients.tStat;
    cur_ps=cur_mdl.Coefficients.pValue;
    all_betas=[all_betas;cur_betas'];
    all_tstats=[all_tstats;cur_tstats'];
    all_ps=[all_ps;cur_ps'];
end
end