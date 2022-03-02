% diffeq_forecast.m
% Jesse Brown
% 01/2022
% jesse.brown@ucsf.edu

% gradient timeseries simulation
% this script is intended to run after running diffeq_setup.m with a
% specific task specified ('wm','motor','language', or 'emotion')
% other required variables from diffeq_setup.m: 
% 'grad_slopes_rest', 'grad_slope_deltas_rest', 'block1_inds_all', 'block2_inds_all', 
% 'components_pca', 'components_pca_val', 'components_pca_task',
% 'grad_means_rest', 'grad_means_cond1', 'grad_means_cond2'

% notes:
% In matlab, the standard syntax for differential equation solutions uses 'y' as the output variable
% for second order equations. Here we use 'y_pos'/'yp' as the output variable for y position 
% and 'yv' for y velocity


%% simulation
% simulate timeseries starting from resting intial conditions: rest/condition_1, OR rest/condition_2
load('rand10k.mat'); % set of 10000 random resting state timepoints used in paper
n_sims=1000;
keep_comps=[1:6]; % keep_comps is the subset of latent dimensions/gradients to use

% load coupling parameters
switch task
    case 'wm'
        load('coupling_parameters_wm_block1_discovery.mat'); % load 'betas_cond1'
        load('coupling_parameters_wm_block2_discovery.mat'); % load 'betas_cond2'
    case 'motor'
        load('coupling_parameters_motor_block1_discovery.mat');
        load('coupling_parameters_motor_block2_discovery.mat');
    case 'language'
        load('coupling_parameters_language_block1_discovery.mat');
        load('coupling_parameters_language_block2_discovery.mat');
    case 'emotion'
        load('coupling_parameters_emotion_block1_discovery.mat');
        load('coupling_parameters_emotion_block2_discovery.mat');
end
load('coupling_parameters_rest_discovery.mat'); % load 'betas_rest'  

% make simulations using gradient_ode.m function
% 'pred_interval' is set to 200 timepoints; and tspan from 0-200 with stepsize of 1;
pred_interval=200;
for i=1:n_sims
    cur_ind=rand10k(i);
    
    % simulate 200 timepoints with gradient slope/velocity
    % initial conditions randomly selected from resting state data sample, using
    % task-free (resting state) coupling parameters
    [yp0 yv0]=gradient_ode(grad_slopes_rest(cur_ind,:),grad_slope_deltas_rest(cur_ind,:),betas_rest,pred_interval);
    if demean_grads
        yp0=yp0+grad_means_rest;
    end
    yp0_all(:,:,i)=yp0;

    % simulate 200 timepoints with gradient slope/velocity
    % initial conditions randomly selected from resting state data sample, using
    % task conditon 1 coupling parameters
    [yp1 yv1]=gradient_ode(grad_slopes_rest(cur_ind,:),grad_slope_deltas_rest(cur_ind,:),betas_cond1,pred_interval);
    if demean_grads
        yp1=yp1+grad_means_cond1;
    end
    yp1_all(:,:,i)=yp1;
    
    % simulate 200 timepoints with gradient slope/velocity
    % initial conditions randomly selected from resting state data sample, using
    % task conditon 2 coupling parameters
    [yp2 yv2]=gradient_ode(grad_slopes_rest(cur_ind,:),grad_slope_deltas_rest(cur_ind,:),betas_cond2,pred_interval);
    if demean_grads
        yp2=yp2+grad_means_cond2;
    end
    yp2_all(:,:,i)=yp2;

    disp(i)
end

% compute simulated data FC matrices
fcmats_0=zeros(273,273,n_sims);
fcmats_1=zeros(273,273,n_sims);
fcmats_2=zeros(273,273,n_sims);
% for each simulation, determine 273 region timeseries from n gradient timeseries
% then calculate 273 x 273 FC matrix

% load in region gradient weights
load('roi_comp_weights_disc.mat'); % loads 'roi_comp_weights_disc' matrix
roi_comp_weights=roi_comp_weights_disc;

sim_region_ts_0_mean=zeros(273,n_sims);
sim_region_ts_1_mean=zeros(273,n_sims);
sim_region_ts_2_mean=zeros(273,n_sims);
sim_region_ts_1_all=zeros(n_sims*(pred_interval+1),273);
sim_region_ts_2_all=zeros(n_sims*(pred_interval+1),273);
for i=1:n_sims
    cur=squeeze(yp0_all(:,:,i));
    sim_region_ts_0=cur*roi_comp_weights(:,keep_comps)';
    fcmats_0(:,:,i)=corr(sim_region_ts_0);
end
for i=1:n_sims
    cur=squeeze(yp1_all(:,:,i));
    sim_region_ts_1=cur*roi_comp_weights(:,keep_comps)';
    fcmats_1(:,:,i)=corr(sim_region_ts_1);
end
for i=1:n_sims
    cur=squeeze(yp2_all(:,:,i));
    sim_region_ts_2=cur*roi_comp_weights(:,keep_comps)';
    fcmats_2(:,:,i)=corr(sim_region_ts_2);
end
fcmat_0_sim_mean=mean(fcmats_0,3);
fcmat_1_sim_mean=mean(fcmats_1,3);
fcmat_2_sim_mean=mean(fcmats_2,3);

% compute actual data FC matrices
if false
    % get matrix based on all discovery + validation subjects
    fcmat_0_real=corr(components_pca(:,keep_comps)*roi_comp_weights(:,keep_comps)');
    fcmat_1_real=corr(components_pca_task(block1_inds_all,keep_comps)*roi_comp_weights(:,keep_comps)');
    fcmat_2_real=corr(components_pca_task(block2_inds_all,keep_comps)*roi_comp_weights(:,keep_comps)');
else
    % get matrix based only on validation subjects
    fcmat_0_real=corr(components_pca_val(:,keep_comps)*roi_comp_weights(:,keep_comps)');
    fcmat_1_real=corr(components_pca_task(block1_inds_all(find(block1_inds_all_subj>100)),keep_comps)*roi_comp_weights(:,keep_comps)');
    fcmat_2_real=corr(components_pca_task(block2_inds_all(find(block2_inds_all_subj>100)),keep_comps)*roi_comp_weights(:,keep_comps)');
end




