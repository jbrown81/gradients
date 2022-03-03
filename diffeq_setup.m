% diffeq_setup.m
% Jesse Brown
% 01/2022
% jesse.brown@ucsf.edu

% this script will:
% 1) load latent data (ie gradient slope timeseries), which is fmri voxelwise timeseries projected into PCA space
% 2) load task regressors
% 3) compute gradient timeseries first and second derivatives for rest and task data
% 4) compute rest and task condition coupling parameters

%% load latent data
maindir='/Users/jessebrown/Dropbox/latent_fmri/final_paper_data/';
cd(maindir)
load_data=false;
if load_data
    % resting gradient timeseries data for 100 discovery subjects + 100 validation subjects
    components_pca=readNPY('grad_ts_pca.npy');
    components_pca_val=readNPY('grad_ts_pca_val_proj_disc.npy');
    % task gradient timeseries data for all 200 discovery + validation subjects
    components_pca_task_disc=readNPY('grad_ts_pca_task_proj_disc.npy');
    components_pca_task_val=readNPY('grad_ts_pca_task_val_proj_disc.npy');
end

%% get fMRI volume indices for specific phases of task
task='wm'; % choose one of: 'wm', 'motor', 'language', 'emotion'

% original scan length/original number of volumes/final number after preprocessing
%Working Memory 5:01/405/400
%Motor 3:34/284/274
%Language  3:57/316/311
%Emotion Processing 2:16/176/171
n_subjs_task=200;
switch task
    case 'wm'
        n_vols_per_scan_task=400;
        components_pca_task=[components_pca_task_disc(1:40000,:);components_pca_task_val(1:40000,:)];
    case 'motor'
        n_vols_per_scan_task=274;
        components_pca_task_pre=[components_pca_task_disc(40001:67900,:);components_pca_task_val(40001:67900,:)];
        % for motor task, need to chop off first 5 volumes to align with task regressors file
        keep_inds=[];
        for i=1:200
            keep_inds=[keep_inds;((6:279)+(279*(i-1)))'];
            %keep_inds=[keep_inds;((1:274)+(279*(i-1)))'];
        end
        components_pca_task=components_pca_task_pre(keep_inds,:);
    case 'language'
        n_vols_per_scan_task=311;
        components_pca_task=[components_pca_task_disc(67901:99000,:);components_pca_task_val(67901:99000,:)];
    case 'emotion'
        n_vols_per_scan_task=171;
        components_pca_task=[components_pca_task_disc(99001:116100,:);components_pca_task_val(99001:116100,:)];
end

% specify current task regressors (blocks and HRF-convolved blocks)
pidns=[000001 000002 000003 000004 000005 000006 000007 000008 000009 000010 000011 000012 000013 000014 000015 000016 000017 000018 000019 000020 000021 000022 000023 000024 000025 000026 000027 000028 000029 000030 000031 000032 000033 000034 000035 000036 000037 000038 000039 000040 000041 000042 000043 000044 000045 000046 000047 000048 000049 000050 000051 000052 000053 000054 000055 000056 000057 000058 000059 000060 000061 000062 000063 000064 000065 000066 000067 000068 000069 000070 000071 000072 000073 000074 000075 000076 000077 000078 000079 000080 000081 000082 000083 000084 000085 000086 000087 000088 000089 000090 000091 000092 000093 000094 000095 000096 000097 000098 000099 000100 000101 000102 000103 000104 000105 000106 000107 000108 000109 000110 000111 000112 000113 000114 000115 000116 000117 000118 000119 000120 000121 000122 000123 000124 000125 000126 000127 000128 000129 000130 000131 000132 000133 000134 000135 000136 000137 000138 000139 000140 000141 000142 000143 000144 000145 000146 000147 000148 000149 000150 000151 000152 000153 000154 000155 000156 000157 000158 000159 000160 000161 000162 000163 000164 000165 000166 000167 000168 000169 000170 000171 000172 000173 000174 000175 000176 000177 000178 000179 000180 000181 000182 000183 000184 000185 000186 000187 000188 000189 000190 000191 000192 000193 000194 000195 000196 000197 000198 000199 000200];
task_regressors_all=zeros(n_vols_per_scan_task,4,n_subjs_task);
for i=1:n_subjs_task
    task_regressors=load(sprintf('hcp_task_regressors/%s/regressors_%06s.txt',task,int2str(pidns(i))));
    task_regressors_all(:,:,i) = task_regressors;
end

% find indices (and associated measures) for timepoints that are during task condition 1 or 2
block1_inds_all=[];
block2_inds_all=[];
block1_inds_all_subj=[];
block2_inds_all_subj=[];
active_thr=.5;
for i=1:n_subjs_task
    k1=find(task_regressors_all(:,2,i)>active_thr); % find timepoints where HRF-convolved task regressor 1 is 'active'
    k2=find(task_regressors_all(:,4,i)>active_thr); % find timepoints where HRF-convolved task regressor 2 is 'active'
    cur_offset=n_vols_per_scan_task*(i-1);
    cur_inds=(1:n_vols_per_scan_task)+cur_offset;
    block1_inds=cur_inds(k1);
    block1_inds_all=[block1_inds_all;block1_inds'];
    block1_inds_all_subj=[block1_inds_all_subj;ones(length(block1_inds),1)*i];
    block2_inds=cur_inds(k2);
    block2_inds_all=[block2_inds_all;block2_inds'];
    block2_inds_all_subj=[block2_inds_all_subj;ones(length(block2_inds),1)*i];
end

block1_seqs=zeros(length(block1_inds_all),1);
k=find([1 diff(block1_inds_all)']>1);
count=1;
for i=1:length(block1_inds_all)
    if ismember(i,k)
        count=count+1;
    end
    block1_seqs(i)=count;
end
block1_block_count=count;

block2_seqs=zeros(length(block2_inds_all),1);
k=find([1 diff(block2_inds_all)']>1);
count=1;
for i=1:length(block2_inds_all)
    if ismember(i,k)
        count=count+1;
    end
    block2_seqs(i)=count;
end
block2_block_count=count;


%% get gradient slopes, derivatives (velocities), second derivatives (accelerations)
% determined using 'gradient()' function in matlab, not to be confused with
% the term 'gradient' used throughout the paper
keep_comps=[1:6]; % keep_comps is the subset of latent dimensions/gradients to use

% task data
grad_slopes_task=components_pca_task(:,keep_comps);
[grad_slope_deltas_task grad_slope_2deltas_task]=latent_derivatives(grad_slopes_task,n_subjs_task,n_vols_per_scan_task);

% set up state-specific gradient slope/derivative/second derivative arrays
grad_slopes_task_block1=grad_slopes_task(block1_inds_all,:);
%grad_slopes_task_block1=grad_slopes_task_block1+10;
%grad_slopes_task_block1=grad_slopes_task_block1-mean(grad_slopes_task_block1);
grad_slope_deltas_task_block1=grad_slope_deltas_task(block1_inds_all,:);
grad_slope_2deltas_task_block1=grad_slope_2deltas_task(block1_inds_all,:);
grad_slopes_task_block2=grad_slopes_task(block2_inds_all,:);
%grad_slopes_task_block2=grad_slopes_task_block2-mean(grad_slopes_task_block2);
grad_slope_deltas_task_block2=grad_slope_deltas_task(block2_inds_all,:);
grad_slope_2deltas_task_block2=grad_slope_2deltas_task(block2_inds_all,:);
% if wanting to estimate coupling parameters only from discovery dataset
if true
    grad_slopes_task_block1=grad_slopes_task_block1(find(block1_inds_all_subj<=100),:);
    grad_slopes_task_block2=grad_slopes_task_block2(find(block2_inds_all_subj<=100),:);
    grad_slope_deltas_task_block1=grad_slope_deltas_task_block1(find(block1_inds_all_subj<=100),:);
    grad_slope_deltas_task_block2=grad_slope_deltas_task_block2(find(block2_inds_all_subj<=100),:);
    grad_slope_2deltas_task_block1=grad_slope_2deltas_task_block1(find(block1_inds_all_subj<=100),:);
    grad_slope_2deltas_task_block2=grad_slope_2deltas_task_block2(find(block2_inds_all_subj<=100),:);
end


%% estimate coupling parameters of differential equation
% gradient acceleration = gradient position + gradient velocity
demean_grads=true;
[betas_cond1 ts_cond1 ps_cond1]=coupling_parameters(grad_slopes_task_block1,grad_slope_deltas_task_block1,grad_slope_2deltas_task_block1,demean_grads);
[betas_cond2 ts_cond2 ps_cond2]=coupling_parameters(grad_slopes_task_block2,grad_slope_deltas_task_block2,grad_slope_2deltas_task_block2,demean_grads);
grad_means_cond1=mean(grad_slopes_task_block1);
grad_means_cond2=mean(grad_slopes_task_block2);

process_rest=true;
if process_rest
    % resting data
    n_subjs_rest=100;
    n_vols_per_scan_rest=1195;
    grad_slopes_rest=components_pca(:,keep_comps);
    [grad_slope_deltas_rest grad_slope_2deltas_rest]=latent_derivatives(components_pca(:,keep_comps),n_subjs_rest,n_vols_per_scan_rest);
    [betas_rest ts_rest ps_rest]=coupling_parameters(grad_slopes_rest,grad_slope_deltas_rest,grad_slope_2deltas_rest,demean_grads);
    grad_means_rest=mean(grad_slopes_rest);
end



