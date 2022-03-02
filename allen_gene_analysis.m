%% Allen gene expression/gradient correlations
maindir='/Users/jessebrown/Dropbox/latent_fmri/final_paper_data/';
cd(maindir)
% load region gradient weights
load('roi_comp_weights_disc.mat'); % loads 'roi_comp_weights_disc' matrix
%load('roi_comp_weights_val.mat'); % loads 'roi_comp_weights_val' matrix
load('roi_comp_weights_val_aligned.mat'); % loads 'roi_comp_weights_val_aligned' matrix
roi_comp_weights_val=roi_comp_weights_val_aligned;

% load region allen gene expression values from abagen
load('allen_expression_brainnetome_261_15655.mat'); % loads 'exp_data' matrix
n_genes_total=size(exp_data,2);
load('all_allen_genes.mat'); % loads 'gene_names' cell array
gene_names=vertcat(gene_names{:});

% exclude brainnetome/SUIT nodes with no allen gene expression from abagen
load('keep_nodes_allen_brainnetome_261.mat'); % loads 'keep_nodes' array
out_nodes=setdiff(1:273,keep_nodes);
keep_nodes_cort=find(keep_nodes<=210); % subset of indices out of 261 that are cortical (have original 273 index <=210)
keep_nodes_subcort=intersect(find(keep_nodes>210),find(keep_nodes<247));

% exclude nodes based on data-driven outlier detection
% based on 261 indexing, not 273
out_nodes_kmeans=[111 113 116 158 159 162 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 244 253 257 258 259 261];
%out_nodes_kmeans=[];
keep_nodes_kmeans=setdiff(1:261,out_nodes_kmeans);
% keep final set of 196 nodes
% indexed out of 261
keep_nodes_final=intersect(keep_nodes_cort,keep_nodes_kmeans);
out_nodes_final=setdiff(1:261,keep_nodes_final);
% indexed out of 273
keep_nodes_final_273=keep_nodes(keep_nodes_final);
out_nodes_final_273=setdiff(1:273,keep_nodes_final_273);


%% gradient-by-gene expression spatial correlations
n_grads=6;
use_spearman=false;
if use_spearman
    [rs_indep_discovery ps_indep_discovery]=corr(exp_data(keep_nodes_final,:),roi_comp_weights_disc(keep_nodes_final_273,1:n_grads),'type','spearman');
    [rs_indep_validation ps_indep_validation]=corr(exp_data(keep_nodes_final,:),roi_comp_weights_val(keep_nodes_final_273,1:n_grads),'type','spearman');
else
    [rs_indep_discovery, ps_indep_discovery]=corr(exp_data(keep_nodes_final,:),roi_comp_weights_disc(keep_nodes_final_273,1:n_grads));
    [rs_indep_validation ps_indep_validation]=corr(exp_data(keep_nodes_final,:),roi_comp_weights_val(keep_nodes_final_273,1:n_grads));
end

%% creating csv files with all r and p-values
% gene names and column headers added posthoc in excel
out1=zeros(15655,n_grads);
out2=zeros(15655,n_grads);
count=1;
for i=1:n_grads
    out1(:,count)=rs_indep_discovery(:,i);
    out1(:,count+1)=ps_indep_discovery(:,i);
    out2(:,count)=rs_indep_validation(:,i);
    out2(:,count+1)=ps_indep_validation(:,i);
    count=count+2;
end

csvwrite('gradients_genes_correlation_discovery.csv',out1)
csvwrite('gradients_genes_correlation_validation.csv',out2)