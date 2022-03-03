Data and code associated with:  
**A dynamic gradient architecture generates brain activity states**  
Jesse A. Brown, Alex J. Lee, Lorenzo Pasquini, William W. Seeley  
currently under review. Preprint is available at: https://www.biorxiv.org/content/10.1101/2020.08.12.248112v3

## **Notes:**
- Gradient map nii files are in **pca_grad_maps/**
- To run the dynamical systems model, the main scripts are **diffeq_setup.m** and **diffeq_forecast.m**
- Voxelwise PCA code is in **fmri_voxelwise_pca.ipynb**
- Allen gene expression analysis is in **allen_gene_analysis.m**


## **Scripts and Functions:**
- **diffeq_setup.m:** script to load gradient slope timeseries, HCP task regressors, and compute gradient coupling parameters
- **diffeq_forecast.m:** scripts to use gradient coupling parameters and simulate gradient timeseries 
- **coupling_parameters.m:** function to compute coupling parameters
- **gradient_ode.m:** function to solve coupled differential equations
- **allen_gene_analysis.m:** script to perform gradient/gene expression spatial correlations
- **fmri_voxelwise_pca.ipynb:** notebook to perform voxelwise PCA and project task data or validation data into discovery latent space

## **Data:**
- **grad_ts_pca.npy:** 119500x100 gradient slope timeseries for task-free discovery dataset
- **grad_ts_pca_val_proj_disc.npy:** 119500x100 gradient slope timeseries for task-free validation dataset, projected into task-free discovery dataset latent space
- **grad_ts_pca_task_proj_disc.npy:** 116100x100 gradient slope timeseries for task discovery dataset, projected into task-free discovery dataset latent space
- **grad_ts_pca_task_val_proj_disc.npy:** 116100x100 gradient slope timeseries for task validation dataset, projected into task-free discovery dataset latent space
- **roi_comp_weights_disc.mat, roi_comp_weights_val.mat, roi_comp_weights_task.mat:** region-wise gradient weight matrices for the first 100 PCA components from the task-free discovery dataset, task-free validation dataset, and task discovery dataset, each based on independent PCA estimation
- **roi_comp_weights_val_aligned.mat:** region-wise gradient weight matrices for the first 12 PCA components from the task-free validation dataset, Procrustes aligned to the first 12 PCA components from the task-free discovery dataset
- **pca_grad_maps/pca_grad_disc_001.nii.gz, ... pca_grad_disc_100.nii.gz, pca_grad_val_001.nii.gz, ... pca_grad_val_100.nii.gz, pca_grad_disc_task_001.nii.gz, ... pca_grad_disc_task_100.nii.gz:** voxel-wise gradient weight maps for the first 100 PCA components from the task-free discovery dataset, task-free validation dataset, and task discovery dataset, each based on independent PCA estimation
- **coupling_parameters_emotion_block1_discovery.mat, coupling_parameters_motor_block2_discovery.mat, coupling_parameters_emotion_block2_discovery.mat, coupling_parameters_rest_discovery.mat, coupling_parameters_language_block1_discovery.mat coupling_parameters_wm_block1_discovery.mat, coupling_parameters_language_block2_discovery.mat coupling_parameters_wm_block2_discovery.mat, coupling_parameters_motor_block1_discovery.mat:** the gradient coupling parameter matrices for task-free discovery dataset and each task condition from the discovery dataset
- **273_rois/vol_1.nii.gz, ... vol_273.nii.gz:** 273 regions of interest using a parcellation of 246 cortical and subcortical regions from the Brainnetome atlas (http://www.brainnetome.org/) and 27 cerebellar regions from the SUIT atlas (http://www.diedrichsenlab.org/imaging/suit.htm); region names are listed in brainnetome_region_abbrevs_names_273.txt
- **allen_expression_brainnetome_261_15655.mat:** Allen Human Brain spatial gene expression patterns for 261 out 273 regions for 15655 genes; 261/273 indices are listed in keep_nodes_allen_brainnetome_261.mat, gene names are listed in all_allen_genes.mat; expression values were derived using abagen (https://abagen.readthedocs.io/en/latest/usage.html)
- **gradients_genes_correlation_discovery.csv:** Gradient/gene expression spatial correlation coefficients for 15655 genes for the first six gradients in the task-free discovery dataset
- **gradients_genes_correlation_validation.csv:** first six gradients in the task-free validation dataset (after Procrustes alignment to task-free discovery dataset)


## **Helper functions:**
- readNPY.m, readNPYheader.m
