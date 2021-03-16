%% Plot dipole/centroid locations (eNeuro Fig. 2)
% Start eeglab
eeglab_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/';
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
clusters_to_plot=3:10;
colors2use={[1 0 0], [0 1 0], [0 0 1], [1 0 1], [0 1 1], ...
            [.8549 .6471 .1255], [1 .4902 0], [1 .0784 .5765], ...
            [.8 0 1], [.6 0 0], [0 0 0]};

% Load cluster file
load('cluster_data.mat');

% Load ICA_STRUCT files
ica_struct_files = dir('data_release_ICA_structs/*ICA_struct.mat');
dipfits = PLTFUNCS.load_ica_structs(ica_struct_files);

% Aggregate component info
[compsall,centsall,colorsall] = PLTFUNCS.aggregate_comps(clusters_to_plot,...
                                                         dipfits, cluster,...
                                                         colors2use);

% Plot dipoles
PLTFUNCS.plot_dipoles(compsall,centsall,dipfits,eeglab_pth,...
                      colors2use,colorsall)
