% Plot dipole/centroid locations (Fig. 2 of eNeuro paper)
%Start eeglab
if ~exist('ALLCOM')
    curr_dir = pwd;
    run('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/eeglab.m');
    cd(curr_dir); close all;
    rmpath('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/plugins/groupSIFT_03/');
end

%%TO DO: check EEG files vs STUDY EEG files
% [STUDY ALLEEG] = pop_loadstudy('filename', 'all_1222018.study',...
%                                'filepath', '/Volumes/Local_Data/VR_connectivity/Data/');

clusters_to_plot=3:10;
colors2use={[1 0 0], [0 1 0], [0 0 1], [1 0 1], [0 1 1], [.8549 .6471 .1255], [1 .4902 0], ... %[1 1 0], [1 .4902 0], ...
    [1 .0784 .5765], [.8 0 1], [.6 0 0], [0 0 0]};

dipfit_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/plugins/dipfit2.3/'; % path to dipfit plugin
hdmfile = [dipfit_pth 'standard_BEM/standard_vol.mat'];
mrifile = [dipfit_pth 'standard_BEM/standard_mri.mat'];


% Load cluster file
load('cluster_data.mat');

% Load ICA_STRUCT files
ica_struct_files = dir('data_release_ICA_structs/*ICA_struct.mat');
dipfits = {};
for i=1:length(ica_struct_files)
    load(['data_release_ICA_structs/' ica_struct_files(i).name]);
    good_comps = ICA_STRUCT.good_comps.brain;
    ICA_STRUCT.dipfit.model = ICA_STRUCT.dipfit.model(good_comps);
    dipfits{i} = ICA_STRUCT.dipfit;
end



k1=1;
k2=1;
compsall=[]; % structure with all of the component posxyz, momxyz, rv
centsall=[]; % structure with all of the cluster centroid posxyz, momxyz, rv
dipsizes = [];

ct = 1;
for i=1:length(clusters_to_plot)
    k2=k2+1;
    for j=1:length(cluster(clusters_to_plot(i)).comps)
        cls_set_i = cluster(clusters_to_plot(i)).sets(1,j);
        comp = cluster(clusters_to_plot(i)).comps(j);
        cluster_dip_models(1,j).posxyz = dipfits{cls_set_i}.model(comp).posxyz;
        cluster_dip_models(1,j).momxyz = dipfits{cls_set_i}.model(comp).momxyz;
        cluster_dip_models(1,j).rv = dipfits{cls_set_i}.model(comp).rv;
        colorsall{ct} = colors2use{i};
        ct= ct+1; k1=k1+1;
    end
    compsall=[compsall,cluster_dip_models(1,:)];
    centsall = [centsall, computecentroid(cluster_dip_models)];
    dipsizes = [dipsizes 25*ones(size(1,length(cluster_dip_models(1,:)))) 40];
    
    clear cluster_dip_models colors1cls
end

%% Plot results
% Plot dipoles
views = [0 0 1; 1 0 0; 0 -1 0];
figure;
for i=1:size(views,1)
    subplot(2,3,i);
    dipplot(compsall,'spheres','on','dipolelength',0,'dipolesize',20,...
        'mri',mrifile,'meshdata',hdmfile,...
        'coordformat',dipfits{1}.coordformat,'color',colorsall,...
        'gui','off');
    view(views(i,:));
end

% Plot centroids
for i=1:size(views,1)
    subplot(2,3,i+size(views,1));
    dipplot(centsall,'spheres','on','dipolelength',0,'dipolesize',50,...
        'mri',mrifile,'meshdata',hdmfile,...
        'coordformat',dipfits{1}.coordformat,'color',colors2use,...
        'gui','off');
    view(views(i,:));
end
