%% Plot EEG spectra, ERSP's, and ERSP timing (eNeuro Figs. 7-10)
tic;
% Start eeglab
eeglab_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/';
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
cluster2plt = 'Supplementary Motor Area'; % select cluster to plot ERSPs from

eeg_files = dir('data_release_final/*.set');
evs = {'M_on_SVZ','M_on_WVZ','pull_Stn','pull_Wlk'};
ev_labels = {'Stand Rotate','Walk Rotate','Stand Pull','Walk Pull'};
ep_add_time = [-1 1]*0.6; % extra epoch times that will be cut off during ERSP computation
epochTimes=[-.5 1.5]+ep_add_time;
alpha = 0.05; icachansind = 1:128;
cluster_labels = {'Left Occipital','Right Occipital','Left Sensorimotor','Anterior Cingulate',...
                  'Right Sensorimotor','Posterior Parietal','Supplementary Motor Area','Anterior Parietal'};

% Identify EEG and ICA_struct files
ica_struct_files = dir('data_release_ICA_structs/*ICA_struct.mat');
eeg_files = dir('data_release_final/*.set');

% Load cluster info
load('cluster_data.mat');        

% Compute time-frequency information
cluster_ind = find(strcmp(cluster2plt,cluster_labels)) + 2;
n_comps = length(cluster(cluster_ind).comps);
erps_all = zeros(n_comps,length(evs),1,100,200); % sbj, event, chan, freq, time
for i=1:n_comps
    cls_set_i = cluster(cluster_ind).sets(1,i);
    comp = cluster(cluster_ind).comps(i);
    
    % Load EEG set
    EEG = pop_loadset('filename', eeg_files(cls_set_i).name,...
                      'filepath','data_release_final/');
    EEG = EEGFUNCS.rem_ev_cwlabels(EEG);
    EEG = EEGFUNCS.rem_ev_cwlabels_pulls(EEG);
    
    % Retain only EEG electrodes
    EEG = pop_select(EEG,'channel',icachansind);
    
    % Load ICA_STRUCT and project components from ICA
    load(['data_release_ICA_structs/' ica_struct_files(cls_set_i).name]);
    [EEG, comp_ind] = EEGFUNCS.apply_ica_to_eeg(EEG, ICA_STRUCT, comp,...
                                                icachansind);
    
    for k=1:length(evs)
        % Epoch data
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname','Epc_tmp','epochinfo','yes');

        % Retain only one component
        EEG_tmp = pop_select(EEG_tmp,'channel',comp_ind);

        % Compute ERSP's
        for p=1:size(erps_all,3)
            [erps_all(i,k,p,:,:),times,freqs] = EEGFUNCS.compute_ersps(EEG_tmp,p);
        end
    end
end

% Plot power spectra
specidx = find(times>=-500 & times<=1500);
PLTFUNCS.plot_eeg_specs(erps_all, specidx, freqs, cluster2plt)

% Combine and baseline-subtract ERSPs
baseidx = find(times>=-500 & times<=0);
erps_all_masked = EEGFUNCS.rem_eeg_ersp_baseline(erps_all, baseidx,...
                                                 times, alpha);

% Plot ERSPs
tlim = [-0.5 1.5]*1000; flim = [4 100]; clim = [-1 1]; nrows = 1;
PLTFUNCS.plot_ersps(erps_all, erps_all_masked, times, freqs,...
                    ev_labels, tlim, flim, clim, nrows)

% Compute ERSP onsets/offsets
tLims=[-0.2 0.5];
fLims_sync=[4 13]; fLims_desync=[8 30];
threshVal=1; %set values below this to zero
[t_onsets_all,fInds] = EEGFUNCS.compute_ersp_onsets(erps_all, tLims, fLims_sync,...
                                                    fLims_desync, threshVal,...
                                                    baseidx, times, freqs,...
                                                    n_comps);

% Plot ERSP onset timings
PLTFUNCS.plot_eeg_ersp_onsets(t_onsets_all, fInds)

toc;