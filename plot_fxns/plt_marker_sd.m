%% Marker SD plot (eNeuro Fig. 4)
% Start eeglab
eeglab_pth = '.../eeglab13_5_4b/'; % EEGLAB directory
root_pth = 'BIDS/'; % top-level data directory
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
num_conds = 4;
eeg_files = dir([root_pth '*/*/*/sub*_ses-01_task*.set']);
n_subjs = length(eeg_files); clear eeg_files;

% Compute marker SD for each condition
sds = zeros(n_subjs, num_conds, 2);
for i=1:n_subjs
    eeg_files = dir([root_pth 'sub-' num2str(i,'%03.f') '/*/*/sub-' ...
                     num2str(i,'%03.f') '*_ses*_task*.set']);
    
    for j=1:length(eeg_files)
        EEG_all(j) = pop_loadset('filename', eeg_files(j).name,...
                                 'filepath',eeg_files(j).folder);
    end
    EEG = pop_mergeset(EEG_all, 1:length(EEG_all));
    good_inds = EEGFUNCS.find_good_nonEEG_inds(EEG);
    
    % Get boundary events
    boundary_latencies2 = EEGFUNCS.get_boundary_evs(EEG);
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    marker_inds(1) = EEGFUNCS.find_chan_ind(EEG, 'Mocap_sacrum_mediolateral');
    marker_inds(2) = EEGFUNCS.find_chan_ind(EEG, 'Mocap_head_mediolateral');
    
    % Compute SD
    sds = EEGFUNCS.compute_marker_sd(marker_inds, num_conds, good_inds,...
                                     boundary_latencies2, sds, EEG, i);
end

% Plot result
PLTFUNCS.plot_marker_sd(marker_inds,sds,n_subjs)
