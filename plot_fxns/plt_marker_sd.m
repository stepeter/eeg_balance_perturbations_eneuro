%% Marker SD plot (eNeuro Fig. 4)
% Start eeglab
eeglab_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/';
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
num_conds = 4;
eeg_files = dir('data_release_final/*.set');
n_subjs = length(eeg_files);

% Compute marker SD for each condition
sds = zeros(n_subjs, num_conds, 2);
for i=1:n_subjs
    EEG = pop_loadset('filename', eeg_files(i).name,...
                      'filepath','data_release_final/');
    good_inds = EEGFUNCS.find_good_nonEEG_inds(EEG);
    
    % Get boundary events
    boundary_latencies2 = EEGFUNCS.get_boundary_evs(EEG);
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    marker_inds(1) = EEGFUNCS.find_chan_ind(EEG, 'SACRx');
    marker_inds(2) = EEGFUNCS.find_chan_ind(EEG, 'HEADx');
    
    % Compute SD
    sds = EEGFUNCS.compute_marker_sd(marker_inds, num_conds, good_inds,...
                                     boundary_latencies2, sds, EEG, i);
end

% Plot result
PLTFUNCS.plot_marker_sd(marker_inds,sds,n_subjs)
