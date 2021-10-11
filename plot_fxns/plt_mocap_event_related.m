%% Marker event-related activity plot (eNeuro Fig. 3)
% Start eeglab
eeglab_pth = '.../eeglab13_5_4b/'; % EEGLAB directory
root_pth = 'BIDS/'; % top-level data directory
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
num_conds = 4; n_markers = 2; lo_cutoff = 6;
epochTimes=[-.5 1.5]; fs =256; % Hz
eeg_files = dir([root_pth '*/*/*/sub*_ses-01_task*.set']);
n_subjs = length(eeg_files); clear eeg_files;

% Compute marker SD
ev_traces = zeros(n_markers,n_subjs,num_conds,fs*(sum(abs(epochTimes))));
evs = {'pull_Stn','pull_Wlk','M_on_SVZ','M_on_WVZ'};
for i=1:n_subjs
    eeg_files = dir([root_pth 'sub-' num2str(i,'%03.f') '/*/*/sub-' ...
                     num2str(i,'%03.f') '*_ses*_task*.set']);
    
    for j=1:length(eeg_files)
        EEG_all(j) = pop_loadset('filename', eeg_files(j).name,...
                                 'filepath',eeg_files(j).folder);
    end
    EEG = pop_mergeset(EEG_all, 1:length(EEG_all));
    EEG = EEGFUNCS.rem_ev_cwlabels(EEG);
    EEG = EEGFUNCS.rem_ev_cwlabels_pulls(EEG);
    
    % Get boundary events
    boundary_latencies2 = EEGFUNCS.get_boundary_evs(EEG);
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    marker_inds(1) = EEGFUNCS.find_chan_ind(EEG, 'Mocap_sacrum_mediolateral');
    marker_inds(2) = EEGFUNCS.find_chan_ind(EEG, 'Mocap_head_mediolateral');
    
    % 6 Hz low-pass filter mocap data
    EEG = EEGFUNCS.filt_data(EEG, marker_inds, lo_cutoff, 'low');
    
    % Compute event-related data
    for k=1:num_conds
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname',['Epc_' num2str(i)],'epochinfo','yes');
        good_inds = EEGFUNCS.find_good_nonEEG_inds(EEG_tmp);
        EEG_tmp = pop_select(EEG_tmp,'trial',good_inds);
        EEG_tmp.data = abs(EEG_tmp.data);
        
        for j=1:length(marker_inds)
            % Save epoched traces to common file
            dat_tmp = EEG_tmp.data(marker_inds(j),:,:);
            for p=1:size(EEG_tmp.data,3)
                dat_tmp(:,:,p) = dat_tmp(:,:,p) - mean(dat_tmp(:,1:(fs*0.5),p));
            end
            ev_traces(j,i,k,:) = mean(dat_tmp,3);
            
        end
    end
end

% Plot mocap traces
PLTFUNCS.plot_mocap_evs(ev_traces,marker_inds,EEG_tmp,n_subjs)