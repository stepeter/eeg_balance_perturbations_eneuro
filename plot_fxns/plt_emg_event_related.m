%% EMG event-related activity plot (eNeuro Fig. 5)
% Start eeglab
eeglab_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/';
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
num_conds = 4; epochTimes=[-.5 1.5]; 
hi_cutoff = 20; fs =256; % Hz
eeg_files = dir('data_release_final/*.set');
emg_marks = {'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};

% Compute marker SD
n_markers = length(emg_marks); n_subjs = length(eeg_files);
ev_traces = zeros(n_markers,n_subjs,num_conds,fs*(sum(abs(epochTimes))));
evs = {'pull_Stn','pull_Wlk','M_on_SVZ','M_on_WVZ'};
for i=1:n_subjs
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    EEG = EEGFUNCS.rem_ev_cwlabels(EEG);
    EEG = EEGFUNCS.rem_ev_cwlabels_pulls(EEG);
    
    % Get boundary events
    boundary_latencies2 = EEGFUNCS.get_boundary_evs(EEG);
    
    % Find EMG electrodes
    marker_inds = [];
    for j=1:length(emg_marks)
        marker_inds(j) = EEGFUNCS.find_chan_ind(EEG, emg_marks(j));
    end
    
    % 20 Hz high-pass filter emg data
    EEG = EEGFUNCS.filt_data(EEG, marker_inds, hi_cutoff, 'high');
    
    % Compute event-related data
    for k=1:num_conds
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname',['Epc_' num2str(i)],'epochinfo','yes');
        good_inds = EEGFUNCS.find_good_nonEEG_inds(EEG_tmp);
        EEG_tmp = pop_select(EEG_tmp,'trial',good_inds);
        EEG_tmp.data = abs(EEG_tmp.data); % full-wave rectification
        
        % Save epoched traces to common file
        for j=1:length(marker_inds)
            dat_tmp = EEG_tmp.data(marker_inds(j),:,:);
            for p=1:size(EEG_tmp.data,3)
                dat_tmp(:,:,p) = dat_tmp(:,:,p) - mean(dat_tmp(:,1:(fs*0.5),p));
            end
            ev_traces(j,i,k,:) = mean(dat_tmp,3);
            
        end
    end
end

% Divide by baseline peak values (obtained baseline walking condition)
load('baseEMGPeakVals.mat');
ev_traces_norm = ev_traces;
for i=1:size(ev_traces,1)
    for j=1:size(ev_traces,2)
        ev_traces_norm(i,j,:,:) = ev_traces(i,j,:,:)/baseEMGPeakVals(j,i)*100;
    end
end

% Plot EMG traces
PLTFUNCS.plot_emg_evs(ev_traces_norm,EEG_tmp,emg_marks,n_subjs)