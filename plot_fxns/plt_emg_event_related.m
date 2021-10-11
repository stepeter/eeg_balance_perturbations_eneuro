%% EMG event-related activity plot (eNeuro Fig. 5)
% Start eeglab
eeglab_pth = '.../eeglab13_5_4b/'; % EEGLAB directory
root_pth = 'BIDS/'; % top-level data directory
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
num_conds = 4; epochTimes=[-.5 1.5]; 
hi_cutoff = 20; fs =256; % Hz
eeg_files = dir([root_pth '*/*/*/sub*_ses-01_task*.set']);
emg_marks = {'L_leg_emg_tibialis_anterior','L_leg_emg_soleus',...
             'L_leg_emg_medial_gastrocnemius','L_leg_emg_peroneus_longus',...
             'R_leg_emg_tibialis_anterior','R_leg_emg_soleus',...
             'R_leg_emg_medial_gastrocnemius','R_leg_emg_peroneus_longus'};

% Compute marker SD
n_markers = length(emg_marks);
n_subjs = length(eeg_files); clear eeg_files;
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
emg_marks_plt = {'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};
PLTFUNCS.plot_emg_evs(ev_traces_norm,EEG_tmp,emg_marks_plt,n_subjs)