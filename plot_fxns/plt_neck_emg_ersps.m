%% Neck EMG ERSP's (eNeuro Fig. 11)
tic;
% Start eeglab
eeglab_pth = '.../eeglab13_5_4b/'; % EEGLAB directory
root_pth = 'BIDS/'; % top-level data directory
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
eeg_files = dir([root_pth '*/*/*/sub*_ses-01_task*.set']);
n_parts = length(eeg_files); clear eeg_files;
evs = {'M_on_SVZ','M_on_WVZ','pull_Stn','pull_Wlk'};
ev_labels = {'Stand Rotate','Walk Rotate','Stand Pull','Walk Pull'};
ch_labels = {'L_neck_emg','R_neck_emg'};
ep_add_time = [-1 1]*0.6; % extra epoch times that will be cut off during ERSP computation
epochTimes=[-.5 1.5]+ep_add_time;
alpha = 0.05;

% Compute time-frequency information
erps_all = zeros(n_parts,length(evs),2,100,200); % sbj, event, chan, freq, time
for i=1:n_parts
    % Load EEG files
    eeg_files = dir([root_pth 'sub-' num2str(i,'%03.f') '/*/*/sub-' ...
                     num2str(i,'%03.f') '*_ses*_task*.set']);
    
    for j=1:length(eeg_files)
        EEG_all(j) = pop_loadset('filename', eeg_files(j).name,...
                                 'filepath',eeg_files(j).folder);
    end
    EEG = pop_mergeset(EEG_all, 1:length(EEG_all));
    EEG = EEGFUNCS.rem_ev_cwlabels(EEG);
    EEG = EEGFUNCS.rem_ev_cwlabels_pulls(EEG);
    
    for k=1:length(evs)
        % Epoch data
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname','Epc_tmp','epochinfo','yes');

        % Retain only neck EMG electrodes
        EEG_tmp = pop_select(EEG_tmp,'channel',[129 130]);

        % Compute ERSP's
        for p=1:size(erps_all,3)
            [erps_all(i,k,p,:,:),times,freqs] = EEGFUNCS.compute_ersps(EEG_tmp,p);
        end
    end
end

% Plot power spectra
specidx = find(times>=-500 & times<=1500);
PLTFUNCS.plot_neck_specs(erps_all, specidx, freqs, ch_labels)

% Combine and baseline-subtract ERSPs
baseidx = find(times>=-500 & times<=0);
erps_all_masked = EEGFUNCS.rem_ersp_baseline(erps_all, baseidx, times, alpha);

% Plot ERSPs
tlim = [-0.5 1.5]*1000; flim = [4 100]; clim = [-1 1]; nrows = 2;
PLTFUNCS.plot_ersps(erps_all, erps_all_masked, times, freqs, ev_labels,...
                    tlim, flim, clim, nrows)
toc;