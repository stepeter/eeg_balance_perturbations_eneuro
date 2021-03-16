%% Early/late behavior plots (eNeuro Fig. 6)
% Start eeglab
eeglab_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/';
if ~exist('ALLCOM')
    PLTFUNCS.start_eeglab(eeglab_pth)
end

% Set parameters
t_lims_ave = [0.2 0.5]*1000;
emg_marks = {'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};
lcell_marks = {'LoadCell_l','LoadCell_r'};
evs = {'pull_Stn','pull_Wlk','M_on_SVZ','M_on_WVZ'};
evs_lcell = {{'L_pull_Stn','R_pull_Stn'},{'L_pull_Wlk','R_pull_Wlk'}};
num_conds = 4;
eeg_files = dir('data_release_final/*.set');
load('baseEMGPeakVals.mat');
n_sbjs = length(eeg_files);
epochTimes=[-.5 1.5];
hi_cutoff = 20; % Hz

% Compute behavioral measures
all_vals = zeros(4,2,n_sbjs,length(evs)); % plots x early/late x sbjs x condition
for i=1:n_sbjs
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    EEG = EEGFUNCS.rem_ev_cwlabels(EEG);
    EEG = EEGFUNCS.rem_ev_cwlabels_pulls(EEG);
    
    % Get boundary events
    boundary_latencies2 = EEGFUNCS.get_boundary_evs(EEG);
    boundary_latencies2 = double(boundary_latencies2);
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    marker_inds(1) = EEGFUNCS.find_chan_ind(EEG, 'SACRx');
    marker_inds(2) = EEGFUNCS.find_chan_ind(EEG, 'HEADx');
    
    % Find EMG electrodes
    emg_inds = [];
    for j=1:length(emg_marks)
        emg_inds(j) = EEGFUNCS.find_chan_ind(EEG, emg_marks(j));
    end
    
    % Find load cell electrodes
    lcell_inds = [];
    for j=1:length(lcell_marks)
        lcell_inds(j) = EEGFUNCS.find_chan_ind(EEG, lcell_marks(j));
    end
    
    % 20 Hz high-pass filter emg data and full-wave rectify
    EEG = EEGFUNCS.filt_data(EEG, emg_inds, hi_cutoff, 'high');
    EEG.data(emg_inds,:) = abs(EEG.data(emg_inds,:));
    
    for k=1:num_conds
        % Create early/late EEG variables
        [EEG_early,EEG_late] = EEGFUNCS.early_late_eeg(EEG,boundary_latencies2,...
                                                       emg_inds,marker_inds,...
                                                       k, i);
        
        % Compute mocap SD
        plt_num = 2;
        for j=1:length(marker_inds)
            all_vals(plt_num+j,1,i,k) = std(EEG_early.data(marker_inds(j),:));
            all_vals(plt_num+j,2,i,k) = std(EEG_late.data(marker_inds(j),:));
        end
        
        % Epoch data by event
        EEG_ep_early = pop_epoch(EEG_early,{evs{k}},epochTimes,...
                                 'newname',['Epc_' num2str(i)],'epochinfo','yes');
        EEG_ep_late = pop_epoch(EEG_late,{evs{k}},epochTimes,...
                                 'newname',['Epc_' num2str(i)],'epochinfo','yes');
        
        % Find max pull force for pull perturbations
        if k<3
            % Compute peak pull force (max of 0.5 second interval after event 
            % onset)
            plt_num = 1;
            sel_times = find((EEG_ep_early.times>=t_lims_ave(1)) & ...
                             (EEG_ep_early.times<=t_lims_ave(2)));
            max_vals_early = squeeze(max(max(EEG_ep_early.data(lcell_inds,...
                                             sel_times,:),[],2),[],1));
            all_vals(plt_num,1,i,k) = mean(max_vals_early);
            max_vals_late = squeeze(max(max(EEG_ep_late.data(lcell_inds,...
                                            sel_times,:),[],2),[],1));
            all_vals(plt_num,2,i,k) = mean(max_vals_late);
        end

        % Baseline subtraction
        fs = EEG_ep_early.srate;
        for p=1:size(EEG_ep_early.data,3)
            EEG_ep_early.data(:,:,p) = EEG_ep_early.data(:,:,p) - ...
                repmat(mean(EEG_ep_early.data(:,1:(fs*0.5),p),2),1,size(EEG_ep_early.data,2));
        end
        for p=1:size(EEG_ep_late.data,3)
            EEG_ep_late.data(:,:,p) = EEG_ep_late.data(:,:,p) - ...
                repmat(mean(EEG_ep_late.data(:,1:(fs*0.5),p),2),1,size(EEG_ep_early.data,2));
        end
        
        % Compute peak EMG (average between t_lims_ave, after baseline
        % subtraction)
        plt_num = 2; sum_emg_early = []; sum_emg_late = [];
        for j=1:length(emg_inds)
            sel_times = find((EEG_ep_early.times>=t_lims_ave(1)) & ...
                             (EEG_ep_early.times<=t_lims_ave(2)));
            sum_emg_early = [sum_emg_early mean(mean(EEG_ep_early.data(emg_inds(j),...
                                sel_times,:),2),3)/baseEMGPeakVals(i,j)];
            sum_emg_late = [sum_emg_late mean(mean(EEG_ep_late.data(emg_inds(j),...
                                sel_times,:),2),3)/baseEMGPeakVals(i,j)];
        end
        all_vals(plt_num,1,i,k) = 100*mean(sum_emg_early);
        all_vals(plt_num,2,i,k) = 100*mean(sum_emg_late);
    end
end

% Plot result
PLTFUNCS.plot_early_late_behavior(all_vals)
