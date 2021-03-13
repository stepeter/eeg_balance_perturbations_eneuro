%% Code to recreate marker SD plot (eNeuro Fig. 4)
% clear; close all; clc;
%Start eeglab
if ~exist('ALLCOM')
    curr_dir = pwd;
    run('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/eeglab.m');
    cd(curr_dir); close all;
    rmpath('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/plugins/groupSIFT_03/');
end
close all;

%% Compute marker SD
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
all_vals = zeros(4,2,n_sbjs,length(evs)); % plots x early/late x sbjs x condition
for i=1:n_sbjs
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    EEG = rem_ev_cwlabels(EEG);
    EEG = rem_ev_cwlabels_pulls(EEG);
    
    % Get boundary events
    boundary_latencies=[];
    for j=1:length(EEG.event)
        if strcmp(EEG.event(j).type,'boundary')
            boundary_latencies=[boundary_latencies EEG.event(j).latency];
        end
    end
    boundary_latencies2=[1 boundary_latencies EEG.pnts]; %int32();
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    for j=1:length(EEG.chanlocs)
        if strcmp(EEG.chanlocs(j).labels,'HEADx')
            marker_inds(2) = j;
        elseif strcmp(EEG.chanlocs(j).labels,'SACRx')
            marker_inds(1) = j;
        end
    end
    
    % Find EMG electrodes
    emg_inds = [];
    for j=1:length(EEG.chanlocs)
        if any(strcmp(EEG.chanlocs(j).labels,emg_marks))
            emg_inds = [emg_inds j];
        end
    end
    
    % Find load cell electrodes
    lcell_inds = [];
    for j=1:length(EEG.chanlocs)
        if any(strcmpi(EEG.chanlocs(j).labels,lcell_marks))
            lcell_inds = [lcell_inds j];
        end
    end
    
    % 20 Hz high-pass filter emg data and full-wave rectify
    lo_cutoff = 20;
    [b,a]=butter(4,lo_cutoff/(EEG.srate/2),'high');
    for k=1:length(emg_inds)
        EEG.data(emg_inds(k),:) = abs(filtfilt(b,a,double(EEG.data(emg_inds(k),:))));
    end
    
    for k=1:num_conds
        % Create early/late EEG variables
        EEG_tmp = pop_select(EEG,'time',[boundary_latencies2(k) boundary_latencies2(k+1)]/EEG.srate);
        if i==25
            EEG_tmp2 = pop_select(EEG_tmp,'channel',emg_inds(1));
        else
            EEG_tmp2 = pop_select(EEG_tmp,'channel',marker_inds(1));
        end
        good_inds_tmp = find_good_nonEEG_inds(EEG_tmp2,1);
        EEG_early = pop_select(EEG_tmp,'time',[good_inds_tmp(1) good_inds_tmp(1)+EEG.srate*60]/EEG.srate);
        EEG_late = pop_select(EEG_tmp,'time',[good_inds_tmp(end)-EEG.srate*60 good_inds_tmp(end)]/EEG.srate);

        % Compute mocap SD
        plt_num = 2;
        for j=1:length(marker_inds)
            all_vals(plt_num+j,1,i,k) = std(EEG_early.data(marker_inds(j),:));
            all_vals(plt_num+j,2,i,k) = std(EEG_late.data(marker_inds(j),:));
        end
        
        % Epoch data by event
        EEG_ep_early = pop_epoch(EEG_early,{evs{k}},epochTimes,...
                                 'newname',['Epc_' num2str(i)],'epochinfo','yes');
%         good_inds = find_good_nonEEG_inds(EEG_ep_early);
%         EEG_ep_early = pop_select(EEG_ep_early,'trial',good_inds);

        EEG_ep_late = pop_epoch(EEG_late,{evs{k}},epochTimes,...
                                 'newname',['Epc_' num2str(i)],'epochinfo','yes');
%         good_inds = find_good_nonEEG_inds(EEG_ep_late);
%         EEG_ep_late = pop_select(EEG_ep_late,'trial',good_inds);
        
        % Find max pull force for pull perturbations
        if k<3
            % Compute peak pull force (max of 0.5 second interval after event 
            % onset)
            plt_num = 1;
            sel_times = find((EEG_ep_early.times>=t_lims_ave(1)) & (EEG_ep_early.times<=t_lims_ave(2)));
            max_vals_early = squeeze(max(max(EEG_ep_early.data(lcell_inds,sel_times,:),[],2),[],1));
            all_vals(plt_num,1,i,k) = mean(max_vals_early);
            max_vals_late = squeeze(max(max(EEG_ep_late.data(lcell_inds,sel_times,:),[],2),[],1));
            all_vals(plt_num,2,i,k) = mean(max_vals_late);
        end

        % Baseline subtraction
        fs = EEG_ep_early.srate;
        for p=1:size(EEG_ep_early.data,3)
            EEG_ep_early.data(:,:,p) = EEG_ep_early.data(:,:,p) - repmat(mean(EEG_ep_early.data(:,1:(fs*0.5),p),2),1,size(EEG_ep_early.data,2));
        end
        for p=1:size(EEG_ep_late.data,3)
            EEG_ep_late.data(:,:,p) = EEG_ep_late.data(:,:,p) - repmat(mean(EEG_ep_late.data(:,1:(fs*0.5),p),2),1,size(EEG_ep_early.data,2));
        end
        
        % Compute peak EMG (average between t_lims_ave, after baseline
        % subtraction)
        plt_num = 2; sum_emg_early = []; sum_emg_late = [];
        for j=1:length(emg_inds)
            sel_times = find((EEG_ep_early.times>=t_lims_ave(1)) & (EEG_ep_early.times<=t_lims_ave(2)));
            sum_emg_early = [sum_emg_early mean(mean(EEG_ep_early.data(emg_inds(j),sel_times,:),2),3)/baseEMGPeakVals(i,j)];
            sum_emg_late = [sum_emg_late mean(mean(EEG_ep_late.data(emg_inds(j),sel_times,:),2),3)/baseEMGPeakVals(i,j)];
        end
        all_vals(plt_num,1,i,k) = 100*mean(sum_emg_early);
        all_vals(plt_num,2,i,k) = 100*mean(sum_emg_late);
    end
end

%% Plot result
figure;
titles = {'Pull force','Peak EMG','Sacrum SD','Head SD'};
ylab = {'Force (N)','Normalized EMG (%)','Marker SD (mm)','Marker SD (mm)'};
cols = {'r','r','m','m','b','b','c','c'};
ylims = {[0 20],[-2 25],[0 80],[0 80]};
           
t_inds = [2 3 5 6 8 9 11 12];
for k=1:size(all_vals,1)
    subplot(2,2,k);
    q = 0; hold on;
    for i=1:size(all_vals,4)
        for j=1:size(all_vals,2)
            q = q + 1;
            if (k>1) || (i<3)
                bar(t_inds(q),mean(all_vals(k,j,:,i)),cols{q});
                errorbar(t_inds(q),mean(all_vals(k,j,:,i)),std(all_vals(k,j,:,i))/sqrt(size(all_vals,3)),'k');
            end
        end
    end
    hold off;
    ylim(ylims{k});
    title(titles{k},'fontsize',14);
    xticks([2.5,5.5,8.5,11.5]);
    xlim([.7,13.3]);
    xticklabels({'St Pull','Wlk Pull','St Rotate','Wlk Rotate'});
    ylabel(ylab{k});
end
