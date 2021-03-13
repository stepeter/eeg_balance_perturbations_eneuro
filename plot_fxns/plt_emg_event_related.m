%% Code to recreate EMG event-related activity plot (eNeuro Fig. 5)
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
num_conds = 4;
epochTimes=[-.5 1.5]; fs =256; % Hz
eeg_files = dir('data_release_final/*.set');
emg_marks = {'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};

n_markers = length(emg_marks);
ev_traces = zeros(n_markers,length(eeg_files),num_conds,fs*(sum(abs(epochTimes))));
evs = {'pull_Stn','pull_Wlk','M_on_SVZ','M_on_WVZ'};
for i=1:length(eeg_files)
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
    boundary_latencies2=int32([1 boundary_latencies EEG.pnts]);
    
    % Find EMG electrodes
    marker_inds = [];
    for j=1:length(EEG.chanlocs)
        if any(strcmp(EEG.chanlocs(j).labels,emg_marks))
            marker_inds = [marker_inds j];
        end
    end
    
    % 20 Hz high-pass filter emg data
    lo_cutoff = 20;
    [b,a]=butter(4,lo_cutoff/(EEG.srate/2),'high');
    for k=1:length(marker_inds)
        EEG.data(marker_inds(k),:) = filtfilt(b,a,double(EEG.data(marker_inds(k),:)));
    end
    
    % Compute SD (just remove any zeros for quick computation)
    for k=1:num_conds
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname',['Epc_' num2str(i)],'epochinfo','yes');
        good_inds = find_good_nonEEG_inds(EEG_tmp);
        EEG_tmp = pop_select(EEG_tmp,'trial',good_inds);
        EEG_tmp.data = abs(EEG_tmp.data); % full-wave rectification
        
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

%% Divide by baseline peak values (obtained baseline walking condition)
load('baseEMGPeakVals.mat');
ev_traces_norm = ev_traces;
for i=1:size(ev_traces,1)
    for j=1:size(ev_traces,2)
        ev_traces_norm(i,j,:,:) = ev_traces(i,j,:,:)/baseEMGPeakVals(j,i)*100;
    end
end

%% Plot mocap traces
figure;
cols = {'r','m','b','c'};
plt_order = [1,3,5,7,2,4,6,8];
for j=1:length(marker_inds)
    subplot(4,2,plt_order(j)); hold on;
    for k=1:num_conds
        errorbar(EEG_tmp.times/1000,mean(squeeze(ev_traces_norm(j,:,k,:)),1),...
                 squeeze(std(ev_traces_norm(j,:,k,:),0,2))/sqrt(length(eeg_files)),cols{k});
    end
%     plot([EEG_tmp.times(1)/1000 EEG_tmp.times(end)/1000],[100 100],':g','linewidth',1.5);
    hold off;
    xlim([-0.5 1.5]);
    ylim([0,50]);
    title(emg_marks{j});
    if j==1
        legend('Stand Pull','Walk Pull',...
               'Stand Rotate','Walk Rotate');
    end
    if j<5
        ylabel([sprintf('Normalized\nEMG ') '(%)']);
    end
    ylp = get(gca,'ylabel');
    ext = get(ylp,'Extent');
    set(ylp, 'Rotation',0,...
        'Position',get(ylp,'Position')-[ext(3) 0 0]);
end