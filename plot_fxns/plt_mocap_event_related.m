%% Code to recreate marker event-related activity plot (eNeuro Fig. 3)
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
num_conds = 4; n_markers = 2;
epochTimes=[-.5 1.5]; fs =256; % Hz
eeg_files = dir('data_release_final/*.set');

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
    
    % Find SACR and HEAD electrodes
    marker_inds = [];
    for j=1:length(EEG.chanlocs)
        if strcmp(EEG.chanlocs(j).labels,'HEADx')
            marker_inds(2) = j;
        elseif strcmp(EEG.chanlocs(j).labels,'SACRx')
            marker_inds(1) = j;
        end
    end
    
    % 6 Hz low-pass filter mocap data
    lo_cutoff = 6;
    [b,a]=butter(4,lo_cutoff/(EEG.srate/2),'low');
    for k=1:length(marker_inds)
        EEG.data(marker_inds(k),:) = filtfilt(b,a,double(EEG.data(marker_inds(k),:)));
    end
    
    % Compute SD (just remove any zeros for quick computation)
    for k=1:num_conds
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname',['Epc_' num2str(i)],'epochinfo','yes');
        good_inds = find_good_nonEEG_inds(EEG_tmp);
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

%% Plot mocap traces
figure;
title_labs = {'Sacrum','Head'};
cols = {'r','m','b','c'};
for j=1:length(marker_inds)
    subplot(2,1,j); hold on;
    for k=1:num_conds
        errorbar(EEG_tmp.times/1000,mean(squeeze(ev_traces(j,:,k,:)),1),...
                 squeeze(std(ev_traces(j,:,k,:),0,2))/sqrt(length(eeg_files)),cols{k});
    end
    hold off;
    xlim([-0.5 1.5]);
    title(title_labs{j});
    legend('Stand Pull','Walk Pull',...
           'Stand Rotate','Walk Rotate');
    ylabel(sprintf('Marker\nSD (mm)'));
    ylp = get(gca,'ylabel');
    ext = get(ylp,'Extent');
    set(ylp, 'Rotation',0,...
        'Position',get(ylp,'Position')-[ext(3) 0 0]);
end