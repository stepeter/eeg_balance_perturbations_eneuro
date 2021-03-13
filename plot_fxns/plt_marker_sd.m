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
num_conds = 4;
eeg_files = dir('data_release_final/*.set');

sds = zeros(length(eeg_files),num_conds,2);
for i=1:length(eeg_files)
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    good_inds = find_good_nonEEG_inds(EEG);
    
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
    
    % Compute SD (just remove any zeros for quick computation)
    for j=1:length(marker_inds)
        for k=1:num_conds
            inds_use = good_inds((good_inds>=boundary_latencies2(k)) & (good_inds<=boundary_latencies2(k+1)));
%             dat_val = EEG.data(marker_inds(j),boundary_latencies2(k):boundary_latencies2(k+1));
            dat_val = EEG.data(marker_inds(j),inds_use);
%             dat_val = dat_val(dat_val~=0);
            dat_val = dat_val(~isnan(dat_val));
            if ~isempty(dat_val)
                sds(i,k,j) = std(dat_val);
            end
        end
    end
end

%% Plot result
figure;
title_labs = {'Sacrum','Head'};
cols = {'r','m','b','c'};
for j=1:length(marker_inds)
    subplot(2,1,j); hold on;
    for k=1:num_conds
        bar(k,mean(sds(:,k,j),1),cols{k});
        errorbar(k,mean(sds(:,k,j),1),std(sds(:,k,j),1)/sqrt(length(eeg_files)),'k');
    end
    hold off;
    xlim([0.5 4.5]);
    xticks([1,2,3,4]);
    title(title_labs{j});
    xticklabels({'Stand Pull','Walk Pull',...
                 'Stand Rotate','Walk Rotate'});
    ylabel(sprintf('Marker\nSD (mm)'));
    ylp = get(gca,'ylabel');
    ext = get(ylp,'Extent');
    set(ylp, 'Rotation',0,...
        'Position',get(ylp,'Position')-[ext(3) 0 0]);
end
