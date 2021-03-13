%% Plot EMG ERSP's (eNeuro Fig. 11)
% clear; close all; clc;
tic;
%Start eeglab
if ~exist('ALLCOM')
    curr_dir = pwd;
    run('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/eeglab.m');
    cd(curr_dir); close all;
    rmpath('/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/plugins/groupSIFT_03/');
end

%% Compute time-frequency information
eeg_files = dir('data_release_final/*.set');

n_parts = length(eeg_files);
evs = {'M_on_SVZ','M_on_WVZ','pull_Stn','pull_Wlk'};
ev_labels = {'Stand Rotate','Walk Rotate','Stand Pull','Walk Pull'};
ch_labels = {'Left Neck','Right Neck'};
ep_add_time = [-1 1]*0.6; % extra epoch times that will be cut off during ERSP computation
epochTimes=[-.5 1.5]+ep_add_time;
alpha = 0.05;
erps_all = zeros(n_parts,length(evs),2,100,200); % sbj, event, chan, freq, time
for i=1:n_parts
    % Load EEG file
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    EEG = rem_ev_cwlabels(EEG);
    EEG = rem_ev_cwlabels_pulls(EEG);
    
    for k=1:length(evs)
        % Epoch data
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname','Epc_tmp','epochinfo','yes');

        % Retain only neck EMG electrodes
        EEG_tmp = pop_select(EEG_tmp,'channel',[129 130]);

        % Compute ERSP's
        for p=1:size(erps_all,3)
            tlimits = [EEG_tmp.xmin, EEG_tmp.xmax]*1000;
            pointrange1 = round(max((tlimits(1)/1000-EEG_tmp.xmin)*EEG_tmp.srate, 1));
            pointrange2 = round(min((tlimits(2)/1000-EEG_tmp.xmin)*EEG_tmp.srate, EEG_tmp.pnts));
            pointrange = [pointrange1:pointrange2];
            tmpsig = EEG_tmp.data(p,pointrange,:);
            tmpsig = reshape(tmpsig, 1, length(pointrange)*size(EEG_tmp.data,3));
            [erps_all(i,k,p,:,:),~,~,times,freqs,~,~] = mod_newtimef(tmpsig, length(pointrange), [tlimits(1) tlimits(2)], EEG_tmp.srate,[3 0.9],...
                                                      'baseline',nan,'freqs',[3 256],'nfreqs',100,'freqscale','log',...
                                                      'plotersp','off','plotitc','off');
        end
    end
end

%% Plot specs
specidx = find(times>=-500 & times<=1500);
specs_all = zeros(size(squeeze(erps_all(:,:,:,:,1))));
for ch=1:size(erps_all,3)
    for ev=1:size(erps_all,2)
        specs_all(:,ev,ch,:) = squeeze(mean(erps_all(:,ev,ch,:,specidx),5));
    end
end

figure; cols = {'b','c','r','m'};
for ch=1:size(specs_all,3)
    subplot(1,size(specs_all,3),ch);
    hold on;
    for ev=1:size(specs_all,2)
        plot(freqs,squeeze(mean(specs_all(:,ev,ch,:),1)),cols{ev},'linewidth',1.5);
    end
    hold off;
    title(ch_labels{ch},'FontSize',16);
    if ch==1
        ylabel(sprintf('Power\n(dB)'),'Rotation',0,'fontsize',16,'fontweight','bold');
    end
    xlabel('Frequency (Hz)','Fontsize',16,'fontweight','bold');
    xlim([4,100]); ylim([35 65]);
    set(gca,'Fontsize',14,'fontweight','bold');
    set(gca,'box','off');
    set(gca,'Xscale','log');
    set(gca,'XMinorTick','Off')
    xticks([4,8,13,30,50,100]);
end

%% Combine and baseline-subtract ERSPs
baseidx = find(times>=-500 & times<=0);
erps_all_masked = zeros(size(squeeze(erps_all(1,:,:,:,:))));
for ch=1:size(erps_all,3)
    for ev=1:size(erps_all,2)
        curr_ersp = permute(squeeze(erps_all(:,ev,ch,:,:)),[2,3,1]);
        
        baseline = mean(curr_ersp(:,baseidx,:),2);
        curr_ersp = curr_ersp-repmat(baseline,1,length(times));
        
        pboot = bootstat(curr_ersp,'mean(arg1,3);','boottype','shuffle',...
                        'label','ERSP','bootside','both','naccu',200,...
                        'basevect',baseidx,'alpha',alpha,'dimaccu',2);         
        curr_ersp = mean(curr_ersp,3);
        curr_maskedersp = curr_ersp;
        curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
        erps_all_masked(ev,ch,:,:) = curr_maskedersp;
    end
end

%% Plot ERSPs
tlim = [-0.5 1.5]*1000; flim = [4 100]; clim = [-1 1];

figure; q = 0;
set(gcf,'Position',[1 1 1102 317]);
for ch=1:size(erps_all,3)
    for ev=1:size(erps_all,2)
        q = q + 1;
        subplot(2,4,q);
        if ev<3
            verts = [0 500];
        else
            verts = [0 1000];
        end
        tftopo(squeeze(erps_all_masked(ev,ch,:,:)),times,freqs,'vert',verts,... 
               'limits',[tlim flim clim],'logfreq','native');
        
        if ch == 1
            title(ev_labels{ev},'FontSize',16);
        else
            title('');
        end
        ylimits = ylim; %(gca);
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843]));%100-.05])); %4.009255933
        if ev==1
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency\n(Hz)'),'Rotation',0,'fontsize',16,'fontweight','bold');
        else
            set(gca,'YTickLabel',{'','','','','','',''});
            ylabel('');
        end
        if ch==1
            xlabel('');
            set(gca,'XTickLabel',{'','','',''});
        else
            xlabel('Time (sec)','Fontsize',16,'fontweight','bold');
            set(gca,'XTickLabel',{'-0.5','0','0.5','1'},'Fontsize',12);
        end
        ylhand = get(gca,'ylabel');
        %set(ylhand,'Rotation',0,'fontsize',16,'fontweight','bold');
        set(gca,'Fontsize',14,'fontweight','bold');
        set(gca,'box','off');
    end
end
toc;