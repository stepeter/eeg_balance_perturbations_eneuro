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
cluster2plt = 'Supplementary Motor Area'; % select cluster to plot ERSPs from

eeg_files = dir('data_release_final/*.set');
evs = {'M_on_SVZ','M_on_WVZ','pull_Stn','pull_Wlk'};
ev_labels = {'Stand Rotate','Walk Rotate','Stand Pull','Walk Pull'};
ch_labels = {'Left Neck','Right Neck'};
ep_add_time = [-1 1]*0.6; % extra epoch times that will be cut off during ERSP computation
epochTimes=[-.5 1.5]+ep_add_time;
alpha = 0.05;
cluster_labels = {'Left Occipital','Right Occipital','Left Sensorimotor','Anterior Cingulate',...
                  'Right Sensorimotor','Posterior Parietal','Supplementary Motor Area','Anterior Parietal'};

              
dipfit_pth = '/Users/stepeter/Documents/BruntonLab/eeglab13_5_4b/plugins/dipfit2.3/'; % path to dipfit plugin
hdmfile = [dipfit_pth 'standard_BEM/standard_vol.mat'];
mrifile = [dipfit_pth 'standard_BEM/standard_mri.mat'];
ica_struct_files = dir('data_release_ICA_structs/*ICA_struct.mat');
eeg_files = dir('data_release_final/*.set');

% Load cluster info
load('cluster_data.mat');        

cluster_ind = find(strcmp(cluster2plt,cluster_labels)) + 2;
n_comps = length(cluster(cluster_ind).comps);
erps_all = zeros(n_comps,length(evs),1,100,200); % sbj, event, chan, freq, time
for i=1:n_comps
    cls_set_i = cluster(cluster_ind).sets(1,i);
    comp = cluster(cluster_ind).comps(i);
    
    % Load EEG set
    EEG = pop_loadset('filename', eeg_files(cls_set_i).name, 'filepath','data_release_final/');
    EEG = rem_ev_cwlabels(EEG);
    EEG = rem_ev_cwlabels_pulls(EEG);
    
    % Retain only EEG electrodes
    EEG = pop_select(EEG,'channel',1:128);
    
    % Load ICA_STRUCT and project components from ICA
    load(['data_release_ICA_structs/' ica_struct_files(cls_set_i).name]);
    good_comps = ICA_STRUCT.good_comps.brain;
    comp_ind = good_comps(comp);
    EEG.icaweights = ICA_STRUCT.weights;
    EEG.icasphere = ICA_STRUCT.sphere;
    EEG.icachansind = 1:128;
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); %compute ICA activations
    EEG.data(1:size(EEG.icaact,1),:) = EEG.icaact;
    
    for k=1:length(evs)
        % Epoch data
        EEG_tmp = pop_epoch(EEG,{evs{k}},epochTimes,...
                            'newname','Epc_tmp','epochinfo','yes');

        % Retain only one component
        EEG_tmp = pop_select(EEG_tmp,'channel',comp_ind);

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
ersp_size = size(erps_all);
specs_all = zeros(ersp_size(1:4));
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
    title(cluster2plt,'FontSize',16);
    if ch==1
        ylabel(sprintf('Power\n(dB)'),'Rotation',0,'fontsize',16,'fontweight','bold');
    end
    xlabel('Frequency (Hz)','Fontsize',16,'fontweight','bold');
    xlim([4,100]); ylim([30 55]);
    set(gca,'Fontsize',14,'fontweight','bold');
    set(gca,'box','off');
    set(gca,'Xscale','log');
    set(gca,'XMinorTick','Off')
    xticks([4,8,13,30,50,100]);
end

%% Combine and baseline-subtract ERSPs
baseidx = find(times>=-500 & times<=0);
ersp_size = size(erps_all);
erps_all_masked = zeros(ersp_size(2:end));
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
set(gcf,'Position',[1 1 1102 250]);
for ch=1:size(erps_all,3)
    for ev=1:size(erps_all,2)
        q = q + 1;
        subplot(1,4,q);
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
        xlabel('Time (sec)','Fontsize',16,'fontweight','bold');
        set(gca,'XTickLabel',{'-0.5','0','0.5','1'},'Fontsize',12);
        ylhand = get(gca,'ylabel');
        %set(ylhand,'Rotation',0,'fontsize',16,'fontweight','bold');
        set(gca,'Fontsize',14,'fontweight','bold');
        set(gca,'box','off');
    end
end

%% Compute ERSP onsets/offsets
tLims=[-0.2 0.5];
fLims_sync=[4 13]; fLims_desync=[8 30];
threshVal=1; %set values below this to zero

% Baseline subtract ERSP again
erps_all_base = zeros(size(erps_all));
for ch=1:size(erps_all,3)
    for ev=1:size(erps_all,2)
        curr_ersp = permute(squeeze(erps_all(:,ev,ch,:,:)),[2,3,1]);
        
        baseline = mean(curr_ersp(:,baseidx,:),2);
        curr_ersp = curr_ersp-repmat(baseline,1,length(times));
        erps_all_base(:,ev,ch,:,:) = permute(curr_ersp,[3,1,2]);
    end
end

tInds=find((tLims(1)*1000)<=times & (tLims(2)*1000)>=times);
fInds{1}=find(fLims_sync(1)<=freqs & fLims_sync(2)>=freqs);
fInds{2}=find(fLims_desync(1)<=freqs & fLims_desync(2)>=freqs);
n_conds = size(erps_all_base,2);
t_onsets_all = zeros(length(fInds),n_conds,n_comps);
for cond=1:n_conds
    for plt_num=1:length(fInds)
        cdata=squeeze(erps_all_base(:,cond,1,fInds{plt_num},tInds));
        
%         maskedersp = squeeze(erps_all_masked(cond,1,fInds{plt_num},tInds));
%         for dip=1:n_comps
%             tmp_val = cdata(dip,:,:);
%             tmp_val(maskedersp==0)=0;
%             cdata(dip,:,:) = tmp_val;
%         end
        
        if plt_num==1
            cdata(cdata<0)=0;
        else
            cdata(cdata>0)=0;
            cdata=abs(cdata);
        end
        
        cdata(cdata>=threshVal)=1; cdata(cdata<threshVal)=0;
        for dip=1:n_comps
            dip_data = squeeze(cdata(dip,:,:));
            if all(dip_data==0)
                timeVal=nan;
            else
                % Search ERSP 
                CC = bwconncomp(dip_data);
                num_pixels = cellfun(@numel,CC.PixelIdxList); % number of pixels in each region
                pix_list = regionprops(dip_data,'PixelIdxList','PixelList'); % x y pixel coordinates
                
                % Find largest connected pixels
                sizes=[];
                for i = 1:length(num_pixels)
                    sizes(i) = length(CC.PixelIdxList{i});
                end
                [largest,largest_loc] = max(sizes);
                pix_locs = cell2mat(CC.PixelIdxList(largest_loc)); % pixel locations
                
                % Find x y coordinates of largest region
                pix_cells=[];
                for i = 1:length(pix_locs)
                    pix_cells(i,1) = find(pix_list.PixelIdxList == pix_locs(i));
                end 
                pixel_xy = pix_list.PixelList(pix_cells,1:2);
                onset = min(pixel_xy(:,1)); % earliest frame of largest region
                t_onsets_all(plt_num,cond,dip)=times(tInds(onset)); %in msec
            end
        end
    end
end

%% Plot ERSP onset timings
figure;
titles = {'Synchronization','Desynchronization'};
cols = {'r','m','b','c'};
ylims = {[0 20],[-2 25],[0 80],[0 80]};
           
t_inds = [1 2 3 4];
plt_ord = [3 4 1 2];
for k=1:length(fInds)
    subplot(2,1,k); hold on;
    for i=1:n_conds
        bar(t_inds(i),mean(t_onsets_all(k,plt_ord(i),:)),cols{i});
        errorbar(t_inds(i),mean(t_onsets_all(k,plt_ord(i),:)),std(t_onsets_all(k,plt_ord(i),:))/sqrt(size(t_onsets_all,3)),'k');
    end
    hold off;
    ylim([0,300]);
    ylabel(sprintf('Onset\nLatency (ms)'));
    title(titles{k},'fontsize',14);
    xticks(t_inds);
    xticklabels({'Stand Pull','Walk Pull','Stand Rotate','Walk Rotate'});
end

toc;