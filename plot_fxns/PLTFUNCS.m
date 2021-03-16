classdef PLTFUNCS
    % Plotting functions for recreating eNeuro figures
    methods(Static)
        function start_eeglab(eeglab_pth)
            % Run EEGLAB and close all windows
            curr_dir = pwd;
            run([eeglab_pth '/eeglab.m']);
            cd(curr_dir); close all;
        end
        
        function dipfits = load_ica_structs(ica_struct_files)
            % Load ICA_struct files
            dipfits = {};
            for i=1:length(ica_struct_files)
                load(['data_release_ICA_structs/' ica_struct_files(i).name]);
                good_comps = ICA_STRUCT.good_comps.brain;
                ICA_STRUCT.dipfit.model = ICA_STRUCT.dipfit.model(good_comps);
                dipfits{i} = ICA_STRUCT.dipfit;
            end
        end
        
        function [compsall,centsall,colorsall] = aggregate_comps(clusters_to_plot,...
                                                                 dipfits,...
                                                                 cluster,...
                                                                 colors2use)
            % Aggregate IC dipole information based on DIPFIT dipole
            % fitting
            k1=1; k2=1; ct=1;
            compsall=[]; centsall=[]; dipsizes = [];
            for i=1:length(clusters_to_plot)
                k2=k2+1;
                for j=1:length(cluster(clusters_to_plot(i)).comps)
                    cls_set_i = cluster(clusters_to_plot(i)).sets(1,j);
                    comp = cluster(clusters_to_plot(i)).comps(j);
                    cluster_dip_models(1,j).posxyz = dipfits{cls_set_i}.model(comp).posxyz;
                    cluster_dip_models(1,j).momxyz = dipfits{cls_set_i}.model(comp).momxyz;
                    cluster_dip_models(1,j).rv = dipfits{cls_set_i}.model(comp).rv;
                    colorsall{ct} = colors2use{i};
                    ct= ct+1; k1=k1+1;
                end
                compsall=[compsall,cluster_dip_models(1,:)];
                centsall = [centsall, EEGFUNCS.computecentroid(cluster_dip_models)];
                dipsizes = [dipsizes 25*ones(size(1,length(cluster_dip_models(1,:)))) 40];

                clear cluster_dip_models colors1cls
            end
        end
        
        function plot_dipoles(compsall,centsall,dipfits,eeglab_pth,...
                              colors2use,colorsall)
            % Plot IC dipoles, using different colors for each cluster
            views = [0 0 1; 1 0 0; 0 -1 0];
            dipfit_pth = [eeglab_pth '/plugins/dipfit2.3/'];
            hdmfile = [dipfit_pth 'standard_BEM/standard_vol.mat'];
            mrifile = [dipfit_pth 'standard_BEM/standard_mri.mat'];
            figure;
            for i=1:size(views,1)
                subplot(2,3,i);
                dipplot(compsall,'spheres','on','dipolelength',0,...
                    'dipolesize',20,'mri',mrifile,'meshdata',hdmfile,...
                    'coordformat',dipfits{1}.coordformat,...
                    'color',colorsall,'gui','off');
                view(views(i,:));
            end

            % Plot centroids
            for i=1:size(views,1)
                subplot(2,3,i+size(views,1));
                dipplot(centsall,'spheres','on','dipolelength',0,...
                    'dipolesize',50,'mri',mrifile,'meshdata',hdmfile,...
                    'coordformat',dipfits{1}.coordformat,...
                    'color',colors2use,'gui','off');
                view(views(i,:));
            end
        end
        
        function plot_marker_sd(marker_inds,sds,n_subjs)
            % Plot head/sacrum standard deviation for each condition
            figure; title_labs = {'Sacrum','Head'};
            cols = {'r','m','b','c'};
            for j=1:length(marker_inds)
                subplot(2,1,j); hold on;
                for k=1:length(cols)
                    bar(k,mean(sds(:,k,j),1),cols{k});
                    errorbar(k,mean(sds(:,k,j),1),...
                             std(sds(:,k,j),1)/sqrt(n_subjs),'k',...
                             'linewidth',1.5);
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
        end
        
        function plot_emg_evs(ev_traces_norm,EEG_tmp,emg_marks,n_subjs)
            % Plot lower leg EMG event-related activity
            figure; cols = {'r','m','b','c'};
            plt_order = [1,3,5,7,2,4,6,8];
            for j=1:length(emg_marks)
                subplot(4,2,plt_order(j)); hold on;
                for k=1:length(cols)
                    errorbar(EEG_tmp.times/1000,...
                             mean(squeeze(ev_traces_norm(j,:,k,:)),1),...
                             squeeze(std(ev_traces_norm(j,:,k,:),0,2))/sqrt(n_subjs),...
                             cols{k});
                end
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
        end
        
        function plot_mocap_evs(ev_traces,marker_inds,EEG_tmp,n_subjs)
            % Plot head/sacrum event-related activity
            figure; title_labs = {'Sacrum','Head'};
            cols = {'r','m','b','c'};
            for j=1:length(marker_inds)
                subplot(2,1,j); hold on;
                for k=1:length(cols)
                    errorbar(EEG_tmp.times/1000,mean(squeeze(ev_traces(j,:,k,:)),1),...
                             squeeze(std(ev_traces(j,:,k,:),0,2))/sqrt(n_subjs),...
                             cols{k},'linewidth',1.5);
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
        end
        
        function plot_early_late_behavior(all_vals)
            % Plot pull force, peak EMG, and sacrum/head SD for the
            % first/last minute of each condition
            titles = {'Pull force','Peak EMG','Sacrum SD','Head SD'};
            ylab = {'Force (N)','Normalized EMG (%)','Marker SD (mm)',...
                    'Marker SD (mm)'};
            cols = {'r','r','m','m','b','b','c','c'};
            ylims = {[0 20],[-2 25],[0 80],[0 80]};

            figure; t_inds = [2 3 5 6 8 9 11 12];
            for k=1:size(all_vals,1)
                subplot(2,2,k);
                q = 0; hold on;
                for i=1:size(all_vals,4)
                    for j=1:size(all_vals,2)
                        q = q + 1;
                        if (k>1) || (i<3)
                            bar(t_inds(q),mean(all_vals(k,j,:,i)),cols{q});
                            errorbar(t_inds(q),mean(all_vals(k,j,:,i)),...
                                     std(all_vals(k,j,:,i))/sqrt(size(all_vals,3)),...
                                     'k','linewidth',1.5);
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
        end
        
        function plot_neck_specs(erps_all, specidx, freqs, ch_labels)
            % Plot neck EMG power spectra
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
                    plot(freqs,squeeze(mean(specs_all(:,ev,ch,:),1)),...
                         cols{ev},'linewidth',1.5);
                end
                hold off;
                title(ch_labels{ch},'FontSize',16);
                if ch==1
                    ylabel(sprintf('Power\n(dB)'),'Rotation',0,...
                           'fontsize',16,'fontweight','bold');
                end
                xlabel('Frequency (Hz)','Fontsize',16,'fontweight','bold');
                xlim([4,100]); ylim([35 65]);
                set(gca,'Fontsize',14,'fontweight','bold');
                set(gca,'box','off');
                set(gca,'Xscale','log');
                set(gca,'XMinorTick','Off')
                xticks([4,8,13,30,50,100]);
            end
        end
        
        function plot_ersps(erps_all, erps_all_masked, times, freqs,...
                            ev_labels, tlim, flim, clim, nrows)
            % Plot EEG or neck EMG time-frequency activity (ERSPs)
            figure; q = 0;
            if nrows==1
                set(gcf,'Position',[1 1 1102 250]);
            elseif nrows==2
                set(gcf,'Position',[1 1 1102 317]);
            end
            for ch=1:size(erps_all,3)
                for ev=1:size(erps_all,2)
                    q = q + 1;
                    subplot(nrows,4,q);
                    if ev<3
                        verts = [0 500];
                    else
                        verts = [0 1000];
                    end
                    tftopo(squeeze(erps_all_masked(ev,ch,:,:)),times,freqs,...
                           'vert',verts,'limits',[tlim flim clim],...
                           'logfreq','native');

                    if ch == 1
                        title(ev_labels{ev},'FontSize',16);
                    else
                        title('');
                    end
                    set(gca,'YTick',log([4.01,8,13,30,50,99.4843]));
                    if ev==1
                        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
                        ylabel(sprintf('Frequency\n(Hz)'),'Rotation',0,...
                                       'fontsize',16,'fontweight','bold');
                    else
                        set(gca,'YTickLabel',{'','','','','','',''});
                        ylabel('');
                    end
                    if ch==nrows
                        xlabel('Time (sec)','Fontsize',16,'fontweight','bold');
                        set(gca,'XTickLabel',{'-0.5','0','0.5','1'},'Fontsize',12);
                    else
                        xlabel('');
                        set(gca,'XTickLabel',{'','','',''});
                    end
                    set(gca,'Fontsize',14,'fontweight','bold');
                    set(gca,'box','off');
                end
            end
        end
        
        function plot_eeg_specs(erps_all, specidx, freqs, cluster2plt)
            % Plot EEG power spectra
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
        end
        
        function plot_eeg_ersp_onsets(t_onsets_all, fInds)
            % Plot EEG ERSP synchronization and descynchronization onset
            % timings for each condition
            figure; titles = {'Synchronization','Desynchronization'};
            cols = {'r','m','b','c'};

            t_inds = [1 2 3 4];
            plt_ord = [3 4 1 2];
            for k=1:length(fInds)
                subplot(2,1,k); hold on;
                for i=1:length(cols)
                    bar(t_inds(i),mean(t_onsets_all(k,plt_ord(i),:)),cols{i});
                    errorbar(t_inds(i),mean(t_onsets_all(k,plt_ord(i),:)),...
                             std(t_onsets_all(k,plt_ord(i),:))/sqrt(size(t_onsets_all,3)),'k');
                end
                hold off;
                ylim([0,300]);
                ylabel(sprintf('Onset\nLatency (ms)'));
                title(titles{k},'fontsize',14);
                xticks(t_inds);
                xticklabels({'Stand Pull','Walk Pull','Stand Rotate','Walk Rotate'});
            end
        end
        
    end
end