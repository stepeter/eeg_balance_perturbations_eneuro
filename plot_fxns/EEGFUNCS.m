classdef EEGFUNCS
    % Analysis functions for recreating eNeuro figures
    methods(Static)
        function dipole = computecentroid(alldipoles)
            % Compute dipole cluster centroid position
            len = length(alldipoles);
            dipole.posxyz = [ 0 0 0 ];
            dipole.momxyz = [ 0 0 0 ];
            dipole.rv = 0;
            count = 0;
            warningon = 1;
            for k = 1:len
                if size(alldipoles(k).posxyz,1) == 2
                    if all(alldipoles(k).posxyz(2,:) == [ 0 0 0 ])
                        alldipoles(k).posxyz(2,:) = [];
                        alldipoles(k).momxyz(2,:) = [];
                    end
                end
                if ~isempty(alldipoles(k).posxyz)
                    dipole.posxyz = dipole.posxyz + mean(alldipoles(k).posxyz,1);
                    dipole.momxyz = dipole.momxyz + mean(alldipoles(k).momxyz,1);
                    dipole.rv     = dipole.rv     + alldipoles(k).rv;
                    count = count+1;
                elseif warningon
                    disp('Some components do not have dipole information');
                    warningon = 0;
                end
            end
            dipole.posxyz = dipole.posxyz/count;
            dipole.momxyz = dipole.momxyz/count;
            dipole.rv     = dipole.rv/count;
            if isfield(alldipoles, 'maxr')
                dipole.maxr = alldipoles(1).max_r;
            end
        end

        function boundary_latencies2 = get_boundary_evs(EEG)
            % Identify start/end times of each condition (identified by
            % boundary events)
            boundary_latencies=[];
            for j=1:length(EEG.event)
                if strcmp(EEG.event(j).type,'boundary')
                    boundary_latencies=[boundary_latencies EEG.event(j).latency];
                end
            end
            boundary_latencies2=int32([1 boundary_latencies EEG.pnts]);
        end
        
        function marker_ind = find_chan_ind(EEG, mark_label)
            % Find the index of a channel in EEG variable based on its
            % label
            for j=1:length(EEG.chanlocs)
                if strcmp(EEG.chanlocs(j).labels,mark_label)
                    marker_ind = j;
                end
            end
        end
        
        function sds = compute_marker_sd(marker_inds,num_conds,good_inds,...
                                         boundary_latencies2,sds,EEG,sbj_ind)
            % Compute marker standard deviation for each condition
            % (removing zeroes for faster computation)
            for j=1:length(marker_inds)
                for k=1:num_conds
                    inds_use = good_inds((good_inds>=boundary_latencies2(k)) ...
                                    & (good_inds<=boundary_latencies2(k+1)));
                    dat_val = EEG.data(marker_inds(j),inds_use);
                    dat_val = dat_val(dat_val~=0);
                    dat_val = dat_val(~isnan(dat_val));
                    if ~isempty(dat_val)
                        sds(sbj_ind,k,j) = std(dat_val);
                    end
                end
            end
        end
        
        function good_inds = find_good_nonEEG_inds(EEG, nb_vicon_chans)
            % Find nonzero timepoints for EMG/mocap/load cell channels
            if nargin == 1
                nb_vicon_chans = 16;
            end

            if ndims(EEG.data)==2
                bad_inds = find(max(abs(EEG.data((end-nb_vicon_chans+1):end,:)),[],1)==0);
                good_inds = setdiff(1:size(EEG.data,2),bad_inds);
            elseif ndims(EEG.data)==3
                bad_inds = [];
                for i=1:EEG.trials
                    tmp_val = max(abs(EEG.data((end-nb_vicon_chans+1):end,:,i)));
                    if any(tmp_val == 0)
                        bad_inds = [bad_inds i];
                    end
                end
                good_inds = setdiff(1:size(EEG.data,3),bad_inds);
            end
        end
        
        function EEG = filt_data(EEG, chan_inds, freq_cutoff, filter_type)
            % Filter data using 4th order butterworth filter (filter_type
            % specifies either 'low' or 'high' pass filtering)
            [b,a]=butter(4,freq_cutoff/(EEG.srate/2),filter_type);
            for k=1:length(chan_inds)
                EEG.data(chan_inds(k),:) = filtfilt(b,a,double(EEG.data(chan_inds(k),:)));
            end
        end
        
        function EEG = rem_ev_cwlabels(EEG)
            % Remove CW/CCW suffixes for rotation perturbation event labels 
            trigPrefix='M_on_';
            trigPrefix1='M_on_CW_';
            trigPrefix2='M_on_CCW_';
            for i=1:length({EEG.event(1:end).type})
                ccwEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigPrefix2 '*']),'once');
                cwEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigPrefix1 '*']),'once');
                if ~isempty(ccwEv)
                    EEG.event(1,i).type=[trigPrefix EEG.event(1,i).type((length(trigPrefix)+5):end)];
                elseif ~isempty(cwEv)
                    EEG.event(1,i).type=[trigPrefix EEG.event(1,i).type((length(trigPrefix)+4):end)];
                end
            end
        end
        
        function EEG = rem_ev_cwlabels_pulls(EEG)
            % Remove L/R prefixes for pull perturbation event labels 
            trigSuffix='_pull_';
            for i=1:length({EEG.event(1:end).type})
                pullEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigSuffix '*']),'once');
                if ~isempty(pullEv)
                    EEG.event(1,i).type=EEG.event(1,i).type(3:end);
                end
            end
        end
        
        function [EEG_early,EEG_late] = early_late_eeg(EEG,boundary_latencies2,...
                                                       emg_inds,marker_inds,...
                                                       cond_ind, sbj_ind)
            % Select the first and last minute of data for each condition,
            % after removing timepoints with no EMG/mocap/load cell data
            EEG_tmp = pop_select(EEG,'time',[boundary_latencies2(cond_ind) ...
                                 boundary_latencies2(cond_ind+1)]/EEG.srate);
            if sbj_ind==25
                EEG_tmp2 = pop_select(EEG_tmp,'channel',emg_inds(1));
            else
                EEG_tmp2 = pop_select(EEG_tmp,'channel',marker_inds(1));
            end
            good_inds_tmp = EEGFUNCS.find_good_nonEEG_inds(EEG_tmp2,1);
            EEG_early = pop_select(EEG_tmp,'time',[good_inds_tmp(1) ...
                                   good_inds_tmp(1)+EEG.srate*60]/EEG.srate);
            EEG_late = pop_select(EEG_tmp,'time',[good_inds_tmp(end)-EEG.srate*60 ...
                                  good_inds_tmp(end)]/EEG.srate);
        end
        
        function [dat_out,times,freqs] = compute_ersps(EEG_tmp, chan_ind)
            % Compute ERSP without generating a plot or saving an
            % intermediate file
            tlimits = [EEG_tmp.xmin, EEG_tmp.xmax]*1000;
            pointrange1 = round(max((tlimits(1)/1000-EEG_tmp.xmin)*EEG_tmp.srate, 1));
            pointrange2 = round(min((tlimits(2)/1000-EEG_tmp.xmin)*EEG_tmp.srate, EEG_tmp.pnts));
            pointrange = [pointrange1:pointrange2];
            tmpsig = EEG_tmp.data(chan_ind,pointrange,:);
            tmpsig = reshape(tmpsig, 1, length(pointrange)*size(EEG_tmp.data,3));
            [dat_out,~,~,times,freqs,~,~] = mod_newtimef(tmpsig, length(pointrange),...
                                                [tlimits(1) tlimits(2)],...
                                                EEG_tmp.srate,[3 0.9],...
                                                'baseline',nan,'freqs',[3 256],...
                                                'nfreqs',100,'freqscale','log',...
                                                'plotersp','off','plotitc','off');
        end
        
        function erps_all_masked = rem_ersp_baseline(erps_all, baseidx,...
                                                     times, alpha)
            % Subtract ERSP baseline and perform statistical masking
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
                    curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & ...
                                    curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                    erps_all_masked(ev,ch,:,:) = curr_maskedersp;
                end
            end
        end
        
        function [EEG, comp_ind] = apply_ica_to_eeg(EEG, comp, icachansind)
            % Apply weights from ICA_STRUCT to the EEG dataset
            good_comps = EEG.etc.good_comps;
            comp_ind = good_comps(comp);
            EEG.icaweights = EEG.etc.icaweights;
            EEG.icasphere = EEG.etc.icasphere;
            EEG.icachansind = icachansind;
            
            % Compute ICA activations
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.data(1:size(EEG.icaact,1),:) = EEG.icaact;
        end
        
        function erps_all_masked = rem_eeg_ersp_baseline(erps_all, baseidx,...
                                                         times, alpha)
            % Subtract baseline for EEG ERSP data (different shape than
            % neck EMG ERSP data)
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
                    curr_maskedersp(curr_ersp > repmat(pboot(:,1),...
                        [1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),...
                        [1 size(curr_ersp,2)])) = 0;
                    erps_all_masked(ev,ch,:,:) = curr_maskedersp;
                end
            end
        end
        
        function [t_onsets_all,fInds] = compute_ersp_onsets(erps_all, tLims,...
                                                            fLims_sync,...
                                                            fLims_desync,...
                                                            threshVal,...
                                                            baseidx, times,...
                                                            freqs,...
                                                            n_comps)
            % Compute ERSP onset timings
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
                            num_pixels = cellfun(@numel,CC.PixelIdxList); % num pixels per region
                            pix_list = regionprops(dip_data,'PixelIdxList',...
                                            'PixelList'); % x y pixel coordinates

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
        end
        
    end
end