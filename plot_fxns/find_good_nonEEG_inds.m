function good_inds = find_good_nonEEG_inds(EEG, nb_vicon_chans)
    % Finds nonzero indices for EMG/mocap/load cell channels
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