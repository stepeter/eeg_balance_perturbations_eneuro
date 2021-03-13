function EEG = rem_ev_cwlabels_pulls(EEG)
    % Remove L/R prefixes for pull event labels 
    trigSuffix='_pull_';
    %Remove CW/CCW from event labels
    for i=1:length({EEG.event(1:end).type})
        pullEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigSuffix '*']),'once');
        if ~isempty(pullEv)
            EEG.event(1,i).type=EEG.event(1,i).type(3:end);
        end
    end

end