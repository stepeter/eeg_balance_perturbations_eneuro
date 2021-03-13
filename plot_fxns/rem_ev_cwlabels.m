function EEG = rem_ev_cwlabels(EEG)
    % Remove CW/CCW suffixes for event labels 
    trigPrefix='M_on_';
    trigPrefix1='M_on_CW_';
    trigPrefix2='M_on_CCW_';
    %Remove CW/CCW from event labels
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