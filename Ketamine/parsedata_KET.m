function parsedata_KET(basename, filepath)
% Parse data into doses - important for later.
EEG = pop_loadset('filename',[basename '_epochs.set'],'filepath',filepath);

EEG.epoch(1).dose =[];
% Assign "dose" to epochs
% index of the event of the dose i.e. 0.15
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};

for dosei = 1:length(dose_l)
    dose_s = dose_l{dosei};
    index_dose = find(strcmp({EEG.event.type}, dose_s));
    if isempty(index_dose) % Not in the recording for any reasons - maybe should never happen
        warning('Following concentration missing: %s', dose_s)
        continue;
    end
    % index of the event of the event when dose beceome i.e. 0.15
    index_trial = find(cellfun(@(x) sum(x==index_dose), {EEG.epoch.event}, 'UniformOutput', 1));
    
    EEG.epoch(index_trial).dose =[];
    temp = repmat({dose_s}, 1, EEG.trials-index_trial); % Put conc from that trial onward
    [EEG.epoch(index_trial+1:end).dose] =deal(temp{:});
end

% NOW FOR SECONDS
sec_l = {{'SEC30B', 'SEC30E'}, {'SEC60B', 'SEC60E'}, {'SEC90B', 'SEC90E'}};

for seci = 1:length(sec_l)
    [sec_beg, sec_end] = sec_l{seci}{:};
    
    %Beginning
    index_sec_b = find(strcmp({EEG.event.type}, sec_beg));
    if isempty(index_sec_b) % Not in the recording for any reasons - maybe should never happen
        warning('%s not found', sec_beg)
        continue;
    end
    index_trial_b = find(cellfun(@(x) sum(x==index_sec_b), {EEG.epoch.event}, 'UniformOutput', 1));
    
    
    %End
    index_sec_e = find(strcmp({EEG.event.type}, sec_end));
    if isempty(index_sec_e) % after the end of the recording
        index_trial_e = EEG.trials;
    else
        index_trial_e = find(cellfun(@(x) sum(x==index_sec_e), {EEG.epoch.event}, 'UniformOutput', 1));
    end
    
    temp = repmat({'SEC'}, 1, index_trial_e-index_trial_b+1);
    [EEG.epoch(index_trial_b:index_trial_e).dose] =deal(temp{:});
end

% Eliminate everything after SEC90B
index_sec_b = find(strcmp({EEG.event.type}, 'SEC90B'));
if not(isempty(index_sec_b)) % Not in the recording for any reasons - maybe should never happen
    index_trial_b = find(cellfun(@(x) sum(x==index_sec_b), {EEG.epoch.event}, 'UniformOutput', 1));
    EEG = pop_select(EEG,'notrial', index_trial_b:EEG.trials);
end

% Eliminate everything before 0.15
index_dose = find(strcmp({EEG.event.type}, '0.15'));
if not(isempty(index_dose)) 
    index_trial = find(cellfun(@(x) sum(x==index_dose), {EEG.epoch.event}, 'UniformOutput', 1));
    if index_trial > 1
        EEG = pop_select(EEG,'notrial', 1:index_trial-1);
    end
end

pop_saveset(EEG,'filename',[basename '_epochs.set'],'filepath',filepath);