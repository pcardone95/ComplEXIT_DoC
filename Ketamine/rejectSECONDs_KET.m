function rejectSECONDs_KET(basename, filepath)
EEG = pop_loadset('filename',basename,'filepath',filepath);
% Find bad trials
badtrials = find(strcmp({EEG.epoch.dose}, 'SEC'));

EEG = pop_select(EEG, 'notrial', badtrials);
if isfield(EEG,'rejepoch')
    EEG.rejepoch = [EEG.rejepoch badtrials];
else
    EEG.rejepoch = badtrials;
end
EEG.SECONDs_rejected = 1;
pop_saveset(EEG,'filename',basename,'filepath',filepath);
