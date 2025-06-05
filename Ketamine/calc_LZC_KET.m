function calc_LZC_KET(basename, filepath)
%% script for calculating LZC in KET data

EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);

results = struct();

% LZC calculation all over
data_binar = zeros(size(EEG.data(:,:,:)));
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};

% Preallocate
strings = cell(size(EEG.data,1),size( EEG.data,3)); % electrode x trial
C = zeros(size(EEG.data,1),size( EEG.data,3));
H = cell(size(EEG.data,1),size( EEG.data,3));

% Calculate
for electi = 1:size( EEG.data,1)
    for triali = 1:size( EEG.data,3)
        data_binar(electi,:,triali) = abs(hilbert(EEG.data(electi,:,triali)));
        data_binar(electi,:,triali) = ...
            data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
        strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
        [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
    end
end

results(1).C = C;
results(1).H = H;
results(1).C_mean = mean(C,2);
results(1).concentration = 'all';
results(1).ntrials = EEG.trials;

for dosei = 1:length(dose_l) 
    dose_s = dose_l{dosei};
    dose_trials = find(strcmp({EEG.epoch.dose},dose_s));
    
    % Take only the data with the same concentration
    EEG_temp = pop_select(EEG, 'trial', dose_trials);
    pop_saveset(EEG_temp,'filename',sprintf('%s_%s.set', basename, dose_s),...
        'filepath',filepath);
  
    results(dosei+1).C = C(:,dose_trials);
    results(dosei+1).H = H(:,dose_trials);
    results(dosei+1).C_mean = mean(C(:,dose_trials),2);
    results(dosei+1).concentration = dose_s;
    results(dosei+1).ntrials = EEG_temp.trials;
end


% Save file
cd(filepath)
save([basename '_LZC.mat'], 'results')

% Demographics
labels = strsplit(basename ,'_');

% Write results in table
results_table_PerConcentration = table(repmat(sprintf('KET0%s',labels{1}(3)), 127*5,1),...
    repmat(str2double(labels{2}(2)),127*5,1),...
    sort(repmat([0.15 0.30, 0.45 0.60 0.75]',127,1)), repmat({EEG.chanlocs.labels}', 5, 1),...
    [results(2).C_mean; results(3).C_mean; results(4).C_mean; results(5).C_mean;...
    results(6).C_mean], 'VariableNames', {'Subject', 'Session', 'Concentration', 'Electrode', 'LZC'});
writetable(results_table_PerConcentration, [basename '_LZC_PerConcentration'])

results_table_all = table(repmat(sprintf('KET0%s',labels{1}(3)), 127,1),...
    repmat(str2double(labels{2}(2)),127,1), repmat('all',127,1), {EEG.chanlocs.labels}',...
    results(1).C_mean, 'VariableNames', {'Subject', 'Session','Concentration', 'Electrode', 'LZC'});
writetable(results_table_all, [basename '_LZC_all'])
