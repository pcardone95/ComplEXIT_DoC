%% Compute LZC over time for hypnosis
% Go to folder

pathwork = 'E:\BioWin';
cd(fullfile(pathwork)) % normally separate data and results

%% Compute LZC
results = struct();
name_fold_HYP = ls("*HYP.set");
name_fold_EO = ls("*EO1.set"); % If you have multiple sessions

if size(name_fold_HYP,1) ~= size(name_fold_EO,1)
    % to understand how to do it more efficiently
    name_fold_EO(16,:) = []; % for sub 27 no 
end

% sub to take out
take_out = [3, 9, 10, 17];
name_fold_EO(take_out,:) = [];
name_fold_HYP(take_out,:) = [];

if strcmp(name_fold_EO(:,1:6), name_fold_HYP(:,1:6))
    % to think
    disp("All good!")
end

condition = {'EO', 'HYP'};


%% LZC
results_table = table();

% Should be a loop, but no time right now
for foldi = 1:size(name_fold_HYP,1) % PRoblem with S3 num 14
    % EO
    EEG = pop_loadset(erase(name_fold_EO(foldi,:),' '), cd); %import data
    
    data = EEG.data(:,:,:); % channel x time x trial
    data_binar = zeros(size(data));
    
    % Preallocate
    strings = cell(size(data,1),size(data,3)); % electrode x trial
    C = zeros(size(data,1),size(data,3));
    H = cell(size(data,1),size(data,3));
    
    % Calcule - from Figure 3 (Paper 2)
    
    for electi = 1:size(data,1)
        for triali = 1:size(data,3)
            data_binar(electi,:,triali) = abs(hilbert(data(electi,:,triali)));
            data_binar(electi,:,triali) = ...
                data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
            strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
            [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
        end
    end
  
    results(foldi).C_EO = C; % results(folderi).H = H; %results(folderi).strings = strings;
    results(foldi).label = erase(name_fold_EO(foldi,4:6),'_');
    % Write in table
    temp = table(repmat({erase(name_fold_EO(foldi,4:6),'_')},127,1),... % Sub
        repmat({'EO'},127,1),... % Condition
        {EEG.chanlocs.labels}', ... % channels
        mean(C,2)); % LZC
    results_table = [results_table; temp];
   
    % Hypnosis
    EEG = pop_loadset(erase(name_fold_HYP(foldi,:), ' '), cd); %import data
    
    data = EEG.data(:,:,:); % channel x time x trial
    data_binar = zeros(size(data));
    
    % Preallocate
    strings = cell(size(data,1),size(data,3)); % electrode x trial
    C = zeros(size(data,1),size(data,3));
    H = cell(size(data,1),size(data,3));
    
    % Calcule - from Figure 3 (Paper 2)
    
    for electi = 1:size(data,1)
        for triali = 1:size(data,3)
            data_binar(electi,:,triali) = abs(hilbert(data(electi,:,triali)));
            data_binar(electi,:,triali) = ...
                data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
            strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
            [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
        end
    end
  
    results(foldi).C_HYP = C; % results(folderi).H = H; %results(folderi).strings = strings;
    % Write in table
    temp = table(repmat({erase(name_fold_HYP(foldi,4:6),'_')},127,1),... % Sub
        repmat({'HYP'},127,1),... % Condition
        {EEG.chanlocs.labels}', ... % channels
        mean(C,2)); % LZC
    results_table = [results_table; temp];
end
%% Save file

cd('C:\Users\Paolo\Documents\GitHub\ComplEXIT_DoC')
save('Hypnosis_LZC_complete_v2.mat', 'results')

results_table.Properties.VariableNames = {'Subject', 'Condition', 'Electrode' 'LZC'};
writetable(results_table, 'results_LZC_hypnosis_complete_v2')