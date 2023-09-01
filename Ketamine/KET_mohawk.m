%% Ketamine pipeline
% Most of the scripts will be taken from the Mohawk, with some small
% modification to take out the bad epochs

% %% Import data
% % Organizing paths
% mohawkpath = ...
% 'C:\Users\Paolo\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\MOHAWK\Volumes\bigdisk\Work\MOHAWK';
% 
% filepath = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';
% cd(filepath)
% % Select file
% [filename,filepath] = uigetfile('*.vhdr','MOHAWK - Select file to process',filepath);
% if filename == 0
%     return
% end
% [filename,ext] = strtok(filename,'.');
% 
% 
% % Specify a name for the final file
% answers = inputdlg2({'Specify dataset base name:'},'MOHAWK Dataset',1,{filename});
% 
% if isempty(answers)
%     return
% elseif isempty(answers{1})
%     error('Basename cannot be empty.');
% else
%     basename = answers{1};
% end
% 
% % fid = fopen(loadpathsloc,'w');
% % fprintf(fid,'filepath=''%s'';',filepath);
% % fclose(fid);
% 
% if ~exist(sprintf('%s%sfigures',filepath,filesep),'dir')
%     mkdir(sprintf('%s%sfigures',filepath,filesep));
% end
% 
% cur_wd = pwd;
% cd(mohawkpath);

%% Data import
% Already done by subject with hom_event_K

%% Epoch data
fprintf('\n*** EPOCHING DATA ***\n');
epochdata(basename, filepath);

%% Parse data
fprintf('\n*** PARSE DATA IN DIFFERENT CONCENTRATION ***\n');
parsedata_KET(basename, filepath);
% % Parse data into doses - important for later.
% EEG = pop_loadset('filename',[basename '_epochs.set'],'filepath',filepath);
% 
% EEG.epoch(1).dose =[];
% % Assign "dose" to epochs
% % index of the event of the dose i.e. 0.15
% dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};
% 
% for dosei = 1:length(dose_l)
%     dose_s = dose_l{dosei};
%     index_dose = find(strcmp({EEG.event.type}, dose_s));
%     
%     % index of the event of the event when dose beceome i.e. 0.15
%     index_trial = find(cellfun(@(x) sum(x==index_dose), {EEG.epoch.event}, 'UniformOutput', 1));
%     
%     EEG.epoch(index_trial).dose =[];
%     temp = repmat({dose_s}, 1, EEG.trials-index_trial);
%     [EEG.epoch(index_trial+1:end).dose] =deal(temp{:});
% end
% 
% % NOW FOR SECONDS
% sec_l = {{'SEC30B', 'SEC30E'}, {'SEC60B', 'SEC60E'}, {'SEC90B', 'SEC90E'}};
% 
% for seci = 1:length(sec_l)
%     [sec_beg, sec_end] = sec_l{seci}{:};
%     
%     %Beginning
%     index_sec_b = find(strcmp({EEG.event.type}, sec_beg));
%     if isempty(index_sec_b) % Not in the recording for any reasons - maybe should never happen
%         warning('%s not found', sec_beg)
%         continue;
%     end
%     index_trial_b = find(cellfun(@(x) sum(x==index_sec_b), {EEG.epoch.event}, 'UniformOutput', 1));
%     
%     
%     %End
%     index_sec_e = find(strcmp({EEG.event.type}, sec_end));
%     if isempty(index_sec_e) % after the end of the recording
%         index_trial_e = EEG.trials;
%     else
%         index_trial_e = find(cellfun(@(x) sum(x==index_sec_e), {EEG.epoch.event}, 'UniformOutput', 1));
%     end
%     
%     temp = repmat({'SEC'}, 1, index_trial_e-index_trial_b+1);
%     [EEG.epoch(index_trial_b:index_trial_e).dose] =deal(temp{:});
% end
% 
% % Eliminate everything after SEC90B
% index_sec_b = find(strcmp({EEG.event.type}, 'SEC90B'));
% if not(isempty(index_sec_b)) % Not in the recording for any reasons - maybe should never happen
%     index_trial_b = find(cellfun(@(x) sum(x==index_sec_b), {EEG.epoch.event}, 'UniformOutput', 1));
%     EEG = pop_select(EEG,'notrial', index_trial_b:EEG.trials);
% end
% 
% % Eliminate everything before 0.15
% index_dose = find(strcmp({EEG.event.type}, '0.15'));
% if not(isempty(index_dose)) % Not in the recording for any reasons - maybe should never happen
%     index_trial = find(cellfun(@(x) sum(x==index_dose), {EEG.epoch.event}, 'UniformOutput', 1));
%     if index_trial > 1
%         EEG = pop_select(EEG,'notrial', 1:index_trial-1);
%     end
% end
% 
% pop_saveset(EEG,'filename',[basename '_epochs.set'],'filepath',filepath);

%% REJECT DATA DURING SECONDs
fprintf('\n DELETE TRIALS DURING SECONDs...\n');
rejectSECONDs_KET([basename '_epochs.set'], filepath)

%% CHANNEL AND TRIAL REJECTION
% First pass of quasi-automatic rejection of noisy channels and epochs based on variance
% thresholding
fprintf('\n*** SELECT BAD CHANNELS AND TRIALS ***\n');
rejartifacts([basename '_epochs'], 1, 4, 1, [], [], [], filepath);

%% ICA
fprintf('\n*** COMPUTING IC DECOMPOSITION ***\n');
% Run ICA decomposition with optional PCA pre-processing
computeic([basename '_epochs'], [], [], filepath)

%% MANUAL
fprintf('\n*** SELECT BAD ICs ***\n');
% Visually identify and reject noisy ICs, e.g., eye movements, muscle
% activity, etc.
rejectic(basename, 'prompt', 'off', 'filepath', filepath)


%% MANUAL
% Second and final pass of quasi-automatic rejection of noisy channels and epochs based on variance
% thresholding, to remove any remaining noisy data.
fprintf('\n*** SELECT ANY REMAINING BAD CHANNELS AND TRIALS AND INTERPOLATE ***\n');
rejartifacts([basename '_clean'], 2, 4, 0, [], 500, 250, filepath);
%% Manual
fprintf('\n*** REFERENCING DATA TO COMMON AVERAGE ***\n');
% re-reference data to common average for connectivity estimation.
rereference(basename, 1, [],[], filepath);

%% MEASURING LZC AND DIVIDE PER SESSION
fprintf('\n*** Measuring complexity and dividing in sessions ***\n');
calc_LZC_KET(basename, filepath) 
% % Only for here, make things easier
% EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);
% 
% results = struct();
% 
% % LZC calculation all over
% data_binar = zeros(size(EEG.data(:,:,:)));
% 
% % Preallocate
% strings = cell(size( EEG.data,1),size( EEG.data,3)); % electrode x trial
% C = zeros(size(EEG.data,1),size( EEG.data,3));
% H = cell(size(EEG.data,1),size( EEG.data,3));
% 
% % Calcule - from Figure 3 (Paper 2)
% for electi = 1:size( EEG.data,1)
%     for triali = 1:size( EEG.data,3)
%         data_binar(electi,:,triali) = abs(hilbert( EEG.data(electi,:,triali)));
%         data_binar(electi,:,triali) = ...
%             data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
%         strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
%         [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
%     end
% end
% 
% results(1).C = C;
% results(1).H = H;
% results(1).C_mean = mean(C,2);
% results(1).concentration = 'all';
% results(1).ntrials = EEG.trials;
% 
% for dosei = 1:length(dose_l) % -1 because last index is the end
%     dose_s = dose_l{dosei};
%     dose_trials = find(strcmp({EEG.epoch.dose},dose_s));
%     
%     % Take only the data with the same concentration
%     EEG_temp = pop_select(EEG, 'trial', dose_trials);
%     pop_saveset(EEG_temp,'filename',sprintf('%s_%s.set', basename, dose_s),...
%         'filepath',filepath);
%   
%     results(dosei+1).C = C(:,dose_trials);
%     results(dosei+1).H = H(:,dose_trials);
%     results(dosei+1).C_mean = mean(C(:,dose_trials),2);
%     results(dosei+1).concentration = dose_s;
%     results(dosei+1).ntrials = EEG_temp.trials;
% end
% 
% 
% % Save file
% cd(filepath)
% save([basename '_LZC.mat'], 'results')
% 
% % Write results in table
% results_table = table(sort(repmat([0.15, 0.30, 0.45,0.60,0.75]',127,1)), repmat({EEG.chanlocs.labels}', 4, 1),...
%     [results(1).C_mean; results(2).C_mean; results(3).C_mean; results(4).C_mean;...
%     results(5).C_mean], 'VariableNames', {'Concentration', 'Electrode', 'LZC'});
% writetable(results_table, [basename '_results_LZC'])

%% Divide per concentration
% fprintf('\n*** RETAINING 10 MINUTES (60 EPOCHS) OF DATA ***\n');
% % optionally fix number of epochs contributing to connectivity estimation.
% % 60 epochs below will effectively use 10 minutes of clean data.
% checktrials(basename,60,''); % Prob not need for this
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};
for dosei = 1:length(dose_l)
    dose_s = dose_l{dosei};
    
    fprintf('\n*** CALCULATING MULTITAPER SPECTRUM ***\n');
    % calculate power spectrum using the multi-taper method
    calcftspec(sprintf('%s_%s', basename, dose_s),filepath);
    
    fprintf('\n*** PLOTTING SPECTRUM ***\n');
    % visualise and save the power spectrum of all channels
    plotftspec(sprintf('%s_%s', basename, dose_s),[], filepath);
    
    fprintf('\n*** CALCULATING CONNECTIVITY ***\n');
    % estimate dwPLI connectivity between pairs of channels
    ftcoherence(sprintf('%s_%s', basename, dose_s),filepath);
    
    fprintf('\n*** CALCULATING GRAPH-THEORETIC NETWORK METRICS ***\n');
    % calculate graph theory metrics
    calcgraph(sprintf('%s_%s', basename, dose_s), 'filepath', filepath);
    
    % %% plot 3D connectivity topographs in the delta, theta and alpha bands.
    plothead(sprintf('%s_%s', basename, dose_s),1, 'filepath',filepath);
    plothead(sprintf('%s_%s', basename, dose_s),2, 'filepath',filepath);
    plothead(sprintf('%s_%s', basename, dose_s),3, 'filepath',filepath);
    %
    % %%
    % fprintf('\n*** PLOTTING METRICS ***\n');
    % plotmetric(basename,'participation coefficient',3,'ylabel','Network centrality')
    % plotbands(basename,'participation coefficient','title','Network centrality');
    %
end

%% Parse data for ECI
parse_data_4ECI_temp(basename, filepath)